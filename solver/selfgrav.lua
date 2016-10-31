local table = require 'ext.table'
local ffi = require 'ffi'
local class = require 'ext.class'
local SelfGravitationBehavior = function(parent)
	local template = class(parent)

	function template:createBuffers()
		template.super.createBuffers(self)
		self:clalloc('potentialBuf', self.volume * ffi.sizeof(self.app.real))
	end

	function template:addConvertToTexs()
		template.super.addConvertToTexs(self)
		self:addConvertToTex{
			name = 'potential',
			vars = {'0'},
			displayBodyCode = [[
		value = buf[index];
]],
		}
	end

	function template:getSolverCode()
		return table{
			template.super.getSolverCode(self),
			[[
__kernel void initPotential(
	__global real* potentialBuf,
	const __global cons_t* UBuf)
{
	SETBOUNDS(0,0);
	potentialBuf[index] = -UBuf[index].rho;
}
]],
			'#define gravitationalConstant 1',		-- 6.67384e-11 m^3 / (kg s^2)
			'#include "solver/selfgrav.cl"',
		}:concat'\n'
	end

	function template:refreshSolverProgram()
		template.super.refreshSolverProgram(self)

		self.initPotentialKernel = self.solverProgram:kernel('initPotential', self.potentialBuf, self.UBuf)

		self.solvePoissonKernel = self.solverProgram:kernel('solvePoisson', self.potentialBuf, self.UBuf)
		
		self.calcGravityDerivKernel = self.solverProgram:kernel'calcGravityDeriv'
		self.calcGravityDerivKernel:setArg(1, self.UBuf)
		self.calcGravityDerivKernel:setArg(2, self.potentialBuf)	
	end

	function template:refreshBoundaryProgram()
		template.super.refreshBoundaryProgram(self)

		self.potentialBoundaryProgram, self.potentialBoundaryKernel =
			self:createBoundaryProgramAndKernel{
				type = 'real',
				methods = table.map(self.boundaryMethods, function(v,k)
					return self.app.boundaryMethods[1+v], k
				end),
			}
		self.potentialBoundaryKernel:setArg(0, self.potentialBuf)
	end

	function template:resetState()
		template.super.resetState(self)

		self.app.cmds:enqueueNDRangeKernel{kernel=self.initPotentialKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		self:potentialBoundary()
		for i=1,20 do
			self.app.cmds:enqueueNDRangeKernel{kernel=self.solvePoissonKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
			self:potentialBoundary()
		end
		
		-- TODO
		-- add potential energy into total energy
		-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
	end

	function template:step(dt)
		template.super.step(self, dt)

		self.integrator:integrate(dt, function(derivBuf)
			for i=1,20 do
				self:potentialBoundary()
				self.app.cmds:enqueueNDRangeKernel{kernel=self.solvePoissonKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
			end
			
			self.calcGravityDerivKernel:setArg(0, derivBuf)
			self.app.cmds:enqueueNDRangeKernel{kernel=self.calcGravityDerivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		end)
	end

	function template:potentialBoundary()
		self:applyBoundaryToBuffer(self.potentialBoundaryKernel)
	end

	return template
end

return SelfGravitationBehavior 
