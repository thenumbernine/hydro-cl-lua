--[[
This is a Solver that uses a fluxBuffer
and that uses calcDerivFromFlux to update it
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local Solver = require 'solver.solver'

local FiniteVolumeSolver = class(Solver)

function FiniteVolumeSolver:createBuffers()
	FiniteVolumeSolver.super.createBuffers(self)
	self:clalloc('fluxBuf', self.volume * self.dim * ffi.sizeof(self.eqn.cons_t))
end

function FiniteVolumeSolver:getSolverCode()
	return table{
		FiniteVolumeSolver.super.getSolverCode(self),
		template(file['solver/calcDeriv.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function FiniteVolumeSolver:refreshSolverProgram()
	FiniteVolumeSolver.super.refreshSolverProgram(self)

	-- 'calcFlux' is usually provided, but the args vary, so I'll leave it to the subclass
	
	self.calcDerivFromFluxKernel = self.solverProgram:kernel'calcDerivFromFlux'
	self.calcDerivFromFluxKernel:setArg(1, self.fluxBuf)
end

function FiniteVolumeSolver:addConvertToTexs()
	FiniteVolumeSolver.super.addConvertToTexs(self)

	-- TODO add kernels for each side
	self:addConvertToTex{
		name = 'flux', 
		type = self.eqn.cons_t,
		varCodePrefix = [[
	const global ]]..self.eqn.cons_t..[[* flux = buf + indexInt;
]],
		vars = range(0,self.eqn.numIntStates-1):map(function(i)
			return {[tostring(i)] = '*value = flux->ptr['..i..'];'}
		end),
	}
end

return FiniteVolumeSolver
