local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local PoissonSolver = require 'solver.poisson'


local RemoveDivergence = class(PoissonSolver)

function RemoveDivergence:getCodeParams()
	return {
		args = 'const __global cons_t* UBuf',
		calcRho = [[
	//TODO make this modular
	//4 pi G rho for gravity
	//div(E) for electromagnetism
	const __global cons_t* U = UBuf + index;
	real divE = .5 * ((U[stepSize.x].epsE.x - U[-stepSize.x].epsE.x) / grid_dx0
					+ (U[stepSize.y].epsE.y - U[-stepSize.y].epsE.y) / grid_dx1,
					+ (U[stepSize.z].epsE.z - U[-stepSize.z].epsE.z) / grid_dx2);
	rho = divE;	//times 4 pi?
]],
	}
end


local GravityPotential = class(PoissonSolver)

GravityPotential.gravityConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

function GravityPotential:getCodeParams()
	return {
		args = 'const __global cons_t* UBuf',
		calcRho = '#define gravitationalConstant '..require 'clnumber'(self.gravityConstant)..'\n'..[[
	//TODO make this modular
	//4 pi G rho for gravity
	//div(E) for electromagnetism
	const __global cons_t* U = UBuf + index;
	rho = 4. * M_PI * gravitationalConstant * U->rho;
]],
	}
end

local selfGravBehavior = function(field, poissonClass)
	return function(parent)
		local template = class(parent)

		function template:init(args)
			self.useGravity = not not args.useGravity
			
			-- TODO in refreshGrid
			if not self.useGravity then
				self[field] = poissonClass(self)
			end

			-- init is gonna call
			template.super.init(self, args)
		end

		function template:createBuffers()
			template.super.createBuffers(self)
			self[field]:createBuffers()
		end

		function template:addConvertToTexs()
			template.super.addConvertToTexs(self)
			self[field]:addConvertToTexs()
		end

		function template:getSolverCode()
			return table{
				template.super.getSolverCode(self),
				self[field]:getSolverCode(),
			}:concat'\n'
		end

		function template:refreshSolverProgram()
			template.super.refreshSolverProgram(self)
			self[field]:refreshSolverProgram()
		end

		function template:refreshBoundaryProgram()
			template.super.refreshBoundaryProgram(self)
			self[field]:refreshBoundaryProgram()
		end

		function template:resetState()
			template.super.resetState(self)
			self[field]:resetState()
		end

		function template:step(dt)
			template.super.step(self, dt)
			self[field]:step(dt)	
		end

		function template:potentialBoundary()
			self:applyBoundaryToBuffer(self[field].potentialBoundaryKernel)
		end

		return template
	end
end

return selfGravBehavior('gravityPoisson', GravityPotential)
