local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local real = require 'hydro.real'
local FiniteVolumeSolver = require 'hydro.solver.fvsolver'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local WENO = class(FiniteVolumeSolver)
WENO.name = 'WENO'

WENO.stencilSize = 3

WENO.wenoMethod = '1996 Jiang Shu'	-- (WENO-JS)
--WENO.wenoMethod = '2008 Borges'	-- (WENO-Z)
--WENO.wenoMethod = '2010 Shen Zha'	-- (WENO-BS?)

WENO.fluxMethod = 'Lax-Friedrichs'
--WENO.fluxMethod = 'Marquina'
--WENO.fluxMethod = 'Roe'	-- isn't as accurate 

--[[
args:
	wenoMethod
	fluxMethod
	order (5 or 7 only)
--]]
function WENO:init(args)
	self.wenoMethod = args.wenoMethod
	self.fluxMethod = args.fluxMethod
	
	local order = args.order
	if order then
		self.stencilSize = (order+1)/2
	end
	
	self.numGhost = self.stencilSize
	WENO.super.init(self, args)
end

-- don't let the parent create a flux object or require a flux arg
function WENO:createFlux()
	self.flux = {
		getSolverCode = function() return '' end,
	}
end

-- TODO find what intermediate values to buffer for perf increase
function WENO:createBuffers()
	WENO.super.createBuffers(self)

	-- flux of cell-centered state values
--	self:clalloc('fluxCellBuf', self.eqn.cons_t, self.numCells * self.dim)
end

function WENO:initCodeModules()
	WENO.super.initCodeModules(self)
	
	self.modules:add{
		name = 'WENO.calcFlux',
		code = template(file['hydro/solver/weno.cl'], {
			solver = self,
			eqn = self.eqn,
			clnumber = require 'cl.obj.number',
		}),
	}
	self.solverModulesEnabled['WENO.calcFlux'] = true
end

-- all these are found eqn's cl code
function WENO:refreshSolverProgram()
	WENO.super.refreshSolverProgram(self)
	
--	self.calcCellFluxKernelObj = self.solverProgramObj:kernel'calcCellFlux'

	self.calcFluxKernelObj = self.solverProgramObj:kernel'calcFlux'
	self.calcFluxKernelObj.obj:setArg(2, self.fluxBuf)
end

-- NOTICE this adds the contents of derivBufObj and does not clear it
function WENO:calcDeriv(derivBufObj, dt)
	local dtArg = real(dt)
	
	if self.usePLM then
		self.calcLRKernelObj(self.solverBuf, self.cellBuf, self.UBuf, self.UBuf, dtArg)
	end

--[[
	self.calcCellFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcCellFluxKernelObj.obj:setArg(1, self.fluxCellBuf)
	self.calcCellFluxKernelObj.obj:setArg(2, self.UBuf)
	self.calcCellFluxKernelObj()
--]]
	
	self.calcFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcFluxKernelObj.obj:setArg(1, self.cellBuf)
self.calcFluxKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcFluxKernelObj.obj:setArg(3, self.UBuf)
--	self.calcFluxKernelObj.obj:setArg(4, self.fluxCellBuf)
	self.calcFluxKernelObj()

-- [=[ this is from the 2017 Zingale book
	if self.useCTU then
		self:boundary()
		-- if we're using CTU then ...
		-- 1) calc fluxes based on a slope-limiter method (PLM, etc)
		-- 2) at each interface, integrate each dimension's LR states by all other dimensions' fluxes with a timestep of -dt/2
		--	( don't use the deriv buf because it already has the sum of all dimensions' flux differences)
		self.updateCTUKernelObj(self.solverBuf, self.cellBuf, self.UBuf, self.fluxBuf, dtArg)

--		self.calcCellFluxKernelObj()
	
		-- now we need to calcBounds on the ULR
		-- TODO this will break for mirror conditions
		-- because I haven't got the boundary code flexible enough to operate on specific fields within the L & R fields of the ULRBuf
		self:boundaryLR()

		-- 3) use the final LR states to calculate the flux ...

		-- the rest of this matches above
		-- maybe use 'repeat'?
		
		self.calcFluxKernelObj()
	end
--]=]
	
	self:boundary()
	self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
self.calcDerivFromFluxKernelObj.obj:setArg(0, self.solverBuf)
self.calcDerivFromFluxKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFromFluxKernelObj()
end

return WENO
