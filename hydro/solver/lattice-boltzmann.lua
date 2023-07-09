local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local vec3sz = require 'vec-ffi.vec3sz'
local GridSolver = require 'hydro.solver.gridsolver'

local LatticeBoltzmann = class(GridSolver)
LatticeBoltzmann.name = 'LatticeBoltzmann'

function LatticeBoltzmann:init(...)
	LatticeBoltzmann.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function LatticeBoltzmann:getSymbolFields()
	return LatticeBoltzmann.super.getSymbolFields(self):append{
		'macros',
		'advect',
		'calcPrims',
		'applyCollision',
	}
end

function LatticeBoltzmann:initMeshVars(args)
	LatticeBoltzmann.super.initMeshVars(self, args)
	self.solverStruct.vars:append{
		{name='ofsmax', type='int4'},
		{name='ofsvol', type='int'},
		{name='weights', type='real[27]'},
	}

	self.ofsmax = vec3sz(
		self.dim >= 1 and 3 or 1,
		self.dim >= 2 and 3 or 1,
		self.dim >= 3 and 3 or 1
	)

	self.ofsvol = self.ofsmax:volume()
end

function LatticeBoltzmann:createSolverBuf()
	LatticeBoltzmann.super.createSolverBuf(self)
	self.solverPtr.ofsmax.x = self.ofsmax.x
	self.solverPtr.ofsmax.y = self.ofsmax.y
	self.solverPtr.ofsmax.z = self.ofsmax.z
	self.solverPtr.ofsmax.w = 0
	
	-- but super also calls refreshSolverBuf...
	self.solverPtr.ofsvol = self.ofsvol
	
	if self.dim == 1 then
		self.solverPtr.weights[0] = 1/4
		self.solverPtr.weights[1] = 1/2
		self.solverPtr.weights[2] = 1/4
	elseif self.dim == 2 then
		self.solverPtr.weights[0] = 1/36
		self.solverPtr.weights[1] = 4/36
		self.solverPtr.weights[2] = 1/36
		self.solverPtr.weights[4] = 4/36
		self.solverPtr.weights[5] = 16/36
		self.solverPtr.weights[6] = 4/36
		self.solverPtr.weights[7] = 1/36
		self.solverPtr.weights[8] = 4/36
		self.solverPtr.weights[9] = 1/36
	elseif self.dim == 3 then
	else
		error("unknown dim "..tostring(self.dim))
	end

	self:refreshSolverBuf()
end

-- TODO put this in its own eqn file? eqn/lattice-boltzmann.lua?
-- or eqn/density-only?
-- or ... ??
-- or just trust/warn that only Euler equation is used
--[[
fields we need:
real rho
real3 v
bool solid
real U[3^dim] ... 'f' in LBM literature
--]]

local LatticeBoltzmannEqn = class(require 'hydro.eqn.eqn')
LatticeBoltzmannEqn.name = 'LatticeBoltzmann'
LatticeBoltzmannEqn.initConds = require 'hydro.init.lattice-boltzmann':getList()
function LatticeBoltzmannEqn:buildVars()
	self.lbmNbhdVol = 3 ^ self.solver.dim
	self.consVars = self.consVars or table()
	self.consVars:append{
		{name='rho', type='real'},
		{name='v', type='real3'},
		{name='solid', type='real'},
		{name='nbhd', type='real['..self.lbmNbhdVol..']'},
	}
end
function LatticeBoltzmannEqn:initCodeModule_calcDTCell() end
function LatticeBoltzmannEqn:initCodeModule_calcDT() end
LatticeBoltzmannEqn.predefinedDisplayVars = {
	'U rho',
	'U v x',
	'U v y',
	'U solid',
}
function LatticeBoltzmann:createEqn()
	self.eqn = LatticeBoltzmannEqn{solver=self}
end

function LatticeBoltzmann:createBuffers()
	LatticeBoltzmann.super.createBuffers(self)
	self:clalloc('UNextBuf', self.eqn.symbols.cons_t, self.numCells)
end

function LatticeBoltzmann:refreshCalcDTKernel() end

function LatticeBoltzmann:initCodeModules()
	LatticeBoltzmann.super.initCodeModules(self)
	self.modules:addFromMarkup(self.eqn:template(file'hydro/solver/lattice-boltzmann.cl':read()))
	self.solverModulesEnabled[self.symbols.advect] = true
	self.solverModulesEnabled[self.symbols.calcPrims] = true
	self.solverModulesEnabled[self.symbols.applyCollision] = true
end

function LatticeBoltzmann:refreshSolverProgram()
	LatticeBoltzmann.super.refreshSolverProgram(self)
	self.advectKernelObj = self.solverProgramObj:kernel(self.symbols.advect)
	self.calcPrimsKernelObj = self.solverProgramObj:kernel(self.symbols.calcPrims)
	self.applyCollisionKernelObj = self.solverProgramObj:kernel(self.symbols.applyCollision)
end

function LatticeBoltzmann:step(dt)
	if cmdline.printBufs then
		print()
		print('LBM step:')
		self:printBuf(self.UBuf)
	end
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	self.advectKernelObj(self.solverBuf, self.UNextBuf, self.UBuf)
	local bufferSize = self.numCells * ffi.sizeof(self.eqn.symbols.cons_t)
	self.cmds:enqueueCopyBuffer{dst=self.UBuf, src=self.UNextBuf, size=bufferSize}
	
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	self:boundary()
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	if cmdline.printBufs then
		print()
		print('LBM advect:')
		self:printBuf(self.UBuf)
	end

	self.calcPrimsKernelObj(self.solverBuf, self.UBuf)
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	self:boundary()
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	if cmdline.printBufs then
		print()
		print('LBM calcPrims:')
		self:printBuf(self.UBuf)
	end

	self.applyCollisionKernelObj(self.solverBuf, self.UBuf)
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	self:boundary()
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	if cmdline.printBufs then
		print()
		print('LBM applyCollisions:')
		self:printBuf(self.UBuf)
	end
end

return LatticeBoltzmann 
