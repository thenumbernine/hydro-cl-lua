local ffi = require 'ffi'
local table = require 'ext.table'
local path = require 'ext.path'
local vec3i = require 'vec-ffi.vec3i'
local GridSolver = require 'hydro.solver.gridsolver'
local real = require 'hydro.real'

local LatticeBoltzmann = GridSolver:subclass()

LatticeBoltzmann.name = 'LatticeBoltzmann'

--LatticeBoltzmann.numGhost = 1 -- I think is enough for a Lattice Boltzmann offset range of +-1

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

function LatticeBoltzmann:initLBOffsets()

	-- when I first wrote this, I was trying to minimize the amount of Lua code, for quick port to clcpp/spirv toolchain separate of Lua.
	-- but meh, too much weird offset math ... so when I do port this to clcpp then I'll move this to template code.

	-- i is 1-based, from [1,dim]
	-- [-1,1] for valid dimensions, 0 otherwise

	self.ofssize = vec3i(
		self.dim >= 1 and 3 or 1,
		self.dim >= 2 and 3 or 1,
		self.dim >= 3 and 3 or 1
	)

	self.ofsmin = vec3i(
		self.dim >= 1 and -1 or 0,
		self.dim >= 2 and -1 or 0,
		self.dim >= 3 and -1 or 0
	)

	self.offsets = table()
	for oz=0,tonumber(self.ofssize.z)-1 do
		local noz = self.ofssize.z - 1 - oz
		for oy=0,tonumber(self.ofssize.y)-1 do
			local noy = self.ofssize.y - 1 - oy
			for ox=0,tonumber(self.ofssize.x)-1 do
				local nox = self.ofssize.x - 1 - ox
				self.offsets:insert{
					c = vec3i(ox,oy,oz) + self.ofsmin,
					oppositeOffsetIndex = nox + self.ofssize.x * (noy + self.ofssize.y * noz),
				}
			end
		end
	end

	-- TODO merge with self.offsets
	local lbWeights
	if self.dim == 1 then
		lbWeights = {1/4, 1/2, 1/4}
	elseif self.dim == 2 then
		lbWeights = {
			1/36, 4/36, 1/36,
			4/36, 16/36, 4/36,
			1/36, 4/36, 1/36,
		}
	elseif self.dim == 3 then
		error'TODO'
	else
		error("unknown dim "..tostring(self.dim))
	end

	assert(#self.offsets == #lbWeights)
	for i,ofs in ipairs(self.offsets) do
		ofs.weight = lbWeights[i]
	end
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

local LatticeBoltzmannEqn = require 'hydro.eqn.eqn':subclass()
LatticeBoltzmannEqn.name = 'LatticeBoltzmann'
LatticeBoltzmannEqn.initConds = require 'hydro.init.lattice-boltzmann':getList()
function LatticeBoltzmannEqn:buildVars()
	self.consVars = self.consVars or table()

	self.solver:initLBOffsets()
	
	self.consVars:append{
		{name='rho', type='real'},
		{name='v', type='real3'},
		{name='solid', type='real'},
	}:append(self.solver.offsets:mapi(function(ofs, i)
		return {name='F'..(i-1), type='real'}
	end))
end
function LatticeBoltzmannEqn:initCodeModule_calcDTCell() end
function LatticeBoltzmannEqn:initCodeModule_calcDT() end
LatticeBoltzmannEqn.predefinedDisplayVars = {
	'U rho',
	'U v x',
	'U v y',
	'U solid',
	'U curl v',
	'U div v',
}
function LatticeBoltzmannEqn:getDisplayVars()
	local vars = LatticeBoltzmannEqn.super.getDisplayVars(self)
	vars:append{
		{name='rhoSum', code=self:template[[
value.vreal = 0;
<? for i=0,#solver.offsets-1 do
?>value.vreal += U->F<?=i?>;
<? end
?>
]], type='real'},
		{name='rho-rhoSum', code=self:template[[
value.vreal = -U->rho;
<? for i=0,#solver.offsets-1 do
?>value.vreal += U->F<?=i?>;
<? end
?>
]], type='real'},
	}

	-- TODO this for all vectors?
	vars:insert(self:createDivDisplayVar{
		field = 'v',
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		units = '1/s',
	} or nil)

	return vars
end
function LatticeBoltzmann:createEqn()
	self.eqn = LatticeBoltzmannEqn{solver=self}
end

function LatticeBoltzmann:refreshIntegrator()
	self.integrator = {}
end

function LatticeBoltzmann:createBuffers()
	LatticeBoltzmann.super.createBuffers(self)
	self:clalloc('UNextBuf', self.eqn.symbols.cons_t, self.numCells)
end

function LatticeBoltzmann:refreshCalcDTKernel() end

function LatticeBoltzmann:initCodeModules()
	LatticeBoltzmann.super.initCodeModules(self)
	self.modules:addFromMarkup(self.eqn:template(path'hydro/solver/lattice-boltzmann.cl':read()))
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
	self:boundary()
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

	self.applyCollisionKernelObj(self.solverBuf, self.UBuf, real(dt))
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
