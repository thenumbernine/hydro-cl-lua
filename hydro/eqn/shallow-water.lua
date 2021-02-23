--[[
based off of:
	2011 Berger et al - The GeoClaw software for depth-averaged flows with adaptive refinement
with help from:
	2016 Titov et al - Development of MOST for Real-Time Tsunami Forecasting
	https://en.wikipedia.org/wiki/Shallow_water_equations
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local ShallowWater = class(Equation)
ShallowWater.name = 'shallow_water'

ShallowWater.roeUseFluxFromCons = true

ShallowWater.initConds = require 'hydro.init.euler':getList()

ShallowWater.numWaves = 4
ShallowWater.numIntStates = 4	-- don't integrate depth.  instead read that from an image and keep it.
-- TODO store it in cellBuf instead of UBuf, so derivs don't mess with it
-- TODO same with all other non-integratable state-vars in all other equations.

-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
ShallowWater.primVars = {
	{name='h', type='real', units='m'},						-- total height
	{name='v', type='real3', units='m/s', variance='u'},	-- contravariant
	{name='depth', type='real', units='m'},					-- bathymetry depth / resting height.  TODO move to cellBuf
}

ShallowWater.consVars = {
	{name='h', type='real', units='m'},
	{name='m', type='real3', units='m^2/s', variance='u'},	-- contravariant
	{name='depth', type='real', units='m'},					-- TODO move to cellBuf
}

function ShallowWater:resetState()
	local solver = self.solver
	
	-- if it's a grid, read from image
	-- if it's a mesh then ... ? 
	if require 'hydro.solver.meshsolver':isa(solver) 
	or solver.dim ~= 2
	then
		print("ShallowWater only does depth images for 2D -- skipping")
		return
	end
	
	local Image = require 'image'

	-- 1-channel, 16-bit signed, negative = below sea level
	local image = Image'bathymetry/world - pacific.tif'
		:resize(
			tonumber(solver.sizeWithoutBorder.x),
			tonumber(solver.sizeWithoutBorder.y))

	local ptr = solver.UBufObj:toCPU()
	for j=0,tonumber(solver.sizeWithoutBorder.y)-1 do
		for i=0,tonumber(solver.sizeWithoutBorder.x)-1 do
			local index = i + solver.numGhost + solver.gridSize.x * (j + solver.numGhost)
			local d = tonumber(image.buffer[i + image.width * (image.height - 1 - j)])
			-- in the image, bathymetry values are negative
			-- throw away altitude values (for now?)
			ptr[index].depth = d	-- bathymetry values are negative
			-- initialize to steady-state
			ptr[index].h = math.max(0, -d)		-- total height is positive
			-- cells are 'wet' where h > 0 -- and there they should be evolved
		end
	end
	solver.UBufObj:fromCPU(ptr)
end

function ShallowWater:createInitState()
	ShallowWater.super.createInitState(self)
	self:addGuiVars{	
		{name='gravity', value=9.8, units='m/s^2'},
		{name='Manning', value=0.025},	-- 2011 Berger eqn 3, Manning coefficient
	}
end

-- don't use default
function ShallowWater:initCodeModule_fluxFromCons() end
function ShallowWater:initCodeModule_consFromPrim_primFromCons() end

function ShallowWater:getModuleDepends_waveCode()
	return {
		self.symbols.eqn_common,
		self.symbols.primFromCons,
	}
end

function ShallowWater:getModuleDepends_displayCode()
	return {
		self.symbols.eqn_common,
	}
end

ShallowWater.solverCodeFile = 'hydro/eqn/shallow-water.cl'

ShallowWater.displayVarCodeUsesPrims = true

-- [=[
ShallowWater.predefinedDisplayVars = {
	'U h',
	'U v',
	'U v mag',
	'U depth',
}
--]=]

function ShallowWater:getDisplayVars()
	local vars = ShallowWater.super.getDisplayVars(self)
	
	vars:append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='wavespeed', code='value.vreal = calc_C(solver, U);', units='m/s'},
	}

	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->h'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->h'
		end,
		units = '1/s',
	} or nil)

	return vars
end

ShallowWater.eigenVars = table{
	-- Roe-averaged vars
	{name='h', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	-- derived vars
	{name='C', type='real', units='m/s'},
}

function ShallowWater:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
real C_nLen = <?=eig?>->C * normal_len(<?=n?>);
real v_n = normal_vecDotN1(<?=n?>, <?=eig?>->v);
]], {
		eig = '('..eig..')',
		x = x,
		n = n,
	})
end

function ShallowWater:consWaveCodePrefix(n, U, x, W)
	return self:template([[
real C_nLen = calc_C(solver, <?=U?>) * normal_len(<?=n?>);
<? if not W then 
	W = 'W'
?><?=prim_t?> W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=x?>);
<? end 
?>real v_n = normal_vecDotN1(n, <?=W?>.v);
]], {
		n = n,
		U = '('..U..')',
		x = x,
		W = W and '('..W..')' or nil,
	})
end

function ShallowWater:consWaveCode(n, U, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - C_nLen)'
	elseif waveIndex >= 1 and waveIndex <= 2 then
		return 'v_n'
	elseif waveIndex == 3 then
		return '(v_n + C_nLen)'
	end
	error'got a bad waveIndex'
end

ShallowWater.eigenWaveCode = ShallowWater.consWaveCode

--[[
try to do the wet/try interface stuff here
if the height of a cell is negative then set the left/right flux to zero, thus skipping integration on this cell
... but there is no UBuf argument of calcDerivFromFlux
this is another reason to put the static / non-integrated / shallow water depth values into cellBuf
all static values should go into cellBuf
--]]
function ShallowWater:postComputeFluxCode()
	return [[
	//if (U->
]]
end

return ShallowWater
