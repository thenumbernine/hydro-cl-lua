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

--[[
This replaces the h^2 with (h^2 - 2 h H) in the flux.
In turn it also causes the wave to change from c = sqrt(g h) to c = sqrt(g (h - H)).
So setting this to 'true' causes 'depth' to be used within the wave code, which I can't calculate well,
so that means if this is true then you must set depthInCell to false - for now.
--]]
ShallowWater.extraTermInFlux = false

--[[
Whether to insert the 'depth' field into the cell_t or the cons_t structure.
Inserting it into cons_t means more wasted memory and computations.
Inserting it into cell_t means making the code a bit more flexible.

Ok so I changed the flux to go from h^2 to (h^2 - 2 h H).
When you do that, the waves go from c = sqrt(g h) to c = sqrt(g (h - H))
And when you have 'depth' in the wavefunction, you need to add the cell avg to 'calcFluxForInterface' and 'eigen_forInterface'.
So until you do that, or remove depth from the waves, this must be false.
--]]
ShallowWater.depthInCell = true

if ShallowWater.extraTermInFlux and ShallowWater.depthInCell then
	error("some of the eigen calcs don't have access to cell, so this won't work")
end

--[[
TODO
when Roe uses fluxFromCons then *only* the mv flux is 1/2 what it should be
also HLL (which uses fluxFromCons) *only* mv flux is 1/2 what it should be
So what's wrong with fluxFromCons?
Sure enough, if I remove it and replace it with the R*Lambda*L transform, things work fine.
And it looks like it works fine when I remove the 1/2 from the F->m term (BUT IT SHOULD BE THERE).
So what's wrong?

So keep this at 'true' for correct flux values (right?)
But doesn't ==true's math depend on the incorrect assumption that dF/dU * U = F?
Why does ==false produce bad values here?
--]]
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
}

ShallowWater.consVars = {
	{name='h', type='real', units='m'},
	{name='m', type='real3', units='m^2/s', variance='u'},	-- contravariant
}

if not ShallowWater.depthInCell then
	-- 'depth' in UBuf:
	table.insert(ShallowWater.primVars, {name='depth', type='real', units='m'})
	table.insert(ShallowWater.consVars, {name='depth', type='real', units='m'})
end
if ShallowWater.depthInCell then
	-- 'depth' in cellBuf:
	-- add our cell depth value here.  no more storing it in non-integrated UBuf (wasting memory and slowing performance).
	-- 'depth' = height-above-seafloor, in m, so sealevel values are negative
	-- but in order to complete this, I need to pass cellL and cellR into the eigen_forInterface function ...
	function ShallowWater:createCellStruct()
		self.solver.coord.cellStruct.vars:insert{name='depth', type='real'}
	end
end

function ShallowWater:getEnv()
	local env = table(ShallowWater.super.getEnv(self))
	env.getDepthSource = function(U, cell)
		U = U or 'U'
		cell = cell or 'cell'
		return self.depthInCell and cell or U
	end
	env.getDepth = function(U, cell)
		return '('..env.getDepthSource(U, cell)..')->depth'
	end
	return env
end

function ShallowWater:resetState()
	local solver = self.solver

	-- TODO in absense of 'readFromImage', how about providing this info in the init/euler? or init/shallow-water?
	-- TODO and for that, if eqn can (now) add custom fields to cell_t, shouldn't applyInitCond also be allowed to write to cellBuf?

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
	-- 4320x2160
	local filename = 'bathymetry/world - pacific.tif'
	local image = Image(filename)
	image = image:resize(
		tonumber(solver.sizeWithoutBorder.x),
		tonumber(solver.sizeWithoutBorder.y))

	local cpuU = solver.UBufObj:toCPU()
	-- cellCpuBuf should already be allocated in applyInitCond
	assert(self.cellCpuBuf)
	if ShallowWater.depthInCell then
		solver.cellBufObj:toCPU(self.cellCpuBuf)
	end
	for j=0,tonumber(solver.sizeWithoutBorder.y)-1 do
		for i=0,tonumber(solver.sizeWithoutBorder.x)-1 do
			local index = i + solver.numGhost + solver.gridSize.x * (j + solver.numGhost)
			
			local d = -tonumber(image.buffer[i + image.width * (image.height - 1 - j)])
			
			-- in the image, bathymetry values are negative
			-- throw away altitude values (for now?)
			if ShallowWater.depthInCell then
				self.cellCpuBuf[index].depth = d	-- bathymetry values are negative
			else
				cpuU[index].depth = d	-- bathymetry values are negative
			end
			-- initialize to steady-state
			-- this happens *after* applyInitCond, so we have to modify the UBuf 'h' values here
			cpuU[index].h = math.max(0, d)		-- total height is positive
			-- cells are 'wet' where h > depth -- and there they should be evolved
		end
	end
	solver.UBufObj:fromCPU(cpuU)
	if ShallowWater.depthInCell then
		solver.cellBufObj:fromCPU(self.cellCpuBuf)
	end
end

function ShallowWater:createInitState()
	ShallowWater.super.createInitState(self)
	self:addGuiVars{	
		{name='water_D', value=0, units='m'},	-- maximum depth.  completely trivial, only influence height values.
		
		--{name='gravity', value=9.8, units='m/s^2'},
		{name='gravity', value=1, units='m/s^2'},
		
		{name='Manning', value=0.025},	-- 2011 Berger eqn 3, Manning coefficient, used for velocity drag in source term
	}
end

-- don't use default
function ShallowWater:initCodeModule_consFromPrim_primFromCons() end

function ShallowWater:getModuleDepends_waveCode()
	return {
		self.symbols.eqn_common,
		self.symbols.primFromCons,
	}
end

function ShallowWater:getModuleDepends_displayCode()
	return table(ShallowWater.super.getModuleDepends_displayCode(self)):append{
		self.symbols.eqn_common,
	}
end

ShallowWater.solverCodeFile = 'hydro/eqn/shallow-water.cl'

ShallowWater.displayVarCodeUsesPrims = true

--[[
ShallowWater.predefinedDisplayVars = {
	'U h',
	--'U h+B',
	--'U v',
	--'U v mag',
	--ShallowWater.depthInCell and 'cell depth' or 'U depth',
	--'U B',	-- this will show the sea floor height above the maximum depth
	'U m x',

	'deriv h',
	'deriv m x',
	'flux x 0',
	'flux x 1',
}
--]]
-- [[
ShallowWater.predefinedDisplayVars = {
	'U h',
	'U h+B',
	'U B',
	'U v x',
}
--]]

function ShallowWater:getDisplayVars()
	local vars = ShallowWater.super.getDisplayVars(self)

	vars:append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='wavespeed', code='value.vreal = calc_C(solver, U);', units='m/s'},
		-- D(x) = maximum depth = constant
		-- H(x) = cell->depth = displacement of seafloor below resting depth.
		-- D(x) = H(x) + B(x)
		-- B(x) + h(x) = wave height above maximum depth
		-- h(x) - H(x) = wave height above resting depth
		-- h(x) - H(x) + D(x) = wave height above maximum depth
		{name='h+B', code = self:template'value.vreal = solver->water_D + U->h - <?=getDepth()?>;', units='m'},
		-- so everything is at rest when h(x) == H(x)
		
		-- B(x) = D(x) - H(x) = height of seafloor above maximum depth
		{name='B', code = self:template'value.vreal = solver->water_D - <?=getDepth()?>;', units='m'},
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
if not ShallowWater.depthInCell then
	-- we don't need this here if it is in cellBuf already
	-- but we do if it's not in cellBuf
	ShallowWater.eigenVars:insert{name='depth', type='real', units='m'}
end

ShallowWater.heightEpsilon = 0

function ShallowWater:eigenWaveCodePrefix(args)
	return self:template([[
real const CSq = solver->gravity * ((<?=eig?>)->h 
<? if eqn.extraTermInFlux then ?>
	- <?=getDepth(eig, 'cell')?>
<? end ?>
);
real const C = sqrt(CSq);
real C_nLen = C * normal_len(<?=n?>);
real v_n = normal_vecDotN1(<?=n?>, (<?=eig?>)->v);

if (C <= <?=clnumber(eqn.heightEpsilon)?>) {
	C_nLen = 0.;
	v_n = 0.;
}
]], args)
end

function ShallowWater:eigenWaveCode(args)
	if args.waveIndex == 0 then
		return 'v_n - C_nLen'
	elseif args.waveIndex >= 1 and args.waveIndex <= 2 then
		return 'v_n'
	elseif args.waveIndex == 3 then
		return 'v_n + C_nLen'
	end
	error'got a bad waveIndex'
end

function ShallowWater:consWaveCodePrefix(args)
	return self:template([[
real const CSq = solver->gravity * ((<?=U?>)->h 
<? if eqn.extraTermInFlux then ?>
	- <?=getDepth(U, 'cell')?>
<? end ?>
);
real const C = sqrt(CSq);
real C_nLen = C * normal_len(<?=n?>);
<?=prim_t?> W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=pt?>);
real v_n = normal_vecDotN1(<?=n?>, W.v);

if (C <= <?=clnumber(eqn.heightEpsilon)?>) {
	C_nLen = 0.;
	v_n = 0.;
}
]], args)
end

ShallowWater.consWaveCode = ShallowWater.eigenWaveCode

--ShallowWater.eigenWaveCodeMinMax uses default
--ShallowWater.consWaveCodeMinMax uses default

function ShallowWater:consWaveCodeMinMaxAllSidesPrefix(args)
	return self:template([[
real const CSq = solver->gravity * ((<?=U?>)->h 
<? if eqn.extraTermInFlux then ?>
	- <?=getDepth(U, 'cell')?>
<? end ?>
);
real const C = sqrt(CSq);
<?=prim_t?> W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=pt?>);
]], args)
end

function ShallowWater:consWaveCodeMinMaxAllSides(args)
	return self:template([[
real C_nLen = C * normal_len(<?=n?>);
real v_n = normal_vecDotN1(<?=n?>, W.v);

if (C <= <?=clnumber(eqn.heightEpsilon)?>) {
	C_nLen = 0.;
	v_n = 0.;
}

<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'v_n - C_nLen',
	'v_n + C_nLen'
)?>
]], args)
end

return ShallowWater
