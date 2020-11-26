--[[
This is a Solver that uses a fluxBuffer
and that uses calcDerivFromFlux to update it
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local math = require 'ext.math'
local file = require 'ext.file'
local template = require 'template'
local real = require 'hydro.real'
local GridSolver = require 'hydro.solver.gridsolver'

local common = require 'hydro.common'
local xNames = common.xNames


local FiniteVolumeSolver = class(GridSolver)

FiniteVolumeSolver.name = 'fvsolver'

function FiniteVolumeSolver:initObjs(args)
	FiniteVolumeSolver.super.initObjs(self, args)	
	self:createFlux(args.flux, args.fluxArgs)
end

function FiniteVolumeSolver:initCodeModules()
	FiniteVolumeSolver.super.initCodeModules(self)

	self.flux:initCodeModules()

	-- the calcFlux kernel is in fvsolver.cl for gridsolvers and in meshsolver.lua for meshsolvers
	-- both are dependent on calcFluxForInterface
	-- which is set up in hydro/flux/flux.lua
	self.solverModulesEnabled['calcFlux'] = true

	self.modules:addFromMarkup(
		template(file['hydro/solver/fvsolver.cl'], {
			solver = self,
			eqn = self.eqn,
			flux = self.flux,
		})
	)
	
	self.solverModulesEnabled['calcDerivFromFlux'] = true

	self:initCodeModule_calcFlux()
end
	
function FiniteVolumeSolver:initCodeModule_calcFlux()
	self.modules:addFromMarkup{
		code = template([[
//// MODULE_NAME: calcFlux
//// MODULE_DEPENDS: calcFluxForInterface cons_parallelPropagate normal_t
// used by all gridsolvers.  the meshsolver alternative is in solver/meshsolver.lua

<? 
local useFlux = solver.fluxLimiter > 1 
	and flux.usesFluxLimiter -- just flux/roe.lua right now
?>

kernel void calcFlux(
	constant solver_t const * const solver,
	global cons_t * const fluxBuf,
	const global <?=solver.getULRArg?>,
	realparam const dt,	//not used by HLL, just making this match Roe / other FV solvers
	global cell_t const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	int const indexR = index;
	real3 const xR = cellBuf[index].pos;
	
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;	

		real const dx = solver->grid_dx.s<?=side?>;

		int const indexL = index - solver->stepsize.s<?=side?>;
		real3 xL = xR;
		xL.s<?=side?> -= dx;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * dx;

		int const indexInt = side + dim * index;
		global cons_t * const flux = fluxBuf + indexInt;


<? if solver.coord.vectorComponent == 'cartesian' 
	or solver.coord.vectorComponent == 'anholonomic'
then ?>
//// MODULE_DEPENDS: cell_area#
		real area = cell_area<?=side?>(xInt);
<? else ?>
		real area = 1.<?
	for i=0,solver.dim-1 do
		if i ~= side then
			?> * solver->grid_dx.s<?=i?><?
		end
	end
?>;
<? end ?>
		if (area <= 1e-7) {
			for (int j = 0; j < numStates; ++j) {
				flux->ptr[j] = 0;
			}
		} else {
		
			<?=solver:getULRCode():gsub('\n', '\n\t\t\t')?>

			//the single act of removing the copy of the U's from global to local memory
			// increases the framerate from 78 to 127
			cons_parallelPropagate<?=side?>(ppUL, UL, xL, .5 * dx);
			cons_parallelPropagate<?=side?>(ppUR, UR, xR, -.5 * dx);

			normal_t const n = normal_forSide<?=side?>(xInt);

<? 
if useFlux then 
?>			//this is used for the flux limiter
			//should it be using the coordinate dx or the grid dx?
			//real dt_dx = dt / cell_dx<?=side?>(xInt);
<? 
	if solver.coord.vectorComponent == 'cartesian' 
	and not require 'hydro.coord.cartesian'.is(solver.coord)
	then 
?>			real const dt_dx = dt / cell_dx<?=side?>(xInt);
<? 	else 
?>			real const dt_dx = dt / dx;
<? 	end
?>
			real3 xIntL = xInt;
			xIntL.s<?=side?> -= dx;
			
			real3 xIntR = xInt;
			xIntR.s<?=side?> += dx;
			
			int const indexR2 = indexR + solver->stepsize.s<?=side?>;
			int const indexL2 = indexL - solver->stepsize.s<?=side?>;
			<?=solver:getULRCode{indexL = 'indexL2', indexR = 'indexL', suffix='_L'}:gsub('\n', '\n\t\t\t')?>
			<?=solver:getULRCode{indexL = 'indexR', indexR = 'indexR2', suffix='_R'}:gsub('\n', '\n\t\t\t')?>

			cons_parallelPropagate<?=side?>(ppUL_L, UL_L, xIntL, 1.5 * dx);		//xIntL2?
			cons_parallelPropagate<?=side?>(ppUL_R, UL_R, xIntL, .5 * dx);
			cons_parallelPropagate<?=side?>(ppUR_L, UR_L, xIntR, -.5 * dx);
			cons_parallelPropagate<?=side?>(ppUR_R, UR_R, xIntR, -1.5 * dx);	// xIntR2?

<?
end
?>
			calcFluxForInterface(
				flux,
				solver,
				UL,//ppULptr,
				UR,//ppURptr,
				xInt,
				n<? if useFlux then ?>,
				dt_dx,
				UL_L,//ppUL_Lptr,
				UL_R,//ppUL_Rptr,
				UR_L,//ppUR_Lptr,
				UR_R,//ppUR_Rptr,
				xIntL,
				xIntR<? end ?>
			);
		}
	}<? end ?>
}

]], 	{
			solver = self,
			eqn = self.eqn,
			flux = self.flux,
		}),
	}
end

function FiniteVolumeSolver:createFlux(fluxName, fluxArgs)
	assert(fluxName, "expected flux")
	local fluxClass = require('hydro.flux.'..fluxName)
	fluxArgs = table(fluxArgs, {solver=self})
	self.flux = fluxClass(fluxArgs)
end

function FiniteVolumeSolver:createBuffers()
	FiniteVolumeSolver.super.createBuffers(self)
	self:clalloc('fluxBuf', self.eqn.cons_t, self.numCells * self.dim)
end

function FiniteVolumeSolver:refreshSolverProgram()
	FiniteVolumeSolver.super.refreshSolverProgram(self)

	self.calcFluxKernelObj = self.solverProgramObj:kernel'calcFlux'
	
	self.calcDerivFromFluxKernelObj = self.solverProgramObj:kernel{name='calcDerivFromFlux', domain=self.domainWithoutBorder}
	self.calcDerivFromFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(3, self.cellBuf)
end

-- NOTICE this adds the contents of derivBufObj and does not clear it
function FiniteVolumeSolver:calcDeriv(derivBufObj, dt)
if self.checkNaNs then assert(math.isfinite(dt)) end
	local dtArg = real(dt)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
	
	if self.usePLM then
		self.calcLRKernelObj(self.solverBuf, self.cellBuf, self:getULRBuf(), self.UBuf, dtArg)
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

	self.calcFluxKernelObj(self.solverBuf, self.fluxBuf, self:getULRBuf(), dtArg, self.cellBuf)

if self.checkNaNs then assert(self:checkFinite(self.fluxBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

-- [=[ this is from the 2017 Zingale book
	if self.useCTU then
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		-- if we're using CTU then ...
		-- 1) calc fluxes based on a slope-limiter method (PLM, etc)
		-- 2) at each interface, integrate each dimension's LR states by all other dimensions' fluxes with a timestep of -dt/2
		--	( don't use the deriv buf because it already has the sum of all dimensions' flux differences)
		self.updateCTUKernelObj(self.solverBuf, self.cellBuf, self:getULRBuf(), self.fluxBuf, dtArg)
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		
		-- now we need to calcBounds on the ULR
		-- TODO this will break for mirror conditions
		-- because I haven't got the boundary code flexible enough to operate on specific fields within the L & R fields of the ULRBuf
		self:boundaryLR()
		
		-- 3) use the final LR states to calculate the flux ...
		
		-- the rest of this matches above
		-- maybe use 'repeat'?
		
--if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
--		self.calcEigenBasisKernelObj()
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		self.calcFluxKernelObj()
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.fluxBufObj)) end
	end
--]=]

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.fluxBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
	
	self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFromFluxKernelObj()

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

end


function FiniteVolumeSolver:getModuleDepends_displayCode()
	return table(FiniteVolumeSolver.super.getModuleDepends_displayCode(self)):append{
		-- wave #
		'normal_t',
		'eigen_t',
		'eqn.waveCode',
		'eigen_forInterface',
	}
end

function FiniteVolumeSolver:addDisplayVars()
	FiniteVolumeSolver.super.addDisplayVars(self)

	for side=0,self.dim-1 do
		local xj = xNames[side+1]
		self:addDisplayVarGroup{
			name = 'flux '..xj, 
			bufferField = 'fluxBuf',
			bufferType = self.eqn.cons_t,
			codePrefix = template([[
	int indexInt = <?=side?> + dim * index;
	const global cons_t* flux = buf + indexInt;
]],			{
				eqn = self.eqn,
				side = side,
			}),
			vars = range(0,self.eqn.numIntStates-1):map(function(i)
				return {name=tostring(i), code='value.vreal = flux->ptr['..i..'];'}
			end),
		}
	end

	-- code for getting the interface eigensystem variables
	local function getEigenCode(args)
		return template([[
	int indexR = index;
	int indexL = index - solver->stepsize.s<?=side?>;
	real3 xInt = x;
	xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	<?=solver:getULRCode{bufName='buf', side=side}:gsub('\n', '\n\t')?>
	normal_t n = normal_forSide<?=side?>(xInt);
	eigen_t eig;
	eigen_forInterface(&eig, solver, UL, UR, xInt, n);
]], 	{
			solver = self,
			eqn = self.eqn,
			side = args.side,
		})
	end
	
	for side=0,self.dim-1 do
		local xj = xNames[side+1]
		self:addDisplayVarGroup{
			name = 'wave '..xj,
			bufferField = self.getULRBufName,
			bufferType = self.getULRBufType,
			codePrefix = table{
				getEigenCode{side=side},
				template([[
	normal_t n<?=side?> = normal_forSide<?=side?>(xInt);
<?=eqn:eigenWaveCodePrefix('n', '&eig', 'xInt')?>
]], 			{
					eqn = self.eqn,
					side = side,
				}),
			}:concat'\n',
			vars = range(0, self.eqn.numWaves-1):map(function(i)
				return {name=tostring(i), code=template([[
	value.vreal = <?=eqn:eigenWaveCode('n'..side, '&eig', 'xInt', i)?>;
]], 			{
					eqn = self.eqn,
					side = side,
					i = i,
				})}
			end),
		}
	end

	-- TODO rename to 'getEigenDisplayVarDescs()'
	local eigenDisplayVars = self.eqn:getEigenDisplayVars()
	if eigenDisplayVars and #eigenDisplayVars > 0 then
		for side=0,self.dim-1 do
			local xj = xNames[side+1]
			self:addDisplayVarGroup{
				name = 'eigen '..xj,
				bufferField = self.getULRBufName,
				bufferType = self.getULRBufType,
				codePrefix = getEigenCode{side=side},
				vars = eigenDisplayVars,
			}
		end
	end

	if self:isModuleUsed'eigen_left/rightTransform' then
		-- ortho
		-- TODO why is the x error getting 'nans' after a few iterations?
		-- TODO WARNING - this only does ForSide, which doesn't match non-cartesian grids w/cartesian components
		for side=0,self.dim-1 do
			self:addDisplayVarGroup{
				name = 'ortho error '..xNames[side+1],
				bufferField = self.getULRBufName,
				bufferType = self.getULRBufType,
				codePrefix = '',
				useLog = true,
				vars = {
					{name='0', code=table{
						getEigenCode{side=side},
						template([[
	value.vreal = 0;
	//the flux transform is F v = R Lambda L v, I = R L
	//but if numWaves < numIntStates then certain v will map to the nullspace 
	//so to test orthogonality for only numWaves dimensions, I will verify that Qinv Q v = v 
	//I = L R
	//Also note (courtesy of Trangenstein) consider summing across outer products of basis vectors to fulfill rank
	for (int k = 0; k < numWaves; ++k) {
		cons_t basis;
		for (int j = 0; j < numStates; ++j) {
			basis.ptr[j] = k == j ? 1 : 0;
		}
		
		normal_t n = normal_forSide<?=side?>(xInt);
		<?=eqn.waves_t?> chars;
		eigen_leftTransform(&chars, solver, &eig, &basis, xInt, n);
		cons_t newbasis;
		eigen_rightTransform(&newbasis, solver, &eig, &chars, xInt, n);
	
		for (int j = 0; j < numStates; ++j) {
			value.vreal += fabs(newbasis.ptr[j] - basis.ptr[j]);
		}
	}
]], 						{
								eqn = self.eqn,
								side = side,
							}),
						}:concat'\n',
					},
				}
			}
		end
	end

	if self:isModuleUsed'eigen_fluxTransform' 
	and self:isModuleUsed'eigen_left/rightTransform'
	then
		-- flux
		-- TODO same as above, why is the x error getting 'nans' after a few iterations?
		-- TODO WARNING - this only does ForSide, which doesn't match non-cartesian grids w/cartesian components
		for side=0,self.dim-1 do
			self:addDisplayVarGroup{
				name = 'flux error '..xNames[side+1],
				bufferField = self.getULRBufName,
				bufferType = self.getULRBufType,
				codePrefix = '',
				useLog = true,
				vars = {
					{name='0', code=table{
						getEigenCode{side=side},
						template([[
	normal_t n<?=side?> = normal_forSide<?=side?>(x);
	<?=eqn:eigenWaveCodePrefix('n'..side, '&eig', 'xInt'):gsub('\n', '\n\t')?>
	
	value.vreal = 0;
	for (int k = 0; k < numIntStates; ++k) {
		//This only needs to be numIntStates in size, but just in case the left/right transforms are reaching past that memory boundary ...
		//Then again, how do I know what the non-integrated states should be?  Defaulting to zero is a bad idea.
		cons_t basis;
		for (int j = 0; j < numStates; ++j) {
			basis.ptr[j] = k == j ? 1 : 0;
		}

		normal_t n = normal_forSide<?=side?>(xInt);
		<?=eqn.waves_t?> chars;
		eigen_leftTransform(&chars, solver, &eig, &basis, xInt, n);

		<?=eqn.waves_t?> charScaled;
		<? for j=0,eqn.numWaves-1 do ?>{
			real lambda_j = <?=eqn:eigenWaveCode('n'..side, '&eig', 'xInt', j)?>;
			charScaled.ptr[<?=j?>] = chars.ptr[<?=j?>] * lambda_j;
		}<? end ?>
	
		//once again, only needs to be numIntStates
		cons_t newtransformed;
		eigen_rightTransform(&newtransformed, solver, &eig, &charScaled, xInt, n);

#if 1
		//this shouldn't need to be reset here
		// but it will if leftTransform does anything destructive
		for (int j = 0; j < numStates; ++j) {
			basis.ptr[j] = k == j ? 1 : 0;
		}
#endif

		//once again, only needs to be numIntStates
		cons_t transformed;
		eigen_fluxTransform(&transformed, solver, &eig, &basis, xInt, n);
		
		for (int j = 0; j < numIntStates; ++j) {
			value.vreal += fabs(newtransformed.ptr[j] - transformed.ptr[j]);
		}
	}
]], 						{
								solver = self,
								eqn = self.eqn,
								side = side,
							}),
						}:concat'\n',
					},
				}
			}
		end
	end
end

return FiniteVolumeSolver
