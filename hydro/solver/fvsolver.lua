--[[
This is a Solver that uses a fluxBuffer
and that uses calcDerivFromFlux to update it
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local real = require 'hydro.real'
local GridSolver = require 'hydro.solver.gridsolver'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local FiniteVolumeSolver = class(GridSolver)

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
end

function FiniteVolumeSolver:getSolverCode()
	return table{
		FiniteVolumeSolver.super.getSolverCode(self),
		template(file['hydro/solver/calcDerivFV.cl'], {solver=self}),
	}:concat'\n'
end



-- NOTICE this adds the contents of derivBufObj and does not clear it
function FiniteVolumeSolver:calcDeriv(derivBufObj, dt)
if self.checkNaNs then assert(math.isfinite(dt)) end
	local dtArg = real(dt)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
	
	if self.usePLM then
		self.calcLRKernelObj(self.solverBuf, self:getULRBuf(), self.UBuf, dtArg)
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

	self.calcFluxKernelObj(self.solverBuf, self.fluxBuf, self:getULRBuf(), dtArg)

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
		self.updateCTUKernelObj(self.solverBuf, self:getULRBuf(), self.fluxBuf, dtArg)
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
	const global <?=eqn.cons_t?>* flux = buf + indexInt;
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
	<?=solver:getULRCode{bufName='buf'}:gsub('\n', '\n\t')?>
	normalInfo_t n = normalInfo_forSide<?=side?>(xInt);
	<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);
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
	normalInfo_t n<?=side?> = normalInfo_forSide<?=side?>(xInt);
	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'xInt'):gsub('\n', '\n\t')?>
]], 			{
					eqn = self.eqn,
					side = side,
				}),
			}:concat'\n',
			vars = range(0, self.eqn.numWaves-1):map(function(i)
				return {name=tostring(i), code=template([[
	value.vreal = <?=eqn:eigenWaveCode('n'..side, 'eig', 'xInt', i)?>;
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
		<?=eqn.cons_t?> basis;
		for (int j = 0; j < numStates; ++j) {
			basis.ptr[j] = k == j ? 1 : 0;
		}
		
		normalInfo_t n = normalInfo_forSide<?=side?>(xInt);
		<?=eqn.waves_t?> chars = eigen_leftTransform(solver, eig, basis, xInt, n);
		<?=eqn.cons_t?> newbasis = eigen_rightTransform(solver, eig, chars, xInt, n);
	
		for (int j = 0; j < numStates; ++j) {
			value.vreal += fabs(newbasis.ptr[j] - basis.ptr[j]);
		}
	}
]], 					{
							eqn = self.eqn,
							side = side,
						}),
					}:concat'\n',
				},
			}
		}
	end

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
	normalInfo_t n<?=side?> = normalInfo_forSide<?=side?>(x);
	<?=eqn:eigenWaveCodePrefix('n'..side, 'eig', 'xInt'):gsub('\n', '\n\t')?>
	
	value.vreal = 0;
	for (int k = 0; k < numIntStates; ++k) {
		//This only needs to be numIntStates in size, but just in case the left/right transforms are reaching past that memory boundary ...
		//Then again, how do I know what the non-integrated states should be?  Defaulting to zero is a bad idea.
		<?=eqn.cons_t?> basis;
		for (int j = 0; j < numStates; ++j) {
			basis.ptr[j] = k == j ? 1 : 0;
		}

		normalInfo_t n = normalInfo_forSide<?=side?>(xInt);
		<?=eqn.waves_t?> chars = eigen_leftTransform(solver, eig, basis, xInt, n);

		<?=eqn.waves_t?> charScaled;
		<? for j=0,eqn.numWaves-1 do ?>{
			real lambda_j = <?=eqn:eigenWaveCode('n'..side, 'eig', 'xInt', j)?>;
			charScaled.ptr[<?=j?>] = chars.ptr[<?=j?>] * lambda_j;
		}<? end ?>
	
		//once again, only needs to be numIntStates
		<?=eqn.cons_t?> newtransformed = eigen_rightTransform(solver, eig, charScaled, xInt, n);

#if 1
		//this shouldn't need to be reset here
		// but it will if leftTransform does anything destructive
		for (int j = 0; j < numStates; ++j) {
			basis.ptr[j] = k == j ? 1 : 0;
		}
#endif

		//once again, only needs to be numIntStates
		<?=eqn.cons_t?> transformed = eigen_fluxTransform(solver, eig, basis, xInt, n);
		
		for (int j = 0; j < numIntStates; ++j) {
			value.vreal += fabs(newtransformed.ptr[j] - transformed.ptr[j]);
		}
	}
]], 					{
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

return FiniteVolumeSolver
