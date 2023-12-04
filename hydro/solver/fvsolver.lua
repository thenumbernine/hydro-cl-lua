--[[
This is a Solver that uses a fluxBuffer
and that uses calcDerivFromFlux to update it
--]]
local table = require 'ext.table'
local range = require 'ext.range'
local math = require 'ext.math'
local path = require 'ext.path'
local real = require 'hydro.real'
local GridSolver = require 'hydro.solver.gridsolver'

local common = require 'hydro.common'
local xNames = common.xNames

local FiniteVolumeSolver = GridSolver:subclass()

FiniteVolumeSolver.name = 'fvsolver'

function FiniteVolumeSolver:getSymbolFields()
	return FiniteVolumeSolver.super.getSymbolFields(self):append{
		'FVSolver',
		'calcFlux',
		'calcDerivFromFlux',
	}
end

function FiniteVolumeSolver:initObjs(args)
	-- used to go after but I need it before so I can manipulate solver_t
	self:createFlux(args.flux, args.fluxArgs)

	FiniteVolumeSolver.super.initObjs(self, args)
end

function FiniteVolumeSolver:initCodeModules()
	FiniteVolumeSolver.super.initCodeModules(self)

	self.flux:initCodeModules()

	-- the calcFlux kernel is in fvsolver.clcpp for gridsolvers and in meshsolver.lua for meshsolvers
	-- both are dependent on calcFluxForInterface
	-- which is set up in hydro/flux/flux.lua
	self.solverModulesEnabled[self.symbols.calcFlux] = true

	self.modules:addFromMarkup(
		self.eqn:template(path'hydro/solver/fvsolver.clcpp':read(), {
			flux = self.flux,
		})
	)

	self.solverModulesEnabled[self.symbols.calcDerivFromFlux] = true

	self:initCodeModule_calcFlux()
end

-- overridden by weno subclass
function FiniteVolumeSolver:initCodeModule_calcFlux()
	self.modules:addFromMarkup{
		code = self.eqn:template([[
//// MODULE_NAME: <?=calcFlux?>
// used by all gridsolvers.  the meshsolver alternative is in solver/meshsolver.lua

// I'm putting the function in the FVSolver module.
// I want uniquely-generated names for calcFlux and calcDerivFromFlux because they are each associated with kernels.
// but I don't want to separate out the code since I'm merging everything into cpp files
// I do want this here and overrideable since WENO etc use dif args than this
//// MODULE_DEPENDS: <?=FVSolver?>

kernel void <?=calcFlux?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const fluxBuf,
	global const <?=solver.getULRArg?>,
	realparam const dt,	//not used by HLL, just making this match Roe / other FV solvers
	global <?=cell_t?> const * const cellBuf
) {
	using namespace <?=Solver?>;
	auto const & solver = *psolver;
	FVSolver::calcFlux<
		// TODO here put the lua solver.flux.cppClassName
		// but I"m basically doing the same thing with the 'using Flux =' inside each flux's C++ file
		Flux
	>(
		solver,
		fluxBuf,
		<?=solver.getULRBufName?>,
		dt,
		cellBuf
	);
}

]], 	{
			flux = self.flux,
		}),
	}
end

function FiniteVolumeSolver:createFlux(fluxName, fluxArgs)
	assert(fluxName, "expected flux")
	local fluxClass = require('hydro.flux.'..fluxName)
	fluxArgs = table(fluxArgs, {solver=self}):setmetatable(nil)
	self.flux = fluxClass(fluxArgs)
end

function FiniteVolumeSolver:createBuffers()
	FiniteVolumeSolver.super.createBuffers(self)
	self:clalloc('fluxBuf', self.eqn.symbols.cons_t, self.numCells * self.dim)
end

function FiniteVolumeSolver:refreshSolverProgram()
	FiniteVolumeSolver.super.refreshSolverProgram(self)

	self.calcFluxKernelObj = self.solverProgramObj:kernel(self.symbols.calcFlux)

	self.calcDerivFromFluxKernelObj = self.solverProgramObj:kernel{name=self.symbols.calcDerivFromFlux, domain=self.domainWithoutBorder}
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

-- [=[ this is from the 2017 Zingale "Introduction to Computational Astrophysics"
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
--[[
eqn = euler
solver = fvsolver
flux = roe
real = half,
gridSize = {512,512} or {512,511}
we get our first nan here at derivBufObj cell position {2,2}
seems intermediate now (thanks to replacing a few 0's with toreal(0)'s?)

// within fvsolver.cl calcDerivFromFlux():
grid_dx.x = 2/512 = 0.00390625
grid_dx.y = 2/512 = 0.00390625
volume = (2/512)^2 = 1.52587890625e-05
invVolume = 65536.0
areaL = grid_dx.y * invVolume = 256.0
areaR = grid_dx.x * invVolume = 256.0
volume threshold is 1e-7
--]]
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

end

function FiniteVolumeSolver:addDisplayVars()
	FiniteVolumeSolver.super.addDisplayVars(self)

	for side=0,self.dim-1 do
		local xj = xNames[side+1]
		self:addDisplayVarGroup{
			name = 'flux '..xj,
			bufferField = 'fluxBuf',
			bufferType = self.eqn.symbols.cons_t,
			codePrefix = self.eqn:template([[
int const indexInt = <?=side?> + dim * index;
global <?=cons_t?> const & flux = buf[indexInt];
]],			{
				side = side,
			}),
			vars = range(0,self.eqn.numIntStates-1):mapi(function(i)
				return {name=tostring(i), code='value.vreal = flux['..i..'];'}
			end),
		}
	end

	-- code for getting the interface eigensystem variables
	local function getEigenCode(args)
		return self.eqn:template([[
int indexR = index;
int indexL = index - solver.stepsize[<?=side?>];
global <?=cell_t?> const & cellL = cellBuf[indexL];
global <?=cell_t?> const & cellR = cellBuf[indexR];

/* TODO this isn't always used by normal_forSide or eigen_forInterface ... */
real3 xInt = x;
xInt[<?=side?>] -= .5 * solver.grid_dx[<?=side?>];

<?=solver:getULRCode{bufName="buf", side=side}:gsub("\n", "\n\t")?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=Solver?>
auto n = <?=Solver?>::Normal::forSide<<?=side?>>(xInt);
auto eig = <?=Solver?>::Eqn::eigen_forInterface(solver, UL, UR, cellL, cellR, xInt, n);
]], 	{
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
				self.eqn:template([[
//// MODULE_DEPENDS: <?=Solver?>
auto n<?=side?> = <?=Solver?>::Normal::forSide<<?=side?>>(xInt);
auto calcWaves = <?=Solver?>::Eqn::EigenWaveCode(solver, eig, n<?=side?>, xInt);
]], 			{
					side = side,
				}),
			}:concat'\n',
			vars = range(0, self.eqn.numWaves-1):mapi(function(i)
				return {name=tostring(i), code=self.eqn:template([[
value.vreal = calcWaves(solver, eig, n<?=side?>, xInt, <?=i?>);
]], 			{
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

	--if self:isModuleUsed(self.eqn.symbols.eigen_leftTransform)
	--and self:isModuleUsed(self.eqn.symbols.eigen_rightTransform)
	--then
	do
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
						self.eqn:template([[
value.vreal = 0;
//the flux transform is F v = R Lambda L v, I = R L
//but if numWaves < numIntStates then certain v will map to the nullspace
//so to test orthogonality for only numWaves dimensions, I will verify that Qinv Q v = v
//I = L R
//Also note (courtesy of Trangenstein) consider summing across outer products of basis vectors to fulfill rank
for (int k = 0; k < numWaves; ++k) {
	<?=Solver?>::Cons basis;
	for (int j = 0; j < numStates; ++j) {
		basis[j] = k == j ? 1 : 0;
	}

//// MODULE_DEPENDS: <?=Solver?>
	auto n = <?=Solver?>::Normal::forSide<<?=side?>>(xInt);
	auto chars = <?=Solver?>::Eqn::eigen_leftTransform(solver, eig, basis, xInt, n);
	auto newbasis = <?=Solver?>::Eqn::eigen_rightTransform(solver, eig, chars, xInt, n);

	for (int j = 0; j < numStates; ++j) {
		value.vreal += fabs(newbasis[j] - basis[j]);
	}
}
]], 						{
								side = side,
							}),
						}:concat'\n',
					},
				}
			}
		end
	end

	-- TODO detect if the clcpp has these ... hmm
	-- but that would mean moving 'addDisplayVarGroup' entirely to clcpp ...
	--if self:isModuleUsed(self.eqn.symbols.eigen_fluxTransform)
	--and self:isModuleUsed(self.eqn.symbols.eigen_leftTransform)
	--and self:isModuleUsed(self.eqn.symbols.eigen_rightTransform)
	--then
	do
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
						self.eqn:template([[
//// MODULE_DEPENDS: <?=cell_calcAvg_withPt?>

using namespace <?=Solver?>;

auto n<?=side?> = Normal::forSide<<?=side?>>(x);
auto calcWaves = Eqn::EigenWaveCode(solver, eig, n<?=side?>, xInt);

value.vreal = 0;
for (int k = 0; k < numIntStates; ++k) {
	//This only needs to be numIntStates in size, but just in case the left/right transforms are reaching past that memory boundary ...
	//Then again, how do I know what the non-integrated states should be?  Defaulting to zero is a bad idea.
	Cons basis;
	for (int j = 0; j < numStates; ++j) {
		basis[j] = k == j ? 1 : 0;
	}

	auto n = Normal::forSide<<?=side?>>(xInt);
	auto chars = Eqn::eigen_leftTransform(solver, eig, basis, xInt, n);

	Waves charScaled;
	for (int j = 0; j < numWaves; ++j) {
		charScaled[j] = chars[j] * calcWaves(solver, eig, n<?=side?>, xInt, j);
	}

	//once again, only needs to be numIntStates
	auto newtransformed = Eqn::eigen_rightTransform(solver, eig, charScaled, xInt, n);

#if 1
	//this shouldn't need to be reset here
	// but it will if leftTransform does anything destructive
	for (int j = 0; j < numStates; ++j) {
		basis[j] = k == j ? 1 : 0;
	}
#endif

	<?=cell_t?> cellAvg;
	cell_calcAvg_withPt(&cellAvg, cellL, cellR, xInt);

	//once again, only needs to be numIntStates
	auto transformed = Eqn::eigen_fluxTransform(solver, eig, basis, cellAvg, n);

	for (int j = 0; j < numIntStates; ++j) {
		value.vreal += fabs(newtransformed[j] - transformed[j]);
	}
}
]], 						{
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
