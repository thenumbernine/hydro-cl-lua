<?
--[[
I'm using from Dullemond:
dt = cfl*min(dx / (|lambdaMax| - |lambdaMin| + epsilon)
for cartesian grid: dx = grid_dx
for curvilinear grid: dx = cell_dx() (which includes metric influence, right?)
for mesh: dx = faceArea (should I use cellDist? volume? etc?)

but "I Do Like CFD" uses for mesh:
dt = cfl*min(dx/(.5*wsn))
dx = cell.volume
edu2d_module_ccfv_residual.f90: compute_residual: it mentions "wsn = the sum of (max_wave_speed)*(face length))"
	calls interface_flux:
	in edu2d_module_flux.f90: "wsn = maximum wave speed (eigenvalue)"
		interface_flux code: wsn = abs(qn) + a = |v_n| + cs = |lambdaMax|
		this is then returned as the varaible 'wave_speed'
		and then wsn(cellIndex) = sum of face[i]'s area * face[i] flux max wavespeed
--]]


-- set to false to use cell-centered U for lambda for dt
-- set to true to use face-averaged U for lambda for dt
-- TODO do this option for the gridsolver above as well
local calcDTFromFaceU = true

-- TODO put all the epsilons in one place?
--local faceAreaEpsilon = 1e-7
--local cellVolumeEpsilon = 1e-7
local faceAreaEpsilon = 0
local cellVolumeEpsilon = 0

?>
//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=normal_t?> <?=SETBOUNDS?> <?=cell_dx_i?>

<? if not require "hydro.solver.meshsolver":isa(solver) then ?>

static inline void <?=calcDTCell?>(
	real & dt,
	constant <?=solver_t?> const & solver,
	global <?=cons_t?> const & U,
	global <?=cell_t?> const & cell
) {
	<?=eqn:consWaveCodeMinMaxAllSidesPrefix{
		U = "U",
		pt = "cell.pos",
	}:gsub("\\*\n", "\\\n\t")?>
	<? for side=0,solver.dim-1 do ?>{
<?
if solver.coord.vectorComponent == "holonomic"
or require "hydro.coord.cartesian":isa(solver.coord)
then
?>		real const dx = solver.grid_dx.s<?=side?>;
<? else
?>		real const dx = cell_dx<?=side?>(cell.pos);
<? end
?>		if (dx > <?=clnumber(faceAreaEpsilon)?>) {
			<?=normal_t?> const n = normal_forSide<?=side?>(cell.pos);
			/* use cell-centered eigenvalues */
			<?=eqn:consWaveCodeMinMaxAllSides{
				n = "n",
				U = "U",
				pt = "cell.pos",
				resultMin = "lambdaMin",
				resultMax = "lambdaMax",
				declare = true,
			}:gsub("\\*\n", "\\\n\t\t\t")?>
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
			absLambdaMax = max((real)1e-9, absLambdaMax);
			dt = (real)min(dt, dx / absLambdaMax);
		}
	}<? end ?>
}

<? else -- meshsolver ?>
//// MODULE_DEPENDS: <?=face_t?>

<? if calcDTFromFaceU then ?>
//// MODULE_DEPENDS: <?=getEdgeStates?> <?=eigen_forInterface?>
<? end ?>

static inline void <?=calcDTCell?>(
	real & dt,
	constant <?=solver_t?> const & solver,
	global <?=cons_t?> const & U,
	global <?=cell_t?> const & cell,
	global <?=face_t?> const * const faces,		/* [numFaces] */
	global int const * const cellFaceIndexes	/* [numCellFaceIndexes] */
) {
	if (cell->volume > <?=clnumber(cellVolumeEpsilon)?>) {
<? if not calcDTFromFaceU then -- cell-centered lambdas ?>
		<?=eqn:consWaveCodeMinMaxAllSidesPrefix{
			U = "U",
			pt = "cell.pos",
		}:gsub("\\*\n", "\\\n\t")?>
<? else -- face-centered lambdas ?>
		/* eqn:eigenWaveCodeMinMax doesn't use a prefix (right?) */
<? end ?>
		for (int i = 0; i < cell.faceCount; ++i) {
			global <?=face_t?> const & face = faces[cellFaceIndexes[i + cell.faceOffset]];
			if (face.area > <?=clnumber(faceAreaEpsilon)?>) {
				if (face.cells.x != -1 && face.cells.y != -1) {
					<?=normal_t?> const n = normal_forFace(face);
					/* all sides? or only the most prominent side? */
					/* which should we pick eigenvalues from? */
					/* use cell-centered eigenvalues */
<? if not calcDTFromFaceU then -- cell-centered lambdas ?>
					<?=eqn:consWaveCodeMinMaxAllSides{
						n = "n",
						U = "U",
						pt = "cell.pos",
						resultMin = "lambdaMin",
						resultMax = "lambdaMax",
						declare = true,
					}:gsub("\\*\n", "\\\n\t\t\t\t\t")?>
<? else -- face-centered lambdas ?>
					/* copied from meshsolver.cl's calcFlux() ... */
					<?=cell_t?> cellL, cellR;
					<?=cons_t?> UL, UR;
					<?=getEdgeStates?>(solver, UL, UR, cellL, cellR, face, UBuf);
					<?=eigen_t?> eig = <?=eigen_forInterface?>(solver, UL, UR, cellL, cellR, face.pos, n);
					<?=eqn:eigenWaveCodeMinMax{
						eig = "eig",
						n = "n",
						pt = "face.pos",
						resultMin = "lambdaMin",
						resultMax = "lambdaMax",
						declare = true,
					}:gsub("\\*\n", "\\\n\t\t\t\t\t")?>
<? end ?>
					real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
					absLambdaMax = max((real)1e-9, absLambdaMax);
					dt = (real)min(dt, cell.volume / (face.area * absLambdaMax));
				}
			}
		}
	}
}

<? end -- meshsolver ?>
