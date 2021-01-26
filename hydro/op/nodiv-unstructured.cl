<?
local scalar = op.scalar
local sub = scalar.."_sub"
local real_mul = scalar.."_real_mul"
?>

//// MODULE_NAME: <?=initPotential?>

kernel void <?=initPotential?>(
	constant <?=solver_t?> const * const solver,
	global <?=op:getPotBufType()?> * const UBuf,

	global <?=cell_t?> const * const cell,
	global <?=face_t?> const * const faces,		/* [numFaces] */\
	global int const * const cellFaceIndexes	/* [numCellFaceIndexes] */\

) {
	int const cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	global <?=op:getPotBufType()?> const * const U = UBuf + cellIndex;
	global <?=cell_t?> const * const cell = cells + cellIndex;
	
	<?=scalar?> source = <?=zero?>;
<?=op:getPoissonDivCode() or ""?>
	
<? if cmdline.selfGravInitPotential == "+" then
?>	U-><?=op.potentialField?> = source;
<? else
?>	U-><?=op.potentialField?> = <?=neg?>(source);
<? end
?>
}

//// MODULE_NAME: <?=solveJacobi?>
//// MODULE_DEPENDS: <?=cell_sqrt_det_g?> <?=cell_dx_i?> <?=table.concat(op.codeDepends or {}, ' ')?>

kernel void <?=solveJacobi?>(
	constant <?=solver_t?> const * const solver,
	global real * const writeBuf,
	global <?=op:getPotBufType()?> const * const UBuf,
	
	global <?=cell_t?> const * const cellBuf,
	global <?=face_t?> const * const faces,		/* [numFaces] */\
	global int const * const cellFaceIndexes	/* [numCellFaceIndexes] */\

<?
if op.stopOnEpsilon then ?>,
	global real * const reduceBuf<? 
end ?>
) {
	int const cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	real3 const x = cellBuf[cellIndex].pos;

	global <?=op:getPotBufType()?> const * const U = UBuf + cellIndex;
	global <?=cell_t?> const * const cell = cellBuf + cellIndex;

	real d2phi_dx2_diagCoeff = 0.;	//diagonal elements
	real d2phi_dx2_skewSum = 0.;	//off-diagonal elements
	real count = 0.;		//average all pairs of neighbors to cell.  TODO weight by cell size?

	for (int i = 0; i < cell->faceCount-1; ++i) {
		global <?=face_t?> const * const fi = faces + cellFaceIndexes[i + cell->faceOffset];
		int iA = fi->cells.x;
		int iB = fi->cells.y;
		if (iA != -1 && iB != -1) {
			if (iA != cellIndex) {	//assert iB is then
				iB = iA;
				iA = cellIndex;
			}

			global <?=cell_t?> const * const ciA = cells + iA;
			global <?=cell_t?> const * const ciB = cells + iB;
			global <?=cons_t?> const * const UiA = UBuf + iA;
			global <?=cons_t?> const * const UiB = UBuf + iB;

			real const phi_iA = UiA-><?=op.potentialField?>;
			real const phi_iB = UiB-><?=op.potentialField?>;
			real3 const dx_i = real3_sub(ciB->pos, ciA->pos); 
			real3 const dphi_dx_iB = _real3( phi_iB / dx_i.x,  phi_iB / dx_i.y,  phi_iB / dx_i.z);
			real3 const d_dx_iA = _real3(-1. / dx_i.x, -1. / dx_i.y, -1. / dx_i.z);
			
			for (int j = i+1; j < cell->faceCount-1; ++j) {
				global <?=face_t?> const * const fj = faces + cellFaceIndexes[j + cell->faceOffset];
				int jA = fj->cells.x;
				int jB = fj->cells.y;
				if (jA != -1 && jB != -1) {
					if (jA != cellIndex) {	//assert jB is then
						jB = jA;
						jA = cellIndex;
					}

					global <?=cell_t?> const * const cjA = cells + jA;
					global <?=cell_t?> const * const cjB = cells + jB;
					global <?=cons_t?> const * const UjA = UBuf + jA;
					global <?=cons_t?> const * const UjB = UBuf + jB;

					real const phi_jA = UjA-><?=op.potentialField?>;
					real const phi_jB = UjB-><?=op.potentialField?>;
					real3 const dx_j = real3_sub(cjB->pos, cjA->pos); 
					real3 const dphi_dx_jB = _real3( phi_jB / dx_j.x,  phi_jB / dx_j.y, phi_jB / dx_j.z);
					real3 const d_dx_jA = _real3(-1. / dx_j.x, -1. / dx_j.y, 1. / dx_j.z);

					real3 const dx2 = real3_sub(fi->pos, fj->pos);
					d2phi_dx2_skewSum += 
						  (dphi_dx_iB.x - dphi_dx_jB.x) / dx2.x
						+ (dphi_dx_iB.y - dphi_dx_jB.y) / dx2.y
						+ (dphi_dx_iB.z - dphi_dx_jB.z) / dx2.z;
					d2phi_dx2_diagCoeff += 
						  (d_dx_iA.x - d_dx_jA.x) / dx2.x
						+ (d_dx_iA.y - d_dx_jA.y) / dx2.y
						+ (d_dx_iA.z - d_dx_jA.z) / dx2.z;

					//TODO weight by cell volume?
					count++;
				}
			}
		}
	}

	real const invCount = 1. / count;
	d2phi_dx2_skewSum *= count;
	d2phi_dx2_diagCoeff *= count;


	//source is 4 pi G rho delta(x) is the laplacian of the gravitational potential field, which is integrated across discretely here
	//in units of 1/s^2
	<?=scalar?> source = <?=zero?>;
<?=op:getPoissonDivCode() or ""?>

	<?=scalar?> oldU = U-><?=op.potentialField?>;
	
	//Jacobi iteration: x_i = sum i!=j of (b_i - A_ij x_j) / A_ii
	<?=scalar?> newU = <?=real_mul?>(<?=sub?>(source, d2phi_dx2_skewSum), 1. / d2phi_dx2_diagCoeff);

	writeBuf[cellIndex] = newU;	
<? if op.stopOnEpsilon then
?>	//residual = b_i - A_ij x_j
	<?=scalar?> residual = <?=sub?>(<?=sub?>(source, d2phi_dx2_skewSum), <?=mul?>(diag, U-><?=op.potentialField?>));
	reduceBuf[cellIndex] = <?=lenSq?>(residual);
<? end
?>
}

//// MODULE_NAME: <?=noDiv?>

kernel void <?=noDiv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf
) {
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?> * const U = UBuf + index;

	real3 dv = real3_zero;
<? for j=0,solver.dim-1 do ?> 
	dv.s<?=j?> = <?=real_mul?>(
		<?=sub?>(
			U[solver->stepsize.s<?=j?>].<?=op.potentialField?>,
			U[-solver->stepsize.s<?=j?>].<?=op.potentialField?>
		), 1. / (2. * solver->grid_dx.s<?=j?>)
	);
<? end ?>
	
	<?=op:writeVectorField'dv'?>
}

//// MODULE_NAME: <?=copyPotentialToReduce?>

//TODO just use the display var kernels
kernel void <?=copyPotentialToReduce?>(
	constant <?=solver_t?> const * const solver,
	global real * const reduceBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}
