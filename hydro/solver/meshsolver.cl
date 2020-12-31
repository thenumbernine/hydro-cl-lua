//// MODULE_NAME: INDEX
//// MODULE_HEADER:

#define INDEX(a,b,c)	((a) + solver->gridSize.x * ((b) + solver->gridSize.y * (c)))

//// MODULE_NAME: INDEXV
//// MODULE_DEPENDS: <?=solver_t?>
//// MODULE_HEADER:

#define INDEXV(i)		indexForInt4ForSize(i, solver->gridSize.x, solver->gridSize.y, solver->gridSize.z)

//// MODULE_NAME: <?=OOB?>
//// MODULE_HEADER:
// this only test for bounds of valid mesh cell

#define <?=OOB?>(lhs,rhs) (index >= get_global_size(0))

//// MODULE_NAME: <?=SETBOUNDS?>
//// MODULE_DEPENDS: <?=OOB?>
//// MODULE_HEADER:

#define <?=SETBOUNDS?>(lhs,rhs)	\
	int index = get_global_id(0); \
	int4 i = (int4)(index,0,0,0);	\
	if (<?=OOB?>(0,0)) return;
	
//// MODULE_NAME: <?=SETBOUNDS_NOGHOST?>
//// MODULE_DEPENDS: <?=OOB?>
//// MODULE_HEADER:
// this uses <?=OOB?> only to test if the cell is valid

#define <?=SETBOUNDS_NOGHOST?>()	\
	int index = get_global_id(0); \
	int4 i = (int4)(index,0,0,0); \
	if (<?=OOB?>(0,0)) return;

//// MODULE_NAME: <?=calcFlux?>
//// MODULE_DEPENDS: face_t <?=normal_t?> <?=calcFluxForInterface?> cell_t
// boundary code, since meshsolver doesn't use gridsolver's boundary: 
<?=cons_t?> reflectCons(
	<?=cons_t?> U,
	real3 n,
	float restitution
) {
<?
-- matches BoundaryMirror:getCode for vectorComponent==cartesian
for _,var in ipairs(eqn.consStruct.vars) do
	if var.type == 'real' 
	or var.type == 'cplx'
	then
		-- do nothing
	elseif var.type == 'real3' 
	or var.type == 'cplx3'
	then
		local field = var.name
		local scalar = var.type == 'cplx3' and 'cplx' or 'real'
		local vec3 = var.type
?>
	U.<?=field?> = <?=vec3?>_sub(
		U.<?=field?>,
		<?=vec3?>_<?=scalar?>_mul(
			<?=vec3?>_from_real3(n),
			<?=scalar?>_real_mul(
				<?=vec3?>_real3_dot(
					U.<?=field?>,
					n
				), 
				restitution + 1.
			)
		)
	);
<?
	else
		error("need to support reflect() for type "..var.type)
	end
end
?>	return U;
}

void getEdgeStates(
	<?=cons_t?> * const UL,
	<?=cons_t?> * const UR,
	global <?=face_t?> const * const e,
	global <?=cons_t?> const * const UBuf,		//[numCells]
	real const restitution
) {
	int const iL = e->cells.s0;
	int const iR = e->cells.s1;
	if (iL != -1 && iR != -1) {
		*UL = UBuf[iL];
		*UR = UBuf[iR];
	} else if (iL != -1) {
		*UL = UBuf[iL];
		//TODO  
		*UR = reflectCons(*UL, e->normal, restitution);
	} else if (iR != -1) {
		*UR = UBuf[iR];
		//TODO  
		*UL = reflectCons(*UR, e->normal, restitution);
	} else {	// both iL and iR are null ...
		//error
		for (int i = 0; i < numStates; ++i) {
			UL->ptr[i] = UR->ptr[i] = 0./0.;
		}
	}
}

kernel void <?=calcFlux?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const fluxBuf,
	global <?=cons_t?> const * const UBuf,
	realparam const dt,
//mesh-specific parameters:	
	global <?=cell_t?> const * const cells,			//[numCells]
	global <?=face_t?> const * const faces,			//[numFaces]
	global int const * const cellFaceIndexes	//[numCellFaceIndexes]
) {
	int faceIndex = get_global_id(0);
	if (faceIndex >= get_global_size(0)) return;
	
	global <?=cons_t?> * const flux = fluxBuf + faceIndex;
	
	global <?=face_t?> const * const face = faces + faceIndex;
	if (face->area <= 1e-7) {
		for (int j = 0; j < numStates; ++j) {
			flux->ptr[j] = 0;
		}
		return;
	}

	real3 const x = face->pos;
	<?=normal_t?> const n = normal_forFace(face);
	
	<?=cons_t?> UL, UR;	
	getEdgeStates(&UL, &UR, face, UBuf, solver->boundaryRestitution);

	//TODO option to rotate to align fluxes?
	// then you'd have to build a new normal_t based on the aligned (x-axis) normal.

	*flux = <?=calcFluxForInterface?>(solver, UL, UR, x, n);
}

//// MODULE_NAME: <?=calcDerivFromFlux?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> cell_t <?=solver_macros?> <?=SETBOUNDS_NOGHOST?>

kernel void <?=calcDerivFromFlux?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const fluxBuf,
//mesh-specific parameters:	
	global <?=cell_t?> const * const cells,			//[numCells]
	global <?=face_t?> const * const faces,			//[numFaces]
	global int const * const cellFaceIndexes	//[numCellFaceIndexes]
) {
	int const cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	global <?=cell_t?> const * const cell = cells + cellIndex;
	global <?=cons_t?> * const deriv = derivBuf + cellIndex;
	
	for (int i = 0; i < cell->faceCount; ++i) {
		int const ei = cellFaceIndexes[i + cell->faceOffset];
		global <?=face_t?> const * const e = faces + ei;
		
		global <?=cons_t?> const * const flux = fluxBuf + ei;
		real const areaOverVolume = e->area / cell->volume;
		
		if (cellIndex == e->cells.s0) {
//std::cout << " ... - " << *e << std::endl;
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] -= flux->ptr[j] * areaOverVolume;
			}
		} else {
//std::cout << " ... + " << *e << std::endl;
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] += flux->ptr[j] * areaOverVolume;
			}
		}
	}
}
