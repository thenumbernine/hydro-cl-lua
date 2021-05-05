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
//// MODULE_DEPENDS: <?=face_t?> <?=normal_t?> <?=cell_t?>
// boundary code, since meshsolver doesn't use gridsolver's boundary: 

#define reflectCons(\
	/*<?=cons_t?> * const */result,\
	/*<?=cons_t?> const * const */U,\
	/*real3 */n,\
	/*real const */restitution\
) {\
	*(result) = *(U);\
<? --\
-- matches BoundaryMirror:getCode for vectorComponent==cartesian --\
for _,var in ipairs(eqn.consStruct.vars) do --\
	if var.type == 'real'  --\
	or var.type == 'cplx' --\
	then --\
		-- do nothing --\
	elseif var.type == 'real3'  --\
	or var.type == 'cplx3' --\
	then --\
		local field = var.name --\
		local scalar = var.type == 'cplx3' and 'cplx' or 'real' --\
		local vec3 = var.type --\
?>\
	(result)-><?=field?> = <?=vec3?>_sub(\
		(result)-><?=field?>,\
		<?=vec3?>_<?=scalar?>_mul(\
			<?=vec3?>_from_real3(n),\
			<?=scalar?>_real_mul(\
				<?=vec3?>_real3_dot(\
					(result)-><?=field?>,\
					n\
				),\
				restitution + 1.\
			)\
		)\
	);\
<?--\
	else--\
		error("need to support reflectCons() for type "..var.type)--\
	end--\
end--\
?>\
}

//TODO how to make this modular?  have the config provide the code?
#define boundaryCons(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> * const */result,\
	/*<?=cons_t?> const * const */U,\
	/*global <?=face_t?> const * const */e,\
	/*real const */restitution\
) {\
<? if false then -- reflect: m dot n = 0 ?>\
	reflectCons(result, U, e->normal, 1.);\
<? end ?>\
<? if true then -- for [-1,1]^2 box with cylinder removed ?>\
	*(result) = *(U);\
	real3 const x = e->pos;\
	if (real3_lenSq(e->pos) > .7*.7) {\
		/*outside = freeflow */\
		/**(result) = *(U); */\
<? if true then ?>\
		real rho = 1.;\
		real3 v = _real3(-0.1, 0, 0);\
		real P = 1.;\
		(result)->rho = rho;\
		(result)->m = real3_real_mul(v, rho);\
		(result)->ETotal = P / (solver->heatCapacityRatio - 1.) + (.5 * coordLenSq(v, x) + (U)->ePot) * rho;\
<? end ?>\
	} else {\
		/* inside = reflect */\
		/*reflectCons(result, U, e->normal, -1);*/\
		/*reflectCons(result, U, e->normal, 0.);*/\
		reflectCons(result, U, e->normal, 1.); /* ghost U momentum is reflected from U's, s the velocity is zero (right?) */\
		/*(result)->m = real3_zero;*/\
	}\
<? end ?>\
<? if false then -- for naca 0012 airfoil ?>\
	if (real3_lenSq(e->pos) > 4.) {\
		/* outside boundary: freeflow */\
		*(result) = *(U);\
	} else {\
		/* inside boundary: v=0 */\
		*(result) = *(U);\
		(result)->m = real3_zero;\
		/* inside boundary: reflect */\
		/*reflectCons(result, U, e->normal, 1.);*/\
	}\
<? end ?>\
}

kernel void <?=calcFlux?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const fluxBuf,
	global <?=cons_t?> const * const UBuf,
	realparam const dt,
//mesh-specific parameters:	
	global <?=cell_t?> const * const cellBuf,			//[numCells]
	global <?=face_t?> const * const faceBuf,			//[numFaces]
	global int const * const cellFaceIndexes	//[numCellFaceIndexes]
) {
	int faceIndex = get_global_id(0);
	if (faceIndex >= get_global_size(0)) return;
	
	global <?=cons_t?> * const flux = fluxBuf + faceIndex;
	
	global <?=face_t?> const * const face = faceBuf + faceIndex;
	if (face->area <= 1e-7) {
		for (int j = 0; j < numStates; ++j) {
			flux->ptr[j] = 0;
		}
		return;
	}

	real3 const x = face->pos;
	<?=normal_t?> const n = normal_forFace(face);

	<?=cell_t?> cellL, cellR;
	<?=cons_t?> UL, UR;	
	//getEdgeStates(solver, &UL, &UR, face, UBuf, solver->boundaryRestitution);
	{
		int const iL = face->cells.s0;
		int const iR = face->cells.s1;
		if (iL != -1 && iR != -1) {
			UL = UBuf[iL];
			UR = UBuf[iR];
			cellL = cellBuf[iL];
			cellR = cellBuf[iR];
		} else if (iL != -1) {
			UL = UBuf[iL];
			cellL = cellR = cellBuf[iL];
			boundaryCons(solver, &UR, &UL, face, solver->boundaryRestitution);
		} else if (iR != -1) {
			UR = UBuf[iR];
			cellL = cellR = cellBuf[iR];
			boundaryCons(solver, &UL, &UR, face, solver->boundaryRestitution);
		} else {	// both iL and iR are null ...
			//error
			for (int i = 0; i < numStates; ++i) {
				UL.ptr[i] = UR.ptr[i] = 0./0.;
			}
		}
	}

	//TODO option to rotate to align fluxes?
	// then you'd have to build a new normal_t based on the aligned (x-axis) normal.

//// MODULE_DEPENDS: <?=calcFluxForInterface?>
	<?=calcFluxForInterface?>(flux, solver, &UL, &UR, &cellL, &cellR, x, n);
}

//// MODULE_NAME: <?=calcDerivFromFlux?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=solver_macros?> <?=SETBOUNDS_NOGHOST?>

kernel void <?=calcDerivFromFlux?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const fluxBuf,
//mesh-specific parameters:	
	global <?=cell_t?> const * const cellBuf,			//[numCells]
	global <?=face_t?> const * const faceBuf,			//[numFaces]
	global int const * const cellFaceIndexes	//[numCellFaceIndexes]
) {
	int const cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	global <?=cell_t?> const * const cell = cellBuf + cellIndex;
	global <?=cons_t?> * const deriv = derivBuf + cellIndex;
	
	for (int i = 0; i < cell->faceCount; ++i) {
		int const faceIndex = cellFaceIndexes[i + cell->faceOffset];
		global <?=face_t?> const * const face = faceBuf + faceIndex;
		
		global <?=cons_t?> const * const flux = fluxBuf + faceIndex;
		real const areaOverVolume = face->area / cell->volume;
		
		if (cellIndex == face->cells.s0) {
//std::cout << " ... - " << *face << std::endl;
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] -= flux->ptr[j] * areaOverVolume;
			}
		} else {
//std::cout << " ... + " << *face << std::endl;
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] += flux->ptr[j] * areaOverVolume;
			}
		}
	}
}
