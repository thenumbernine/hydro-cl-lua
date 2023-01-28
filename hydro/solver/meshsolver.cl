<?
--local faceAreaEpsilon = 1e-7
local faceAreaEpsilon = 0

?>
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

//// MODULE_NAME: <?=boundaryCons?>
//// MODULE_DEPENDS: <?=face_t?> <?=normal_t?> <?=cell_t?>
// boundary code, since meshsolver doesn't use gridsolver's boundary:

void <?=boundaryCons?>(
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> * const result,
	<?=cons_t?> const * const U,
	global <?=face_t?> const * const e
) {
<? for boundaryMethodIndex,boundaryMethod in ipairs(solver.boundaryMethods) do ?>
	<?
if boundaryMethodIndex > 1 then
	?>} else <?
end
	?>if (e->boundaryMethodIndex == <?=boundaryMethodIndex?>) {
<?
	local code, depends = boundaryMethod:getCode{
			solver = solver,
			dst = "result",
			src = "U",
			face = "e",
		}
	if depends then
?>//// MODULE_DEPENDS: <?=table.concat(depends, ' ')?>
<?	end
?>
		<?=code:gsub("\n", "\n\t\t")?>
<? end
if #solver.boundaryMethods > 0 then
?>
	} else <?
end
	?>{
		*(result) = *(U);	//default = freeflow?
	}
}

//// MODULE_NAME: <?=getEdgeStates?>
//// MODULE_DEPENDS: <?=boundaryCons?>

#define <?=getEdgeStates?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> * const */UL,\
	/*<?=cons_t?> * const */UR,\
	/*<?=cell_t?> * const & */cellL,\
	/*<?=cell_t?> * const & */cellR,\
	/*global <?=face_t?> const * const */face,\
	/*global <?=cons_t?> const * const */UBuf\
) {\
	int const iL = (face)->cells.s0;\
	int const iR = (face)->cells.s1;\
	if (iL != -1 && iR != -1) {\
		*(UL) = UBuf[iL];\
		*(UR) = UBuf[iR];\
		cellL = cellBuf[iL];\
		cellR = cellBuf[iR];\
	} else if (iL != -1) {\
		*(UL) = UBuf[iL];\
		cellL = cellR = cellBuf[iL];\
		<?=boundaryCons?>(solver, UR, UL, face);\
	} else if (iR != -1) {\
		*(UR) = UBuf[iR];\
		cellL = cellR = cellBuf[iR];\
		<?=boundaryCons?>(solver, UL, UR, face);\
	} else {	/* both iL and iR are null ... */\
		/*error*/\
		for (int i = 0; i < numStates; ++i) {\
			(UL)->ptr[i] = (UR)->ptr[i] = 0./0.;\
		}\
	}\
}

//// MODULE_NAME: <?=calcFlux?>
//// MODULE_DEPENDS: <?=getEdgeStates?> <?=face_t?> <?=normal_t?> <?=cell_t?>

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

//TODO keep this here, or move it somewhere else?
for (int j = 0; j < numStates; ++j) {
	flux->ptr[j] = 0;
}

	global <?=face_t?> const * const face = faceBuf + faceIndex;
	if (face->area <= <?=clnumber(faceAreaEpsilon)?>) {
		for (int j = 0; j < numStates; ++j) {
			flux->ptr[j] = 0;
		}
		return;
	}

	//TODO: cell_t or cell_t*?
	<?=cell_t?> cellL, cellR;
	<?=cons_t?> UL, UR;
	<?=getEdgeStates?>(solver, &UL, &UR, cellL, cellR, face, UBuf);

	//TODO option to rotate to align fluxes?
	// then you'd have to build a new normal_t based on the aligned (x-axis) normal.

//// MODULE_DEPENDS: <?=calcFluxForInterface?>
	<?=normal_t?> const n = normal_forFace(face);
	*flux = <?=calcFluxForInterface?>(solver, &UL, &UR, &cellL, &cellR, face->pos, n);
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
