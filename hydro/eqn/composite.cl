<?
local function subeqnDepends(symbol)
	return eqn.eqns:mapi(function(subeqn)
		return subeqn.symbols[symbol]
	end):concat' '
end
?>

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=subeqnDepends'eigen_forInterface'?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_forInterface?>(\
		&(resultEig)-><?=subeqn.field?>,\
		solver,\
		&(UL)-><?=subeqn.field?>,\
		&(UR)-><?=subeqn.field?>,\
		cellL,\
		cellR,\
		pt,\
		n)\
<? end --\
?>\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=subeqnDepends'eigen_forCell'?>

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_forCell?>(\
		&(resultEig)-><?=subeqn.field?>,\
		solver,\
		&(U)-><?=subeqn.field?>,\
		cell,\
		n)\
<? end --\
?>\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=subeqnDepends'eigen_leftTransform'?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> */n\
) {\
<? --\
for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_leftTransform?>(\
		&(result)-><?=subeqn.field?>,\
		solver,\
		&(eig)-><?=subeqn.field?>,\
		&(X)-><?=subeqn.field?>,\
		pt,\
		n)\
<? --\
end --\
?>\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=subeqnDepends'eigen_rightTransform'?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
<? --\
for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_rightTransform?>(\
		&(result)-><?=subeqn.field?>,\
		solver,\
		&(eig)-><?=subeqn.field?>,\
		&(X)-><?=subeqn.field?>,\
		pt,\
		n)\
<? --\
end --\
?>\
}

//// MODULE_NAME: <?=cons_parallelPropagate?>
//// MODULE_DEPENDS: <?=subeqnDepends'cons_parallelPropagate'?>

<? for side=0,solver.dim-1 do ?>
#define <?=cons_parallelPropagate?><?=side?>(\
	resultName,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*real const */dx\
)\
<?=cons_t?> resultName##base;\
{\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.cons_parallelPropagate?><?=side?>(\
		<?=subeqn.symbolPrefix?>resultName,\
		&(U)-><?=subeqn.field?>,\
		pt,\
		dx);\
	resultName##base.<?=subeqn.field?> = *<?=subeqn.symbolPrefix?>resultName;\
<? end --\
?>\
}\
<?=cons_t?> * const resultName = &resultName##base;\
<? end ?>


//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=subeqnDepends'fluxFromCons'?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultF,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const*/U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.fluxFromCons?>(\
		&(resultF)-><?=subeqn.field?>,\
		solver,\
		&(U)-><?=subeqn.field?>,\
		cell,\
		n);\
<? end --\
?>\
}


//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: <?=subeqnDepends'primFromCons'?>

#define <?=primFromCons?>(\
	/*<?=cons_t?> * const */W,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const*/U,\
	/*real3 const */pt\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.primFromCons?>(\
		&(W)-><?=subeqn.field?>,\
		solver,\
		&(U)-><?=subeqn.field?>,\
		pt);\
<? end --\
?>\
}


//// MODULE_NAME: <?=consFromPrim?>
//// MODULE_DEPENDS: <?=subeqnDepends'consFromPrim'?>

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */U,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const*/W,\
	/*real3 const */pt\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.consFromPrim?>(\
		&(U)-><?=subeqn.field?>,\
		solver,\
		&(W)-><?=subeqn.field?>,\
		pt);\
<? end --\
?>\
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?> <?=subeqnDepends'applyInitCondCell'?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
<? for i,subeqn in ipairs(eqn.eqns) do
?>	<?=subeqn.symbols.applyInitCondCell?>(
		solver,
		initCond,
		&U-><?=subeqn.field?>,
		cell
	);
<? end
?>
}

//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?> <?=subeqnDepends'calcDTCell'?>

<? if require "hydro.solver.meshsolver":isa(solver) then
?>
//// MODULE_DEPENDS: <?=face_t?>
<? end
?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell<? --\
if require "hydro.solver.meshsolver":isa(solver) then --\
?>,\
	/*global <?=face_t?> const * const */faces,\
	/*global int const * const */cellFaceIndexes<? --\
end --\
?>\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.calcDTCell?>(\
		dt,\
		solver,\
		&(U)-><?=subeqn.field?>,\
		cell<? --\
if require "hydro.solver.meshsolver":isa(solver) then --\
?>,\
		faces,\
		cellFaceIndexes<? --\
end --\
?>\
	);\
<? end --\
?>\
}
