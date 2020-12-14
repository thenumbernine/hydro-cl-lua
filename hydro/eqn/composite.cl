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
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
<? for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_forInterface?>(\
		&(resultEig)->eqn<?=i?>,\
		solver,\
		&(UL)->eqn<?=i?>,\
		&(UR)->eqn<?=i?>,\
		pt,\
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
	/*normal_t */n\
) {\
<? --\
for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_leftTransform?>(\
		&(result)->eqn<?=i?>,\
		solver,\
		&(eig)->eqn<?=i?>,\
		&(X)->eqn<?=i?>,\
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
	/*normal_t const */n\
) {\
<? --\
for i,subeqn in ipairs(eqn.eqns) do --\
?>	<?=subeqn.symbols.eigen_rightTransform?>(\
		&(result)->eqn<?=i?>,\
		solver,\
		&(eig)->eqn<?=i?>,\
		&(X)->eqn<?=i?>,\
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
		&(U)->eqn<?=i?>,\
		pt,\
		dx);\
	resultName##base.eqn<?=i?> = *<?=subeqn.symbolPrefix?>resultName;\
<? end --\
?>\
}\
<?=cons_t?> * const resultName = &resultName##base;\
<? end ?>

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
		&U->eqn<?=i?>,
		cell
	);
<? end
?>
}

//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?> <?=subeqnDepends'calcDTCell'?>

<? if require "hydro.solver.meshsolver".is(solver) then
?>
//// MODULE_DEPENDS: <?=face_t?>
<? end
?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell<? --\
if require "hydro.solver.meshsolver".is(solver) then --\
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
		&(U)->eqn<?=i?>,\
		cell<? --\
if require "hydro.solver.meshsolver".is(solver) then --\
?>,\
		faces,\
		cellFaceIndexes<? --\
end --\
?>\
	);\
<? end --\
?>\
}
