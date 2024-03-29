//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=cartesianToCoord?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->initCondMins, solver->initCondMaxs), .5);
	bool const lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;
	
	<?=initCode()?>

<? if eqn.usePressure then
?>	U->Pi_tt = P;
<? else
?>	U->Pi_tt = rho;
<? end
?>	U->Pi_ti = real3_zero;
	U->Pi_ij = real3s3_zero;
	U->Psi_ttk = cartesianToCoord(v, x);
	U->Psi_tik = real3x3_zero;
	U->Psi_ijk = real3x3s3_zero;
}

//// MODULE_NAME: <?=fluxFromCons?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultF,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const c = solver->wavespeed;\
	real3 const nL = normal_l1(n);\
	real3 const nU = normal_u1(n);\
\
<? --\
for _,fields in ipairs{ --\
	{"Pi_tt", "Psi_ttk", {""}}, --\
	{"Pi_ti", "Psi_tik", xNames:mapi(function(xi) return "."..xi end)}, --\
	{"Pi_ij", "Psi_ijk", symNames:mapi(function(xij) return "."..xij end)}, --\
} do --\
	local Pi, Psi, xabs = table.unpack(fields) --\
	for _,xab in ipairs(xabs) do --\
?>\
	/* F^Pi_ab = -c Psi_ab,k n^k */\
	(resultF)-><?=Pi?><?=xab?> = -c * ((U)-><?=Psi?>.x<?=xab?> * nU.x + (U)-><?=Psi?>.y<?=xab?> * nU.y + (U)-><?=Psi?>.z<?=xab?> * nU.z);\
\
	/* F^{Psi_ab,k} = -c (Pi_ab n_k) */\
<? 		for k,xk in ipairs(xNames) do --\
?>	(resultF)-><?=Psi?>.<?=xk?><?=xab?> = -c * ((U)-><?=Pi?><?=xab?> * nL.<?=xk?>);\
<? 		end --\
	end --\
end --\
?>\
}

//// MODULE_NAME: <?=eigen_forInterface?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
)

//// MODULE_NAME: <?=eigen_forCell?>

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*real3 const */n\
)

//// MODULE_NAME: <?=eigen_leftTransform?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real * Yp = (Y)->ptr;\
\
	real nLen = normal_len(n);\
	real invDenom = 1. / (nLen * nLen);\
\
<? --\
for _,fields in ipairs{ --\
	{"Pi_tt", "Psi_ttk", {""}}, --\
	{"Pi_ti", "Psi_tik", xNames:mapi(function(xi) return "."..xi end)}, --\
	{"Pi_ij", "Psi_ijk", symNames:mapi(function(xij) return "."..xij end)}, --\
} do --\
	local Pi, Psi, xabs = table.unpack(fields) --\
	for _,xab in ipairs(xabs) do --\
?>\
	Yp[0] = .5 * invDenom * (\
		(X)-><?=Pi?><?=xab?> * nLen\
		+ (X)-><?=Psi?>.x<?=xab?> * normal_u1x(n)\
		+ (X)-><?=Psi?>.y<?=xab?> * normal_u1y(n)\
		+ (X)-><?=Psi?>.z<?=xab?> * normal_u1z(n)\
	);\
	Yp[1] = invDenom * (\
		(X)-><?=Psi?>.x<?=xab?> * normal_u2x(n)\
		+ (X)-><?=Psi?>.y<?=xab?> * normal_u2y(n)\
		+ (X)-><?=Psi?>.z<?=xab?> * normal_u2z(n)\
	);\
	Yp[2] = invDenom * (\
		(X)-><?=Psi?>.x<?=xab?> * normal_u3x(n)\
		+ (X)-><?=Psi?>.y<?=xab?> * normal_u3y(n)\
		+ (X)-><?=Psi?>.z<?=xab?> * normal_u3z(n)\
	);\
	Yp[3] = .5 * invDenom * (\
		-(X)-><?=Pi?><?=xab?> * nLen\
		+ (X)-><?=Psi?>.x<?=xab?> * normal_u1x(n)\
		+ (X)-><?=Psi?>.y<?=xab?> * normal_u1y(n)\
		+ (X)-><?=Psi?>.z<?=xab?> * normal_u1z(n)\
	);\
	Yp += 4;\
<? --\
	end --\
end --\
?>\
}

//// MODULE_NAME: <?=eigen_rightTransform?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const * Xp = (X)->ptr;\
	real const nLen = normal_len(n);\
<? --\
for _,fields in ipairs{ --\
	{"Pi_tt", "Psi_ttk", {""}}, --\
	{"Pi_ti", "Psi_tik", xNames:mapi(function(xi) return "."..xi end)}, --\
	{"Pi_ij", "Psi_ijk", symNames:mapi(function(xij) return "."..xij end)}, --\
} do --\
	local Pi, Psi, xabs = table.unpack(fields) --\
	for _,xab in ipairs(xabs) do --\
?>\
	(Y)-><?=Pi?><?=xab?> = (Xp[0] - Xp[3]) * nLen;\
	(Y)-><?=Psi?>.x<?=xab?> = Xp[0] * normal_l1x(n)\
			+ Xp[1] * normal_l2x(n)\
			+ Xp[2] * normal_l3x(n)\
			+ Xp[3] * normal_l1x(n);\
	(Y)-><?=Psi?>.y<?=xab?> = Xp[0] * normal_l1y(n)\
			+ Xp[1] * normal_l2y(n)\
			+ Xp[2] * normal_l3y(n)\
			+ Xp[3] * normal_l1y(n);\
	(Y)-><?=Psi?>.z<?=xab?> = Xp[0] * normal_l1z(n)\
			+ Xp[1] * normal_l2z(n)\
			+ Xp[2] * normal_l3z(n)\
			+ Xp[3] * normal_l1z(n);\
	Xp += 4;\
<? --\
	end --\
end --\
?>\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

// by default in hydro/eqn/eqn.lua, <?=fluxFromCons?> is defined by <?=eigen_fluxTransform?>
// but since eig is empty, we can define <?=eigen_fluxTransform?> with <?=fluxFromCons?>
#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X_,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
)\
	<?=fluxFromCons?>(result, solver, X_, cell, n)

//// MODULE_NAME: <?=addSource?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	real3 const x = cellBuf[index].pos;
	
	/* TODO add stress-energy influence here */
	/* TODO TODO coupled solvers, once again, but it does seem messy */
}
