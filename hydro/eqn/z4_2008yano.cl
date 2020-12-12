//// MODULE_NAME: calc_gamma_ll

#define calc_gamma_ll(U, x)	((U)->gamma_ll)

//// MODULE_NAME: calc_gamma_uu
//// MODULE_DEPENDS: <?=cons_t?>


sym3 calc_gamma_uu(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	return gamma_uu;
}

//// MODULE_NAME: setFlatSpace

void setFlatSpace(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	U->alpha = 1;
	U->gamma_ll = sym3_ident;
	U->a_l = real3_zero;
	U->d_lll.x = sym3_zero;
	U->d_lll.y = sym3_zero;
	U->d_lll.z = sym3_zero;
	U->K_ll = sym3_zero;
	U->Theta = 0;
	U->Z_l = real3_zero;
}

//// MODULE_NAME: applyInitCond
//// MODULE_DEPENDS: SETBOUNDS coordMap setFlatSpace

kernel void applyInitCond(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	real3 const x = cellBuf[index].pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=cons_t?> * const U = UBuf + index;
	setFlatSpace(solver, U, x);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;

	<?=initCode()?>

	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;
	
	//Z_u n^u = 0
	//Theta = alpha n_u Z^u = alpha Z^u
	//for n_a = (-alpha, 0)
	//n^a_l = (1/alpha, -beta^i/alpha)
	//(Z_t - Z_i beta^i) / alpha = Theta ... = ?
	//Z^t n_t + Z^i n_i = -alpha Z^t = Theta
	U->Theta = 0;
	U->Z_l = real3_zero;
}

//// MODULE_NAME: initDerivs
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> SETBOUNDS numGhost

kernel void initDerivs(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=cons_t?> * const U = UBuf + index;

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (solver->grid_dx.s<?=i-1?> * U->alpha);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> - U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>) / solver->grid_dx.s<?=i-1?>;
	<? end ?>
<? end
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>
}

//// MODULE_NAME: calcDT
//// MODULE_DEPENDS: SETBOUNDS initCond.codeprefix

kernel void calcDT(
	constant <?=solver_t?> const * const solver,
	global real * const dtBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	global <?=cons_t?> const * const U = UBuf + index;
	real const det_gamma = sym3_det(U->gamma_ll);
	real const f = calc_f(U->alpha);
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) is only called once
	real const sqrt_f = sqrt(f);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side==0 then ?>
		real const gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
		<? elseif side==1 then ?>
		real const gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
		<? elseif side==2 then ?>
		real const gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
		<? end ?>	
		real const lambdaLight = U->alpha * sqrt(gammaUjj);
		
		real const lambdaGauge = lambdaLight * sqrt_f;
		real const lambda = (real)max(lambdaGauge, lambdaLight);
		
		real const lambdaMin = (real)min((real)0., -lambda);
		real const lambdaMax = (real)max((real)0., lambda);
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt; 
}

//// MODULE_NAME: eigen_forCell
//// MODULE_DEPENDS: initCond.codeprefix

//used by PLM
#define eigen_forCell(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	(eig)->alpha = U.alpha;\
	(eig)->sqrt_f = sqrt(calc_f(U.alpha));\
	(eig)->gamma_ll = U.gamma_ll;\
	real const det_gamma = sym3_det(U.gamma_ll);\
	(eig)->gamma_uu = sym3_inv(U.gamma_ll, det_gamma);\
	(eig)->sqrt_gammaUjj = _real3(sqrt((eig)->gamma_uu.xx), sqrt((eig)->gamma_uu.yy), sqrt((eig)->gamma_uu.zz));\
}

//// MODULE_NAME: calcCellMinMaxEigenvalues
//// MODULE_DEPENDS: initCond.codeprefix

#define calcCellMinMaxEigenvalues(\
	/*range_t * const */result,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const det_gamma = sym3_det((U)->gamma_ll);\
	\
	<? if side==0 then ?>\
	real const gammaUjj = ((U)->gamma_ll.yy * (U)->gamma_ll.zz - (U)->gamma_ll.yz * (U)->gamma_ll.yz) / det_gamma;\
	<? elseif side==1 then ?>\
	real const gammaUjj = ((U)->gamma_ll.xx * (U)->gamma_ll.zz - (U)->gamma_ll.xz * (U)->gamma_ll.xz) / det_gamma;\
	<? elseif side==2 then ?>\
	real const gammaUjj = ((U)->gamma_ll.xx * (U)->gamma_ll.yy - (U)->gamma_ll.xy * (U)->gamma_ll.xy) / det_gamma;\
	<? end ?>\
	\
	real const lambdaLight = (U)->alpha * sqrt(gammaUjj);\
	\
	real const f = calc_f((U)->alpha);\
	real const lambdaGauge = lambdaLight * sqrt(f);\
\
	real const lambdaMax = max(lambdaGauge, lambdaLight);\
	/*= lambdaLight * max(sqrt(f), 1)*/\
\
	(result)->min = -lambdaMax;\
	(result)->max = lambdaMax;\
}

//// MODULE_NAME: eigen_forInterface
//// MODULE_DEPENDS: initCond.codeprefix

#define eigen_forInterface(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const alpha = .5 * ((UL)->alpha + (UR)->alpha);\
	sym3 const avg_gamma = (sym3){\
		.xx = .5 * ((UL)->gamma_ll.xx + (UR)->gamma_ll.xx),\
		.xy = .5 * ((UL)->gamma_ll.xy + (UR)->gamma_ll.xy),\
		.xz = .5 * ((UL)->gamma_ll.xz + (UR)->gamma_ll.xz),\
		.yy = .5 * ((UL)->gamma_ll.yy + (UR)->gamma_ll.yy),\
		.yz = .5 * ((UL)->gamma_ll.yz + (UR)->gamma_ll.yz),\
		.zz = .5 * ((UL)->gamma_ll.zz + (UR)->gamma_ll.zz),\
	};\
	real const det_avg_gamma = sym3_det(avg_gamma);\
	(eig)->alpha = alpha;\
	(eig)->sqrt_f = sqrt(calc_f(alpha));\
	(eig)->gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);\
	(eig)->sqrt_gammaUjj = _real3(sqrt((eig)->gamma_uu.xx), sqrt((eig)->gamma_uu.yy), sqrt((eig)->gamma_uu.zz));\
}

//// MODULE_NAME: eigen_left/rightTransform
//// MODULE_DEPENDS: rotate

#define eigen_leftTransform(\
	/*<?=waves_t?> * const */results,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */inputU,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	/* input */\
	real3 const a_l = real3_swap((inputU)->a_l, n.side);							/* 0-2 */\
	_3sym3 const d_lll = _3sym3_swap((inputU)->d_lll, n.side);					/* 3-20 */\
	sym3 const K_ll = sym3_swap((inputU)->K_ll, n.side);							/* 21-26 */\
	real const Theta = (inputU)->Theta;											/* 27 */\
	real3 const Z_l = real3_swap((inputU)->Z_l, n.side);							/* 28-30 */\
\
	/* eig */\
	real const sqrt_f = (eig)->sqrt_f;\
	real const f = sqrt_f * sqrt_f;	\
	real const _1_sqrt_f = 1. / sqrt_f;\
	real const _1_f = _1_sqrt_f * _1_sqrt_f;\
	real const fSq = f * f;\
	real const f_toThe_3_2 = f * sqrt_f;\
	real const f_toThe_5_2 = f * f * sqrt_f;\
\
	real const f_minus_1 = f - 1.;\
	real const f_minus_1_sq = f_minus_1 * f_minus_1;\
\
	sym3 const gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 const gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
	\
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);\
	\
	/* d_llu = d_ij^k = d_ijl * gamma^lk */\
	real3x3 d_llu[3] = {\
<? for i,xi in ipairs(xNames) do --\
?>		sym3_sym3_mul(d_lll.<?=xi?>, gamma_uu),\
<? end --\
?>	};\
\
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, K_ll);\
\
	real const tr_K = real3x3_trace(K_ul);\
\
	real sqrt_gUxx;\
	if (false) {}\
<? for side=0,2 do --\
?>	else if (n.side == <?=side?>) {\
		sqrt_gUxx = (eig)->sqrt_gammaUjj.s<?=side?>;\
	}\
<? end --\
?>\
	real const _1_sqrt_gUxx = 1. / sqrt_gUxx;\
	real const _1_gammaUxx = _1_sqrt_gUxx * _1_sqrt_gUxx;\
\
	real const a_u_x = a_l.x * gamma_uu.xx + a_l.y * gamma_uu.xy + a_l.z * gamma_uu.xz;\
	real const Z_u_x = Z_l.x * gamma_uu.xx + Z_l.y * gamma_uu.xy + Z_l.z * gamma_uu.xz;\
\
	/* d_i = d_ijk gamma^jk */\
	real3 const d_l = _real3(\
		sym3_dot(d_lll.x, gamma_uu),\
		sym3_dot(d_lll.y, gamma_uu),\
		sym3_dot(d_lll.z, gamma_uu));\
\
	real const d_u_x = d_l.x * gamma_uu.xx + d_l.y * gamma_uu.xy + d_l.z * gamma_uu.xz;\
\
	/* e_i = d_jki gamma^jk */\
	real3 const e_l = (real3){\
<? for i,xi in ipairs(xNames) do --\
?>		.<?=xi?> = 0. <? --\
	for j,xj in ipairs(xNames) do --\
		for k,xk in ipairs(xNames) do --\
?> + d_lll.<?=xj?>.<?=sym(k,i)?> * gamma_uu.<?=sym(j,k)?><? --\
		end --\
	end --\
?>,\
<? end --\
?>	};\
\
	real const e_u_x = e_l.x * gamma_uu.xx + e_l.y * gamma_uu.xy + e_l.z * gamma_uu.xz;\
\
	real const m = solver->m;\
\
	(results)->ptr[0] = (\
		-(2. * gamma_uu.xy * a_l.y * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * a_l.y * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xy * a_l.y * f * _1_sqrt_gUxx \
		+ 8. * gamma_uu.xy * K_ll.xy * f_toThe_3_2 \
		- 4. * gamma_uu.xy * K_ll.xy * sqrt_f \
		- 4. * gamma_uu.xy * K_ll.xy * f_toThe_5_2 \
		+ 8. * gamma_uu.xz * K_ll.xz * f_toThe_3_2 \
		- 4. * gamma_uu.xz * K_ll.xz * f_toThe_5_2 \
		- 4. * gamma_uu.xz * K_ll.xz * sqrt_f \
		+ 2. * gamma_uu.xz * a_l.z * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xz * a_l.z * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xz * a_l.z * f * _1_sqrt_gUxx \
		+ 4. * gamma_uu.yy * K_ll.yy * f_toThe_3_2 \
		- 2. * gamma_uu.yy * K_ll.yy * sqrt_f \
		- 2. * gamma_uu.yy * K_ll.yy * f_toThe_5_2 \
		- 4. * gamma_uu.yz * K_ll.yz * f_toThe_5_2 \
		- 4. * gamma_uu.yz * K_ll.yz * sqrt_f \
		+ 8. * gamma_uu.yz * K_ll.yz * f_toThe_3_2 \
		+ 4. * gamma_uu.zz * K_ll.zz * f_toThe_3_2 \
		- 2. * gamma_uu.zz * K_ll.zz * sqrt_f \
		- 2. * gamma_uu.zz * K_ll.zz * f_toThe_5_2 \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		- 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 8. * f * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 8. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		+ 8. * f * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 8. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * Z_l.x * sqrt_gUxx \
		- 4. * fSq * Z_l.x * sqrt_gUxx \
		+ 2. * f_toThe_5_2 * m * Theta \
		- 2. * f_toThe_3_2 * m * Theta \
		- 2. * f * m * Z_l.x * sqrt_gUxx \
		+ 2. * fSq * m * Z_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx * fSq \
		- 4. * a_l.x * sqrt_gUxx * f \
		+ 4. * K_ll.xx * gamma_uu.xx * f_toThe_3_2 \
		- 2. * K_ll.xx * gamma_uu.xx * sqrt_f \
		- 2. * K_ll.xx * f_toThe_5_2 * gamma_uu.xx \
		+ 4. * Theta * sqrt_f \
		- 4. * Theta * f_toThe_3_2)) / (gamma_uu.xx * sqrt_f * (4. \
		- 8. * f \
		+ 4. * fSq));\
	(results)->ptr[1] = (2. * gamma_uu.xx * sqrt_gUxx * K_ll.xz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.zz \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.z.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yz * sqrt_gUxx \
		+ 2. * gamma_uu.xz * K_ll.zz * sqrt_gUxx \
		+ gamma_uu.yy * gamma_uu.xx * d_lll.y.yz \
		- gamma_uu.yy * gamma_uu.xx * d_lll.z.yy \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.y.zz \
		- gamma_uu.yz * gamma_uu.xx * d_lll.z.yz \
		- a_l.z * gamma_uu.xx \
		+ 2. * Z_l.z * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[2] = (2. * gamma_uu.xx * sqrt_gUxx * K_ll.xy \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yy * sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.yz \
		+ 2. * gamma_uu.xz * K_ll.yz * sqrt_gUxx \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.zz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.zz * gamma_uu.xx * d_lll.z.yz \
		- a_l.y * gamma_uu.xx \
		+ 2. * Z_l.y * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[3] = (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		+ gamma_uu.xy * a_l.y \
		- 2. * gamma_uu.xy * K_ll.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		- gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		- 2. * gamma_uu.xz * K_ll.xz * sqrt_gUxx \
		+ gamma_uu.xz * a_l.z \
		- 2. * gamma_uu.yy * K_ll.yy * sqrt_gUxx \
		- 4. * gamma_uu.yz * K_ll.yz * sqrt_gUxx \
		- 2. * gamma_uu.zz * K_ll.zz * sqrt_gUxx \
		+ 2. * Theta * sqrt_gUxx \
		+ 2. * Z_l.x * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[4] = (gamma_uu.xx * gamma_uu.yy * d_lll.y.xy \
		- gamma_uu.xx * gamma_uu.yy * d_lll.x.yy \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.y.xz \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.z.xy \
		- 2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz \
		+ gamma_uu.xx * gamma_uu.zz * d_lll.z.xz \
		- gamma_uu.xx * gamma_uu.zz * d_lll.x.zz \
		+ gamma_uu.xx * Z_l.x \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.xz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.xy \
		+ 2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		- gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.y.xy \
		+ gamma_uu.xy * Z_l.y \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.x.yy \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		- gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.xz \
		+ gamma_uu.xz * Z_l.z \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.x.zz \
		+ Theta * sqrt_gUxx) / (2. * sqrt_gUxx);\
	(results)->ptr[5] = (gamma_uu.xx * d_lll.z.xz \
		- gamma_uu.xx * d_lll.x.zz \
		- gamma_uu.xy * d_lll.y.zz \
		+ gamma_uu.xy * d_lll.z.yz \
		+ K_ll.zz * sqrt_gUxx) / (2. * sqrt_gUxx);\
	(results)->ptr[6] = (gamma_uu.xx * d_lll.y.xz \
		+ gamma_uu.xx * d_lll.z.xy \
		- 2. * gamma_uu.xx * d_lll.x.yz \
		- gamma_uu.xy * d_lll.y.yz \
		+ gamma_uu.xy * d_lll.z.yy \
		+ gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * d_lll.z.yz \
		+ 2. * K_ll.yz * sqrt_gUxx) / (4. * sqrt_gUxx);\
	(results)->ptr[7] = (\
		-(gamma_uu.xx * gamma_uu.yy * d_lll.y.xy * f \
		- gamma_uu.xx * gamma_uu.yy * d_lll.x.yy * f \
		- gamma_uu.xx * gamma_uu.yy * m * d_lll.y.xy * f \
		+ gamma_uu.xx * gamma_uu.yy * m * d_lll.x.yy * f \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.y.xz * f \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.z.xy * f \
		- 2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz * f \
		- gamma_uu.xx * gamma_uu.yz * m * d_lll.y.xz * f \
		- gamma_uu.xx * gamma_uu.yz * m * d_lll.z.xy * f \
		+ 2. * gamma_uu.xx * gamma_uu.yz * m * d_lll.x.yz * f \
		+ gamma_uu.xx * gamma_uu.zz * d_lll.z.xz * f \
		- gamma_uu.xx * gamma_uu.zz * d_lll.x.zz * f \
		- gamma_uu.xx * gamma_uu.zz * m * d_lll.z.xz * f \
		+ gamma_uu.xx * gamma_uu.zz * m * d_lll.x.zz * f \
		- 2. * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * f \
		- 2. * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * f \
		+ 4. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * f \
		+ gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * f \
		+ gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * f \
		- 2. * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * f \
		- 2. * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * f \
		+ 2. * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * f \
		+ gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * f \
		- gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * f \
		- 2. * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * f \
		+ 2. * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * f \
		+ gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * f \
		- gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * f \
		- 2. * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * f \
		+ gamma_uu.xy * a_l.y \
		+ 2. * gamma_uu.xy * Z_l.y * f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * f \
		+ gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * f \
		- gamma_uu.xy * m * Z_l.y * f \
		- gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * f \
		+ 2. * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * f \
		- 2. * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * f \
		- gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * f \
		+ gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * f \
		+ 2. * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * f \
		- 2. * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * f \
		- gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * f \
		+ gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * f \
		- 2. * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * f \
		+ gamma_uu.xz * a_l.z \
		+ 2. * gamma_uu.xz * Z_l.z * f \
		+ 2. * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * f \
		+ gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * f \
		- gamma_uu.xz * m * Z_l.z * f \
		- gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * f \
		- f * gamma_uu.xy * a_l.y \
		- f * gamma_uu.xz * a_l.z \
		+ a_l.x * gamma_uu.xx \
		- d_lll.x.xx * gamma_uu.xx * gamma_uu.xx * f \
		- m * Z_l.x * gamma_uu.xx * f)) / (gamma_uu.xx * gamma_uu.xx * f);\
	(results)->ptr[8] = (\
		-(2. * gamma_uu.xy * d_lll.y.xy \
		- 2. * gamma_uu.xy * d_lll.x.yy \
		+ gamma_uu.xz * d_lll.y.xz \
		+ gamma_uu.xz * d_lll.z.xy \
		- 2. * gamma_uu.xz * d_lll.x.yz \
		+ gamma_uu.yz * d_lll.y.yz \
		- gamma_uu.yz * d_lll.z.yy \
		+ gamma_uu.zz * d_lll.y.zz \
		- gamma_uu.zz * d_lll.z.yz \
		+ a_l.y \
		- 2. * Z_l.y \
		- 2. * d_lll.x.xy * gamma_uu.xx)) / (2. * gamma_uu.xx);\
	(results)->ptr[9] = (\
		-(gamma_uu.xy * d_lll.y.xz \
		+ gamma_uu.xy * d_lll.z.xy \
		- 2. * gamma_uu.xy * d_lll.x.yz \
		+ 2. * gamma_uu.xz * d_lll.z.xz \
		- 2. * gamma_uu.xz * d_lll.x.zz \
		- gamma_uu.yy * d_lll.y.yz \
		+ gamma_uu.yy * d_lll.z.yy \
		- gamma_uu.yz * d_lll.y.zz \
		+ gamma_uu.yz * d_lll.z.yz \
		+ a_l.z \
		- 2. * Z_l.z \
		- 2. * d_lll.x.xz * gamma_uu.xx)) / (2. * gamma_uu.xx);\
	(results)->ptr[10] = d_lll.y.xx;\
	(results)->ptr[11] = d_lll.y.xy;\
	(results)->ptr[12] = d_lll.y.xz;\
	(results)->ptr[13] = d_lll.y.yy;\
	(results)->ptr[14] = d_lll.y.yz;\
	(results)->ptr[15] = d_lll.y.zz;\
	(results)->ptr[16] = d_lll.z.xx;\
	(results)->ptr[17] = d_lll.z.xy;\
	(results)->ptr[18] = d_lll.z.xz;\
	(results)->ptr[19] = d_lll.z.yy;\
	(results)->ptr[20] = d_lll.z.yz;\
	(results)->ptr[21] = d_lll.z.zz;\
	(results)->ptr[22] = (\
		-(gamma_uu.xy * gamma_uu.xz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yy \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.yz \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.zz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.zz * gamma_uu.xx * d_lll.z.yz \
		- a_l.y * gamma_uu.xx)) / (2. * gamma_uu.xx);\
	(results)->ptr[23] = (gamma_uu.xy * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yz \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.z.yy \
		- gamma_uu.yy * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yy * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yz \
		+ a_l.z * gamma_uu.xx) / (2. * gamma_uu.xx);\
	(results)->ptr[24] = (\
		-(2. * gamma_uu.xx * sqrt_gUxx * K_ll.xz \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yz \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yz * sqrt_gUxx \
		+ 2. * gamma_uu.xz * K_ll.zz * sqrt_gUxx \
		- gamma_uu.yy * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yy * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yz \
		+ a_l.z * gamma_uu.xx \
		- 2. * Z_l.z * gamma_uu.xx)) / (4. * gamma_uu.xx);\
	(results)->ptr[25] = (\
		-(2. * gamma_uu.xx * sqrt_gUxx * K_ll.xy \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.xz * d_lll.y.zz \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.z.yz \
		+ 2. * gamma_uu.xz * K_ll.yz * sqrt_gUxx \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.y.yz \
		- gamma_uu.yz * gamma_uu.xx * d_lll.z.yy \
		+ gamma_uu.zz * gamma_uu.xx * d_lll.y.zz \
		- gamma_uu.zz * gamma_uu.xx * d_lll.z.yz \
		+ a_l.y * gamma_uu.xx \
		- 2. * Z_l.y * gamma_uu.xx)) / (4. * gamma_uu.xx);\
	(results)->ptr[26] = (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		+ gamma_uu.xy * a_l.y \
		+ 2. * gamma_uu.xy * K_ll.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		- gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		+ 2. * gamma_uu.xz * K_ll.xz * sqrt_gUxx \
		+ gamma_uu.xz * a_l.z \
		+ 2. * gamma_uu.yy * K_ll.yy * sqrt_gUxx \
		+ 4. * gamma_uu.yz * K_ll.yz * sqrt_gUxx \
		+ 2. * gamma_uu.zz * K_ll.zz * sqrt_gUxx \
		- 2. * Theta * sqrt_gUxx \
		+ 2. * Z_l.x * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[27] = (\
		-(gamma_uu.xx * gamma_uu.yy * d_lll.y.xy \
		- gamma_uu.xx * gamma_uu.yy * d_lll.x.yy \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.y.xz \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.z.xy \
		- 2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz \
		+ gamma_uu.xx * gamma_uu.zz * d_lll.z.xz \
		- gamma_uu.xx * gamma_uu.zz * d_lll.x.zz \
		+ gamma_uu.xx * Z_l.x \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.xz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.xy \
		+ 2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		- gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.y.xy \
		+ gamma_uu.xy * Z_l.y \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.x.yy \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		- gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.xz \
		+ gamma_uu.xz * Z_l.z \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.x.zz \
		- Theta * sqrt_gUxx)) / (2. * sqrt_gUxx);\
	(results)->ptr[28] = (\
		-(gamma_uu.xx * d_lll.z.xz \
		- gamma_uu.xx * d_lll.x.zz \
		- gamma_uu.xy * d_lll.y.zz \
		+ gamma_uu.xy * d_lll.z.yz \
		- K_ll.zz * sqrt_gUxx)) / (2. * sqrt_gUxx);\
	(results)->ptr[29] = (\
		-(gamma_uu.xx * d_lll.y.xz \
		+ gamma_uu.xx * d_lll.z.xy \
		- 2. * gamma_uu.xx * d_lll.x.yz \
		- gamma_uu.xy * d_lll.y.yz \
		+ gamma_uu.xy * d_lll.z.yy \
		+ gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * d_lll.z.yz \
		- 2. * K_ll.yz * sqrt_gUxx)) / (4. * sqrt_gUxx);\
	(results)->ptr[30] = (2. * gamma_uu.xy * a_l.y * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * a_l.y * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xy * a_l.y * f * _1_sqrt_gUxx \
		- 8. * gamma_uu.xy * K_ll.xy * f_toThe_3_2 \
		+ 4. * gamma_uu.xy * K_ll.xy * sqrt_f \
		+ 4. * gamma_uu.xy * K_ll.xy * f_toThe_5_2 \
		- 8. * gamma_uu.xz * K_ll.xz * f_toThe_3_2 \
		+ 4. * gamma_uu.xz * K_ll.xz * f_toThe_5_2 \
		+ 4. * gamma_uu.xz * K_ll.xz * sqrt_f \
		+ 2. * gamma_uu.xz * a_l.z * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xz * a_l.z * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xz * a_l.z * f * _1_sqrt_gUxx \
		- 4. * gamma_uu.yy * K_ll.yy * f_toThe_3_2 \
		+ 2. * gamma_uu.yy * K_ll.yy * sqrt_f \
		+ 2. * gamma_uu.yy * K_ll.yy * f_toThe_5_2 \
		+ 4. * gamma_uu.yz * K_ll.yz * f_toThe_5_2 \
		+ 4. * gamma_uu.yz * K_ll.yz * sqrt_f \
		- 8. * gamma_uu.yz * K_ll.yz * f_toThe_3_2 \
		- 4. * gamma_uu.zz * K_ll.zz * f_toThe_3_2 \
		+ 2. * gamma_uu.zz * K_ll.zz * sqrt_f \
		+ 2. * gamma_uu.zz * K_ll.zz * f_toThe_5_2 \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		- 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 8. * f * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 8. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		+ 8. * f * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 8. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * Z_l.x * sqrt_gUxx \
		- 4. * fSq * Z_l.x * sqrt_gUxx \
		- 2. * f_toThe_5_2 * m * Theta \
		+ 2. * f_toThe_3_2 * m * Theta \
		- 2. * f * m * Z_l.x * sqrt_gUxx \
		+ 2. * fSq * m * Z_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx * fSq \
		- 4. * a_l.x * sqrt_gUxx * f \
		- 4. * K_ll.xx * gamma_uu.xx * f_toThe_3_2 \
		+ 2. * K_ll.xx * gamma_uu.xx * sqrt_f \
		+ 2. * K_ll.xx * f_toThe_5_2 * gamma_uu.xx \
		- 4. * Theta * sqrt_f \
		+ 4. * Theta * f_toThe_3_2) / (gamma_uu.xx * sqrt_f * (4. \
		- 8. * f \
		+ 4. * fSq));\
}

#define eigen_rightTransform(\
	/*<?=cons_t?> * const */resultU,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */input,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	for (int j = 0; j < numStates; ++j) {\
		(resultU)->ptr[j] = 0;\
	}\
	\
	sym3 gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
\
	real sqrt_f = (eig)->sqrt_f;\
	real f = sqrt_f * sqrt_f;\
	real fSq = f * f;\
	real f_toThe_3_2 = f * sqrt_f;\
\
	real f_minus_1 = f - 1.;\
	real sqrt_gUxx = sqrt(gamma_uu.xx);\
	real _1_sqrt_gUxx = 1. / sqrt_gUxx;\
	real _1_gammaUxx = _1_sqrt_gUxx * _1_sqrt_gUxx;\
	\
	real m = solver->m;\
\
	(resultU)->a_l.x = (gamma_uu.xy * gamma_uu.yz * (input)->ptr[14] * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yz * (input)->ptr[14] * f * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yz * (input)->ptr[19] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.yz * (input)->ptr[19] * f * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.zz * (input)->ptr[15] * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.zz * (input)->ptr[15] * f * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.zz * (input)->ptr[20] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.zz * (input)->ptr[20] * f * _1_sqrt_gUxx \
		- 2. * gamma_uu.xy * (input)->ptr[22] * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * (input)->ptr[22] * f * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[14] * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[14] * f * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[19] * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[19] * f * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yz * (input)->ptr[15] * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yz * (input)->ptr[15] * f * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yz * (input)->ptr[20] * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yz * (input)->ptr[20] * f * _1_sqrt_gUxx \
		- 2. * gamma_uu.xz * (input)->ptr[23] * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xz * (input)->ptr[23] * f * _1_sqrt_gUxx \
		+ f_toThe_3_2 * gamma_uu.xx * (input)->ptr[0] \
		- sqrt_f * gamma_uu.xx * (input)->ptr[0] \
		+ sqrt_f * gamma_uu.xx * (input)->ptr[30] \
		- f_toThe_3_2 * gamma_uu.xx * (input)->ptr[30] \
		+ 2. * f * (input)->ptr[27] \
		- 2. * f * (input)->ptr[4] \
		- f * m * (input)->ptr[27] \
		+ f * m * (input)->ptr[4]) / (sqrt_gUxx * (1. \
		- f));\
	(resultU)->a_l.y = (gamma_uu.xy * gamma_uu.xz * (input)->ptr[14] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[19] \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[15] \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[20] \
		- gamma_uu.yz * gamma_uu.xx * (input)->ptr[14] \
		+ gamma_uu.yz * gamma_uu.xx * (input)->ptr[19] \
		- gamma_uu.zz * gamma_uu.xx * (input)->ptr[15] \
		+ gamma_uu.zz * gamma_uu.xx * (input)->ptr[20] \
		+ 2. * (input)->ptr[22] * gamma_uu.xx) / gamma_uu.xx;\
	(resultU)->a_l.z = (\
		-(gamma_uu.xy * gamma_uu.xz * (input)->ptr[15] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[20] \
		+ gamma_uu.xy * gamma_uu.xy * (input)->ptr[14] \
		- gamma_uu.xy * gamma_uu.xy * (input)->ptr[19] \
		- gamma_uu.yy * gamma_uu.xx * (input)->ptr[14] \
		+ gamma_uu.yy * gamma_uu.xx * (input)->ptr[19] \
		- gamma_uu.yz * gamma_uu.xx * (input)->ptr[15] \
		+ gamma_uu.yz * gamma_uu.xx * (input)->ptr[20] \
		- 2. * (input)->ptr[23] * gamma_uu.xx)) / gamma_uu.xx;\
	(resultU)->d_lll.x.xx = (gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f_toThe_3_2 * (input)->ptr[4] \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f_toThe_3_2 * m * (input)->ptr[4] \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * sqrt_f \
		- sqrt_gUxx * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] * _1_sqrt_gUxx * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] * f_toThe_3_2 * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * f_toThe_3_2 * sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] * _1_sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] * f_toThe_3_2 * _1_sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * f_toThe_3_2 * sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * f_toThe_3_2 \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * sqrt_f * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * f_toThe_3_2 \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * sqrt_f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * _1_sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f_toThe_3_2 * _1_sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * _1_sqrt_gUxx * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f_toThe_3_2 * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * _1_sqrt_gUxx * sqrt_f \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f_toThe_3_2 * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * _1_sqrt_gUxx * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f_toThe_3_2 * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f_toThe_3_2 * (input)->ptr[27] * _1_sqrt_gUxx \
		+ 3. * gamma_uu.xy * gamma_uu.xy * f_toThe_3_2 * (input)->ptr[27] * sqrt_gUxx * gamma_uu.yy \
		- 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.yy * sqrt_f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * f_toThe_3_2 \
		+ 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * f_toThe_3_2 * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		+ 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.yy * f_toThe_3_2 * gamma_uu.xx \
		- 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.yy * sqrt_f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f_toThe_3_2 \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * sqrt_f * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * sqrt_f * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * sqrt_f * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * sqrt_gUxx * sqrt_f \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * f_toThe_3_2 * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * f_toThe_3_2 * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * sqrt_gUxx * sqrt_f \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * f_toThe_3_2 * sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * f_toThe_3_2 * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] * sqrt_f * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] * f_toThe_3_2 * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] * f_toThe_3_2 * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] * sqrt_f * gamma_uu.xx * gamma_uu.xx \
		- f_toThe_3_2 * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] \
		+ 2. * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * _1_sqrt_gUxx \
		- 3. * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * sqrt_gUxx * gamma_uu.yy \
		- f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * m * (input)->ptr[4] * _1_sqrt_gUxx \
		+ 2. * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * m * (input)->ptr[4] * sqrt_gUxx * gamma_uu.yy \
		+ f_toThe_3_2 * m * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] \
		+ f_toThe_3_2 * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * _1_sqrt_gUxx \
		- 2. * f_toThe_3_2 * m * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx * gamma_uu.yy \
		+ (input)->ptr[0] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * (input)->ptr[0] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ 2. * (input)->ptr[0] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * f \
		- (input)->ptr[0] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f \
		+ (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		- (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f \
		- (input)->ptr[30] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * (input)->ptr[30] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		- (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ (input)->ptr[30] * sqrt_gUxx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * (input)->ptr[30] * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.yy * gamma_uu.yy \
		- (input)->ptr[7] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * sqrt_f \
		+ (input)->ptr[7] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * f_toThe_3_2 * gamma_uu.yy * gamma_uu.yy \
		+ 2. * (input)->ptr[7] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * sqrt_f * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[7] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_f * gamma_uu.xx \
		+ (input)->ptr[7] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		- 2. * (input)->ptr[7] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx) / (\
		-gamma_uu.xx * sqrt_f * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ 2. * gamma_uu.xx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.xx * f * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->d_lll.x.xy = (2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] \
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		- 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] \
		+ 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] \
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] \
		+ gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] \
		+ gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] \
		+ gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx \
		+ gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[27] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx \
		+ gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[4] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] \
		- gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] \
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * (input)->ptr[25] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy \
		- (input)->ptr[8] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy \
		- (input)->ptr[8] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ 2. * (input)->ptr[8] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy) / (\
		-gamma_uu.xx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->d_lll.x.xz = (\
		-(gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		- gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		+ (input)->ptr[1] \
		+ (input)->ptr[24] \
		- (input)->ptr[9] * gamma_uu.xx)) / gamma_uu.xx;\
	(resultU)->d_lll.x.yy = (2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		- gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.yy \
		- gamma_uu.xx * (input)->ptr[27] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * (input)->ptr[27] * gamma_uu.yy \
		- gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.yy \
		+ gamma_uu.xx * (input)->ptr[4] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * (input)->ptr[4] * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] \
		+ 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] * gamma_uu.xx * gamma_uu.yy \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] \
		- 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[25] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[25] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[2] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[2] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * (input)->ptr[14] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * (input)->ptr[14] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xz * (input)->ptr[14] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xz * (input)->ptr[19] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xz * (input)->ptr[19] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		- gamma_uu.xz * (input)->ptr[19] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ gamma_uu.xz * (input)->ptr[1] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xz * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx \
		- gamma_uu.xz * (input)->ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * (input)->ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy \
		+ (input)->ptr[11] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * (input)->ptr[11] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ (input)->ptr[11] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy) / (sqrt_gUxx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->d_lll.x.yz = (\
		-(gamma_uu.xy * (input)->ptr[14] * _1_sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[19] * _1_sqrt_gUxx \
		- gamma_uu.xz * (input)->ptr[15] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[20] * _1_sqrt_gUxx \
		- (input)->ptr[12] * sqrt_gUxx \
		- (input)->ptr[17] * sqrt_gUxx \
		- 2. * (input)->ptr[29] \
		+ 2. * (input)->ptr[6])) / (2. * sqrt_gUxx);\
	(resultU)->d_lll.x.zz = (\
		-(gamma_uu.xy * (input)->ptr[15] * _1_sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[20] * _1_sqrt_gUxx \
		- (input)->ptr[18] * sqrt_gUxx \
		- (input)->ptr[28] \
		+ (input)->ptr[5])) / sqrt_gUxx;\
	(resultU)->d_lll.y.xx = (input)->ptr[10];\
	(resultU)->d_lll.y.xy = (input)->ptr[11];\
	(resultU)->d_lll.y.xz = (input)->ptr[12];\
	(resultU)->d_lll.y.yy = (input)->ptr[13];\
	(resultU)->d_lll.y.yz = (input)->ptr[14];\
	(resultU)->d_lll.y.zz = (input)->ptr[15];\
	(resultU)->d_lll.z.xx = (input)->ptr[16];\
	(resultU)->d_lll.z.xy = (input)->ptr[17];\
	(resultU)->d_lll.z.xz = (input)->ptr[18];\
	(resultU)->d_lll.z.yy = (input)->ptr[19];\
	(resultU)->d_lll.z.yz = (input)->ptr[20];\
	(resultU)->d_lll.z.zz = (input)->ptr[21];\
	(resultU)->K_ll.xx = (gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] * f \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] * f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] * f \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * f * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * f * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] * f \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * f * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * sqrt_gUxx * f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * sqrt_gUxx * f \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * f * sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * f * gamma_uu.xx * sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy * gamma_uu.yy * f * gamma_uu.xx * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy * gamma_uu.yy * f * gamma_uu.xx * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[27] \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[27] * gamma_uu.yy * gamma_uu.xx \
		+ 4. * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[4] \
		+ 4. * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[4] * gamma_uu.yy * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * sqrt_gUxx \
		- 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * sqrt_gUxx * f \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * sqrt_gUxx \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f * sqrt_gUxx \
		- 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		- 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		+ 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		+ 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * f \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * f * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * f * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- f * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- f * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] \
		- 3. * f * m * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 3. * f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.xx \
		- f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] \
		- 3. * f * m * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 3. * f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.xx \
		+ f * m * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		+ f * m * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy \
		+ (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f \
		+ (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.xx \
		+ 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy \
		+ (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f \
		+ (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.xx \
		- 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		+ 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx) \
		/ (gamma_uu.xx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 3. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ 3. * gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f \
		+ 3. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy \
		- 3. * gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f));\
	(resultU)->K_ll.xy = (2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] \
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		+ 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] \
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * (input)->ptr[26] * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * (input)->ptr[3] * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] \
		- gamma_uu.xx * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] \
		+ gamma_uu.xx * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] \
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] \
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] \
		+ gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[27] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[4] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] \
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] \
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx) / (\
		-sqrt_gUxx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->K_ll.xz = (\
		-(gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		- (input)->ptr[1] \
		+ (input)->ptr[24])) / sqrt_gUxx;\
	(resultU)->K_ll.yy = (sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[25] \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * (input)->ptr[25] * gamma_uu.yy \
		- sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[2] \
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * (input)->ptr[2] * gamma_uu.yy \
		- sqrt_gUxx * gamma_uu.xz * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.xz * (input)->ptr[1] * gamma_uu.yy \
		+ sqrt_gUxx * gamma_uu.xz * (input)->ptr[24] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.xz * (input)->ptr[24] * gamma_uu.yy \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		+ gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.yy \
		+ gamma_uu.xx * (input)->ptr[27] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * (input)->ptr[27] * gamma_uu.yy \
		- gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.yy \
		+ gamma_uu.xx * (input)->ptr[4] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * (input)->ptr[4] * gamma_uu.yy \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] \
		- 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] * gamma_uu.xx * gamma_uu.yy \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] \
		- 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy) / (\
		-(gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.yy));\
	(resultU)->K_ll.yz = (input)->ptr[29] \
		+ (input)->ptr[6];\
	(resultU)->K_ll.zz = (input)->ptr[28] \
		+ (input)->ptr[5];\
	(resultU)->Theta = (input)->ptr[27] \
		+ (input)->ptr[4];\
	(resultU)->Z_l.x = (\
		-(gamma_uu.xy * (input)->ptr[22] \
		+ gamma_uu.xz * (input)->ptr[23] \
		- (input)->ptr[26] * gamma_uu.xx \
		- (input)->ptr[3] * gamma_uu.xx)) / gamma_uu.xx;\
	(resultU)->Z_l.y = (input)->ptr[22] \
		+ (input)->ptr[25] \
		+ (input)->ptr[2];\
	(resultU)->Z_l.z = (input)->ptr[1] \
		+ (input)->ptr[23] \
		+ (input)->ptr[24];\
\
	(resultU)->a_l = real3_swap((resultU)->a_l, n.side);							/* 0-2 */\
	(resultU)->d_lll = _3sym3_swap((resultU)->d_lll, n.side);						/* 3-20 */\
	(resultU)->K_ll = sym3_swap((resultU)->K_ll, n.side);							/* 21-26 */\
	(resultU)->Z_l = real3_swap((resultU)->Z_l, n.side);							/* 28-30 */\
}

//// MODULE_NAME: eigen_fluxTransform
//// MODULE_DEPENDS: rotate

//notice this paper uses the decomposition alpha A = R Lambda L
// so this computation is for alpha A
#define eigen_fluxTransform(\
	/*<?=cons_t?> * const */results,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */input,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	for (int i = 0; i < numStates; ++i) {\
		(results)->ptr[i] = 0;\
	}\
\
	/* input */\
	real3 a_l = real3_swap((input)->a_l, n.side);							/* 0-2 */\
	_3sym3 d_lll = _3sym3_swap((input)->d_lll, n.side);					/* 3-20 */\
	sym3 K_ll = sym3_swap((input)->K_ll, n.side);							/* 21-26 */\
	real Theta = (input)->Theta;											/* 27 */\
	real3 Z_l = real3_swap((input)->Z_l, n.side);							/* 28-30 */\
\
	sym3 gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
	\
	real sqrt_f = (eig)->sqrt_f;\
	real f = sqrt_f * sqrt_f;\
	\
	real m = solver->m;\
\
	/* TODO this still needs to know what vars are going into it */\
	/* so it can use the correctly swapped vars */\
\
	(results)->ptr[0] = f * (gamma_uu.xx * (input)->ptr[21] \
		+ 2. * gamma_uu.xy * (input)->ptr[22] \
		+ 2. * gamma_uu.xz * (input)->ptr[23] \
		+ gamma_uu.yy * (input)->ptr[24] \
		+ 2. * gamma_uu.yz * (input)->ptr[25] \
		+ gamma_uu.zz * (input)->ptr[26] \
		- m * (input)->ptr[27]);\
	(results)->ptr[1] = 0.;\
	(results)->ptr[2] = 0.;\
	(results)->ptr[3] = (input)->ptr[21];\
	(results)->ptr[4] = (input)->ptr[22];\
	(results)->ptr[5] = (input)->ptr[23];\
	(results)->ptr[6] = (input)->ptr[24];\
	(results)->ptr[7] = (input)->ptr[25];\
	(results)->ptr[8] = (input)->ptr[26];\
	(results)->ptr[9] = 0.;\
	(results)->ptr[10] = 0.;\
	(results)->ptr[11] = 0.;\
	(results)->ptr[12] = 0.;\
	(results)->ptr[13] = 0.;\
	(results)->ptr[14] = 0.;\
	(results)->ptr[15] = 0.;\
	(results)->ptr[16] = 0.;\
	(results)->ptr[17] = 0.;\
	(results)->ptr[18] = 0.;\
	(results)->ptr[19] = 0.;\
	(results)->ptr[20] = 0.;\
	(results)->ptr[21] = \
		-(gamma_uu.yy * (input)->ptr[10] \
		- gamma_uu.yy * (input)->ptr[6] \
		+ gamma_uu.yz * (input)->ptr[11] \
		+ gamma_uu.yz * (input)->ptr[16] \
		- 2. * gamma_uu.yz * (input)->ptr[7] \
		+ gamma_uu.zz * (input)->ptr[17] \
		- gamma_uu.zz * (input)->ptr[8] \
		- (input)->ptr[0] \
		+ 2. * (input)->ptr[28]);\
	(results)->ptr[22] = (2. * gamma_uu.xy * (input)->ptr[10] \
		- 2. * gamma_uu.xy * (input)->ptr[6] \
		+ gamma_uu.xz * (input)->ptr[11] \
		+ gamma_uu.xz * (input)->ptr[16] \
		- 2. * gamma_uu.xz * (input)->ptr[7] \
		+ gamma_uu.yz * (input)->ptr[13] \
		- gamma_uu.yz * (input)->ptr[18] \
		+ gamma_uu.zz * (input)->ptr[14] \
		- gamma_uu.zz * (input)->ptr[19] \
		+ (input)->ptr[1] \
		- 2. * (input)->ptr[29]) / 2.;\
	(results)->ptr[23] = (gamma_uu.xy * (input)->ptr[11] \
		+ gamma_uu.xy * (input)->ptr[16] \
		- 2. * gamma_uu.xy * (input)->ptr[7] \
		+ 2. * gamma_uu.xz * (input)->ptr[17] \
		- 2. * gamma_uu.xz * (input)->ptr[8] \
		- gamma_uu.yy * (input)->ptr[13] \
		+ gamma_uu.yy * (input)->ptr[18] \
		- gamma_uu.yz * (input)->ptr[14] \
		+ gamma_uu.yz * (input)->ptr[19] \
		+ (input)->ptr[2] \
		- 2. * (input)->ptr[30]) / 2.;\
	(results)->ptr[24] = \
		-(gamma_uu.xx * (input)->ptr[10] \
		- gamma_uu.xx * (input)->ptr[6] \
		+ gamma_uu.xz * (input)->ptr[13] \
		- gamma_uu.xz * (input)->ptr[18]);\
	(results)->ptr[25] = (\
		-(gamma_uu.xx * (input)->ptr[11] \
		+ gamma_uu.xx * (input)->ptr[16] \
		- 2. * gamma_uu.xx * (input)->ptr[7] \
		- gamma_uu.xy * (input)->ptr[13] \
		+ gamma_uu.xy * (input)->ptr[18] \
		+ gamma_uu.xz * (input)->ptr[14] \
		- gamma_uu.xz * (input)->ptr[19])) / 2.;\
	(results)->ptr[26] = \
		-(gamma_uu.xx * (input)->ptr[17] \
		- gamma_uu.xx * (input)->ptr[8] \
		- gamma_uu.xy * (input)->ptr[14] \
		+ gamma_uu.xy * (input)->ptr[19]);\
	(results)->ptr[27] = \
		-(gamma_uu.xx * gamma_uu.yy * (input)->ptr[10] \
		- gamma_uu.xx * gamma_uu.yy * (input)->ptr[6] \
		+ gamma_uu.xx * gamma_uu.yz * (input)->ptr[11] \
		+ gamma_uu.xx * gamma_uu.yz * (input)->ptr[16] \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[7] \
		+ gamma_uu.xx * gamma_uu.zz * (input)->ptr[17] \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[8] \
		+ gamma_uu.xx * (input)->ptr[28] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[11] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[16] \
		+ 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[7] \
		- gamma_uu.xy * gamma_uu.yz * (input)->ptr[13] \
		+ gamma_uu.xy * gamma_uu.yz * (input)->ptr[18] \
		- gamma_uu.xy * gamma_uu.zz * (input)->ptr[14] \
		+ gamma_uu.xy * gamma_uu.zz * (input)->ptr[19] \
		- gamma_uu.xy * gamma_uu.xy * (input)->ptr[10] \
		+ gamma_uu.xy * (input)->ptr[29] \
		+ gamma_uu.xy * gamma_uu.xy * (input)->ptr[6] \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[13] \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[18] \
		+ gamma_uu.xz * gamma_uu.yz * (input)->ptr[14] \
		- gamma_uu.xz * gamma_uu.yz * (input)->ptr[19] \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[17] \
		+ gamma_uu.xz * (input)->ptr[30] \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[8]);\
	(results)->ptr[28] = gamma_uu.xy * (input)->ptr[22] \
		+ gamma_uu.xz * (input)->ptr[23] \
		+ gamma_uu.yy * (input)->ptr[24] \
		+ 2. * gamma_uu.yz * (input)->ptr[25] \
		+ gamma_uu.zz * (input)->ptr[26] \
		- (input)->ptr[27];\
	(results)->ptr[29] = \
		-(gamma_uu.xx * (input)->ptr[22] \
		+ gamma_uu.xy * (input)->ptr[24] \
		+ gamma_uu.xz * (input)->ptr[25]);\
	(results)->ptr[30] = \
		-(gamma_uu.xx * (input)->ptr[23] \
		+ gamma_uu.xy * (input)->ptr[25] \
		+ gamma_uu.xz * (input)->ptr[26]);\
}

//// MODULE_NAME: addSource
//// MODULE_DEPENDS: SETBOUNDS_NOGHOST initCond.codeprefix

kernel void addSource(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const deriv = derivBuf + index;

	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real const f = calc_f(U->alpha);	/* could be based on alpha... */

	real3 const S_l = real3_zero;
	sym3 const S_ll = sym3_zero;
	real const S = 0.;
	real const rho = 0.;

	/*  source terms */
	
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			/* K^i_j */
	real const trK = real3x3_trace(K_ul);								/* K^k_k */
	sym3 const KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		/* KSq_ij = K_ik K^k_j */

	/* d_llu = d_ij^k = d_ijl * gamma^lk */
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};

	/* d_ull = d^i_jk = gamma^il d_ljk */
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	/* e_l = d^j_ji */
	real3 const e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end	?>,
<? end
?>	};

	/* conn^k_ij = d_ij^k + d_ji^k - d^k_ij */
	_3sym3 const conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d_lll.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};


	/* d_l = d_ij^j */
	real3 const d_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_trace(d_llu[<?=i-1?>]),
<? end
?>	};
	
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 const Z_u = sym3_real3_mul(gamma_uu, U->Z_l);

	real3 const V_l = real3_sub(d_l, e_l);


	/* d_luu */
	_3sym3 const d_luu = (_3sym3){
<? for i,xi in ipairs(xNames) do		
?>		.<?=xi?> = sym3_real3x3_to_sym3_mul(gamma_uu, d_llu[<?=i-1?>]),
<? end
?>	};

	/* alpha_,t = shift terms - alpha^2 f gamma^ij K_ij */
	deriv->alpha += -U->alpha * U->alpha * f * trK;
	
	/* gamma_ij,t = shift terms - 2 alpha K_ij */
	sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));

	/* 2005 Bona et al A.1 */
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->K_ll.<?=xij?> += U->alpha * (
		.5 * (U->a_l.<?=xi?> * d_l.<?=xj?> + U->a_l.<?=xj?> * d_l.<?=xi?>)
		- .25 * (U->a_l.<?=xj?> * (2. * e_l.<?=xi?> - d_l.<?=xi?>) + U->a_l.<?=xi?> * (2. * e_l.<?=xj?> - d_l.<?=xj?>))
		- U->a_l.<?=xi?> * U->Z_l.<?=xj?> - U->a_l.<?=xj?> * U->Z_l.<?=xi?>
		+ (trK - 2. * U->Theta) * U->K_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- .5 * U->a_l.<?=xk?> * conn_ull.<?=xk?>.<?=xij?>
		+ .5 * U->a_l.<?=xk?> * d_ull.<?=xk?>.<?=xij?>
		- 2. * e_l.<?=xk?> * (d_llu[<?=i-1?>].<?=xj?>.<?=xk?> + d_llu[<?=j-1?>].<?=xi?>.<?=xk?>)
		+ (d_l.<?=xk?> + U->a_l.<?=xk?> - 2. * U->Z_l.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		- 2. * K_ul.<?=xk?>.<?=xi?> * U->K_ll.<?=sym(k,j)?>
<?		for l,xl in ipairs(xNames) do
?>		+ 2. * (d_llu[<?=i-1?>].<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,j)?> + d_llu[<?=j-1?>].<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,i)?>)
		- conn_ull.<?=xk?>.<?=sym(l,j)?> * conn_ull.<?=xl?>.<?=sym(k,i)?>
<?		end
	end
?> - 8. * M_PI * (S_ll.<?=xij?> - .5 * U->gamma_ll.<?=xij?> * (S - rho)));
<? end
?>

	/* 2005 Bona et al A.2 */
<? for i,xi in ipairs(xNames) do
?>	deriv->Z_l.<?=xi?> += U->alpha * (
		U->a_l.<?=xi?> * (trK - 2. * U->Theta)
<?	for k,xk in ipairs(xNames) do
?>		- U->a_l.<?=xk?> * K_ul.<?=xk?>.<?=xi?>
		+ K_ul.<?=xk?>.<?=xi?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?		for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * conn_ull.<?=xr?>.<?=sym(k,i)?>
<?		end
	end
?>
	) - 8 * M_PI * U->alpha * S_l.<?=xi?>;
<? end 
?>
	/* 2005 Bona et al A.3 */
	deriv->Theta += U->alpha * .5 * ( trK * (trK - 2. * U->Theta)
<? 
for k,xk in ipairs(xNames) do 
?>		+ 2. * U->a_l.<?=xk?> * (d_u.<?=xk?> - e_u.<?=xk?> - 2. * Z_u.<?=xk?>) 
		- d_u.<?=xk?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?	for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * K_ul.<?=xr?>.<?=xk?>
<?		for s,xs in ipairs(xNames) do
?>		+ d_luu.<?=xk?>.<?=sym(r,k)?> * conn_ull.<?=xk?>.<?=sym(r,s)?>
<?		end
	end
end?>
	) - 8. * M_PI * U->alpha * rho;

}

//// MODULE_NAME: constrainU

kernel void constrainU(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost);		
	global <?=cons_t?> * const U = UBuf + index;

	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real const f = calc_f(U->alpha);	/* could be based on alpha... */

	real const rho = 0.;


	sym3 const R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>
			+ conn_ull.<?=xk?>.<?=xij?> * (V_l.<?=xk?> - e_l.<?=xk?>)

<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(j,l)?>
			- 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=l-1?>].<?=xj?>.<?=xk?>
			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=j-1?>].<?=xl?>.<?=xk?>
			+ 2. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=k-1?>].<?=xj?>.<?=xl?>
			- 3. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=j-1?>].<?=xk?>.<?=xl?>
<? 		end
	end
?>		,
<? end
?>	};


	/* calculate the Hamiltonian and momentum constraints */
	/* scaled down by 1/8 to match B&S BSSNOK equations ... maybe I'll scale theirs up by 8 ... */
	/* B&S eqn 2.125 ... divded by two */
	/* Alcubierre eqn 2.5.9 */
	/* H = 1/2 (R + K^2 - K_ij K^ij) - 8 pi rho */
	real const R = sym3_dot(R_ll, gamma_uu);
	real const tr_KSq = sym3_dot(KSq_ll, gamma_uu);
	U->H = .5 * (R + trK * trK - tr_KSq) - 8. * M_PI * rho;
}
