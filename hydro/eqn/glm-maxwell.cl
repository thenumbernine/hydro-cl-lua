typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=solver.solver_t?> solver_t;

<? if moduleName == nil then ?>
<? elseif moduleName == "sqrt_1_2" then ?>

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? elseif moduleName == "eqn.common" then ?>

<? if scalar == 'real' then ?>

#define eqn_coordLenSq coordLenSq
#define eqn_cartesianToCoord cartesianToCoord
#define eqn_coord_lower coord_lower

<? elseif scalar == 'cplx' then ?>

real eqn_coordLenSq(cplx3 v, real3 x) {
	return coordLenSq(cplx3_re(v), x)
		+ coordLenSq(cplx3_im(v), x);
}

cplx3 eqn_cartesianToCoord(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		cartesianToCoord(cplx3_re(v), x),
		cartesianToCoord(cplx3_im(v), x));
}

cplx3 eqn_coord_lower(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		coord_lower(cplx3_re(v), x),
		coord_lower(cplx3_im(v), x));
}

<? end -- scalar ?>

<?=vec3?> calc_E(<?=eqn.cons_t?> U) { 
	return <?=vec3?>_<?=susc_t?>_mul(U.D, <?=susc_t?>_mul(U.sqrt_1_eps, U.sqrt_1_eps));
}
<?=vec3?> calc_H(<?=eqn.cons_t?> U) { 
	return <?=vec3?>_<?=susc_t?>_mul(U.B, <?=susc_t?>_mul(U.sqrt_1_mu, U.sqrt_1_mu));
}

<? elseif moduleName == "applyInitCond" then ?>

<? 
local cons_t = eqn.cons_t
local susc_t = eqn.susc_t
local scalar = eqn.scalar
local vec3 = eqn.vec3
local zero = scalar..'_zero'
?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	global <?=cons_t?>* U = UBuf + index;

	//used
	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=scalar?>_from_real(1.);
	<?=susc_t?> permittivity = <?=susc_t?>_from_real(1.);
	<?=susc_t?> permeability = <?=susc_t?>_from_real(1.);
	
	//throw-away
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
	
	<?=initCode()?>
	
	U->D = eqn_cartesianToCoord(D, x);
	U->B = eqn_cartesianToCoord(B, x);
	U->phi = <?=zero?>;
	U->psi = <?=zero?>;
	U->sigma = conductivity;
	U->rhoCharge = <?=zero?>;
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permittivity));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permeability));
}


<? elseif moduleName == "fluxFromCons" then ?>

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=vec3?> E = calc_E(U);
	<?=vec3?> H = calc_H(U);
	
	<?=eqn.cons_t?> F;
	if (n.side == 0) {
		F.D = _<?=vec3?>(<?=real_mul?>(U.phi, solver->divPhiWavespeed / unit_m_per_s),  H.z, <?=neg?>(H.y));
		F.B = _<?=vec3?>(<?=real_mul?>(U.psi, solver->divPsiWavespeed / unit_m_per_s), <?=neg?>(E.z),  E.y);
	} else if (n.side == 1) {
		F.D = _<?=vec3?>(<?=neg?>(H.z), <?=real_mul?>(U.phi, solver->divPhiWavespeed / unit_m_per_s),  H.x);
		F.B = _<?=vec3?>( E.z, <?=real_mul?>(U.psi, solver->divPsiWavespeed / unit_m_per_s), <?=neg?>(E.x));
	} else if (n.side == 2) {
		F.D = _<?=vec3?>( H.y, <?=neg?>(H.x), <?=real_mul?>(U.phi, solver->divPhiWavespeed / unit_m_per_s));
		F.B = _<?=vec3?>(<?=neg?>(E.y),  E.x, <?=real_mul?>(U.psi, solver->divPsiWavespeed / unit_m_per_s));
	}
	real D_n = normal_vecDotN1(n, U.D);
	real B_n = normal_vecDotN1(n, U.B);
	F.phi = <?=real_mul?>(D_n, solver->divPhiWavespeed / unit_m_per_s);
	F.psi = <?=real_mul?>(B_n, solver->divPsiWavespeed / unit_m_per_s);
	F.sigma = <?=zero?>;
	F.rhoCharge = <?=zero?>;
	F.sqrt_1_eps = <?=susc_t?>_zero;
	F.sqrt_1_mu = <?=susc_t?>_zero;
	return F;
}


<? elseif moduleName == "eigen_forInterface" then ?>

typedef <?=eqn.eigen_t?> eigen_t;

eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	normal_t n
) {
	//this will fail with tensor susceptibility
	//but it doesn't belong here -- this is only the scalar case
	//for tensors, I should be eigen-decomposing the levi-civita times the tensor	
	return (eigen_t){
		.sqrt_1_eps = <?=susc_t?>_real_mul(<?=susc_t?>_add(UL.sqrt_1_eps, UR.sqrt_1_eps), .5),
		.sqrt_1_mu = <?=susc_t?>_real_mul(<?=susc_t?>_add(UL.sqrt_1_mu, UR.sqrt_1_mu), .5),
	};
}

<? elseif moduleName == "eigen_forCell" then ?>

typedef <?=eqn.eigen_t?> eigen_t;

//used by PLM
eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	normal_t n
) {
	return (eigen_t){
		.sqrt_1_eps = U.sqrt_1_eps,
		.sqrt_1_mu = U.sqrt_1_mu,
	};
}

<? elseif moduleName == "eigen_left/rightTransform" then ?>

typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/
waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x,
	normal_t n
) {
	waves_t Y;
	<?=scalar?>* Xp = (<?=scalar?>*)X.ptr;
	<?=scalar?>* Yp = (<?=scalar?>*)Y.ptr;

	<?=susc_t?> sqrt_1_eps = eig.sqrt_1_eps;					//(m^3 kg)^.5/(C s)
	<?=susc_t?> sqrt_eps = <?=susc_t?>_inv(sqrt_1_eps);			//(C s)/(m^3 kg)^.5
	<?=susc_t?> sqrt_1_mu = eig.sqrt_1_mu;						//C/(kg m)^.5
	<?=susc_t?> sqrt_mu = <?=susc_t?>_inv(sqrt_1_mu);			//(kg m)^.5/C
	<?=susc_t?> v_p = <?=susc_t?>_mul(sqrt_1_eps, sqrt_1_mu);	//m/s

	if (n.side == 0) {

		Yp[0] = (Xp[6] - Xp[0]) * sqrt_eps * sqrt_1_2;					//(C^2 s)/(m^3.5 kg^.5)
		Yp[1] = (Xp[7] - Xp[3]) * sqrt_eps * sqrt_1_2;					//kg^.5/m^1.5
		Yp[2] = (Xp[4] * sqrt_eps + Xp[2] * sqrt_mu) * v_p * sqrt_1_2;	//kg^.5/(m^.5 s) 
		Yp[3] = (Xp[5] * sqrt_eps - Xp[1] * sqrt_mu) * v_p * sqrt_1_2;	//kg^.5/(m^.5 s) 
		Yp[4] = (Xp[4] * sqrt_eps - Xp[2] * sqrt_mu) * v_p * sqrt_1_2;	//kg^.5/(m^.5 s) 
		Yp[5] = (Xp[5] * sqrt_eps + Xp[1] * sqrt_mu) * v_p * sqrt_1_2;	//kg^.5/(m^.5 s) 
		Yp[6] = (Xp[0] + Xp[6]) * sqrt_eps * sqrt_1_2;					//(C^2 s)/(m^3.5 kg^.5)
		Yp[7] = (Xp[3] + Xp[7]) * sqrt_eps * sqrt_1_2;					//kg^.5/m^1.5

	} else if (n.side == 1) {

		//same units as x dir
		Yp[0] = (Xp[6] - Xp[1]) * sqrt_eps * sqrt_1_2;
		Yp[1] = (Xp[7] - Xp[4]) * sqrt_eps * sqrt_1_2;
		Yp[2] = (Xp[3] * sqrt_eps - Xp[2] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[3] = (Xp[5] * sqrt_eps + Xp[0] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[4] = (Xp[3] * sqrt_eps + Xp[2] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[5] = (Xp[5] * sqrt_eps - Xp[0] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[6] = (Xp[1] + Xp[6]) * sqrt_eps * sqrt_1_2;
		Yp[7] = (Xp[4] + Xp[7]) * sqrt_eps * sqrt_1_2;

	} else if (n.side == 2) {
	
		//same units as x dir
		Yp[0] = (Xp[6] - Xp[2]) * sqrt_eps * sqrt_1_2;
		Yp[1] = (Xp[7] - Xp[5]) * sqrt_eps * sqrt_1_2;
		Yp[2] = (Xp[3] * sqrt_eps + Xp[1] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[3] = (Xp[4] * sqrt_eps - Xp[0] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[4] = (Xp[3] * sqrt_eps - Xp[1] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[5] = (Xp[4] * sqrt_eps + Xp[0] * sqrt_mu) * v_p * sqrt_1_2;
		Yp[6] = (Xp[2] + Xp[6]) * sqrt_eps * sqrt_1_2;
		Yp[7] = (Xp[5] + Xp[7]) * sqrt_eps * sqrt_1_2;

	}

	return Y;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x,
	normal_t n
) {
	cons_t Y;
	<?=scalar?>* Xp = (<?=scalar?>*)X.ptr;
	<?=scalar?>* Yp = (<?=scalar?>*)Y.ptr;

	<?=scalar?> sqrt_1_eps = eig.sqrt_1_eps;
	<?=scalar?> sqrt_1_mu = eig.sqrt_1_mu;
	<?=scalar?> sqrt_eps = <?=inv?>(eig.sqrt_1_eps);
	<?=scalar?> sqrt_mu = <?=inv?>(eig.sqrt_1_mu);
	<?=scalar?> sqrt_2 = <?=inv?>(sqrt_1_2);

	if (n.side == 0) {

		Yp[0] = ((-(Xp[0] - Xp[6])) / (sqrt_2 * sqrt_eps));
		Yp[1] = ((-(sqrt_eps * (Xp[3] - Xp[5]))) / sqrt_2);
		Yp[2] = ((sqrt_eps * (Xp[2] - Xp[4])) / sqrt_2);
		Yp[3] = ((-(Xp[1] - Xp[7])) / (sqrt_2 * sqrt_eps));
		Yp[4] = ((sqrt_mu * (Xp[2] + Xp[4])) / sqrt_2);
		Yp[5] = ((sqrt_mu * (Xp[3] + Xp[5])) / sqrt_2);
		Yp[6] = ((Xp[0] + Xp[6]) / (sqrt_2 * sqrt_eps));
		Yp[7] = ((Xp[1] + Xp[7]) / (sqrt_2 * sqrt_eps));
	
	} else if (n.side == 1) {

		Yp[0] = ((sqrt_eps * (Xp[3] - Xp[5])) / sqrt_2);
		Yp[1] = ((-(Xp[0] - Xp[6])) / (sqrt_2 * sqrt_eps));
		Yp[2] = ((-(sqrt_eps * (Xp[2] - Xp[4]))) / sqrt_2);
		Yp[3] = ((sqrt_mu * (Xp[2] + Xp[4])) / sqrt_2);
		Yp[4] = ((-(Xp[1] - Xp[7])) / (sqrt_2 * sqrt_eps));
		Yp[5] = ((sqrt_mu * (Xp[3] + Xp[5])) / sqrt_2);
		Yp[6] = ((Xp[0] + Xp[6]) / (sqrt_2 * sqrt_eps));
		Yp[7] = ((Xp[1] + Xp[7]) / (sqrt_2 * sqrt_eps));
	
	} else if (n.side == 2) {

		Yp[0] = ((-(sqrt_eps * (Xp[3] - Xp[5]))) / sqrt_2);
		Yp[1] = ((sqrt_eps * (Xp[2] - Xp[4])) / sqrt_2);
		Yp[2] = ((-(Xp[0] - Xp[6])) / (sqrt_2 * sqrt_eps));
		Yp[3] = ((sqrt_mu * (Xp[2] + Xp[4])) / sqrt_2);
		Yp[4] = ((sqrt_mu * (Xp[3] + Xp[5])) / sqrt_2);
		Yp[5] = ((-(Xp[1] - Xp[7])) / (sqrt_2 * sqrt_eps));
		Yp[6] = ((Xp[0] + Xp[6]) / (sqrt_2 * sqrt_eps));
		Yp[7] = ((Xp[1] + Xp[7]) / (sqrt_2 * sqrt_eps));

	}

	for (int i = <?=eqn.numWaves?>; i < <?=eqn.numStates?>; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

<? elseif moduleName == "eigen_fluxTransform" then ?>

#define eigen_fluxTransform(solver, eig, X, x, n) fluxFromCons(solver, X, x, n)

<? elseif moduleName == "addSource" then ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cellBuf[index].pos;
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	//TODO J = J_f + J_b = J_f + J_P + J_M = J_f + dP/dt + curl M
	deriv->D = <?=vec3?>_sub(
		deriv->D, 
		<?=vec3?>_<?=scalar?>_mul(
			U->D, 
			<?=mul?>(<?=mul?>(U->sqrt_1_eps, U->sqrt_1_eps), U->sigma)
		)
	);
	
	//for non-time-varying susceptibilities, here's the source term:
	//D_i,t ... = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_,j B_k - J_i
	//B_i,t ... = 1/sqrt(g) g_il epsBar^ljk (1/eps)_,j B_k
	
	<?=vec3?> grad_1_mu = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[solver->stepsize.<?=xj?>].sqrt_1_mu, U[solver->stepsize.<?=xj?>].sqrt_1_mu),
			<?=mul?>(U[-solver->stepsize.<?=xj?>].sqrt_1_mu, U[-solver->stepsize.<?=xj?>].sqrt_1_mu)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>
	
	<?=vec3?> grad_1_eps = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[solver->stepsize.<?=xj?>].sqrt_1_eps, U[solver->stepsize.<?=xj?>].sqrt_1_eps),
			<?=mul?>(U[-solver->stepsize.<?=xj?>].sqrt_1_eps, U[-solver->stepsize.<?=xj?>].sqrt_1_eps)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>
	
	real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>{
		cons_t flux = fluxFromCons(solver, *U, x, normal_forSide<?=j?>(x));
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = <?=sub?>(deriv->D.<?=xj?>, <?=vec3?>_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = <?=sub?>(deriv->B.<?=xj?>, <?=vec3?>_dot(flux.B, grad_1_eps));
	}<? end ?>
	
	deriv->phi = <?=add?>(deriv->phi, <?=real_mul?>(U->rhoCharge, solver->divPhiWavespeed / unit_m_per_s));
}

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
