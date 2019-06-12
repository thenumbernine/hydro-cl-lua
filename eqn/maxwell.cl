<?
local solver = eqn.solver
local common = require 'common'()
local xNames = common.xNames
local sym = common.sym
?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
cons_t fluxFromCons_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	<?=vec3?> E = calc_E(U);
	<?=vec3?> H = calc_H(U);
	return (cons_t){
<? if side == 0 then 
?>		.D = _<?=vec3?>(<?=zero?>, H.z, <?=neg?>(H.y)),
		.B = _<?=vec3?>(<?=zero?>, <?=neg?>(E.z), E.y),
<? elseif side == 1 then 
?>		.D = _<?=vec3?>(<?=neg?>(H.z), <?=zero?>, H.x),
		.B = _<?=vec3?>(E.z, <?=zero?>, <?=neg?>(E.x)),
<? elseif side == 2 then 
?>		.D = _<?=vec3?>(H.y, <?=neg?>(H.x), <?=zero?>),
		.B = _<?=vec3?>(<?=neg?>(E.y), E.x, <?=zero?>),
<? end 
?>		.divDPot = <?=zero?>,
		.divBPot = <?=zero?>,
		.rhoCharge = <?=zero?>,
		.sigma = <?=zero?>,
		.sqrt_1_eps = <?=susc_t?>_zero,
		.sqrt_1_mu = <?=susc_t?>_zero,
	};
}
<? end ?>

eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3 n
) {
	//this will fail with tensor susceptibility
	//but it doesn't belong here -- this is only the scalar case
	//for tensors, I should be eigen-decomposing the levi-civita times the tensor	
	return (eigen_t){
		.sqrt_1_eps = <?=susc_t?>_mul(
			<?=susc_t?>_add(UL.sqrt_1_eps, UR.sqrt_1_eps),
			<?=susc_t?>_from_real(.5)
		),
		.sqrt_1_mu = <?=susc_t?>_mul(
			<?=susc_t?>_add(UL.sqrt_1_mu, UR.sqrt_1_mu),
			<?=susc_t?>_from_real(.5)
		),
	};
}

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/
<? for side=0,solver.dim-1 do ?>

waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	waves_t Y;
	<?=scalar?> *Yp = (<?=scalar?>*)Y.ptr;
	<?=scalar?> *Xp = (<?=scalar?>*)X.ptr;

	const <?=scalar?> ise = <?=real_mul?>(eig.sqrt_1_eps, sqrt_1_2);
	const <?=scalar?> isu = <?=real_mul?>(eig.sqrt_1_mu, sqrt_1_2);

	<? if side==0 then ?>
	
	Yp[0] = <?=add?>(<?=mul?>(isu, Xp[4]), <?=mul?>(ise, Xp[2]));
	Yp[1] = <?=sub?>(<?=mul?>(isu, Xp[5]), <?=mul?>(ise, Xp[1]));
	Yp[2] = <?=sub?>(<?=mul?>(isu, Xp[3]), <?=mul?>(ise, Xp[0]));
	Yp[3] = <?=add?>(<?=mul?>(isu, Xp[3]), <?=mul?>(ise, Xp[0]));
	Yp[4] = <?=add?>(<?=mul?>(isu, Xp[5]), <?=mul?>(ise, Xp[1]));
	Yp[5] = <?=sub?>(<?=mul?>(isu, Xp[4]), <?=mul?>(ise, Xp[2]));
	
	<? elseif side==1 then ?>
	
	Yp[0] = <?=add?>(<?=mul?>(isu, Xp[5]), <?=mul?>(ise, Xp[0]));
	Yp[1] = <?=sub?>(<?=mul?>(isu, Xp[3]), <?=mul?>(ise, Xp[2]));
	Yp[2] = <?=sub?>(<?=mul?>(isu, Xp[4]), <?=mul?>(ise, Xp[1]));
	Yp[3] = <?=add?>(<?=mul?>(isu, Xp[4]), <?=mul?>(ise, Xp[1]));
	Yp[4] = <?=add?>(<?=mul?>(isu, Xp[3]), <?=mul?>(ise, Xp[2]));
	Yp[5] = <?=sub?>(<?=mul?>(isu, Xp[5]), <?=mul?>(ise, Xp[0]));
	
	<? elseif side==2 then ?>
	
	Yp[0] = <?=add?>(<?=mul?>(isu, Xp[3]), <?=mul?>(ise, Xp[1]));
	Yp[1] = <?=sub?>(<?=mul?>(isu, Xp[4]), <?=mul?>(ise, Xp[0]));
	Yp[2] = <?=sub?>(<?=mul?>(isu, Xp[5]), <?=mul?>(ise, Xp[2]));
	Yp[3] = <?=add?>(<?=mul?>(isu, Xp[5]), <?=mul?>(ise, Xp[2]));
	Yp[4] = <?=add?>(<?=mul?>(isu, Xp[4]), <?=mul?>(ise, Xp[0]));
	Yp[5] = <?=sub?>(<?=mul?>(isu, Xp[3]), <?=mul?>(ise, Xp[1]));
	
	<? end ?>

	return Y;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x
) {
	cons_t Y;
	<?=scalar?> *Yp = (<?=scalar?>*)Y.ptr;
	<?=scalar?> *Xp = (<?=scalar?>*)X.ptr;

	const <?=scalar?> se = <?=real_mul?>(<?=inv?>(eig.sqrt_1_eps), sqrt_1_2);
	const <?=scalar?> su = <?=real_mul?>(<?=inv?>(eig.sqrt_1_mu), sqrt_1_2);

	<? if side==0 then ?>
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
	Yp[0] = <?=mul?>(se, <?=sub?>(Xp[3], Xp[2]));
	Yp[1] = <?=mul?>(se, <?=sub?>(Xp[4], Xp[1]));
	Yp[2] = <?=mul?>(se, <?=sub?>(Xp[0], Xp[5]));
	Yp[3] = <?=mul?>(su, <?=add?>(Xp[2], Xp[3]));
	Yp[4] = <?=mul?>(su, <?=add?>(Xp[0], Xp[5]));
	Yp[5] = <?=mul?>(su, <?=add?>(Xp[1], Xp[4]));
	
	<? elseif side==1 then ?>

/*
x, -z, -y, y, z, -x
z,  x,  y, y, x,  z

1  0  0 0 0 -1
0  0 -1 1 0  0
0 -1  0 0 1  0
0  1  0 0 1  0
0  0  1 1 0  0
1  0  0 0 0  1
*/
	Yp[0] = <?=mul?>(se, <?=sub?>(Xp[0], Xp[5]));
	Yp[1] = <?=mul?>(se, <?=sub?>(Xp[3], Xp[2]));
	Yp[2] = <?=mul?>(se, <?=sub?>(Xp[4], Xp[1]));
	Yp[3] = <?=mul?>(su, <?=add?>(Xp[1], Xp[4]));
	Yp[4] = <?=mul?>(su, <?=add?>(Xp[2], Xp[3]));
	Yp[5] = <?=mul?>(su, <?=add?>(Xp[0], Xp[5]));
	
	<? elseif side==2 then ?>

/*
y, -x, -z, z, x, -y
x,  y,  z, z,  y,  x

0 -1  0 0 1  0
1  0  0 0 0 -1
0  0 -1 1 0  0
1  0  0 0 0  1
0  1  0 0 1  0
0  0  1 1 0  0
*/
	Yp[0] = <?=mul?>(se, <?=sub?>(Xp[4], Xp[1]));
	Yp[1] = <?=mul?>(se, <?=sub?>(Xp[0], Xp[5]));
	Yp[2] = <?=mul?>(se, <?=sub?>(Xp[3], Xp[2]));
	Yp[3] = <?=mul?>(su, <?=add?>(Xp[0], Xp[5]));
	Yp[4] = <?=mul?>(su, <?=add?>(Xp[1], Xp[4]));
	Yp[5] = <?=mul?>(su, <?=add?>(Xp[2], Xp[3]));
	
	<? end ?>
	
	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	cons_t Y;
	<?=scalar?> *Yp = (<?=scalar?>*)Y.ptr;
	<?=scalar?> *Xp = (<?=scalar?>*)X.ptr;

	<?=vec3?> D = X.D;
	<?=vec3?> B = X.B;
	<?=scalar?> _1_eps = <?=mul?>(eig.sqrt_1_eps, eig.sqrt_1_eps);
	<?=scalar?> _1_mu = <?=mul?>(eig.sqrt_1_mu, eig.sqrt_1_mu);

	<? if side==0 then ?>
	
	Yp[0] = <?=zero?>;
	Yp[1] = <?=mul?>(B.z, _1_mu);
	Yp[2] = <?=neg?>(<?=mul?>(B.y, _1_mu));
	Yp[3] = <?=zero?>;
	Yp[4] = <?=neg?>(<?=mul?>(D.z, _1_eps));
	Yp[5] = <?=mul?>(D.y, _1_eps);

	<? elseif side==1 then ?>
		
	Yp[0] = <?=neg?>(<?=mul?>(B.z, _1_mu));
	Yp[1] = <?=zero?>;
	Yp[2] = <?=mul?>(B.x, _1_mu);
	Yp[3] = <?=mul?>(D.z, _1_eps);
	Yp[4] = <?=zero?>;
	Yp[5] = <?=neg?>(<?=mul?>(D.x, _1_eps));
		
	<? elseif side==2 then ?>
		
	Yp[0] = <?=mul?>(B.y, _1_mu);
	Yp[1] = <?=neg?>(<?=mul?>(B.x, _1_mu));
	Yp[2] = <?=zero?>;
	Yp[3] = <?=neg?>(<?=mul?>(D.y, _1_eps));
	Yp[4] = <?=mul?>(D.x, _1_eps);
	Yp[5] = <?=zero?>;
		
	<? end ?>
	
	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}
<? end ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	//TODO J = J_f + J_b = J_f + J_P + J_M = J_f + dP/dt + curl M
<? for i,xi in ipairs(xNames) do 
?>	deriv->D.<?=xi?> = <?=sub?>(
		deriv->D.<?=xi?>,
			<?=mul?>(U->D.<?=xi?>, <?=mul?>(
				<?=susc_t?>_mul(U->sqrt_1_eps, U->sqrt_1_eps),
				U->sigma
			)
		)
	);
<? end
?>
	//for non-time-varying susceptibilities, here's the source term:
	//D_i,t ... = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_,j B_k - J_i
	//B_i,t ... = 1/sqrt(g) g_il epsBar^ljk (1/eps)_,j B_k

	<?=vec3?> grad_1_mu = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(
				U[solver->stepsize.<?=xj?>].sqrt_1_mu,
				U[solver->stepsize.<?=xj?>].sqrt_1_mu),
			<?=mul?>(
				U[-solver->stepsize.<?=xj?>].sqrt_1_mu,
				U[-solver->stepsize.<?=xj?>].sqrt_1_mu)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>
	
	<?=vec3?> grad_1_eps = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1]
?>	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(
				U[solver->stepsize.<?=xj?>].sqrt_1_eps,
				U[solver->stepsize.<?=xj?>].sqrt_1_eps),
			<?=mul?>(
				U[-solver->stepsize.<?=xj?>].sqrt_1_eps,
				U[-solver->stepsize.<?=xj?>].sqrt_1_eps)
		), 1. / solver->grid_dx.s<?=j?>);
<? end
?>
	real _1_sqrt_det_g = 1. / coord_volume(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] 
	?>{
		cons_t flux = fluxFromCons_<?=j?>(solver, *U, x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = <?=sub?>(deriv->D.<?=xj?>, <?=vec3?>_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = <?=sub?>(deriv->B.<?=xj?>, <?=vec3?>_dot(flux.B, grad_1_eps));
	}<? end ?>
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
eigen_t eigen_forCell_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	eigen_t eig;
	eig.sqrt_1_eps = U.sqrt_1_eps;
	eig.sqrt_1_mu = U.sqrt_1_mu;
	return eig;
}
<? end ?>
