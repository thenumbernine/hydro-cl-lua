typedef <?=eqn.cons_t?> cons_t;
typedef <?=solver.solver_t?> solver_t;

<? if moduleName == nil then ?>
<? elseif moduleName == "calc_gamma_ll" then ?>

#define calc_gamma_ll(U, x)	((U)->gamma_ll)

<? elseif moduleName == "calc_gamma_uu" then 
depmod{
	"cons_t",
}
?>

sym3 calc_gamma_uu(const global <?=eqn.cons_t?>* U, real3 x) {
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	return gamma_uu;
}

<? elseif moduleName == "setFlatSpace" then 
depmod{
	"solver_t",
	"cons_t",
}
?>

void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->gamma_ll = sym3_ident;
	U->a_l = real3_zero;
	U->d_lll = _3sym3_zero;
	U->K_ll = sym3_zero;
	U->V_l = real3_zero;
<? if eqn.useShift ~= 'none' then 
?>	U->beta_u = real3_zero;
<? end 
?>	
<? if eqn.useStressEnergyTerms then ?>
	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

<? elseif moduleName == "applyInitCond" then 
depmod{
	"coordMap",
	"coord_g_ll",
	"rescaleFromCoord/rescaleToCoord",
}
?>

<?
-- eqn.einstein compatability hack ...
-- if initCond.initAnalytical is *NOT* set then it will call initDerivs(), which is usually only done for finite-differencing the original state variables.
-- if initCond.initAnalytical is set then it expects the analytical expressions to be written *HERE*
-- ... which adm3d does not do yet.
if eqn.initCond.initAnalytical then
	error("TODO - adm3d can't handle analytical initial conditions yet")
end
?>
<? if eqn.initCond.useBSSNVars then ?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real W = 1.;
	real K = 0.;
	real3 LambdaBar_U = real3_zero;
	real3 beta_U = real3_zero;
	real3 B_U = real3_zero;
	sym3 epsilon_LL = sym3_zero;
	sym3 ABar_LL = sym3_zero;

	//TODO more stress-energy vars 
	real rho = 0.;

	<?=initCode()?>

	U->alpha = alpha;

	// gammaHat_IJ = delta_IJ
	// gamma_ij = e_i^I e_j^J (epsilon_IJ + gammaHat_IJ) / W^2
	sym3 gammaBar_LL = sym3_add(epsilon_LL, sym3_ident);
	sym3 gamma_LL = sym3_real_mul(gammaBar_LL, 1. / (W*W));
	U->gamma_ll = sym3_rescaleToCoord_LL(gamma_LL, x);
	
	// K_ij = e_i^I e_j^J (ABar_IJ + gammaBar_IJ K/3) / W^2
	U->K_ll = sym3_rescaleToCoord_LL(
		sym3_add(
			sym3_real_mul(ABar_LL, 1. / (W*W)),
			sym3_real_mul(gamma_LL, K / 3.)
		), x);

	// TODO maybe derive this from LambdaBar_U ?
	U->V_l = real3_zero;

<? if eqn.useShift ~= 'none' then
?>	U->beta_u = real3_rescaleFromCoord_U(beta_U);
<? end -- TODO support for hyperbolic gamma driver, so we can read B_U
?>

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	global <?=solver.coord.cell_t?> const *cellBuf 
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (solver->grid_dx.s<?=i-1?> * U->alpha);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> - U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>) / solver->grid_dx.s<?=i-1?>;
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>

//V_i = d_ik^k - d^k_ki 
<? for i,xi in ipairs(xNames) do ?>
	U->V_l.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + gamma_uu.<?=sym(j,k)?> * ( U->d_lll.<?=xi?>.<?=sym(j,k)?> - U->d_lll.<?=xj?>.<?=sym(k,i)?> )<?
		end
	end ?>;
<? end ?>
}

<? else	-- not eqn.initCond.useBSSNVars ?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = coord_g_ll(x);
	sym3 K_ll = sym3_zero;

	//TODO more stress-energy vars 
	real rho = 0.;

	<?=initCode()?>

	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;
	U->V_l = real3_zero;
<? if eqn.useShift ~= 'none' then
?>	U->beta_u = beta_u;
<? end
?>
<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	global <?=solver.coord.cell_t?> const *cellBuf 
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (solver->grid_dx.s<?=i-1?> * U->alpha);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> - U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>) / solver->grid_dx.s<?=i-1?>;
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>

//V_i = d_ik^k - d^k_ki 
<? for i,xi in ipairs(xNames) do ?>
	U->V_l.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + gamma_uu.<?=sym(j,k)?> * ( U->d_lll.<?=xi?>.<?=sym(j,k)?> - U->d_lll.<?=xj?>.<?=sym(k,i)?> )<?
		end
	end ?>;
<? end ?>
}

<? end	-- eqn.initCond.useBSSNVars ?>



<? elseif moduleName == "fluxFromCons" then 
depmod{
	"rotate",
	"cons_t",
	"solver_t",
	"normal_t",
	"rotate",	-- real3_swap*
	"initCond.codeprefix",		-- calc_f
}
?>

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	real f_alpha = calc_f_alpha(U.alpha);
	
	real det_gamma = sym3_det(U.gamma_ll);
	sym3 gamma_uu = sym3_inv(U.gamma_ll, det_gamma);
	
	real alpha = U.alpha;
	
	real3 V_l = U.V_l;
	real3 a_l = U.a_l;
	sym3 gamma_ll = U.gamma_ll;
	sym3 K_ll = U.K_ll;
	_3sym3 d_lll = U.d_lll;
	
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		V_l = real3_swap<?=side?>(V_l);
		a_l = real3_swap<?=side?>(a_l);
		gamma_ll = sym3_swap<?=side?>(gamma_ll);
		K_ll = sym3_swap<?=side?>(K_ll);
		d_lll = _3sym3_swap<?=side?>(d_lll);
		gamma_uu = sym3_swap<?=side?>(gamma_uu);
	}
	<? end ?>

	<?=eqn.cons_t?> F = {.ptr={0}};

	// BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/adm_noZeroRows.html
	//(except me replacing alpha * f with f_alpha)
	F.a_l.x = f_alpha * (2. * K_ll.xy * gamma_uu.xy + 2. * K_ll.xz * gamma_uu.xz + 2. * K_ll.yz * gamma_uu.yz + K_ll.xx * gamma_uu.xx + K_ll.yy * gamma_uu.yy + K_ll.zz * gamma_uu.zz);
	F.d_lll.x.xx = K_ll.xx * alpha;
	F.d_lll.x.xy = K_ll.xy * alpha;
	F.d_lll.x.xz = K_ll.xz * alpha;
	F.d_lll.x.yy = K_ll.yy * alpha;
	F.d_lll.x.yz = K_ll.yz * alpha;
	F.d_lll.x.zz = K_ll.zz * alpha;
	F.K_ll.xx = alpha * (a_l.x + d_lll.x.yy * gamma_uu.yy + 2. * d_lll.x.yz * gamma_uu.yz + d_lll.x.zz * gamma_uu.zz - d_lll.y.xx * gamma_uu.xy - 2. * d_lll.y.xy * gamma_uu.yy - 2. * d_lll.y.xz * gamma_uu.yz - d_lll.z.xx * gamma_uu.xz - 2. * d_lll.z.xy * gamma_uu.yz - 2. * d_lll.z.xz * gamma_uu.zz);
	F.K_ll.xy = (alpha * (a_l.y - 2. * d_lll.x.yy * gamma_uu.xy - 2. * d_lll.x.yz * gamma_uu.xz - 2. * d_lll.y.yy * gamma_uu.yy - 2. * d_lll.y.yz * gamma_uu.yz - 2. * d_lll.z.yy * gamma_uu.yz - 2. * d_lll.z.yz * gamma_uu.zz)) / 2.;
	F.K_ll.xz = (alpha * (a_l.z - 2. * d_lll.x.yz * gamma_uu.xy - 2. * d_lll.x.zz * gamma_uu.xz - 2. * d_lll.y.yz * gamma_uu.yy - 2. * d_lll.y.zz * gamma_uu.yz - 2. * d_lll.z.yz * gamma_uu.yz - 2. * d_lll.z.zz * gamma_uu.zz)) / 2.;
	F.K_ll.yy = alpha * (d_lll.x.yy * gamma_uu.xx + d_lll.y.yy * gamma_uu.xy + d_lll.z.yy * gamma_uu.xz);
	F.K_ll.yz = alpha * (d_lll.x.yz * gamma_uu.xx + d_lll.y.yz * gamma_uu.xy + d_lll.z.yz * gamma_uu.xz);
	F.K_ll.zz = alpha * (d_lll.x.zz * gamma_uu.xx + d_lll.y.zz * gamma_uu.xy + d_lll.z.zz * gamma_uu.xz);	
	// END CUT

	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		F.V_l = real3_swap<?=side?>(F.V_l);
		F.a_l = real3_swap<?=side?>(F.a_l);
		F.gamma_ll = sym3_swap<?=side?>(F.gamma_ll);
		F.K_ll = sym3_swap<?=side?>(F.K_ll);
		F.d_lll = _3sym3_swap<?=side?>(F.d_lll);
	}
	<? end ?>

	return F;
}

<? elseif moduleName == "calcDT" then 
depmod{
	"solver_t",
	"SETBOUNDS",
	"eqn.guiVars.compileTime",
	"initCond.codeprefix",		-- calc_f
}
?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma_ij) is only called once
	real f_alphaSq = calc_f_alphaSq(U->alpha);
	real det_gamma = sym3_det(U->gamma_ll);
	real alpha_sqrt_f = sqrt(f_alphaSq);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side==0 then ?>
		real gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
		<? elseif side==1 then ?>
		real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
		<? elseif side==2 then ?>
		real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
		<? end ?>	
		real sqrt_gammaUjj = sqrt(gammaUjj);
		real lambdaLight = sqrt_gammaUjj * U->alpha;
		real lambdaGauge = sqrt_gammaUjj * alpha_sqrt_f;
		
		real lambda = (real)max(lambdaGauge, lambdaLight);

		<? if eqn.useShift ~= 'none' then ?>
		real betaUi = U->beta_u.s<?=side?>;
		<? else ?>
		const real betaUi = 0.;
		<? end ?>
		
		real lambdaMin = (real)min((real)0., -betaUi - lambda);
		real lambdaMax = (real)max((real)0., -betaUi + lambda);
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt; 
}

<? elseif moduleName == "eigen_forCell" then 
depmod{
	"initCond.codeprefix",		-- calc_f
}
?>

typedef <?=eqn.eigen_t?> eigen_t;

//used by PLM
eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	
	//This is interesting, because normal_t varies based on our vector components.
	//However in my GR solvers the components are irrespective of the grid -- instead they are based on the metric of the state variables.
	normal_t n
) {
	eigen_t eig;
	eig.alpha = U.alpha;
	eig.alpha_sqrt_f = sqrt(calc_f_alphaSq(U.alpha));
	real det_gamma = sym3_det(U.gamma_ll);
	eig.gamma_uu = sym3_inv(U.gamma_ll, det_gamma);
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gamma_uu.xx), sqrt(eig.gamma_uu.yy), sqrt(eig.gamma_uu.zz));
	
	<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = U.beta_u;
	<? end ?>

	return eig;
}

<? elseif moduleName == "calcCellMinMaxEigenvalues" then 
depmod{
	"initCond.codeprefix",		-- calc_f
}
?>

range_t calcCellMinMaxEigenvalues(
	const global cons_t* U,
	real3 x,
	normal_t n
) {
	real det_gamma = sym3_det(U->gamma_ll);

	real gammaUjj;
	if (n.side == 0) {
		gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
	} else if (n.side == 1) {
		gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
	} else if (n.side == 2) {
		gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
	}
	
	real lambdaLight = U->alpha * sqrt(gammaUjj);
	
	real f = calc_f(U->alpha);
	real lambdaGauge = lambdaLight * sqrt(f);

	real lambdaMax = max(lambdaGauge, lambdaLight);
	//= lambdaLight * max(sqrt(f), 1)
	real lambdaMin = -lambdaMin;

	<? if eqn.useShift ~= 'none' then ?>
	lambdaMin -= normal_vecDotN1(n, U->beta_u);
	lambdaMax -= normal_vecDotN1(n, U->beta_u);
	<? end ?>

	return (range_t){
		.min = lambdaMin, 
		.max = lambdaMax,
	};
}

<? elseif moduleName == "eigen_forInterface" then 
depmod{
	"initCond.codeprefix",		-- calc_f
}
?>

typedef <?=eqn.eigen_t?> eigen_t;

//used for interface eigen basis
eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	normal_t n
) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	sym3 avg_gamma = (sym3){
		.xx = .5 * (UL.gamma_ll.xx + UR.gamma_ll.xx),
		.xy = .5 * (UL.gamma_ll.xy + UR.gamma_ll.xy),
		.xz = .5 * (UL.gamma_ll.xz + UR.gamma_ll.xz),
		.yy = .5 * (UL.gamma_ll.yy + UR.gamma_ll.yy),
		.yz = .5 * (UL.gamma_ll.yz + UR.gamma_ll.yz),
		.zz = .5 * (UL.gamma_ll.zz + UR.gamma_ll.zz),
	};
	real det_avg_gamma = sym3_det(avg_gamma);

	eigen_t eig;
	eig.alpha = alpha;
	eig.alpha_sqrt_f = sqrt(calc_f_alphaSq(alpha));
	eig.gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);
	eig.sqrt_gammaUjj.x = sqrt(eig.gamma_uu.xx);
	eig.sqrt_gammaUjj.y = sqrt(eig.gamma_uu.yy);
	eig.sqrt_gammaUjj.z = sqrt(eig.gamma_uu.zz);
	
	<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = real3_real_mul(real3_add(UL.beta_u, UR.beta_u), .5);
	<? end ?>

	return eig;
}

<? elseif moduleName == "eigen_left/rightTransform" then ?>

typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;

waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x,
	normal_t n
) {
	waves_t results;
<? if not eqn.noZeroRowsInFlux then ?>

	real f = calc_f(eig.alpha);

	if (n.side == 0) {

		//a_y, a_z
		results.ptr[6] = inputU.a_l.y;
		results.ptr[7] = inputU.a_l.z;

		//d_yij
		results.ptr[8] = inputU.d_lll.y.xx;
		results.ptr[9] = inputU.d_lll.y.xy;
		results.ptr[10] = inputU.d_lll.y.xz;
		results.ptr[11] = inputU.d_lll.y.yy;
		results.ptr[12] = inputU.d_lll.y.yz;
		results.ptr[13] = inputU.d_lll.y.zz;

		//d_zij
		results.ptr[14] = inputU.d_lll.z.xx;
		results.ptr[15] = inputU.d_lll.z.xy;
		results.ptr[16] = inputU.d_lll.z.xz;
		results.ptr[17] = inputU.d_lll.z.yy;
		results.ptr[18] = inputU.d_lll.z.yz;
		results.ptr[19] = inputU.d_lll.z.zz;

		//V_j
		results.ptr[20] = inputU.V_l.x;
		results.ptr[21] = inputU.V_l.y;
		results.ptr[22] = inputU.V_l.z;
		
		sym3 K_sqrt_gammaUxx = sym3_real_mul(inputU.K_ll, eig.sqrt_gammaUjj.x);

		//a^x - f d^xj_j

		real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
		real d_x_input = sym3_dot(eig.gamma_uu, inputU.d_lll.x);
		results.ptr[23] = inputU.a_l.x - f * d_x_input;

		//gauge:
		//sqrt(f gamma^xx) K +- (a^x + 2 V^x)

		real ev0a = sqrt_f * sym3_dot(eig.gamma_uu, K_sqrt_gammaUxx);
		real ev0b = eig.gamma_uu.xx * (inputU.a_l.x + 2. * inputU.V_l.x) 
					+ eig.gamma_uu.xy * (inputU.a_l.y + 2. * inputU.V_l.y)
					+ eig.gamma_uu.xz * (inputU.a_l.z + 2. * inputU.V_l.z);
		results.ptr[0] = ev0a - ev0b;
		results.ptr[29] = ev0a + ev0b;

		//light:
		//sqrt(gamma^xx) K_xy +- (d^x_xy + .5 (a_y - d_yj^j) + V_y)

		real d_y_input = sym3_dot(eig.gamma_uu, inputU.d_lll.y);
		real dUx_xy_input = eig.gamma_uu.xx * inputU.d_lll.x.xy + eig.gamma_uu.xy * inputU.d_lll.y.xy + eig.gamma_uu.xz * inputU.d_lll.z.xy;
		real ev1b = .5 * (inputU.a_l.y - d_y_input) + inputU.V_l.y + dUx_xy_input;
		results.ptr[1] = K_sqrt_gammaUxx.xy - ev1b;
		results.ptr[24] = K_sqrt_gammaUxx.xy + ev1b;

		//light:
		//sqrt(gamma^xx) K_xz +- (d^x_xz + .5 (a_z - d_zj^j) + V_z)

		real d_z_input = sym3_dot(eig.gamma_uu, inputU.d_lll.z);
		real dUx_xz_input = eig.gamma_uu.xx * inputU.d_lll.x.xz + eig.gamma_uu.xy * inputU.d_lll.y.xz + eig.gamma_uu.xz * inputU.d_lll.z.xz;
		real ev2b = .5 * (inputU.a_l.z - d_z_input) + inputU.V_l.z + dUx_xz_input;
		results.ptr[2] = K_sqrt_gammaUxx.xz - ev2b;
		results.ptr[25] = K_sqrt_gammaUxx.xz + ev2b;

		//light:
		//sqrt(gamma^xx) K_yy +- d^x_yy

		real dUx_yy_input = eig.gamma_uu.xx * inputU.d_lll.x.yy + eig.gamma_uu.xy * inputU.d_lll.y.yy + eig.gamma_uu.xz * inputU.d_lll.z.yy;
		results.ptr[3] = K_sqrt_gammaUxx.yy - dUx_yy_input;
		results.ptr[26] = K_sqrt_gammaUxx.yy + dUx_yy_input;

		//light:
		//sqrt(gamma^xx) K_yz +- d^x_yz

		real dUx_yz_input = eig.gamma_uu.xx * inputU.d_lll.x.yz + eig.gamma_uu.xy * inputU.d_lll.y.yz + eig.gamma_uu.xz * inputU.d_lll.z.yz;
		results.ptr[4] = K_sqrt_gammaUxx.yz - dUx_yz_input; 
		results.ptr[27] = K_sqrt_gammaUxx.yz + dUx_yz_input;

		//light:
		//sqrt(gamma^xx) K_zz +- d^x_zz

		real dUx_zz_input = eig.gamma_uu.xx * inputU.d_lll.x.zz + eig.gamma_uu.xy * inputU.d_lll.y.zz + eig.gamma_uu.xz * inputU.d_lll.z.zz;
		results.ptr[5] = K_sqrt_gammaUxx.zz - dUx_zz_input;
		results.ptr[28] = K_sqrt_gammaUxx.zz + dUx_zz_input;

	} else if (n.side == 1) {

		//a_x, a_z
		results.ptr[6] = inputU.a_l.x;
		results.ptr[7] = inputU.a_l.z;

		//d_xij
		results.ptr[8] = inputU.d_lll.x.xx;
		results.ptr[9] = inputU.d_lll.x.xy;
		results.ptr[10] = inputU.d_lll.x.xz;
		results.ptr[11] = inputU.d_lll.x.yy;
		results.ptr[12] = inputU.d_lll.x.yz;
		results.ptr[13] = inputU.d_lll.x.zz;
		
		//d_zij
		results.ptr[14] = inputU.d_lll.z.xx;
		results.ptr[15] = inputU.d_lll.z.xy;
		results.ptr[16] = inputU.d_lll.z.xz;
		results.ptr[17] = inputU.d_lll.z.yy;
		results.ptr[18] = inputU.d_lll.z.yz;
		results.ptr[19] = inputU.d_lll.z.zz;
		
		//V_j
		results.ptr[20] = inputU.V_l.x;
		results.ptr[21] = inputU.V_l.y;
		results.ptr[22] = inputU.V_l.z;
		
		sym3 K_sqrt_gammaUyy = sym3_real_mul(inputU.K_ll, eig.sqrt_gammaUjj.y);

		//a^y - f d^yj_j

		real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
		real f = sqrt_f * sqrt_f;
		real d_y_input = sym3_dot(eig.gamma_uu, inputU.d_lll.y);
		results.ptr[23] = inputU.a_l.y - f * d_y_input;
		
		//gauge:
		//sqrt(f gamma^yy) K +- (a^y + 2 V^y)

		real ev0a = sqrt_f * sym3_dot(eig.gamma_uu, K_sqrt_gammaUyy);
		real ev0b = eig.gamma_uu.xy * (inputU.a_l.x + 2. * inputU.V_l.x)
					+ eig.gamma_uu.yy * (inputU.a_l.y + 2. * inputU.V_l.y)
					+ eig.gamma_uu.yz * (inputU.a_l.z + 2. * inputU.V_l.z);
		results.ptr[0] = ev0a - ev0b;
		results.ptr[29] = ev0a + ev0b;

		//light:
		//sqrt(gamma^yy) K_xx +- d^y_xx

		real dUy_xx_input = eig.gamma_uu.xy * inputU.d_lll.x.xx + eig.gamma_uu.yy * inputU.d_lll.y.xx + eig.gamma_uu.yz * inputU.d_lll.z.xx;
		results.ptr[1] = K_sqrt_gammaUyy.xx - dUy_xx_input;
		results.ptr[24] = K_sqrt_gammaUyy.xx + dUy_xx_input;

		//light:
		//sqrt(gamma^yy) K_xy +- (d^y_xy + .5 (a_x - d_xj^j) + V_x)

		real d_x_input = sym3_dot(eig.gamma_uu, inputU.d_lll.x);
		real dUy_xy_input = eig.gamma_uu.xy * inputU.d_lll.x.xy + eig.gamma_uu.yy * inputU.d_lll.y.xy + eig.gamma_uu.yz * inputU.d_lll.z.xy;
		real ev2b = dUy_xy_input + .5 * (inputU.a_l.x - d_x_input) + inputU.V_l.x;
		results.ptr[2] = K_sqrt_gammaUyy.xy - ev2b;
		results.ptr[25] = K_sqrt_gammaUyy.xy + ev2b;

		//light:
		//sqrt(gamma^yy) K_xz +- d^y_xz

		real dUy_xz_input = eig.gamma_uu.xy * inputU.d_lll.x.xz + eig.gamma_uu.yy * inputU.d_lll.y.xz + eig.gamma_uu.yz * inputU.d_lll.z.xz;
		results.ptr[3] = K_sqrt_gammaUyy.xz - dUy_xz_input;
		results.ptr[26] = K_sqrt_gammaUyy.xz + dUy_xz_input;

		//light:
		//sqrt(gamma^yy) K_yz +- (d^y_yz + .5 (a_z - d_zj^j) + V_z)

		real dUy_yz_input = eig.gamma_uu.xy * inputU.d_lll.x.yz + eig.gamma_uu.yy * inputU.d_lll.y.yz + eig.gamma_uu.yz * inputU.d_lll.z.yz;
		real d_z_input = sym3_dot(eig.gamma_uu, inputU.d_lll.z);
		real ev4b = dUy_yz_input + .5 * (inputU.a_l.z - d_z_input) + inputU.V_l.z;
		results.ptr[4] = K_sqrt_gammaUyy.yz - ev4b;
		results.ptr[27] = K_sqrt_gammaUyy.yz + ev4b;

		//light:
		//sqrt(gamma^yy) K_zz +- d^y_zz

		real dUy_zz_input = eig.gamma_uu.xy * inputU.d_lll.x.zz + eig.gamma_uu.yy * inputU.d_lll.y.zz + eig.gamma_uu.yz * inputU.d_lll.z.zz;
		results.ptr[5] = K_sqrt_gammaUyy.zz - dUy_zz_input;
		results.ptr[28] = K_sqrt_gammaUyy.zz - dUy_zz_input;
	
	} else if (n.side == 2) {

		//a_x, a_y
		results.ptr[6] = inputU.a_l.x;
		results.ptr[7] = inputU.a_l.y;
		
		//d_xij
		results.ptr[8] =  inputU.d_lll.x.xx;
		results.ptr[9] =  inputU.d_lll.x.xy;
		results.ptr[10] = inputU.d_lll.x.xz;
		results.ptr[11] = inputU.d_lll.x.yy;
		results.ptr[12] = inputU.d_lll.x.yz;
		results.ptr[13] = inputU.d_lll.x.zz;
		
		//d_yij
		results.ptr[14] = inputU.d_lll.y.xx;
		results.ptr[15] = inputU.d_lll.y.xy;
		results.ptr[16] = inputU.d_lll.y.xz;
		results.ptr[17] = inputU.d_lll.y.yy;
		results.ptr[18] = inputU.d_lll.y.yz;
		results.ptr[19] = inputU.d_lll.y.zz;
		
		//V_j
		results.ptr[20] = inputU.V_l.x;
		results.ptr[21] = inputU.V_l.y;
		results.ptr[22] = inputU.V_l.z;

		sym3 K_sqrt_gammaUzz = sym3_real_mul(inputU.K_ll, eig.sqrt_gammaUjj.z);

		//a^z - f d^zj_j

		real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
		real f = sqrt_f * sqrt_f;
		real d_z_input = sym3_dot(eig.gamma_uu, inputU.d_lll.z);
		results.ptr[23] = inputU.a_l.z - f * d_z_input;

		//gauge:
		//sqrt(f gamma^zz) K +- (a^z + 2 V^z)

		real ev0a = sqrt_f * sym3_dot(eig.gamma_uu, K_sqrt_gammaUzz);
		real ev0b = eig.gamma_uu.xz * (inputU.a_l.x + 2. * inputU.V_l.x)
					+ eig.gamma_uu.yz * (inputU.a_l.y + 2. * inputU.V_l.y)
					+ eig.gamma_uu.zz * (inputU.a_l.z + 2. * inputU.V_l.z);
		results.ptr[0] = ev0a - ev0b;
		results.ptr[29] = ev0a + ev0b;

		//light:
		//sqrt(gamma^zz) K_xx +- d^z_xx
		
		real dUz_xx_input = eig.gamma_uu.xz * inputU.d_lll.x.xx + eig.gamma_uu.yz * inputU.d_lll.y.xx + eig.gamma_uu.zz * inputU.d_lll.z.xx;
		results.ptr[1] = K_sqrt_gammaUzz.xx - dUz_xx_input;
		results.ptr[24] = K_sqrt_gammaUzz.xx + dUz_xx_input;

		//light:
		//sqrt(gamma^zz) K_xy +- d^z_xy

		real dUz_xy_input = eig.gamma_uu.xz * inputU.d_lll.x.xy + eig.gamma_uu.yz * inputU.d_lll.y.xy + eig.gamma_uu.zz * inputU.d_lll.z.xy;
		results.ptr[2] = K_sqrt_gammaUzz.xy - dUz_xy_input;
		results.ptr[25] = K_sqrt_gammaUzz.xy + dUz_xy_input;

		//light:
		//sqrt(gamma^zz) K_xz +- (d^z_xz + .5 (a_x - d_xj^j) + V_x)
		
		real d_x_input = sym3_dot(eig.gamma_uu, inputU.d_lll.x);
		real dUz_xz_input = eig.gamma_uu.xz * inputU.d_lll.x.xz + eig.gamma_uu.yz * inputU.d_lll.y.xz + eig.gamma_uu.zz * inputU.d_lll.z.xz;
		real ev3b = .5 * (inputU.a_l.x - d_x_input) + inputU.V_l.x + dUz_xz_input;
		results.ptr[3] = K_sqrt_gammaUzz.xz - ev3b;
		results.ptr[26] = K_sqrt_gammaUzz.xz + ev3b;

		//light:
		//sqrt(gamma^zz) K_yy +- d^z_yy

		real dUz_yy_input = eig.gamma_uu.xz * inputU.d_lll.x.yy + eig.gamma_uu.yz * inputU.d_lll.y.yy + eig.gamma_uu.zz * inputU.d_lll.z.yy;
		results.ptr[4] = K_sqrt_gammaUzz.yy - dUz_yy_input;
		results.ptr[27] = K_sqrt_gammaUzz.yy + dUz_yy_input;
		
		//light:
		//sqrt(gamma^zz) K_yz

		real d_y_input = sym3_dot(eig.gamma_uu, inputU.d_lll.y);
		real dUz_yz_input = eig.gamma_uu.xz * inputU.d_lll.x.yz + eig.gamma_uu.yz * inputU.d_lll.y.yz + eig.gamma_uu.zz * inputU.d_lll.z.yz;
		real ev5b = .5 * (inputU.a_l.y - d_y_input) + inputU.V_l.y + dUz_yz_input;
		results.ptr[5] = K_sqrt_gammaUzz.yz - ev5b;
		results.ptr[28] = K_sqrt_gammaUzz.yz + ev5b;

	}

<? else -- eqn.noZeroRowsInFlux ?>

	real _1_sqrt_f = eig.alpha / eig.alpha_sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f; 

	real sqrt_gammaUjj, _1_gammaUjj, a_j;
	sym3 d_lll, K_ll, gamma_uu;
	//now swap x and side on the sym3's
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		sqrt_gammaUjj = eig.sqrt_gammaUjj.s<?=side?>;
		_1_gammaUjj = 1. / eig.gamma_uu.s<?=side?><?=side?>;

		a_j = inputU.a_l.s<?=side?>;
		
		d_lll = sym3_swap<?=side?>(inputU.d_lll.v<?=side?>);
		K_ll = sym3_swap<?=side?>(inputU.K_ll);
		gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);
	}
	<? end ?>

	real K_dot_eig_gamma = sym3_dot(K_ll, gamma_uu);
	real dj_dot_eig_gamma = sym3_dot(d_lll, gamma_uu);

	results.ptr[0] = (a_j * -sqrt_gammaUjj * _1_sqrt_f + K_dot_eig_gamma) * .5 * _1_gammaUjj;

		<? for i=1,5 do ?>
	results.ptr[<?=i?>] = .5 * (-sqrt_gammaUjj * d_lll.s[<?=i?>] + K_ll.s[<?=i?>]);
		<? end ?>

	results.ptr[6] = (-a_j * _1_f + dj_dot_eig_gamma) * _1_gammaUjj;

		<? for i=1,5 do ?>
	results.ptr[<?=6+i?>] = .5 * (sqrt_gammaUjj * d_lll.s[<?=i?>] + K_ll.s[<?=i?>]);
		<? end ?>

	results.ptr[12] = (a_j * sqrt_gammaUjj * _1_sqrt_f + K_dot_eig_gamma) * .5 * _1_gammaUjj;
	
<? end -- eqn.noZeroRowsInFlux ?>
	return results;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t input,
	real3 x,
	normal_t n
) {
	cons_t resultU;
	for (int j = 0; j < numStates; ++j) {
		resultU.ptr[j] = 0;
	}

<? if not eqn.noZeroRowsInFlux then ?>
	
	if (n.side == 0) {

		//write zeros to the alpha and gammaLL terms
		resultU.alpha = 0;
		resultU.gamma_ll = sym3_zero;

		resultU.a_l.y = input.ptr[6];
		resultU.a_l.z = input.ptr[7];
		
		resultU.d_lll.y.xx = input.ptr[8];
		resultU.d_lll.y.xy = input.ptr[9];
		resultU.d_lll.y.xz = input.ptr[10];
		resultU.d_lll.y.yy = input.ptr[11];
		resultU.d_lll.y.yz = input.ptr[12];
		resultU.d_lll.y.zz = input.ptr[13];
		
		resultU.d_lll.z.xx = input.ptr[14];
		resultU.d_lll.z.xy = input.ptr[15];
		resultU.d_lll.z.xz = input.ptr[16];
		resultU.d_lll.z.yy = input.ptr[17];
		resultU.d_lll.z.yz = input.ptr[18];
		resultU.d_lll.z.zz = input.ptr[19];
		
		resultU.V_l.x = input.ptr[20];
		resultU.V_l.y = input.ptr[21];
		resultU.V_l.z = input.ptr[22];

		real sqrt_gammaUxx = eig.sqrt_gammaUjj.x;
		real _1_sqrt_gammaUxx = 1. / sqrt_gammaUxx;
		real _1_gammaUxx = _1_sqrt_gammaUxx * _1_sqrt_gammaUxx;
		real invDenom = .5 * _1_gammaUxx;

		real VUx_input = eig.gamma_uu.xx * input.ptr[20]		//V_x
						+ eig.gamma_uu.xy * input.ptr[21]	//V_y
						+ eig.gamma_uu.xz * input.ptr[22];	//V_z

		real gauge_diff = input.ptr[0] - input.ptr[29];
		real gauge_sum = input.ptr[0] + input.ptr[29];

		resultU.a_l.x = -(
				gauge_diff
				+ 2. * (
					eig.gamma_uu.xy * input.ptr[6]	//a_y
					+ eig.gamma_uu.xz * input.ptr[7]	//a_z
				)	
				+ 4. * VUx_input
			) * invDenom;
		
		real K_input_minus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[1]
				+ eig.gamma_uu.xz * input.ptr[2]
				+ eig.gamma_uu.yz * input.ptr[4]
			) + eig.gamma_uu.yy * input.ptr[3]
			+ eig.gamma_uu.zz * input.ptr[5];

		real K_input_plus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[24]
				+ eig.gamma_uu.xz * input.ptr[25]
				+ eig.gamma_uu.yz * input.ptr[27]
			) + eig.gamma_uu.yy * input.ptr[26]
			+ eig.gamma_uu.zz * input.ptr[28];

		real _1_sqrt_f = eig.alpha / eig.alpha_sqrt_f;
		real _1_f = _1_sqrt_f * _1_sqrt_f;

		resultU.d_lll.x.xx = -(
				- K_input_minus
				+ K_input_plus

				+ 2. * eig.gamma_uu.xy * (
					eig.gamma_uu.xx * input.ptr[8]	//d_yxx
					- input.ptr[6]		//a_y
					- 2. * input.ptr[21]	//V_y
				)
				+ 2. * eig.gamma_uu.xz * (
					eig.gamma_uu.xx * input.ptr[14]	//d_zxx
					- input.ptr[7]		//a_z
					- 2. * input.ptr[22]	//V_z
				)
				//(a^x + 2 V^x) / f - 2 gamma^xx d^xj_j
				+ (
					+ gauge_diff	//-ev0b = -a^x - 2 V^x
					+ 2. * (
						+ eig.gamma_uu.xx * input.ptr[23]	//a_x - f d^xj_j
						+ eig.gamma_uu.xy * input.ptr[6]	//a_y
						+ eig.gamma_uu.xz * input.ptr[7]	//a_z
						+ 2. * VUx_input
					)
				) * _1_f
			) * invDenom * _1_gammaUxx;

		//d_yj^j
		real d_y_input = 
			eig.gamma_uu.xx * input.ptr[8]
			+ eig.gamma_uu.yy * input.ptr[11]
			+ eig.gamma_uu.zz * input.ptr[13]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[9]
				+ eig.gamma_uu.xz * input.ptr[10]
				+ eig.gamma_uu.yz * input.ptr[12]
			);
		
		resultU.d_lll.x.xy = -(
				input.ptr[1]
				+ input.ptr[6]
				- d_y_input
				+ 2. * eig.gamma_uu.xy * input.ptr[9]
				+ 2. * eig.gamma_uu.xz * input.ptr[15]
				+ 2. * input.ptr[21]
				- input.ptr[24]
			) * invDenom;
		
		//d_zj^j
		real d_z_input = 
			eig.gamma_uu.xx * input.ptr[14]
			+ eig.gamma_uu.yy * input.ptr[17]
			+ eig.gamma_uu.zz * input.ptr[19]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[15]
				+ eig.gamma_uu.xz * input.ptr[16]
				+ eig.gamma_uu.yz * input.ptr[18]
			);

		resultU.d_lll.x.xz = -(
				input.ptr[2]
				+ input.ptr[7]
				- d_z_input
				+ 2. * eig.gamma_uu.xy * input.ptr[10]
				+ 2. * eig.gamma_uu.xz * input.ptr[16]
				+ 2. * input.ptr[22]
				- input.ptr[25]
			) * invDenom;
		resultU.d_lll.x.yy = -(
				input.ptr[3]
				+ 2. * eig.gamma_uu.xy * input.ptr[11]
				+ 2. * eig.gamma_uu.xz * input.ptr[17]
				- input.ptr[26]
			) * invDenom;
		resultU.d_lll.x.yz = -(
				input.ptr[4]
				+ 2. * eig.gamma_uu.xy * input.ptr[12]
				+ 2. * eig.gamma_uu.xz * input.ptr[18]
				- input.ptr[27]
			) * invDenom;
		resultU.d_lll.x.zz = -(
				input.ptr[5]
				+ 2. * eig.gamma_uu.xy * input.ptr[13]
				+ 2. * eig.gamma_uu.xz * input.ptr[19]
				- input.ptr[28]
			) * invDenom;

		resultU.K_ll.xx = (
				- K_input_minus
				- K_input_plus	
				+ gauge_sum * _1_sqrt_f
			) * invDenom * _1_sqrt_gammaUxx;

		real tmp = .5 * _1_sqrt_gammaUxx;
		resultU.K_ll.xy = (input.ptr[1] + input.ptr[24]) * tmp;
		resultU.K_ll.xz = (input.ptr[2] + input.ptr[25]) * tmp;
		resultU.K_ll.yy = (input.ptr[3] + input.ptr[26]) * tmp;
		resultU.K_ll.yz = (input.ptr[4] + input.ptr[27]) * tmp;
		resultU.K_ll.zz = (input.ptr[5] + input.ptr[28]) * tmp;

	} else if (n.side == 1) {
	
		//write zeros to the alpha and gammaLL terms
		resultU.alpha = 0;
		resultU.gamma_ll = sym3_zero;
		
		resultU.a_l.x = input.ptr[6];
		resultU.a_l.z = input.ptr[7];
		
		resultU.d_lll.x.xx = input.ptr[8];
		resultU.d_lll.x.xy = input.ptr[9];
		resultU.d_lll.x.xz = input.ptr[10];
		resultU.d_lll.x.yy = input.ptr[11];
		resultU.d_lll.x.yz = input.ptr[12];
		resultU.d_lll.x.zz = input.ptr[13];
		
		resultU.d_lll.z.xx = input.ptr[14];
		resultU.d_lll.z.xy = input.ptr[15];
		resultU.d_lll.z.xz = input.ptr[16];
		resultU.d_lll.z.yy = input.ptr[17];
		resultU.d_lll.z.yz = input.ptr[18];
		resultU.d_lll.z.zz = input.ptr[19];
		
		resultU.V_l.x = input.ptr[20];
		resultU.V_l.y = input.ptr[21];
		resultU.V_l.z = input.ptr[22];

		real sqrt_gammaUyy = eig.sqrt_gammaUjj.y;
		real _1_sqrt_gammaUyy = 1. / sqrt_gammaUyy;
		real inv_gammaUyy = _1_sqrt_gammaUyy * _1_sqrt_gammaUyy;
		real invDenom = .5 * inv_gammaUyy;

		real VUy_input = eig.gamma_uu.xy * input.ptr[20]
						+ eig.gamma_uu.yy * input.ptr[21]
						+ eig.gamma_uu.yz * input.ptr[22];
		
		real gauge_diff = input.ptr[0] - input.ptr[29];
		real gauge_sum = input.ptr[0] + input.ptr[29];

		resultU.a_l.y = -(
				gauge_diff
				+ 2. * (
					eig.gamma_uu.xy * input.ptr[6]
					+ eig.gamma_uu.yz * input.ptr[7]
				)
				+ 4. * VUy_input
			) * invDenom;

		resultU.d_lll.y.xx = -(
				+ input.ptr[1]
				- input.ptr[24]
				+ 2. * eig.gamma_uu.xy * input.ptr[8]
				+ 2. * eig.gamma_uu.yz * input.ptr[14]
			) * invDenom;

		//d_xj^j
		real d_x_input = 
			eig.gamma_uu.xx * input.ptr[8]
			+ eig.gamma_uu.yy * input.ptr[11] 
			+ eig.gamma_uu.zz * input.ptr[13]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[9]
				+ eig.gamma_uu.xz * input.ptr[10]
				+ eig.gamma_uu.yz * input.ptr[12]
			);
		
		resultU.d_lll.y.xy = -(
				+ input.ptr[2]
				+ input.ptr[6]
				- d_x_input			
				+ 2. * eig.gamma_uu.xy * input.ptr[9]
				+ 2. * eig.gamma_uu.yz * input.ptr[15]
				+ 2. * input.ptr[20]
				- input.ptr[25]
			) * invDenom;
		
		resultU.d_lll.y.xz = -(
				+ input.ptr[3]
				+ 2. * eig.gamma_uu.xy * input.ptr[10]
				+ 2. * eig.gamma_uu.yz * input.ptr[16]
				- input.ptr[26]
			) * invDenom;
		
		real K_input_minus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[2]
				+ eig.gamma_uu.xz * input.ptr[3]
				+ eig.gamma_uu.yz * input.ptr[4]
			) + eig.gamma_uu.xx * input.ptr[1]
			+ eig.gamma_uu.zz * input.ptr[5];

		real K_input_plus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[25]
				+ eig.gamma_uu.xz * input.ptr[26]
				+ eig.gamma_uu.yz * input.ptr[27]
			) + eig.gamma_uu.xx * input.ptr[24]
			+ eig.gamma_uu.zz * input.ptr[28];

		real _1_sqrt_f = eig.alpha / eig.alpha_sqrt_f;
		real _1_f = _1_sqrt_f  * _1_sqrt_f;

		resultU.d_lll.y.yy = -(
				- K_input_minus	
				+ K_input_plus	
				
				+ 2. * eig.gamma_uu.xy * (
					input.ptr[11] * eig.gamma_uu.yy
					- input.ptr[6]
					- 2. * input.ptr[20]
				)	
				+ 2. * eig.gamma_uu.yz * (
					input.ptr[17] * eig.gamma_uu.yy
					- input.ptr[7]
					- 2. * input.ptr[22]
				)
				+ (
					+ gauge_diff
					+ 2. * (
						+ eig.gamma_uu.yy * input.ptr[23]
						+ eig.gamma_uu.xy * input.ptr[6]
						+ eig.gamma_uu.yz * input.ptr[7]
						+ 2. * VUy_input
					)
				) * _1_f
			) * invDenom * inv_gammaUyy;

		//gamma_zj^j
		real d_z_input = eig.gamma_uu.xx * input.ptr[14]
			+ eig.gamma_uu.yy * input.ptr[17]
			+ eig.gamma_uu.zz * input.ptr[19]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[15]
				+ eig.gamma_uu.xz * input.ptr[16]
				+ eig.gamma_uu.yz * input.ptr[18]);

		resultU.d_lll.y.yz = -(
				+ input.ptr[4]
				+ input.ptr[7]
				- d_z_input	
				+ 2. * eig.gamma_uu.xy * input.ptr[12]
				+ 2. * eig.gamma_uu.yz * input.ptr[18]
				+ 2. * input.ptr[22]
				- input.ptr[27]
			) * invDenom;
		
		resultU.d_lll.y.zz = -(
				+ input.ptr[5]
				- input.ptr[28]
				+ 2. * eig.gamma_uu.xy * input.ptr[13]
				+ 2. * eig.gamma_uu.yz * input.ptr[19]
			) * invDenom;

		resultU.K_ll.yy = (
				- K_input_minus
				- K_input_plus
				+ gauge_sum * _1_sqrt_f
			) * invDenom * _1_sqrt_gammaUyy;

		real tmp = .5 * _1_sqrt_gammaUyy;
		resultU.K_ll.xx = (input.ptr[1] + input.ptr[24]) * tmp;
		resultU.K_ll.xy = (input.ptr[2] + input.ptr[25]) * tmp;		//once this gets enabled, compiling crashes
		resultU.K_ll.xz = (input.ptr[3] + input.ptr[26]) * tmp;
		resultU.K_ll.yz = (input.ptr[4] + input.ptr[27]) * tmp;
		resultU.K_ll.zz = (input.ptr[5] + input.ptr[28]) * tmp;
	
	} else if (n.side == 2) {

		//write zeros to the alpha and gammaLL terms
		resultU.alpha = 0;
		resultU.gamma_ll = sym3_zero;
		
		resultU.a_l.x = input.ptr[6];
		resultU.a_l.y = input.ptr[7];
		
		resultU.d_lll.x.xx = input.ptr[8];
		resultU.d_lll.x.xy = input.ptr[9];
		resultU.d_lll.x.xz = input.ptr[10];
		resultU.d_lll.x.yy = input.ptr[11];
		resultU.d_lll.x.yz = input.ptr[12];
		resultU.d_lll.x.zz = input.ptr[13];
		
		resultU.d_lll.y.xx = input.ptr[14];
		resultU.d_lll.y.xy = input.ptr[15];
		resultU.d_lll.y.xz = input.ptr[16];
		resultU.d_lll.y.yy = input.ptr[17];
		resultU.d_lll.y.yz = input.ptr[18];
		resultU.d_lll.y.zz = input.ptr[19];

		resultU.V_l.x = input.ptr[20];
		resultU.V_l.y = input.ptr[21];
		resultU.V_l.z = input.ptr[22];

		real sqrt_gammaUzz = eig.sqrt_gammaUjj.z;
		real _1_sqrt_gammaUzz = 1. / sqrt_gammaUzz;
		real inv_gammaUzz = _1_sqrt_gammaUzz * _1_sqrt_gammaUzz;
		real invDenom = .5 * inv_gammaUzz;

		real VUz_input = eig.gamma_uu.xz * input.ptr[20]
						+ eig.gamma_uu.yz * input.ptr[21]
						+ eig.gamma_uu.zz * input.ptr[22];

		real gauge_diff = input.ptr[0] - input.ptr[29];
		real gauge_sum = input.ptr[0] + input.ptr[29];
		
		resultU.a_l.z = -(
				gauge_diff
				+ 2. * (
					eig.gamma_uu.xz * input.ptr[6]
					+ eig.gamma_uu.yz * input.ptr[7]
				)
				+ 4. * VUz_input
			) * invDenom;
		
		resultU.d_lll.z.xx = -(
				+ input.ptr[1]
				- input.ptr[24]
				+ 2. * eig.gamma_uu.xz * input.ptr[8]
				+ 2. * eig.gamma_uu.yz * input.ptr[14]
			) * invDenom;
		
		resultU.d_lll.z.xy = -(
				+ input.ptr[2]
				- input.ptr[25]
				+ 2. * eig.gamma_uu.xz * input.ptr[9]
				+ 2. * eig.gamma_uu.yz * input.ptr[15]
			) * invDenom;

		//d_xj^j
		real d_x_input = 
			eig.gamma_uu.xx * input.ptr[8]
			+ eig.gamma_uu.yy * input.ptr[11]
			+ eig.gamma_uu.zz * input.ptr[13]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[9]
				+ eig.gamma_uu.xz * input.ptr[10]
				+ eig.gamma_uu.yz * input.ptr[12]);

		resultU.d_lll.z.xz = -(
				+ input.ptr[3]
				+ input.ptr[6]
				- d_x_input
				+ 2. * eig.gamma_uu.xz * input.ptr[10]
				+ 2. * eig.gamma_uu.yz * input.ptr[16]
				+ 2. * input.ptr[20]
				- input.ptr[26]
			) * invDenom;
		
		resultU.d_lll.z.yy = -(
				+ input.ptr[4]
				- input.ptr[27]
				+ 2. * eig.gamma_uu.xz * input.ptr[11]
				+ 2. * eig.gamma_uu.yz * input.ptr[17]
			) * invDenom;

		//d_yj^j
		real d_y_input = 
			eig.gamma_uu.xx * input.ptr[14]
			+ eig.gamma_uu.yy * input.ptr[17]
			+ eig.gamma_uu.zz * input.ptr[19]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[15]
				+ eig.gamma_uu.xz * input.ptr[16]
				+ eig.gamma_uu.yz * input.ptr[18]);

		resultU.d_lll.z.yz = -(
				+ input.ptr[5]
				+ input.ptr[7]
				- d_y_input	
				+ 2. * eig.gamma_uu.xz * input.ptr[12]
				+ 2. * eig.gamma_uu.yz * input.ptr[18]
				+ 2. * input.ptr[21]
				- input.ptr[28]
			) * invDenom;

		real K_input_minus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[2]
				+ eig.gamma_uu.xz * input.ptr[3]
				+ eig.gamma_uu.yz * input.ptr[5]
			) + eig.gamma_uu.xx * input.ptr[1]
			+ eig.gamma_uu.yy * input.ptr[4];

		real K_input_plus = 
			2. * (eig.gamma_uu.xy * input.ptr[25]
				+ eig.gamma_uu.xz * input.ptr[26]
				+ eig.gamma_uu.yz * input.ptr[28]
			) + eig.gamma_uu.xx * input.ptr[24]
			+ eig.gamma_uu.yy * input.ptr[27];

		real _1_sqrt_f = eig.alpha / eig.alpha_sqrt_f;
		real _1_f = _1_sqrt_f * _1_sqrt_f;

		resultU.d_lll.z.zz = -(
				- K_input_minus
				+ K_input_plus
				
				+ 2. * eig.gamma_uu.xz * (
					input.ptr[13] * eig.gamma_uu.zz
					- input.ptr[6]
					- 2. * input.ptr[20]
				)
				+ 2. * eig.gamma_uu.yz * (
					eig.gamma_uu.zz * input.ptr[19]
					- input.ptr[7]
					- 2. * input.ptr[21]
				)	
				+ (	
					+ gauge_diff
					+ 2. * (
						eig.gamma_uu.zz * input.ptr[23]
						+ eig.gamma_uu.xz * input.ptr[6]
						+ eig.gamma_uu.yz * input.ptr[7]
						+ 2. * VUz_input
					)
				) * _1_f
			) * invDenom * inv_gammaUzz;

		resultU.K_ll.zz = (
				+ gauge_sum * _1_sqrt_f
				- K_input_minus
				- K_input_plus
			) * invDenom * _1_sqrt_gammaUzz;
		
		real tmp = .5 * _1_sqrt_gammaUzz;
		resultU.K_ll.xx = (input.ptr[1] + input.ptr[24]) * tmp;
		resultU.K_ll.xy = (input.ptr[2] + input.ptr[25]) * tmp;
		resultU.K_ll.xz = (input.ptr[3] + input.ptr[26]) * tmp;
		resultU.K_ll.yy = (input.ptr[4] + input.ptr[27]) * tmp;
		resultU.K_ll.yz = (input.ptr[5] + input.ptr[28]) * tmp;
		
	}

<? else -- eqn.noZeroRowsInFlux ?>
	
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		//TODO swap size inside eigen_t structure
		//instead of doing it here
		sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);

		real input1_dot_gammaU = input.ptr[1] * 2. * gamma_uu.xy
			+ input.ptr[2] * 2. * gamma_uu.xz
			+ input.ptr[3] * gamma_uu.yy
			+ input.ptr[4] * 2. * gamma_uu.yz
			+ input.ptr[5] * gamma_uu.zz;
		real input7_dot_gammaU = input.ptr[7] * 2. * gamma_uu.xy
			+ input.ptr[8] * 2. * gamma_uu.xz
			+ input.ptr[9] * gamma_uu.yy
			+ input.ptr[10] * 2. * gamma_uu.yz
			+ input.ptr[11] * gamma_uu.zz;

		real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
		real _1_sqrt_f = 1. / sqrt_f;
		real sqrt_gammaUjj = eig.sqrt_gammaUjj.s<?=side?>;
		real _1_sqrt_gammaUjj = 1. / sqrt_gammaUjj;
		real _1_gammaUjj = _1_sqrt_gammaUjj * _1_sqrt_gammaUjj; 

		resultU.a_l.s<?=side?> = sqrt_f * sqrt_gammaUjj * (input.ptr[12] - input.ptr[0]);

		sym3 d_lll, K_ll;
		d_lll.xx = (
			(input.ptr[12] - input.ptr[0]) * _1_sqrt_f
			+ (input1_dot_gammaU - input7_dot_gammaU) * _1_gammaUjj
		) * _1_sqrt_gammaUjj + input.ptr[6];

		<? for i=1,5 do ?>
		d_lll.s[<?=i?>] = (input.ptr[<?=i+6?>] - input.ptr[<?=i?>]) * _1_sqrt_gammaUjj;
		<? end ?>

		K_ll.xx = input.ptr[0] + input.ptr[12] - (input1_dot_gammaU + input7_dot_gammaU) * _1_gammaUjj;

		<? for i=1,5 do ?>
		K_ll.s[<?=i?>] = input.ptr[<?=i?>] + input.ptr[<?=i+6?>];
		<? end ?>

		//now swap x and side on the sym3's
		resultU.d_lll.v<?=side?> = sym3_swap<?=side?>(d_lll);
		resultU.K_ll = sym3_swap<?=side?>(K_ll);
	}
	<? end ?>

<? end 	-- eqn.noZeroRowsInFlux ?>
	return resultU;
}

<? elseif moduleName == "eigen_fluxTransform" then 
depmod{
	"eigen_t",
	"waves_t",
	"rotate",
}
?>

typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;

cons_t eigen_fluxTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x,
	normal_t n
) {
<? if not eqn.noZeroRowsInFlux then ?>

	// TODO make this a default implementation somewhere.
	// If no one has provided one then just fall back on left / wave / right transform.
	// TODO use that static function for the calc waves as well
	

	waves_t waves = eigen_leftTransform(solver, eig, inputU, x);

	<?=eqn:eigenWaveCodePrefix(n, 'eig', 'x')?>

<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode(n, 'eig', 'x', j)?>;
<? end 
?>
	return eigen_rightTransform(solver, eig, waves, x);

<? else -- noZeroRowsInFlux ?>
<? if false then 	-- by-hand ?>

	cons_t resultU;
	for (int i = 0; i < numStates; ++i) {
		resultU.ptr[i] = 0;
	}

	real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
	real f = sqrt_f * sqrt_f;

	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {

		//now swap x and side on the sym3's
		sym3 input_d = sym3_swap<?=side?>(inputU.d_lll.v<?=side?>);
		sym3 input_K = sym3_swap<?=side?>(inputU.K_ll);
		sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);

		resultU.a_l.s<?=side?> = sym3_dot(input_K, gamma_uu) * eig.alpha * f;
		sym3 result_d = sym3_real_mul(input_K, eig.alpha);
		sym3 result_K = sym3_real_mul(input_d, eig.alpha * gamma_uu.xx);
		result_K.xx += (inputU.a_l.s<?=side?> - sym3_dot(input_d, gamma_uu)) * eig.alpha;

		//now swap x and side on the sym3's
		resultU.d_lll.v<?=side?> = sym3_swap<?=side?>(result_d);
		resultU.K_ll = sym3_swap<?=side?>(result_K);

	}
	<? end ?>

	return resultU;
<? else	-- codegen ?>
	real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
	real f = sqrt_f * sqrt_f;
	
	sym3 gamma_uu = eig.gamma_uu;
	
	real alpha = inputU.alpha;
	real3 V_l = inputU.V_l;
	real3 a_l = inputU.a_l;
	sym3 gamma_ll = inputU.gamma_ll;
	sym3 K_ll = inputU.K_ll;
	_3sym3 d_lll = inputU.d_lll;
	
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		V_l = real3_swap<?=side?>(V_l);
		a_l = real3_swap<?=side?>(a_l);
		gamma_ll = sym3_swap<?=side?>(gamma_ll);
		K_ll = sym3_swap<?=side?>(K_ll);
		d_lll = _3sym3_swap<?=side?>(d_lll);
		gamma_uu = sym3_swap<?=side?>(gamma_uu);
	}
	<? end ?>

	<?=eqn.cons_t?> resultU = {.ptr={0}};

	// BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/adm_noZeroRows.html
	resultU.a_l.x = eig.alpha * f * (2. * K_ll.xy * gamma_uu.xy + 2. * K_ll.xz * gamma_uu.xz + 2. * K_ll.yz * gamma_uu.yz + K_ll.xx * gamma_uu.xx + K_ll.yy * gamma_uu.yy + K_ll.zz * gamma_uu.zz);
	resultU.d_lll.x.xx = K_ll.xx * eig.alpha;
	resultU.d_lll.x.xy = K_ll.xy * eig.alpha;
	resultU.d_lll.x.xz = K_ll.xz * eig.alpha;
	resultU.d_lll.x.yy = K_ll.yy * eig.alpha;
	resultU.d_lll.x.yz = K_ll.yz * eig.alpha;
	resultU.d_lll.x.zz = K_ll.zz * eig.alpha;
	resultU.K_ll.xx = eig.alpha * (a_l.x + d_lll.x.yy * gamma_uu.yy + 2. * d_lll.x.yz * gamma_uu.yz + d_lll.x.zz * gamma_uu.zz - d_lll.y.xx * gamma_uu.xy - 2. * d_lll.y.xy * gamma_uu.yy - 2. * d_lll.y.xz * gamma_uu.yz - d_lll.z.xx * gamma_uu.xz - 2. * d_lll.z.xy * gamma_uu.yz - 2. * d_lll.z.xz * gamma_uu.zz);
	resultU.K_ll.xy = (eig.alpha * (a_l.y - 2. * d_lll.x.yy * gamma_uu.xy - 2. * d_lll.x.yz * gamma_uu.xz - 2. * d_lll.y.yy * gamma_uu.yy - 2. * d_lll.y.yz * gamma_uu.yz - 2. * d_lll.z.yy * gamma_uu.yz - 2. * d_lll.z.yz * gamma_uu.zz)) / 2.;
	resultU.K_ll.xz = (eig.alpha * (a_l.z - 2. * d_lll.x.yz * gamma_uu.xy - 2. * d_lll.x.zz * gamma_uu.xz - 2. * d_lll.y.yz * gamma_uu.yy - 2. * d_lll.y.zz * gamma_uu.yz - 2. * d_lll.z.yz * gamma_uu.yz - 2. * d_lll.z.zz * gamma_uu.zz)) / 2.;
	resultU.K_ll.yy = eig.alpha * (d_lll.x.yy * gamma_uu.xx + d_lll.y.yy * gamma_uu.xy + d_lll.z.yy * gamma_uu.xz);
	resultU.K_ll.yz = eig.alpha * (d_lll.x.yz * gamma_uu.xx + d_lll.y.yz * gamma_uu.xy + d_lll.z.yz * gamma_uu.xz);
	resultU.K_ll.zz = eig.alpha * (d_lll.x.zz * gamma_uu.xx + d_lll.y.zz * gamma_uu.xy + d_lll.z.zz * gamma_uu.xz);
	// END CUT

	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		resultU.V_l = real3_swap<?=side?>(resultU.V_l);
		resultU.a_l = real3_swap<?=side?>(resultU.a_l);
		resultU.gamma_ll = sym3_swap<?=side?>(resultU.gamma_ll);
		resultU.K_ll = sym3_swap<?=side?>(resultU.K_ll);
		resultU.d_lll = _3sym3_swap<?=side?>(resultU.d_lll);
	}
	<? end ?>

	return resultU;
<? end ?>
<? end -- noZeroRowsInFlux ?>
}

<? elseif moduleName == "addSource" then 
depmod{
	"initCond.codeprefix",		-- calc_f
	"real3x3x3",
}
?>

/*
this should just be 

alpha_,t + F^i^alpha_,i = -f alpha gamma^ij K_ij
*/


//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? if eqn.useStressEnergyTerms then ?> 
	real rho = U->rho;
	real3 S_l = sym3_real3_mul(U->gamma_ll, U->S_u);
	sym3 S_ll = U->S_ll;
	real S = sym3_dot(S_ll, gamma_uu);
<? else ?>
	real rho = 0.;
	real3 S_l = real3_zero;
	sym3 S_ll = sym3_zero;
	real S = 0.;
<? end ?>


	// source terms

	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real tr_K = real3x3_trace(K_ul);							//K^k_k
	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		//KSq_ij = K_ik K^k_j

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3x3 d_llu = _3sym3_sym3_mul(U->d_lll, gamma_uu);

	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//e_i = d^j_ji
	real3 e_l = _3sym3_tr12(d_ull);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>			.<?=xij?> = d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?> - d_ull.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	// Rsrc_ll.ij := Rsrc_ij
	//	= 	+ conn^k_ij (V_k - e_k)
	//		+ 2 d^l_ki d_lj^k 
	//		- 2 d^l_ki d^k_lj 
	//		+ d_il^k d_jk^l
	sym3 Rsrc_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * (U->V_l.<?=xk?> - e_l.<?=xk?>)
<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_ull.<?=xl?>.<?=sym(k,i)?> * (d_llu.<?=xl?>.<?=xj?>.<?=xk?> - d_ull.<?=xk?>.<?=sym(l,j)?>)
			+ d_llu.<?=xi?>.<?=xl?>.<?=xk?> * d_llu.<?=xj?>.<?=xk?>.<?=xl?>
<? 		end
	end
?>		,
<? end
?>	};


	//srcK_ij = 
	
	//lapse:
	//	(-a_i a_j 
	//		+ conn^k_ij a_k

	//Riemann curvature:
	//		+ conn^k_ij (V_k - d^l_lk) 
	//		+ 2 d^l_ki d_lj^k 
	//		- 2 d^l_ki d^k_lj 
	//		+ d_il^k d_jk^l
	
	//extrinsic curvature:
	//		+ K K_ij 
	//		- 2 K_ik K^k_j

	//matter terms:
	//		- 8 pi S_ij 
	//		+ 4 pi  gamma_ij (S - rho)
	
	//...and the shift comes later ...

	sym3 stressConstraint_ll = (sym3){
<? for ij,xij in ipairs(symNames) do	
?>		.<?=xij?> = 
			Rsrc_ll.<?=xij?> 
			+ tr_K * U->K_ll.<?=xij?> 
			- KSq_ll.<?=xij?>
			- 8. * M_PI * S_ll.<?=xij?> 
			+ 4. * M_PI * U->gamma_ll.<?=xij?> * (S - rho)
		,
<? end
?>	};

	real HamiltonianConstraint = sym3_dot(stressConstraint_ll, gamma_uu);

#if 0	//hand-rolled

	sym3 srcK_ll_over_alpha = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 
			- U->a_l.<?=xi?> * U->a_l.<?=xj?>
<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * U->a_l.<?=xk?>
<?	end
?>			- KSq_ll.<?=xij?>
		,
<? end
?>	};


	//d_i = d_ij^j
	real3 d_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_trace(d_llu.<?=xi?>),
<? end
?>	};

	//srcV_k = -a_k K 
	//			+ 8 pi j_k
	//			+ a^j K_jk 
	//			- 4 K d^j_jk 
	//			+ 2 K^i_k d^j_ji 
	//			+ K^i_k d_ij^j 
	//			- K^i_j d_ki^j
	//			+ 2 K^i_j d_ik^j 
	real3 srcV_l = (real3){
<? for k,xk in ipairs(xNames) do
?>		.<?=xk?> = -tr_K * U->a_l.<?=xk?>
			- 8. * M_PI * S_l.<?=xk?>
<?	for j,xj in ipairs(xNames) do
?>			+ U->a_l.<?=xj?> * K_ul.<?=xj?>.<?=xk?>
			- 4. * tr_K * d_ull.<?=xj?>.<?=sym(j,k)?>
			+ 2. * K_ul.<?=xj?>.<?=xk?> * e_l.<?=xj?>
			+ 2. * K_ul.<?=xj?>.<?=xk?> * d_l.<?=xj?>
<?		for i,xi in ipairs(xNames) do
?>			- K_ul.<?=xi?>.<?=xj?> * d_llu.<?=xk?>.<?=xi?>.<?=xj?>
			+ 2. * K_ul.<?=xi?>.<?=xj?> * d_llu.<?=xi?>.<?=xk?>.<?=xj?>
<? 		end
	end ?>,
<? end
?>	};

	real f = calc_f(U->alpha);	//could be based on alpha...

	//alpha_,t = shift terms - alpha^2 f gamma^ij K_ij
	deriv->alpha += -U->alpha * U->alpha * f * tr_K;
	
	//gamma_ij,t = shift terms - 2 alpha K_ij
<? for ij,xij in ipairs(symNames) do
?>	deriv->gamma_ll.<?=xij?> += -2. * U->alpha * U->K_ll.<?=xij?>;
<? end
?>

	//K_ij,t = shift terms + alpha srcK_ij
<? for ij,xij in ipairs(symNames) do
?>	deriv->K_ll.<?=xij?> += U->alpha * srcK_ll_over_alpha.<?=xij?>;
<? end
?>
	//V_k,t = shift terms + alpha srcV_k
<? for i,xi in ipairs(xNames) do
?>	deriv->V_l.<?=xi?> += U->alpha * srcV_l.<?=xi?>;
<? end
?>
#else	//code-generated

	real alpha = U->alpha;
	sym3 gamma_ll = U->gamma_ll;
	real3 a_l = U->a_l;
	_3sym3 d_lll = U->d_lll;
	sym3 K_ll = U->K_ll;
	
	real f = calc_f(alpha);									//order 1/alpha
	real f_alpha = calc_f_alpha(alpha);						//order 1
	real f_alphaSq = calc_f_alphaSq(alpha);					//order alpha
	real dalpha_f = calc_dalpha_f(alpha);					//order 1/alpha^2
	real alphaSq_dalpha_f = calc_alphaSq_dalpha_f(alpha);	//order 1

	real3 d_l = _3sym3_sym3_dot23(d_lll, gamma_uu);
	
	// BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/adm_noZeroRows.html
	real tmp1 = K_ll.xy * gamma_uu.xy;
	real tmp2 = K_ll.xz * gamma_uu.xz;
	real tmp3 = K_ll.yz * gamma_uu.yz;
	real tmp5 = K_ll.zz * gamma_uu.zz;
	real tmp6 = K_ll.yy * gamma_uu.yy;
	real tmp23 = gamma_uu.xx * gamma_uu.xx;
	real tmp29 = gamma_uu.xx * gamma_uu.xy;
	real tmp35 = gamma_uu.xx * gamma_uu.xz;
	real tmp41 = gamma_uu.xy * gamma_uu.xy;
	real tmp47 = gamma_uu.xy * gamma_uu.xz;
	real tmp53 = gamma_uu.xz * gamma_uu.xz;
	real tmp72 = gamma_uu.xx * gamma_uu.yy;
	real tmp84 = gamma_uu.xx * gamma_uu.yz;
	real tmp96 = gamma_uu.xy * gamma_uu.yy;
	real tmp102 = gamma_uu.xy * gamma_uu.yz;
	real tmp108 = gamma_uu.xz * gamma_uu.yy;
	real tmp114 = gamma_uu.xz * gamma_uu.yz;
	real tmp145 = gamma_uu.xx * gamma_uu.zz;
	real tmp163 = gamma_uu.xy * gamma_uu.zz;
	real tmp175 = gamma_uu.xz * gamma_uu.zz;
	real tmp204 = gamma_uu.yy * gamma_uu.yy;
	real tmp210 = gamma_uu.yy * gamma_uu.yz;
	real tmp216 = gamma_uu.yz * gamma_uu.yz;
	real tmp265 = gamma_uu.yy * gamma_uu.zz;
	real tmp277 = gamma_uu.yz * gamma_uu.zz;
	real tmp318 = gamma_uu.zz * gamma_uu.zz;
	real tmp1020 = a_l.x * alpha;
	real tmp1032 = a_l.y * alpha;
	real tmp1044 = a_l.z * alpha;
	real tmp1090 = d_lll.y.xx * gamma_uu.xy;
	real tmp1091 = d_lll.z.xx * gamma_uu.xz;
	real tmp1123 = d_lll.x.yz * tmp84;
	real tmp1128 = d_lll.x.zz * tmp145;
	real tmp1132 = d_lll.y.xy * tmp72;
	real tmp1135 = d_lll.y.xz * tmp84;
	real tmp1138 = d_lll.y.yy * tmp96;
	real tmp1140 = d_lll.y.yz * tmp108;
	real tmp1143 = d_lll.y.zz * tmp163;
	real tmp1147 = d_lll.y.zz * tmp114;
	real tmp1150 = d_lll.z.xy * tmp84;
	real tmp1153 = d_lll.z.xz * tmp145;
	real tmp1156 = d_lll.z.yy * tmp102;
	real tmp1159 = d_lll.z.yy * tmp108;
	real tmp1163 = d_lll.z.yz * tmp163;
	real tmp1166 = d_lll.z.zz * tmp175;
	real tmp1173 = d_lll.x.yz * tmp102;
	real tmp1178 = d_lll.x.zz * tmp163;
	real tmp1183 = d_lll.y.xy * tmp96;
	real tmp1186 = d_lll.y.xz * tmp102;
	real tmp1189 = d_lll.y.yy * tmp204;
	real tmp1192 = d_lll.y.yz * tmp210;
	real tmp1195 = d_lll.y.zz * tmp265;
	real tmp1200 = d_lll.y.zz * tmp216;
	real tmp1203 = d_lll.z.xy * tmp102;
	real tmp1206 = d_lll.z.xz * tmp163;
	real tmp1209 = d_lll.z.yy * tmp210;
	real tmp1212 = d_lll.z.yz * tmp265;
	real tmp1215 = d_lll.z.zz * tmp277;
	real tmp1228 = d_lll.x.zz * tmp175;
	real tmp1233 = d_lll.y.xy * tmp108;
	real tmp1239 = d_lll.y.yy * tmp210;
	real tmp1242 = d_lll.y.yz * tmp265;
	real tmp1245 = d_lll.y.zz * tmp277;
	real tmp1248 = d_lll.z.xy * tmp114;
	real tmp1251 = d_lll.z.xz * tmp175;
	real tmp1254 = d_lll.z.yy * tmp265;
	real tmp1259 = d_lll.z.yy * tmp216;
	real tmp1262 = d_lll.z.yz * tmp277;
	real tmp1265 = d_lll.z.zz * tmp318;
	real tmp1287 = d_lll.z.xx * tmp108;
	real tmp1288 = d_lll.x.yy * d_lll.x.yy;
	real tmp1303 = d_lll.z.xx * tmp163;
	real tmp1307 = d_lll.x.yz * d_lll.x.yz;
	real tmp1323 = d_lll.x.zz * d_lll.x.zz;
	real tmp1333 = d_lll.y.xz * tmp108;
	real tmp1354 = d_lll.z.xx * tmp84;
	real tmp1364 = d_lll.z.xy * tmp108;
	real tmp1393 = d_lll.z.xx * tmp114;
	real tmp1398 = d_lll.z.xy * tmp265;
	real tmp1401 = d_lll.z.xy * tmp216;
	real tmp1405 = d_lll.y.xz * d_lll.y.xz;
	real tmp1422 = d_lll.z.xy * tmp163;
	real tmp1447 = d_lll.z.xy * d_lll.z.xy;
	real tmp1451 = d_lll.y.xx * d_lll.y.xx;
	real tmp1455 = d_lll.z.xx * d_lll.z.xx;
	real tmp1459 = K_ll.xy * K_ll.xy;
	real tmp1462 = gamma_uu.zz * tmp1307;
	real tmp1467 = gamma_uu.zz * tmp1405;
	real tmp1472 = gamma_uu.zz * tmp1447;
	real tmp1476 = K_ll.xz * K_ll.xz;
	real tmp1587 = K_ll.yy * gamma_uu.xy;
	real tmp1589 = K_ll.yz * gamma_uu.xz;
	real tmp1606 = d_lll.x.yy * gamma_uu.xy;
	real tmp1609 = d_lll.x.yz * gamma_uu.xz;
	real tmp1615 = d_lll.y.xz * gamma_uu.xz;
	real tmp1618 = d_lll.z.xy * gamma_uu.xz;
	real tmp1622 = d_lll.x.yz * gamma_uu.yz;
	real tmp1628 = d_lll.y.xz * gamma_uu.yz;
	real tmp1631 = d_lll.z.xy * gamma_uu.yz;
	real tmp1646 = d_lll.x.yy * tmp29;
	real tmp1648 = d_lll.x.yz * tmp35;
	real tmp1650 = d_lll.y.xy * tmp29;
	real tmp1655 = d_lll.y.xz * tmp35;
	real tmp1668 = d_lll.y.zz * tmp53;
	real tmp1672 = d_lll.z.xy * tmp35;
	real tmp1679 = d_lll.x.yz * tmp47;
	real tmp1682 = d_lll.y.xy * tmp41;
	real tmp1687 = d_lll.y.xz * tmp47;
	real tmp1697 = d_lll.y.yz * tmp102;
	real tmp1707 = d_lll.z.xy * tmp47;
	real tmp1726 = d_lll.y.xy * tmp47;
	real tmp1731 = d_lll.y.xz * tmp53;
	real tmp1741 = d_lll.y.yz * tmp163;
	real tmp1746 = d_lll.y.zz * tmp175;
	real tmp1751 = d_lll.z.xy * tmp53;
	real tmp1756 = d_lll.z.yy * tmp163;
	real tmp1759 = d_lll.z.yy * tmp114;
	real tmp1773 = d_lll.x.zz * tmp114;
	real tmp1807 = d_lll.y.xx * tmp84;
	real tmp1823 = d_lll.y.yz * tmp216;
	real tmp1832 = d_lll.z.xx * tmp145;
	real tmp1909 = d_lll.z.xx * tmp35;
	real tmp1931 = d_lll.z.xx * tmp47;
	real tmp1950 = d_lll.z.xx * tmp53;
	real tmp1971 = d_lll.z.xy * tmp210;
	real tmp1992 = d_lll.z.xy * tmp145;
	real tmp2024 = gamma_uu.yy * tmp1288;
	real tmp2032 = gamma_uu.yz * tmp1447;
	real tmp2176 = K_ll.xz * gamma_uu.xx;
	real tmp2177 = K_ll.yz * gamma_uu.xy;
	real tmp2179 = K_ll.zz * gamma_uu.xz;
	real tmp2181 = K_ll.yz * gamma_uu.yy;
	real tmp2183 = K_ll.zz * gamma_uu.yz;
	real tmp2196 = d_lll.x.yz * gamma_uu.xy;
	real tmp2199 = d_lll.x.zz * gamma_uu.xz;
	real tmp2202 = d_lll.y.xz * gamma_uu.xy;
	real tmp2206 = d_lll.z.xy * gamma_uu.xy;
	real tmp2236 = d_lll.x.yz * tmp29;
	real tmp2238 = d_lll.x.zz * tmp35;
	real tmp2240 = d_lll.y.xz * tmp29;
	real tmp2244 = d_lll.z.xy * tmp29;
	real tmp2248 = d_lll.z.xz * tmp35;
	real tmp2253 = d_lll.z.yy * tmp41;
	real tmp2257 = d_lll.z.yz * tmp47;
	real tmp2262 = d_lll.z.zz * tmp53;
	real tmp2269 = d_lll.x.zz * tmp47;
	real tmp2272 = d_lll.y.xz * tmp41;
	real tmp2277 = d_lll.y.zz * tmp102;
	real tmp2282 = d_lll.y.zz * tmp108;
	real tmp2285 = d_lll.z.xy * tmp41;
	real tmp2290 = d_lll.z.xz * tmp47;
	real tmp2295 = d_lll.z.yy * tmp96;
	real tmp2300 = d_lll.z.yz * tmp108;
	real tmp2305 = d_lll.z.zz * tmp114;
	real tmp2313 = d_lll.x.zz * tmp53;
	real tmp2334 = d_lll.z.xz * tmp53;
	real tmp2344 = d_lll.z.yz * tmp114;
	real tmp2363 = d_lll.y.xz * tmp96;
	real tmp2367 = d_lll.z.xx * tmp72;
	real tmp2371 = d_lll.z.xx * tmp41;
	real tmp2373 = d_lll.z.xy * tmp96;
	real tmp2377 = d_lll.z.xz * tmp102;
	real tmp2386 = d_lll.z.yz * tmp210;
	real tmp2391 = d_lll.z.zz * tmp216;
	real tmp2432 = d_lll.z.xz * tmp114;
	real tmp2441 = d_lll.z.yz * tmp216;
	real tmp2479 = d_lll.y.xz * tmp72;
	real tmp2491 = d_lll.z.xx * tmp29;
	real tmp2493 = d_lll.z.xy * tmp72;
	real tmp2500 = d_lll.z.xz * tmp84;
	real tmp2508 = d_lll.z.yz * tmp102;
	real tmp2614 = gamma_uu.yz * tmp1405;
	real tmp2619 = gamma_uu.yy * tmp1307;
	real tmp2622 = gamma_uu.yy * tmp1447;
	real tmp2625 = gamma_uu.zz * tmp1323;
	real tmp2810 = d_lll.z.yy * gamma_uu.yz;
	real tmp2853 = d_lll.y.yy * tmp72;
	real tmp2858 = d_lll.y.yz * tmp84;
	real tmp2863 = d_lll.z.yy * tmp84;
	real tmp2884 = d_lll.z.yy * tmp145;
	real tmp2897 = d_lll.y.xx * tmp29;
	real tmp2935 = d_lll.y.yy * tmp108;
	real tmp2940 = d_lll.y.yz * tmp114;
	real tmp2986 = d_lll.z.yy * tmp175;
	real tmp3006 = d_lll.z.yy * tmp47;
	real tmp3117 = d_lll.y.zz * d_lll.y.zz;
	real tmp3160 = d_lll.z.yy * d_lll.z.yy;
	real tmp3162 = K_ll.yz * K_ll.yz;
	real tmp3340 = d_lll.y.zz * gamma_uu.yz;
	real tmp3377 = d_lll.y.zz * tmp84;
	real tmp3387 = d_lll.z.yy * tmp72;
	real tmp3417 = d_lll.x.zz * tmp84;
	real tmp3630 = d_lll.y.zz * tmp210;
	real tmp3748 = gamma_uu.zz * tmp3117;
	real tmp3966 = d_lll.y.zz * tmp72;
	real tmp3969 = d_lll.z.xz * tmp29;
	real tmp3974 = d_lll.z.yz * tmp72;
	real tmp3979 = d_lll.z.zz * tmp84;
	real tmp3995 = d_lll.z.yz * tmp84;
	real tmp4000 = d_lll.z.zz * tmp145;
	real tmp4012 = d_lll.y.zz * tmp96;
	real tmp4014 = d_lll.z.xz * tmp72;
	real tmp4022 = d_lll.z.yz * tmp96;
	real tmp4032 = d_lll.z.zz * tmp108;
	real tmp4063 = d_lll.z.zz * tmp163;
	real tmp4203 = d_lll.z.zz * tmp265;
	deriv->alpha += -f_alphaSq * (tmp6
		+ tmp5
		+ 2. * tmp3
		+ 2. * tmp2
		+ 2. * tmp1
		+ K_ll.xx * gamma_uu.xx);
	deriv->a_l.x += -alphaSq_dalpha_f * (
		K_ll.zz * a_l.x * gamma_uu.zz
		+ 2. * K_ll.yz * a_l.x * gamma_uu.yz
		+ K_ll.yy * a_l.x * gamma_uu.yy
		+ 2. * K_ll.xz * a_l.x * gamma_uu.xz
		+ 2. * K_ll.xy * a_l.x * gamma_uu.xy
		+ K_ll.xx * a_l.x * gamma_uu.xx
	)
	- f_alpha * (
		- 2. * K_ll.zz * d_lll.x.xx * tmp53
		- 4. * K_ll.zz * d_lll.x.xy * tmp114
		- 4. * K_ll.zz * d_lll.x.xz * tmp175
		- 2. * K_ll.zz * d_lll.x.yy * tmp216
		- 4. * K_ll.zz * d_lll.x.yz * tmp277
		- 2. * K_ll.zz * d_lll.x.zz * tmp318
		+ K_ll.zz * a_l.x * gamma_uu.zz
		- 4. * K_ll.yz * d_lll.x.xx * tmp47
		- 4. * K_ll.yz * d_lll.x.xy * tmp102
		- 4. * K_ll.yz * d_lll.x.xy * tmp108
		- 4. * K_ll.yz * d_lll.x.xz * tmp163
		- 4. * K_ll.yz * d_lll.x.xz * tmp114
		- 4. * K_ll.yz * d_lll.x.yy * tmp210
		- 4. * K_ll.yz * d_lll.x.yz * tmp265
		- 4. * K_ll.yz * d_lll.x.yz * tmp216
		- 4. * K_ll.yz * d_lll.x.zz * tmp277
		+ 2. * K_ll.yz * a_l.x * gamma_uu.yz
		- 2. * K_ll.yy * d_lll.x.xx * tmp41
		- 4. * K_ll.yy * d_lll.x.xy * tmp96
		- 4. * K_ll.yy * d_lll.x.xz * tmp102
		- 2. * K_ll.yy * d_lll.x.yy * tmp204
		- 4. * K_ll.yy * d_lll.x.yz * tmp210
		- 2. * K_ll.yy * d_lll.x.zz * tmp216
		+ K_ll.yy * a_l.x * gamma_uu.yy
		- 4. * K_ll.xz * d_lll.x.xx * tmp35
		- 4. * K_ll.xz * d_lll.x.xy * tmp84
		- 4. * K_ll.xz * d_lll.x.xy * tmp47
		- 4. * K_ll.xz * d_lll.x.xz * tmp145
		- 4. * K_ll.xz * d_lll.x.xz * tmp53
		- 4. * K_ll.xz * d_lll.x.yy * tmp102
		- 4. * K_ll.xz * d_lll.x.yz * tmp163
		- 4. * K_ll.xz * d_lll.x.yz * tmp114
		- 4. * K_ll.xz * d_lll.x.zz * tmp175
		+ 2. * K_ll.xz * a_l.x * gamma_uu.xz
		- 4. * K_ll.xy * d_lll.x.xx * tmp29
		- 4. * K_ll.xy * d_lll.x.xy * tmp72
		- 4. * K_ll.xy * d_lll.x.xy * tmp41
		- 4. * K_ll.xy * d_lll.x.xz * tmp84
		- 4. * K_ll.xy * d_lll.x.xz * tmp47
		- 4. * K_ll.xy * d_lll.x.yy * tmp96
		- 4. * K_ll.xy * d_lll.x.yz * tmp102
		- 4. * K_ll.xy * d_lll.x.yz * tmp108
		- 4. * K_ll.xy * d_lll.x.zz * tmp114
		+ 2. * K_ll.xy * a_l.x * gamma_uu.xy
		- 2. * K_ll.xx * d_lll.x.xx * tmp23
		- 4. * K_ll.xx * d_lll.x.xy * tmp29
		- 4. * K_ll.xx * d_lll.x.xz * tmp35
		- 2. * K_ll.xx * d_lll.x.yy * tmp41
		- 4. * K_ll.xx * d_lll.x.yz * tmp47
		- 2. * K_ll.xx * d_lll.x.zz * tmp53
		+ K_ll.xx * a_l.x * gamma_uu.xx
	);
	deriv->a_l.y += -alphaSq_dalpha_f * (
		K_ll.zz * a_l.y * gamma_uu.zz
		+ 2. * K_ll.yz * a_l.y * gamma_uu.yz
		+ K_ll.yy * a_l.y * gamma_uu.yy
		+ 2. * K_ll.xz * a_l.y * gamma_uu.xz
		+ 2. * K_ll.xy * a_l.y * gamma_uu.xy
		+ K_ll.xx * a_l.y * gamma_uu.xx
	)
	- f_alpha * (
		- 2. * K_ll.zz * d_lll.y.xx * tmp53
		- 4. * K_ll.zz * d_lll.y.xy * tmp114
		- 4. * K_ll.zz * d_lll.y.xz * tmp175
		- 2. * K_ll.zz * d_lll.y.yy * tmp216
		- 4. * K_ll.zz * d_lll.y.yz * tmp277
		- 2. * K_ll.zz * d_lll.y.zz * tmp318
		+ K_ll.zz * a_l.y * gamma_uu.zz
		- 4. * K_ll.yz * d_lll.y.xx * tmp47
		- 4. * K_ll.yz * d_lll.y.xy * tmp102
		- 4. * K_ll.yz * d_lll.y.xy * tmp108
		- 4. * K_ll.yz * d_lll.y.xz * tmp163
		- 4. * K_ll.yz * d_lll.y.xz * tmp114
		- 4. * K_ll.yz * d_lll.y.yy * tmp210
		- 4. * K_ll.yz * d_lll.y.yz * tmp265
		- 4. * K_ll.yz * d_lll.y.yz * tmp216
		- 4. * K_ll.yz * d_lll.y.zz * tmp277
		+ 2. * K_ll.yz * a_l.y * gamma_uu.yz
		- 2. * K_ll.yy * d_lll.y.xx * tmp41
		- 4. * K_ll.yy * d_lll.y.xy * tmp96
		- 4. * K_ll.yy * d_lll.y.xz * tmp102
		- 2. * K_ll.yy * d_lll.y.yy * tmp204
		- 4. * K_ll.yy * d_lll.y.yz * tmp210
		- 2. * K_ll.yy * d_lll.y.zz * tmp216
		+ K_ll.yy * a_l.y * gamma_uu.yy
		- 4. * K_ll.xz * d_lll.y.xx * tmp35
		- 4. * K_ll.xz * d_lll.y.xy * tmp84
		- 4. * K_ll.xz * d_lll.y.xy * tmp47
		- 4. * K_ll.xz * d_lll.y.xz * tmp145
		- 4. * K_ll.xz * d_lll.y.xz * tmp53
		- 4. * K_ll.xz * d_lll.y.yy * tmp102
		- 4. * K_ll.xz * d_lll.y.yz * tmp163
		- 4. * K_ll.xz * d_lll.y.yz * tmp114
		- 4. * K_ll.xz * d_lll.y.zz * tmp175
		+ 2. * K_ll.xz * a_l.y * gamma_uu.xz
		- 4. * K_ll.xy * d_lll.y.xx * tmp29
		- 4. * K_ll.xy * d_lll.y.xy * tmp72
		- 4. * K_ll.xy * d_lll.y.xy * tmp41
		- 4. * K_ll.xy * d_lll.y.xz * tmp84
		- 4. * K_ll.xy * d_lll.y.xz * tmp47
		- 4. * K_ll.xy * d_lll.y.yy * tmp96
		- 4. * K_ll.xy * d_lll.y.yz * tmp102
		- 4. * K_ll.xy * d_lll.y.yz * tmp108
		- 4. * K_ll.xy * d_lll.y.zz * tmp114
		+ 2. * K_ll.xy * a_l.y * gamma_uu.xy
		- 2. * K_ll.xx * d_lll.y.xx * tmp23
		- 4. * K_ll.xx * d_lll.y.xy * tmp29
		- 4. * K_ll.xx * d_lll.y.xz * tmp35
		- 2. * K_ll.xx * d_lll.y.yy * tmp41
		- 4. * K_ll.xx * d_lll.y.yz * tmp47
		- 2. * K_ll.xx * d_lll.y.zz * tmp53
		+ K_ll.xx * a_l.y * gamma_uu.xx
	);
	deriv->a_l.z += -alphaSq_dalpha_f * (
		K_ll.zz * a_l.z * gamma_uu.zz
		+ 2. * K_ll.yz * a_l.z * gamma_uu.yz
		+ K_ll.yy * a_l.z * gamma_uu.yy
		+ 2. * K_ll.xz * a_l.z * gamma_uu.xz
		+ 2. * K_ll.xy * a_l.z * gamma_uu.xy
		+ K_ll.xx * a_l.z * gamma_uu.xx
	)
	- f_alpha * (
		- 2. * K_ll.zz * d_lll.z.xx * tmp53
		- 4. * K_ll.zz * d_lll.z.xy * tmp114
		- 4. * K_ll.zz * d_lll.z.xz * tmp175
		- 2. * K_ll.zz * d_lll.z.yy * tmp216
		- 4. * K_ll.zz * d_lll.z.yz * tmp277
		- 2. * K_ll.zz * d_lll.z.zz * tmp318
		+ K_ll.zz * a_l.z * gamma_uu.zz
		- 4. * K_ll.yz * d_lll.z.xx * tmp47
		- 4. * K_ll.yz * d_lll.z.xy * tmp102
		- 4. * K_ll.yz * d_lll.z.xy * tmp108
		- 4. * K_ll.yz * d_lll.z.xz * tmp163
		- 4. * K_ll.yz * d_lll.z.xz * tmp114
		- 4. * K_ll.yz * d_lll.z.yy * tmp210
		- 4. * K_ll.yz * d_lll.z.yz * tmp265
		- 4. * K_ll.yz * d_lll.z.yz * tmp216
		- 4. * K_ll.yz * d_lll.z.zz * tmp277
		+ 2. * K_ll.yz * a_l.z * gamma_uu.yz
		- 2. * K_ll.yy * d_lll.z.xx * tmp41
		- 4. * K_ll.yy * d_lll.z.xy * tmp96
		- 4. * K_ll.yy * d_lll.z.xz * tmp102
		- 2. * K_ll.yy * d_lll.z.yy * tmp204
		- 4. * K_ll.yy * d_lll.z.yz * tmp210
		- 2. * K_ll.yy * d_lll.z.zz * tmp216
		+ K_ll.yy * a_l.z * gamma_uu.yy
		- 4. * K_ll.xz * d_lll.z.xx * tmp35
		- 4. * K_ll.xz * d_lll.z.xy * tmp84
		- 4. * K_ll.xz * d_lll.z.xy * tmp47
		- 4. * K_ll.xz * d_lll.z.xz * tmp145
		- 4. * K_ll.xz * d_lll.z.xz * tmp53
		- 4. * K_ll.xz * d_lll.z.yy * tmp102
		- 4. * K_ll.xz * d_lll.z.yz * tmp163
		- 4. * K_ll.xz * d_lll.z.yz * tmp114
		- 4. * K_ll.xz * d_lll.z.zz * tmp175
		+ 2. * K_ll.xz * a_l.z * gamma_uu.xz
		- 4. * K_ll.xy * d_lll.z.xx * tmp29
		- 4. * K_ll.xy * d_lll.z.xy * tmp72
		- 4. * K_ll.xy * d_lll.z.xy * tmp41
		- 4. * K_ll.xy * d_lll.z.xz * tmp84
		- 4. * K_ll.xy * d_lll.z.xz * tmp47
		- 4. * K_ll.xy * d_lll.z.yy * tmp96
		- 4. * K_ll.xy * d_lll.z.yz * tmp102
		- 4. * K_ll.xy * d_lll.z.yz * tmp108
		- 4. * K_ll.xy * d_lll.z.zz * tmp114
		+ 2. * K_ll.xy * a_l.z * gamma_uu.xy
		- 2. * K_ll.xx * d_lll.z.xx * tmp23
		- 4. * K_ll.xx * d_lll.z.xy * tmp29
		- 4. * K_ll.xx * d_lll.z.xz * tmp35
		- 2. * K_ll.xx * d_lll.z.yy * tmp41
		- 4. * K_ll.xx * d_lll.z.yz * tmp47
		- 2. * K_ll.xx * d_lll.z.zz * tmp53
		+ K_ll.xx * a_l.z * gamma_uu.xx);
	deriv->gamma_ll.xx += -2. * K_ll.xx * alpha;
	deriv->gamma_ll.xy += -2. * K_ll.xy * alpha;
	deriv->gamma_ll.xz += -2. * K_ll.xz * alpha;
	deriv->gamma_ll.yy += -2. * K_ll.yy * alpha;
	deriv->gamma_ll.yz += -2. * K_ll.yz * alpha;
	deriv->gamma_ll.zz += -2. * K_ll.zz * alpha;
	deriv->d_lll.x.xx += -K_ll.xx * tmp1020;
	deriv->d_lll.x.xy += -K_ll.xy * tmp1020;
	deriv->d_lll.x.xz += -K_ll.xz * tmp1020;
	deriv->d_lll.x.yy += -K_ll.yy * tmp1020;
	deriv->d_lll.x.yz += -K_ll.yz * tmp1020;
	deriv->d_lll.x.zz += -K_ll.zz * tmp1020;
	deriv->d_lll.y.xx += -K_ll.xx * tmp1032;
	deriv->d_lll.y.xy += -K_ll.xy * tmp1032;
	deriv->d_lll.y.xz += -K_ll.xz * tmp1032;
	deriv->d_lll.y.yy += -K_ll.yy * tmp1032;
	deriv->d_lll.y.yz += -K_ll.yz * tmp1032;
	deriv->d_lll.y.zz += -K_ll.zz * tmp1032;
	deriv->d_lll.z.xx += -K_ll.xx * tmp1044;
	deriv->d_lll.z.xy += -K_ll.xy * tmp1044;
	deriv->d_lll.z.xz += -K_ll.xz * tmp1044;
	deriv->d_lll.z.yy += -K_ll.yy * tmp1044;
	deriv->d_lll.z.yz += -K_ll.yz * tmp1044;
	deriv->d_lll.z.zz += -K_ll.zz * tmp1044;
	deriv->K_ll.xx += -alpha * (2. * gamma_uu.zz * tmp1476
		+ 2. * gamma_uu.yy * tmp1459
		- 2. * gamma_uu.yy * tmp1462
		- 2. * gamma_uu.yy * tmp1467
		- 2. * gamma_uu.yy * tmp1472
		+ gamma_uu.xx * K_ll.xx * K_ll.xx
		- gamma_uu.xx * gamma_uu.yy * tmp1451
		- gamma_uu.xx * gamma_uu.zz * tmp1455
		+ 2. * tmp1447 * tmp216
		+ d_lll.z.xx * tmp1254
		- 2. * d_lll.z.xx * tmp1259
		- 2. * d_lll.z.xx * tmp1262
		- d_lll.z.xx * tmp1265
		+ 2. * tmp1405 * tmp216
		- d_lll.y.yy * d_lll.z.xx * tmp210
		- 2. * d_lll.y.yz * d_lll.z.xx * tmp265
		- d_lll.y.zz * d_lll.z.xx * tmp277
		- 2. * d_lll.z.xx * tmp1422
		- 2. * d_lll.z.xx * tmp1251
		+ 4. * d_lll.y.xz * tmp1398
		- 4. * d_lll.y.xz * tmp1401
		+ 2. * d_lll.y.xz * tmp1303
		- 4. * d_lll.y.xz * tmp1393
		+ 2. * d_lll.y.xx * tmp1364
		- 2. * d_lll.y.xx * tmp1206
		- d_lll.y.xx * tmp1209
		- 2. * d_lll.y.xx * tmp1212
		- d_lll.y.xx * tmp1215
		- 2. * d_lll.y.xy * tmp1287
		+ d_lll.y.xx * tmp1195
		- 2. * d_lll.y.xx * tmp1200
		- 2. * d_lll.y.xx * tmp1354
		- 4. * d_lll.y.xx * tmp1203
		+ d_lll.x.zz * d_lll.y.xx * tmp163
		- 2. * d_lll.x.zz * d_lll.y.xx * tmp114
		- d_lll.x.zz * d_lll.z.xx * tmp175
		- tmp1323 * tmp318
		- 2. * d_lll.y.xx * tmp1183
		- 2. * d_lll.y.xx * tmp1333
		- d_lll.y.xx * tmp1189
		- 2. * d_lll.y.xx * tmp1192
		+ d_lll.x.yy * tmp1287
		- tmp1288 * tmp204
		- 4. * d_lll.x.yz * d_lll.x.zz * tmp277
		- 2. * d_lll.x.yz * d_lll.y.xx * tmp108
		- 2. * d_lll.x.yz * tmp1303
		- 2. * tmp1307 * tmp216
		+ 2. * d_lll.x.xz * tmp1265
		- 4. * d_lll.x.yy * d_lll.x.yz * tmp210
		- 2. * d_lll.x.yy * d_lll.x.zz * tmp216
		- d_lll.x.yy * d_lll.y.xx * tmp96
		- 2. * d_lll.x.yy * d_lll.z.xx * tmp102
		+ 4. * d_lll.x.xz * tmp1262
		+ 4. * d_lll.x.xz * tmp1259
		+ 4. * d_lll.x.xz * tmp1251
		- 2. * d_lll.x.xz * tmp1254
		+ 4. * d_lll.x.xz * tmp1248
		+ 2. * d_lll.x.xz * tmp1245
		+ 4. * d_lll.x.xz * tmp1242
		+ 2. * d_lll.x.xz * tmp1239
		+ 4. * d_lll.x.xz * d_lll.y.xz * tmp114
		+ 4. * d_lll.x.xz * tmp1233
		+ 2. * d_lll.x.xy * tmp1215
		- 2. * d_lll.x.xz * d_lll.x.yy * tmp108
		- 4. * d_lll.x.xz * d_lll.x.yz * tmp114
		- 2. * d_lll.x.xz * tmp1228
		+ 4. * d_lll.x.xy * tmp1212
		+ 2. * d_lll.x.xy * tmp1209
		+ 4. * d_lll.x.xy * tmp1206
		+ 4. * d_lll.x.xy * tmp1203
		+ 4. * d_lll.x.xy * tmp1200
		+ 4. * d_lll.x.xy * tmp1192
		- 2. * d_lll.x.xy * tmp1195
		+ 2. * d_lll.x.xy * tmp1189
		+ 4. * d_lll.x.xy * tmp1186
		+ 4. * d_lll.x.xy * tmp1183
		+ d_lll.x.xx * tmp1166
		- 2. * d_lll.x.xy * d_lll.x.yy * tmp96
		- 4. * d_lll.x.xy * tmp1173
		- 2. * d_lll.x.xy * tmp1178
		+ 2. * d_lll.x.xx * tmp1163
		+ 2. * d_lll.x.xx * tmp1156
		- d_lll.x.xx * tmp1159
		+ 2. * d_lll.x.xx * tmp1153
		+ 2. * d_lll.x.xx * tmp1150
		+ 2. * d_lll.x.xx * tmp1147
		+ 2. * d_lll.x.xx * tmp1140
		- d_lll.x.xx * tmp1143
		+ d_lll.x.xx * tmp1138
		+ 2. * d_lll.x.xx * tmp1135
		+ 2. * d_lll.x.xx * tmp1132
		+ a_l.x * a_l.x
		- d_lll.x.xx * d_lll.x.yy * tmp72
		- 2. * d_lll.x.xx * tmp1123
		- d_lll.x.xx * tmp1128
		+ a_l.z * d_lll.z.xx * gamma_uu.zz
		+ a_l.z * d_lll.y.xx * gamma_uu.yz
		+ a_l.y * d_lll.z.xx * gamma_uu.yz
		- a_l.z * d_lll.x.xx * gamma_uu.xz
		- 2. * a_l.z * d_lll.x.xy * gamma_uu.yz
		- 2. * a_l.z * d_lll.x.xz * gamma_uu.zz
		+ a_l.y * d_lll.y.xx * gamma_uu.yy
		+ a_l.x * tmp1091
		- a_l.y * d_lll.x.xx * gamma_uu.xy
		- 2. * a_l.y * d_lll.x.xy * gamma_uu.yy
		- 2. * a_l.y * d_lll.x.xz * gamma_uu.yz
		+ a_l.x * tmp1090
		+ 4. * M_PI * gamma_ll.xx * rho
		- a_l.x * d_lll.x.xx * gamma_uu.xx
		- 2. * a_l.x * d_lll.x.xy * gamma_uu.xy
		- 2. * a_l.x * d_lll.x.xz * gamma_uu.xz
		+ 8. * M_PI * S_ll.xx
		+ 4. * K_ll.xy * K_ll.xz * gamma_uu.yz
		- 4. * M_PI * S * gamma_ll.xx
		+ 2. * K_ll.xx * tmp2
		- K_ll.xx * tmp6
		- 2. * K_ll.xx * tmp3
		- K_ll.xx * tmp5
		+ 2. * K_ll.xx * tmp1);
	deriv->K_ll.xy += -alpha * (2. * gamma_uu.xy * tmp1467
		- 2. * gamma_uu.xz * tmp2032
		+ 2. * gamma_uu.xy * tmp1462
		+ gamma_uu.xy * tmp2024
		+ gamma_uu.xx * gamma_uu.xy * tmp1451
		+ 2. * d_lll.z.xx * tmp1759
		- 2. * d_lll.z.xy * tmp1251
		- d_lll.z.xy * tmp1254
		- 2. * d_lll.z.xy * tmp1262
		- d_lll.z.xy * tmp1265
		+ 2. * d_lll.y.yz * tmp1303
		- 2. * d_lll.y.yz * tmp1393
		- 2. * d_lll.y.yz * tmp1401
		- d_lll.y.zz * d_lll.z.xy * tmp277
		- d_lll.z.xx * tmp1992
		- 2. * d_lll.z.xx * tmp1756
		+ d_lll.y.xz * tmp1265
		- d_lll.y.yy * tmp1971
		+ 2. * d_lll.y.xz * tmp1262
		+ d_lll.y.xz * tmp1254
		+ 2. * d_lll.y.xz * tmp1251
		+ 2. * d_lll.y.xz * tmp1248
		+ 2. * d_lll.y.xz * tmp1950
		- 2. * d_lll.y.xz * tmp1422
		+ d_lll.y.xz * tmp1245
		- d_lll.y.xz * tmp1832
		+ 2. * d_lll.y.xz * tmp1823
		+ d_lll.y.xz * tmp1239
		+ 2. * d_lll.y.xy * tmp1931
		- 2. * d_lll.y.xy * tmp1203
		+ 2. * d_lll.y.xy * tmp1186
		- 2. * d_lll.y.xy * tmp1354
		+ d_lll.y.xx * tmp1166
		+ 2. * d_lll.y.xx * tmp1163
		+ d_lll.y.xx * tmp1159
		+ 2. * d_lll.y.xx * tmp1153
		+ d_lll.y.xx * tmp1150
		+ d_lll.y.xx * tmp1909
		+ 2. * d_lll.y.xx * tmp1147
		+ 2. * d_lll.y.xx * tmp1697
		- d_lll.y.xx * tmp1143
		+ d_lll.y.xx * tmp1138
		+ 2. * d_lll.y.xx * tmp1687
		+ d_lll.y.xx * tmp1135
		+ 2. * d_lll.y.xx * tmp1682
		+ d_lll.x.zz * d_lll.y.xx * tmp53
		- 2. * d_lll.x.zz * d_lll.y.xy * tmp114
		- d_lll.x.zz * d_lll.y.xz * tmp175
		- d_lll.x.zz * d_lll.y.yy * tmp216
		- 2. * d_lll.x.zz * d_lll.y.yz * tmp277
		- d_lll.x.zz * d_lll.y.zz * tmp318
		- d_lll.x.zz * d_lll.z.xy * tmp175
		+ d_lll.x.yz * tmp1265
		- d_lll.x.zz * d_lll.y.xx * tmp145
		+ 2. * d_lll.x.yz * tmp1262
		+ 2. * d_lll.x.yz * tmp1259
		+ 2. * d_lll.x.yz * tmp1251
		- d_lll.x.yz * tmp1254
		+ 2. * d_lll.x.yz * tmp1248
		+ d_lll.x.yz * tmp1832
		- 2. * d_lll.x.yz * tmp1422
		+ 2. * d_lll.x.yz * d_lll.y.xx * tmp47
		- 2. * d_lll.x.yz * d_lll.y.xy * tmp102
		- d_lll.x.yz * tmp1239
		- 2. * d_lll.x.yz * tmp1823
		- d_lll.x.yz * tmp1245
		+ d_lll.x.yz * tmp1228
		- d_lll.x.yz * tmp1807
		+ d_lll.x.yy * tmp1215
		+ 2. * d_lll.x.yy * tmp1212
		+ d_lll.x.yy * tmp1209
		+ 2. * d_lll.x.yy * tmp1206
		+ d_lll.x.yy * tmp1364
		+ d_lll.x.yy * tmp1354
		+ d_lll.x.yy * tmp1200
		+ 2. * d_lll.x.yy * tmp1186
		- d_lll.x.yy * tmp1333
		- d_lll.x.yy * tmp1195
		+ d_lll.x.yy * d_lll.y.xx * tmp41
		+ 2. * d_lll.x.yy * tmp1773
		+ d_lll.x.yy * d_lll.x.yz * tmp108
		- d_lll.x.yy * tmp1178
		+ 2. * d_lll.x.yy * tmp1173
		+ 2. * d_lll.x.xz * tmp1756
		- 2. * d_lll.x.xz * tmp1759
		+ 2. * d_lll.x.xz * d_lll.x.yz * tmp53
		- 4. * d_lll.x.xz * tmp1726
		- 2. * d_lll.x.xz * tmp1731
		- 2. * d_lll.x.xz * d_lll.y.yy * tmp102
		- 4. * d_lll.x.xz * tmp1741
		- 2. * d_lll.x.xz * tmp1746
		- 2. * d_lll.x.xz * tmp1751
		+ 2. * d_lll.x.xz * d_lll.x.yy * tmp47
		+ 2. * d_lll.x.xy * tmp1156
		- 2. * d_lll.x.xy * tmp1159
		+ 2. * d_lll.x.xy * tmp1679
		- 4. * d_lll.x.xy * tmp1682
		- 2. * d_lll.x.xy * tmp1687
		- 2. * d_lll.x.xy * tmp1138
		- 4. * d_lll.x.xy * tmp1697
		- 2. * d_lll.x.xy * tmp1147
		- 2. * d_lll.x.xy * tmp1707
		+ 2. * d_lll.x.xy * d_lll.x.yy * tmp41
		+ d_lll.x.xx * tmp1648
		- 2. * d_lll.x.xx * tmp1650
		- d_lll.x.xx * tmp1655
		- d_lll.x.xx * d_lll.y.yy * tmp41
		- 2. * d_lll.x.xx * d_lll.y.yz * tmp47
		- d_lll.x.xx * tmp1668
		- d_lll.x.xx * tmp1672
		+ d_lll.x.xx * tmp1646
		+ a_l.z * d_lll.z.xy * gamma_uu.zz
		+ a_l.y * tmp1631
		- a_l.z * d_lll.x.yy * gamma_uu.yz
		- a_l.z * d_lll.x.yz * gamma_uu.zz
		- a_l.z * d_lll.y.xx * gamma_uu.xz
		- a_l.z * d_lll.y.xz * gamma_uu.zz
		+ a_l.x * tmp1618
		- a_l.y * d_lll.x.yy * gamma_uu.yy
		- a_l.y * tmp1622
		- a_l.y * tmp1090
		- a_l.y * tmp1628
		+ a_l.x * a_l.y
		- a_l.x * tmp1606
		- a_l.x * tmp1609
		- a_l.x * d_lll.y.xx * gamma_uu.xx
		- a_l.x * tmp1615
		+ 4. * M_PI * gamma_ll.xy * rho
		+ 8. * M_PI * S_ll.xy
		+ 2. * K_ll.xz * K_ll.yz * gamma_uu.zz
		- 4. * M_PI * S * gamma_ll.xy
		+ 2. * K_ll.xz * K_ll.yy * gamma_uu.yz
		+ K_ll.xy * tmp6
		- K_ll.xy * tmp5
		+ 2. * K_ll.xx * tmp1589
		+ 2. * K_ll.xx * tmp1587
		+ K_ll.xx * K_ll.xy * gamma_uu.xx);
	deriv->K_ll.xz += -alpha * (gamma_uu.xz * tmp2625
		+ 2. * gamma_uu.xz * tmp2622
		+ 2. * gamma_uu.xz * tmp2619
		+ gamma_uu.xx * gamma_uu.xz * tmp1455
		- 2. * gamma_uu.xy * tmp2614
		+ d_lll.z.xy * tmp1215
		+ 2. * d_lll.z.xy * tmp2441
		+ d_lll.z.xy * tmp1209
		+ 2. * d_lll.z.xy * tmp2432
		+ d_lll.z.xx * tmp1166
		+ 2. * d_lll.z.xx * tmp2344
		+ 2. * d_lll.z.xx * tmp1156
		- d_lll.z.xx * tmp1159
		+ 2. * d_lll.z.xx * tmp2334
		+ 2. * d_lll.z.xx * tmp1707
		+ d_lll.z.xx * tmp1150
		+ d_lll.y.zz * tmp1398
		+ d_lll.y.zz * tmp1303
		+ 2. * d_lll.y.yz * tmp1971
		+ 2. * d_lll.y.yz * tmp1287
		+ d_lll.y.yy * d_lll.z.xy * tmp204
		+ d_lll.y.yy * d_lll.z.xx * tmp96
		+ 2. * d_lll.y.xz * tmp1203
		- 2. * d_lll.y.xz * tmp1364
		- 2. * d_lll.y.xz * tmp2432
		- d_lll.y.xz * tmp1209
		- 2. * d_lll.y.xz * tmp2441
		- d_lll.y.xz * tmp1215
		+ d_lll.y.xz * tmp1354
		+ 2. * d_lll.y.xy * tmp2373
		- d_lll.y.xz * tmp1189
		- 2. * d_lll.y.xz * tmp1192
		- d_lll.y.xz * tmp1195
		+ 2. * d_lll.y.xy * tmp2367
		+ 2. * d_lll.y.xx * tmp2300
		- 2. * d_lll.y.xy * tmp2363
		+ 2. * d_lll.y.xx * tmp2290
		- 2. * d_lll.y.xx * tmp2508
		+ 2. * d_lll.y.xx * tmp2285
		- 2. * d_lll.y.xx * tmp2500
		+ d_lll.y.xx * tmp2491
		- d_lll.y.xx * tmp2493
		+ 2. * d_lll.y.xx * tmp2277
		- 2. * d_lll.y.xx * tmp2282
		+ d_lll.x.zz * tmp1259
		- d_lll.y.xx * tmp2479
		+ 2. * d_lll.x.zz * tmp1248
		- d_lll.x.zz * tmp1254
		+ d_lll.x.zz * tmp1950
		- d_lll.x.zz * tmp1422
		+ d_lll.x.zz * tmp1245
		+ 2. * d_lll.x.zz * tmp1242
		+ d_lll.x.zz * tmp1239
		+ d_lll.x.zz * d_lll.y.xz * tmp163
		+ 2. * d_lll.x.zz * tmp1233
		+ d_lll.x.zz * tmp1807
		+ 2. * d_lll.x.yz * tmp1931
		- 2. * d_lll.x.yz * tmp2432
		- d_lll.x.yz * tmp1209
		- 2. * d_lll.x.yz * tmp2441
		- d_lll.x.yz * tmp1215
		+ 2. * d_lll.x.yz * tmp1200
		- d_lll.x.yz * tmp1354
		+ 2. * d_lll.x.yz * tmp1192
		- d_lll.x.yz * tmp1195
		+ d_lll.x.yz * tmp1189
		+ 2. * d_lll.x.yz * tmp1186
		- 2. * d_lll.x.yz * tmp1333
		+ 2. * d_lll.x.yz * tmp1183
		+ d_lll.x.yz * d_lll.y.xx * tmp72
		+ 2. * d_lll.x.yz * tmp1773
		+ d_lll.x.yz * tmp1178
		+ d_lll.x.yy * tmp2371
		- d_lll.x.yy * tmp2373
		- 2. * d_lll.x.yy * tmp2377
		- d_lll.x.yy * d_lll.z.yy * tmp204
		- 2. * d_lll.x.yy * tmp2386
		- d_lll.x.yy * tmp2391
		+ 2. * d_lll.x.yy * d_lll.x.zz * tmp102
		- d_lll.x.yy * d_lll.x.zz * tmp108
		- d_lll.x.yy * tmp2363
		- d_lll.x.yy * tmp2367
		+ d_lll.x.yy * d_lll.x.yz * tmp96
		+ 2. * d_lll.x.xz * tmp1147
		- 2. * d_lll.x.xz * tmp1707
		- 4. * d_lll.x.xz * tmp2334
		- 2. * d_lll.x.xz * tmp1156
		- 4. * d_lll.x.xz * tmp2344
		- 2. * d_lll.x.xz * tmp1166
		+ 2. * d_lll.x.xz * tmp2313
		- 2. * d_lll.x.xz * tmp1687
		- 2. * d_lll.x.xz * tmp1143
		+ 2. * d_lll.x.xz * tmp1679
		+ 2. * d_lll.x.xy * tmp2282
		- 2. * d_lll.x.xy * tmp2285
		- 4. * d_lll.x.xy * tmp2290
		- 2. * d_lll.x.xy * tmp2295
		- 4. * d_lll.x.xy * tmp2300
		- 2. * d_lll.x.xy * tmp2305
		+ 2. * d_lll.x.xy * tmp2269
		- 2. * d_lll.x.xy * tmp2272
		- 2. * d_lll.x.xy * tmp2277
		+ 2. * d_lll.x.xy * d_lll.x.yz * tmp41
		+ d_lll.x.xx * tmp2238
		- d_lll.x.xx * tmp2240
		- d_lll.x.xx * tmp2244
		- 2. * d_lll.x.xx * tmp2248
		- d_lll.x.xx * tmp2253
		- 2. * d_lll.x.xx * tmp2257
		- d_lll.x.xx * tmp2262
		+ d_lll.x.xx * tmp2236
		+ a_l.z * tmp1628
		- a_l.z * tmp1091
		- a_l.z * tmp1631
		+ a_l.y * d_lll.y.xz * gamma_uu.yy
		- a_l.y * d_lll.z.xx * gamma_uu.xy
		- a_l.y * d_lll.z.xy * gamma_uu.yy
		- a_l.z * tmp1622
		- a_l.z * d_lll.x.zz * gamma_uu.zz
		+ a_l.x * tmp2202
		- a_l.x * d_lll.z.xx * gamma_uu.xx
		- a_l.x * tmp2206
		- a_l.y * d_lll.x.yz * gamma_uu.yy
		- a_l.y * d_lll.x.zz * gamma_uu.yz
		+ a_l.x * a_l.z
		- a_l.x * tmp2196
		- a_l.x * tmp2199
		+ 4. * M_PI * gamma_ll.xz * rho
		+ 8. * M_PI * S_ll.xz
		+ K_ll.xz * tmp5
		- 4. * M_PI * S * gamma_ll.xz
		+ 2. * K_ll.xy * tmp2183
		- K_ll.xz * tmp6
		+ 2. * K_ll.xy * tmp2181
		+ 2. * K_ll.xx * tmp2179
		+ 2. * K_ll.xx * tmp2177
		+ K_ll.xx * tmp2176);
	deriv->K_ll.yy += alpha * (gamma_uu.yy * gamma_uu.zz * tmp3160
		- 2. * gamma_uu.zz * tmp3162
		+ 2. * gamma_uu.xx * tmp1472
		- gamma_uu.yy * K_ll.yy * K_ll.yy
		+ 2. * gamma_uu.xx * tmp1467
		+ 2. * gamma_uu.xx * tmp1462
		+ gamma_uu.xx * tmp2024
		+ d_lll.z.yy * tmp1265
		- 2. * gamma_uu.xx * tmp1459
		+ 2. * d_lll.z.yy * tmp1262
		+ 2. * d_lll.z.xz * tmp2986
		+ 2. * d_lll.z.xy * tmp1756
		- 2. * tmp1447 * tmp53
		+ 2. * d_lll.z.xx * d_lll.z.yy * tmp53
		+ tmp3117 * tmp318
		- d_lll.z.xx * tmp2884
		+ d_lll.y.zz * d_lll.z.yy * tmp277
		+ 2. * d_lll.y.yz * tmp1832
		- 4. * d_lll.y.yz * tmp1950
		- 4. * d_lll.y.yz * tmp1248
		- 4. * d_lll.y.yz * tmp1251
		- 4. * d_lll.y.yz * tmp1262
		- 2. * d_lll.y.yz * tmp1265
		+ 2. * d_lll.y.yz * tmp1245
		+ d_lll.y.yy * tmp1354
		- 2. * d_lll.y.yy * tmp1931
		- 2. * d_lll.y.yy * tmp1364
		- 2. * d_lll.y.yy * tmp1206
		- 2. * d_lll.y.yy * tmp1212
		- d_lll.y.yy * tmp1215
		+ d_lll.y.yy * tmp1195
		+ 2. * tmp1405 * tmp53
		+ 2. * d_lll.y.xz * tmp1756
		+ 4. * d_lll.y.xz * tmp1746
		+ 4. * d_lll.y.xz * tmp2940
		+ 2. * d_lll.y.xz * tmp2935
		+ 2. * d_lll.y.xy * tmp1143
		- 2. * d_lll.y.xy * tmp1909
		- 4. * d_lll.y.xy * tmp1707
		- 4. * d_lll.y.xy * tmp1153
		- 4. * d_lll.y.xy * tmp1163
		- 2. * d_lll.y.xy * tmp1166
		+ 4. * d_lll.y.xy * tmp1687
		+ tmp1451 * tmp23
		+ 2. * d_lll.y.xx * tmp3006
		+ 2. * d_lll.y.xx * tmp1668
		- d_lll.y.xx * tmp2863
		+ 2. * d_lll.y.xx * tmp2858
		+ d_lll.y.xx * tmp2853
		+ 4. * d_lll.y.xx * tmp1655
		+ 2. * d_lll.y.xx * tmp1650
		+ d_lll.x.zz * tmp2986
		+ d_lll.x.zz * d_lll.y.yy * tmp163
		- 2. * d_lll.x.zz * d_lll.y.yy * tmp114
		- 2. * d_lll.x.zz * d_lll.y.yz * tmp175
		+ 2. * d_lll.x.zz * d_lll.y.xy * tmp145
		- 4. * d_lll.x.zz * d_lll.y.xy * tmp53
		+ 4. * d_lll.x.yz * tmp1759
		- 2. * tmp1307 * tmp53
		+ 4. * d_lll.x.yz * tmp1751
		- 2. * d_lll.x.yz * tmp1756
		+ d_lll.x.yy * tmp1166
		- 4. * d_lll.x.yz * tmp1726
		- 2. * d_lll.x.yz * tmp2935
		- 4. * d_lll.x.yz * tmp2940
		- 4. * d_lll.x.yz * tmp1992
		+ 2. * d_lll.x.yy * tmp1163
		+ 2. * d_lll.x.yy * tmp1159
		+ 2. * d_lll.x.yy * tmp1153
		+ 4. * d_lll.x.yy * tmp1707
		+ d_lll.x.yy * tmp1909
		- 2. * d_lll.x.yy * tmp1150
		+ 2. * d_lll.x.yy * tmp1147
		+ 2. * d_lll.x.yy * tmp1135
		- d_lll.x.yy * tmp1143
		+ d_lll.x.yy * tmp2897
		+ 2. * d_lll.x.yy * tmp2313
		+ 2. * d_lll.x.yy * tmp1123
		- d_lll.x.yy * tmp1128
		+ 2. * d_lll.x.xz * tmp2884
		+ 2. * d_lll.x.xz * d_lll.x.yy * tmp35
		- 4. * d_lll.x.xz * d_lll.y.xy * tmp35
		- 2. * d_lll.x.xz * d_lll.y.yy * tmp84
		- 4. * d_lll.x.xz * d_lll.y.yz * tmp145
		+ 2. * d_lll.x.xy * tmp2863
		+ 2. * d_lll.x.xy * tmp1646
		- 4. * d_lll.x.xy * tmp1650
		- 2. * d_lll.x.xy * tmp2853
		- 4. * d_lll.x.xy * tmp2858
		+ d_lll.x.xx * d_lll.z.yy * tmp35
		+ d_lll.x.xx * d_lll.x.yy * tmp23
		- 2. * d_lll.x.xx * d_lll.y.xy * tmp23
		- d_lll.x.xx * d_lll.y.yy * tmp29
		- 2. * d_lll.x.xx * d_lll.y.yz * tmp35
		+ 2. * a_l.z * d_lll.y.yz * gamma_uu.zz
		- a_l.z * d_lll.z.yy * gamma_uu.zz
		- a_l.y * a_l.y
		+ a_l.z * d_lll.y.yy * gamma_uu.yz
		+ 2. * a_l.z * d_lll.y.xy * gamma_uu.xz
		+ 2. * a_l.y * d_lll.y.yz * gamma_uu.yz
		- a_l.y * tmp2810
		- a_l.z * d_lll.x.yy * gamma_uu.xz
		+ a_l.y * d_lll.y.yy * gamma_uu.yy
		+ 2. * a_l.y * d_lll.y.xy * gamma_uu.xy
		+ 2. * a_l.x * d_lll.y.yz * gamma_uu.xz
		- a_l.x * d_lll.z.yy * gamma_uu.xz
		- a_l.y * tmp1606
		+ a_l.x * d_lll.y.yy * gamma_uu.xy
		+ 2. * a_l.x * d_lll.y.xy * gamma_uu.xx
		+ 4. * M_PI * S * gamma_ll.yy
		- 8. * M_PI * S_ll.yy
		- 4. * M_PI * gamma_ll.yy * rho
		- a_l.x * d_lll.x.yy * gamma_uu.xx
		+ K_ll.yy * tmp5
		+ 2. * K_ll.xz * K_ll.yy * gamma_uu.xz
		- 2. * K_ll.yy * tmp3
		+ K_ll.xx * K_ll.yy * gamma_uu.xx
		- 2. * K_ll.xy * tmp1587
		- 4. * K_ll.xy * tmp1589);
	deriv->K_ll.yz += alpha * (2. * gamma_uu.xy * gamma_uu.xz * tmp1307
		- gamma_uu.yy * gamma_uu.yz * tmp3160
		- gamma_uu.yz * tmp3748
		+ d_lll.z.xx * tmp2863
		- 2. * d_lll.z.xx * tmp3006
		- 2. * d_lll.z.xy * tmp2334
		- 2. * d_lll.z.xy * tmp1156
		- d_lll.z.xy * tmp1159
		- 2. * d_lll.z.xy * tmp2344
		- d_lll.z.xy * tmp1166
		- 2. * d_lll.z.xz * tmp1759
		- 2. * d_lll.z.yy * tmp2441
		- d_lll.z.yy * tmp1215
		- 2. * gamma_uu.xx * tmp2614
		- 2. * gamma_uu.xx * tmp2032
		+ d_lll.y.zz * tmp1422
		- 2. * d_lll.y.zz * tmp1248
		- d_lll.y.zz * tmp1259
		- d_lll.z.xx * tmp1672
		+ d_lll.y.zz * tmp1832
		- d_lll.y.zz * tmp1950
		+ 2. * d_lll.y.yz * tmp1215
		+ 4. * d_lll.y.yz * tmp2441
		+ 4. * d_lll.y.yz * tmp2432
		+ 2. * d_lll.y.yz * tmp1203
		+ 2. * d_lll.y.yz * tmp1931
		+ d_lll.y.yy * tmp2391
		- 2. * d_lll.y.yz * tmp1200
		+ 2. * d_lll.y.yy * tmp2386
		+ 2. * d_lll.y.yy * tmp2377
		+ d_lll.y.yy * tmp2373
		+ d_lll.y.yy * tmp2371
		+ d_lll.y.xz * tmp1166
		- d_lll.y.yy * tmp3630
		+ 2. * d_lll.y.xz * tmp2344
		+ d_lll.y.xz * tmp1159
		+ 2. * d_lll.y.xz * tmp2334
		- 2. * d_lll.y.xz * tmp1156
		+ d_lll.y.xz * tmp1909
		+ 2. * d_lll.y.xy * tmp2305
		- d_lll.y.xz * tmp1138
		- 2. * d_lll.y.xz * tmp1697
		- d_lll.y.xz * tmp1143
		- 2. * d_lll.y.xz * tmp1147
		+ 4. * d_lll.y.xy * tmp2508
		+ 4. * d_lll.y.xy * tmp2500
		+ 2. * d_lll.y.xy * tmp2285
		+ 2. * d_lll.y.xy * tmp2491
		+ d_lll.y.xx * tmp2262
		- 2. * d_lll.y.xy * tmp2272
		- 2. * d_lll.y.xy * tmp2277
		+ 2. * d_lll.y.xx * tmp2257
		+ d_lll.y.xx * tmp3387
		- d_lll.y.xx * tmp2253
		+ 2. * d_lll.y.xx * tmp2248
		+ d_lll.y.xx * tmp2244
		+ d_lll.y.xx * d_lll.z.xx * tmp23
		+ d_lll.y.xx * tmp3377
		- 2. * d_lll.y.xx * d_lll.y.zz * tmp47
		+ 2. * d_lll.x.zz * tmp1741
		- 2. * d_lll.x.zz * tmp2940
		- d_lll.x.zz * tmp1746
		- d_lll.x.zz * tmp1992
		- d_lll.x.zz * tmp1756
		- d_lll.y.xx * tmp2240
		+ d_lll.x.zz * d_lll.y.xz * tmp145
		- 2. * d_lll.x.zz * tmp1731
		+ 2. * d_lll.x.zz * tmp1726
		+ d_lll.x.yz * tmp1166
		- 2. * d_lll.x.zz * d_lll.y.xy * tmp84
		+ 2. * d_lll.x.yz * tmp2344
		+ 2. * d_lll.x.yz * tmp2334
		- d_lll.x.yz * tmp1159
		+ 2. * d_lll.x.yz * tmp1150
		- 2. * d_lll.x.yz * tmp1707
		+ d_lll.x.yz * tmp1909
		+ 2. * d_lll.x.yz * tmp1697
		- d_lll.x.yz * tmp1143
		+ d_lll.x.yz * tmp1138
		+ 2. * d_lll.x.yz * tmp1135
		- 2. * d_lll.x.yz * tmp1687
		+ 2. * d_lll.x.yz * tmp1682
		+ d_lll.x.yz * tmp2897
		+ d_lll.x.yz * tmp1128
		+ 2. * d_lll.x.yy * tmp2300
		+ 2. * d_lll.x.yy * tmp2290
		- d_lll.x.yy * tmp2295
		- 2. * d_lll.x.yy * tmp2508
		+ d_lll.x.yy * tmp2493
		- 2. * d_lll.x.yy * tmp2285
		- 2. * d_lll.x.yy * tmp2500
		+ 2. * d_lll.x.yy * tmp3417
		- 2. * d_lll.x.yy * tmp2269
		- d_lll.x.yy * tmp2479
		- d_lll.x.yy * tmp2282
		+ d_lll.x.yy * d_lll.x.yz * tmp72
		+ 2. * d_lll.x.xz * tmp1648
		- 2. * d_lll.x.xz * tmp1655
		- 2. * d_lll.x.xz * d_lll.y.zz * tmp145
		- 2. * d_lll.x.xz * tmp1672
		- 2. * d_lll.x.xz * tmp2863
		+ 2. * d_lll.x.xy * tmp2236
		- 2. * d_lll.x.xy * tmp2240
		- 2. * d_lll.x.xy * tmp3377
		- 2. * d_lll.x.xy * tmp2244
		- 2. * d_lll.x.xy * tmp3387
		+ d_lll.x.xx * d_lll.x.yz * tmp23
		- d_lll.x.xx * d_lll.y.xz * tmp23
		- d_lll.x.xx * d_lll.y.zz * tmp35
		- d_lll.x.xx * d_lll.z.xy * tmp23
		- d_lll.x.xx * d_lll.z.yy * tmp29
		+ a_l.z * tmp2810
		+ a_l.z * tmp1618
		+ a_l.z * d_lll.y.zz * gamma_uu.zz
		+ a_l.z * tmp1615
		+ a_l.y * d_lll.z.yy * gamma_uu.yy
		- a_l.z * tmp1609
		+ a_l.y * tmp2206
		+ a_l.y * tmp3340
		+ a_l.y * tmp2202
		+ a_l.x * d_lll.z.yy * gamma_uu.xy
		- a_l.y * a_l.z
		- a_l.y * tmp2196
		+ a_l.x * d_lll.z.xy * gamma_uu.xx
		+ a_l.x * d_lll.y.zz * gamma_uu.xz
		+ a_l.x * d_lll.y.xz * gamma_uu.xx
		+ 4. * M_PI * S * gamma_ll.yz
		- 8. * M_PI * S_ll.yz
		- 4. * M_PI * gamma_ll.yz * rho
		- a_l.x * d_lll.x.yz * gamma_uu.xx
		+ K_ll.xx * K_ll.yz * gamma_uu.xx
		- 2. * K_ll.xy * tmp2176
		- 2. * K_ll.xy * tmp2179
		- 2. * K_ll.xz * tmp1587
		- K_ll.yy * tmp2181
		- 2. * K_ll.yy * tmp2183
		- K_ll.yz * tmp5);
	deriv->K_ll.zz += alpha * (gamma_uu.yy * tmp3748
		- gamma_uu.zz * K_ll.zz * K_ll.zz
		+ gamma_uu.xx * tmp2625
		- 2. * gamma_uu.yy * tmp3162
		+ 2. * gamma_uu.xx * tmp2622
		+ 2. * gamma_uu.xx * gamma_uu.yy * tmp1405
		+ 2. * gamma_uu.xx * tmp2619
		+ tmp3160 * tmp204
		- 2. * gamma_uu.xx * tmp1476
		+ d_lll.z.yy * tmp4203
		+ 2. * d_lll.z.yy * tmp2386
		+ 2. * d_lll.z.xz * tmp1159
		+ 2. * tmp1447 * tmp41
		+ 2. * d_lll.z.xy * tmp4063
		+ 4. * d_lll.z.xy * tmp2508
		+ 4. * d_lll.z.xy * tmp2295
		+ 4. * d_lll.z.xy * tmp2290
		+ tmp1455 * tmp23
		+ d_lll.z.xx * tmp4000
		+ 2. * d_lll.z.xx * tmp3995
		+ 2. * d_lll.z.xx * tmp2253
		+ 2. * d_lll.z.xx * tmp2248
		+ 4. * d_lll.z.xx * tmp2244
		+ d_lll.y.zz * tmp1209
		+ 2. * d_lll.y.zz * tmp1364
		+ 2. * d_lll.y.zz * tmp1931
		+ 2. * d_lll.y.yz * tmp3630
		- 4. * d_lll.y.yz * d_lll.z.xz * tmp108
		- 4. * d_lll.y.yz * tmp2386
		- 2. * d_lll.y.yz * tmp4203
		- d_lll.y.zz * tmp1354
		+ d_lll.y.yy * d_lll.y.zz * tmp204
		- 2. * d_lll.y.yy * d_lll.z.xz * tmp96
		- 2. * d_lll.y.yy * d_lll.z.yz * tmp204
		- d_lll.y.yy * d_lll.z.zz * tmp210
		+ 2. * d_lll.y.xz * tmp2282
		- 4. * d_lll.y.xz * tmp2290
		- 4. * d_lll.y.xz * tmp2508
		- 2. * d_lll.y.xz * tmp4063
		- 2. * tmp1405 * tmp41
		+ 2. * d_lll.y.xy * tmp4012
		- 4. * d_lll.y.xy * tmp4014
		- 4. * d_lll.y.xy * tmp4022
		- 2. * d_lll.y.xy * tmp4032
		+ d_lll.y.xx * tmp3979
		- 2. * d_lll.y.xx * d_lll.z.zz * tmp47
		+ 2. * d_lll.y.xx * tmp3974
		- 4. * d_lll.y.xx * d_lll.z.yz * tmp41
		+ 2. * d_lll.y.xx * d_lll.y.zz * tmp41
		- 2. * d_lll.y.xx * tmp3969
		+ 2. * d_lll.x.zz * tmp1156
		- d_lll.x.zz * tmp1159
		- d_lll.y.xx * tmp3966
		+ 2. * d_lll.x.zz * tmp1150
		+ d_lll.x.zz * tmp1909
		+ 2. * d_lll.x.zz * tmp1143
		+ 2. * d_lll.x.zz * tmp1140
		+ d_lll.x.zz * tmp1138
		+ 4. * d_lll.x.zz * tmp1687
		+ 2. * d_lll.x.zz * tmp1132
		- 2. * d_lll.x.zz * tmp1135
		+ d_lll.x.zz * tmp2897
		+ 4. * d_lll.x.yz * tmp2277
		- 2. * d_lll.x.yz * tmp2282
		- 4. * d_lll.x.yz * tmp2290
		- 4. * d_lll.x.yz * tmp2508
		- 2. * d_lll.x.yz * tmp4063
		- 2. * tmp1307 * tmp41
		+ 4. * d_lll.x.yz * tmp2272
		+ 2. * d_lll.x.yz * tmp3417
		- 4. * d_lll.x.yz * tmp2479
		+ d_lll.x.yy * tmp4032
		+ 2. * d_lll.x.yy * tmp4014
		- 4. * d_lll.x.yy * d_lll.z.xz * tmp41
		- 2. * d_lll.x.yy * tmp4022
		- 2. * d_lll.x.yy * d_lll.z.zz * tmp102
		+ d_lll.x.yy * tmp4012
		+ 2. * d_lll.x.yy * d_lll.x.zz * tmp41
		+ 2. * d_lll.x.xz * tmp3377
		- 4. * d_lll.x.xz * tmp2248
		- 4. * d_lll.x.xz * tmp3995
		- 2. * d_lll.x.xz * tmp4000
		- d_lll.x.yy * d_lll.x.zz * tmp72
		+ 2. * d_lll.x.xz * tmp2238
		+ 2. * d_lll.x.xy * tmp3966
		- 4. * d_lll.x.xy * tmp3969
		- 4. * d_lll.x.xy * tmp3974
		- 2. * d_lll.x.xy * tmp3979
		+ 2. * d_lll.x.xy * d_lll.x.zz * tmp29
		+ d_lll.x.xx * d_lll.y.zz * tmp29
		- 2. * d_lll.x.xx * d_lll.z.xz * tmp23
		- 2. * d_lll.x.xx * d_lll.z.yz * tmp29
		- d_lll.x.xx * d_lll.z.zz * tmp35
		+ d_lll.x.xx * d_lll.x.zz * tmp23
		+ a_l.z * d_lll.z.zz * gamma_uu.zz
		- a_l.z * a_l.z
		+ 2. * a_l.z * d_lll.z.yz * gamma_uu.yz
		+ 2. * a_l.z * d_lll.z.xz * gamma_uu.xz
		+ a_l.y * d_lll.z.zz * gamma_uu.yz
		- a_l.z * tmp2199
		- a_l.z * tmp3340
		+ 2. * a_l.y * d_lll.z.yz * gamma_uu.yy
		+ 2. * a_l.y * d_lll.z.xz * gamma_uu.xy
		+ a_l.x * d_lll.z.zz * gamma_uu.xz
		- a_l.y * d_lll.x.zz * gamma_uu.xy
		- a_l.y * d_lll.y.zz * gamma_uu.yy
		+ 2. * a_l.x * d_lll.z.yz * gamma_uu.xy
		+ 2. * a_l.x * d_lll.z.xz * gamma_uu.xx
		+ 4. * M_PI * S * gamma_ll.zz
		- 8. * M_PI * S_ll.zz
		- 4. * M_PI * gamma_ll.zz * rho
		- a_l.x * d_lll.x.zz * gamma_uu.xx
		- a_l.x * d_lll.y.zz * gamma_uu.xy
		+ K_ll.yy * K_ll.zz * gamma_uu.yy
		- 2. * K_ll.yz * tmp2183
		+ 2. * K_ll.xy * K_ll.zz * gamma_uu.xy
		- 4. * K_ll.xz * tmp2177
		- 2. * K_ll.xz * tmp2179
		+ K_ll.xx * K_ll.zz * gamma_uu.xx);	
	// END CUT

#endif

<? if eqn.useShift ~= 'none' then ?>


	<?
	if eqn.useShift == 'HarmonicShiftCondition-FiniteDifference'
	or eqn.useShift == 'MinimalDistortionElliptic'
	or eqn.useShift == 'MinimalDistortionEllipticEvolve'
	then 
	?>

	//make each useShift option into an object
	//and give it a method for producing the partial_beta_ul
	// so the finite difference ones can use this
	// and the state ones can just assign a variable
	//then all this can be moved into the "if eqn.useShift ~= 'none' then" block
	//but then again, any FVS-based method will already have the U_,k beta^k terms incorporated into its eigensystem ...
	//this means the useShift will also determine whether the flux has any shift terms
	// for any '-FiniteDifference' shift, the flux won't have any shift terms.

	//partial_beta_ul[j].i := beta^i_,j
<?=eqn:makePartial1'beta_u'?>	

	//= gamma_ik beta^k_,j
	real3x3 partial_beta_u_ll = (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_real3_mul(gamma_uu, partial_beta_ul[<?=i-1?>]),
<? end
?>	};

	//first add the Lie derivative beta terms to the time derivatives

	//alpha_,t = alpha_,i beta^i + ...
	// = alpha a_i beta^i + ..
	deriv->alpha += real3_dot(U->beta_u, U->a_l) * U->alpha;

	//Alcubierre 2.3.11
	//gamma_ij,t = gamma_ij,k beta^k + gamma_kj beta^k_,i + gamma_ik beta^k_,j
	// = 2 d_kij beta^k + gamma_kj beta^k_,i + gamma_ik beta^k_,j
	//
	//Alcubierre eqn 2.3.11 and B&S eqn 2.128 rewrite it as:
	//gamma_ij,t = -2 alpha K_ij + D_i beta_j + D_j beta_i
	//	= -2 alpha K_ij + D_i beta^k gamma_kj + D_j beta^k gamma_ki
	//	= -2 alpha K_ij + (beta^k_,i + Gamma^k_li beta^l) gamma_kj + (beta^k_,j + Gamma^k_lj beta^l) gamma_ki
	//	= -2 alpha K_ij + beta^k_,i gamma_kj + beta^k_,j gamma_ki + (Gamma_jki + Gamma_ikj) beta^k
	//	= -2 alpha K_ij + beta^k_,i gamma_kj + beta^k_,j gamma_ki + (Gamma_jik + Gamma_ijk) beta^k
	//Gamma_ijk = -d_ijk + d_jik + d_kij
	//Gamma_jik = -d_jik + d_ijk + d_kji
	//Gamma_ijk + Gamma_jik = 2 d_kij
	//gamma_ij,t = -2 alpha K_ij + beta^k_,i gamma_kj + beta^k_,j gamma_ki + 2 d_kij beta^k
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	deriv->gamma_ll.<?=xij?> += partial_beta_u_ll.<?=xi?>.<?=xj?>
		+ partial_beta_u_ll.<?=xj?>.<?=xi?>
<?	for k,xk in ipairs(xNames) do
?> 		+ 2. * U->d_lll.<?=xk?>.<?=xij?> * U->beta_u.<?=xk?>
<?	end
?>	;
<? end ?>

#if 0	//the hyperbolic vars should get shifted via the wavespeeds
		//I should only need to manually shift the source-only vars


	//partial_a_l[j].i = a_i,j
<?=eqn:makePartial1'a_l'?>

	//a_i,t = a_i,j beta^j + a_j beta^j_,i
<? for i,xi in ipairs(xNames) do
?>	deriv->a.<?=xi?> += 0.<?
	for j,xj in ipairs(xNames) do
?>		+ partial_a_l[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
		+ U->a.<?=xj?> * partial_beta_ul[<?=i-1?>].<?=xj?>
<?	end
?>	;
<? end
?>

	//technically this can be represented as a sym3sym3 since d_ijk,l is symmetric with jk and il
	//partial_d_llll[k][l].ij = d_kij,l
<?=eqn:makePartial1'd_lll'?>

	//d_kij,t = d_kij,l beta^l + d_lij beta^l_,k + d_klj beta^l_,i + d_kil beta^l_,j
<? for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>	deriv->d_lll.<?=xk?>.<?=xij?> += 0.
<?		for l,xl in ipairs(xNames) do
?>			+ partial_d_llll[<?=l-1?>].<?=xk?>.<?=xij?>
<?		end
?>	;
<?	end
end ?>

	//partial_K_lll[k].ij = K_ij,k
<?=eqn:makePartial1'K_ll'?>

	//K_ij,t = K_ij,k beta^k + K_kj beta^k_,i + K_ik beta^k_,j
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>
	deriv->K.<?=xij?> += 0.
<?	for k,xk in ipairs(xNames) do
?>		+ partial_K_l[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>
		+ U->K.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
		+ U->K.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
<? 	end 
?>	;
<? end ?>

	//partial_V_l[j].i = V_i,j
<?=eqn:makePartial1'V_l'?>

	//V_i,t = V_i,j beta^j + V_j beta^i_,j
<? for i,xi in ipairs(xNames) do
?>	deriv->V_l.<?=xi?> += 0. <?
	for j,xj in ipairs(xNames) do
?>		+ partial_V_l[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
		+ U->V_l.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?	end
?>	;
<? end
?>

#endif


	<? 
	end	-- using partial_beta_ul to integrate alpha, gamma_ll, 
	if eqn.useShift == '2005 Bona / 2008 Yano' then ?>
	
	
	
	<? elseif eqn.useShift == 'HarmonicShiftCondition-FiniteDifference' then 
	?>

	//next add the source for the particular useShift

	// 2008 Alcubierre 4.3.37
	//beta^i_,t = beta^j beta^i_,j - alpha alpha^,i + alpha^2 conn^i + beta^i / alpha (alpha_,t - beta^j alpha_,j + alpha^2 K)
	//using alpha_,t = beta^i alpha_,i - alpha^2 f K
	//= beta^j beta^i_,j + alpha^2 (conn^i - a^i) + beta^i (beta^j a_j - alpha f K - beta^j a_j + alpha K)
	//= beta^j beta^i_,j + alpha^2 (conn^i - a^i) + alpha K beta^i (1 - f)

	real3 dbeta_beta = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + partial_beta_ul[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?><?
	end ?>,
<? end
?>	};
	
	real3 a_u = sym3_real3_mul(gamma_uu, U->a_l);

	//conn^i = conn^i_jk gamma^jk
	real3 conn_u = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(conn_ull.<?=xi?>, gamma_uu),
<? end
?>	};

	deriv->beta_u = real3_add(
		deriv->beta_u,
		real3_add(
			real3_add(
				dbeta_beta,
				real3_real_mul(
					U->beta_u,
					U->alpha * tr_K * (1. - f)
				)
			),
			real3_real_mul(
				real3_sub(
					conn_u,
					a_u
				),
				U->alpha * U->alpha
			)
		)
	);

	<?
	elseif eqn.useShift == 'MinimalDistortionEllipticEvolve' then
	?>

	/*
	beta^i_,t = epsilon (gamma^jk D_j D_k beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j - D_j (2 alpha (K^ij - 1/3 K gamma^ij))
	= epsilon (
		gamma^jk D_j (beta^i_,k + Gamma^i_lk beta^l) 
		+ 1/3 gamma^ik D_k (beta^j_,j + Gamma^j_lj beta^l) 
		+ R^i_j beta^j 
		- (2 alpha_,j (K^ij - 1/3 K gamma^ij)
		- (2 alpha (D_j K^ij - 1/3 K_,j gamma^ij)
	)
	= epsilon (
		- 2 alpha_,j K^ij 
		+ 2/3 alpha_,j K gamma^ij
		+ gamma^jk Gamma^i_lk,j beta^l 
		+ 2 Gamma^ik_j beta^j_,k 
		- Gamma^kj_j beta^i_,k 
		+ 1/3 gamma^ik Gamma^j_lj beta^l_,k 
		+ gamma^jk beta^i_,kj 
		+ 1/3 gamma^ik beta^j_,jk 
		+ 1/3 gamma^ik Gamma^j_lj,k beta^l 
		- 2 alpha K^ij_,j 
		+ 2/3 alpha K_,j gamma^ij
		+ Gamma^i_jm Gamma^mj_l beta^l 
		- Gamma^m_jl Gamma^ij_m beta^l 
		- Gamma^mj_j Gamma^i_lm beta^l
		+ Gamma^ij_l Gamma^l_kj beta^k
		+ 1/3 gamma^ik Gamma^j_lj Gamma^l_mk beta^m
		+ R^i_j beta^j 
		- 2 alpha Gamma^i_kj K^ki 
		- 2 alpha Gamma^j_kj K^ik 
	)
	*/

	<?
	elseif eqn.useShift == 'LagrangianCoordinates' 
	or eqn.useShift == 'MinimalDistortionElliptic' 
	then
		-- do nothing
	else
		error("I don't have any source terms implemented for this particular useShift")
	end ?>

<? end	-- eqn.useShift ?>

	// and now for the first-order constraints

	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	<? for i,xi in ipairs(xNames) do ?>{
		<? if i <= solver.dim then ?>
		real partial_i_log_alpha = (
			log(U[solver->stepsize.<?=xi?>].alpha)
			- log(U[-solver->stepsize.<?=xi?>].alpha)
		) / (2. * solver->grid_dx.s<?=i-1?>);
		<? else ?>
		real partial_i_log_alpha = 0.;
		<? end ?>
		deriv->a_l.<?=xi?> += solver->a_convCoeff * (partial_i_log_alpha - U->a_l.<?=xi?>);
	}<? end ?>	
	
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	<? 
for i,xi in ipairs(xNames) do 
	for jk,xjk in ipairs(symNames) do ?>{
		<? if i <= solver.dim then ?>
		real partial_i_gamma_jk = (
			U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> 
			- U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>
		) / (2. * solver->grid_dx.s<?=i-1?>);
		<? else ?>
		real partial_i_gamma_jk = 0;
		<? end ?>
		deriv->d_lll.<?=xi?>.<?=xjk?> += solver->d_convCoeff * (.5 * partial_i_gamma_jk - U->d_lll.<?=xi?>.<?=xjk?>);
	}<? 
	end
end ?>


	//add stress-energy terms and H damping:
	//TODO I see it's a term in the K_ij,t ADM formalism
	//but where is the idea of adding this to gamma_ij,t from?

<? for ij,xij in ipairs(symNames) do
?>	deriv->gamma_ll.<?=xij?> += U->alpha * (
		+ solver->gamma_ll_srcStressCoeff * stressConstraint_ll.<?=xij?>
		+ solver->gamma_ll_srcHCoeff * U->gamma_ll.<?=xij?> * HamiltonianConstraint
	);
<? end 
?>

<? for ij,xij in ipairs(symNames) do
?>	deriv->K_ll.<?=xij?> += U->alpha * (
			+ solver->K_ll_srcStressCoeff * stressConstraint_ll.<?=xij?>
			+ solver->K_ll_srcHCoeff * U->gamma_ll.<?=xij?> * HamiltonianConstraint
	);
<? end 
?>

	//converge V_i to its constrained value:


	//V_i = d_ik^k - d^k_ki <=> V_i += eta (d_ik^k - d^k_ki - V_i)
	deriv->V_l = real3_add(
		deriv->V_l,
		real3_real_mul(
			real3_sub(real3_sub(d_l, e_l), U->V_l),
			solver->V_convCoeff));

	//Kreiss-Oligar diffusion, for stability's sake?
	//isn't that just a hack to improve the stability of finite-difference, which diverges by nature?
}

<? elseif moduleName == "constrainU" then 
depmod{
	"real3x3x3",
	"sym3sym3",
}
?>

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost);		
	global cons_t* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real tr_K = real3x3_trace(K_ul);							//K^k_k
	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		//KSq_ij = K_ik K^k_j
	sym3 K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);			//K^ij

<?
local constrainVGuiVar = eqn.guiVars['constrain V']
local constrainV = constrainVGuiVar.options[constrainVGuiVar.value]  
if constrainV ~= 'none' then 
?>	//gravitational_wave_sim uses this (for 1D), HydroGPU doesn't (for 2D/3D)

	real3 delta;
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d_lll.<?=xi?>, gamma_uu);
		real d2 = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + U->d_lll.<?=xj?>.<?=sym(k,i)?> * gamma_uu.<?=sym(j,k)?><?
		end
	end ?>;
		delta.<?=xi?> = U->V_l.<?=xi?> - (d1 - d2);
	}<? end ?>

<?
	if constrainV == 'replace V' then
?>
	//directly assign to V
	U->V_l = real3_sub(U->V_l, delta);
<?
	elseif constrainV == 'average' then
?>
	//average between V and d

/*
interpolation between V and d
weight = 1 means only adjust V (gives a 1e-20 error in the constraint eqns)
weight = 0 means only adjust d (hmm, is only at 1e-2 ...)

V_i - (d_im^m - d^m_mi) = Q_i ...
V'_i = V_i - alpha Q_i
(d_im^m - d^m_mi)' = (d_im^m - d^m_mi) + (1-alpha) Q_i
so (V_i - (d_im^m - d^m_mi))' = V'_i - (d_im^m - d^m_mi)' = V_i - alpha Q_i - (d_im^m - d^m_mi) - (1-alpha) Q_i = Q_i - Q_i = 0
works for d'_ijk = d_ijk + (1-alpha)/4 (Q_i gamma_jk - Q_k gamma_ij)
so d'_im^m = d'_ijk gamma^jk = d_im^m + (1-alpha)/4 (Q_i gamma_jk - Q_k gamma_ij) gamma^jk = d_im^m + (1-alpha)/2 Q_i
and d'^m_mi = d'_kji gamma^jk = d^m_mi + (1-alpha)/4 (Q_k gamma_ij - Q_i gamma_jk) gamma^jk = d^m_mi - (1-alpha)/2 Q_i
therefore (d_im^m - d^m_mi)' = d_im^m - d^m_mi + (1-alpha) Q_i

...but how come it's not working?
...probably because neighbors are influencing one another
so to solve this, you must solve a giant linear system of all variables
*/
	const real weight = .5;
	const real v_weight = weight;
	const real d_weight = (1. - weight) / 4.;

	U->V_l = real3_sub(U->V_l, real3_real_mul(delta, v_weight));
<? 
		for i,xi in ipairs(xNames) do 
			for jk,xjk in ipairs(symNames) do
				local j,k = from6to3x3(jk)
				local xk = xNames[k]
?>	U->d_lll.<?=xi?>.<?=xjk?> += (
		delta.<?=xi?> * U->gamma_ll.<?=xjk?> 
		- delta.<?=xk?> * U->gamma_ll.<?=sym(i,j)?>
	) * d_weight;
<?	
			end
		end 
	end
?>

//...or linearly project out the [V_i, U->d_ijk] vector
//...or do a single gradient descent step
<?
end	-- constrain V
?>

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3x3 d_llu = _3sym3_sym3_mul(U->d_lll, gamma_uu);
	
	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//e_i = d^j_ji
	real3 e_l = _3sym3_tr12(d_ull);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?> - d_ull.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	//partial_d_lll.ij.kl = d_kij,l = d_(k|(ij),|l)
	//so this object's indexes are rearranged compared to the papers 
	sym3sym3 partial_d_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
		local k,l,xk,xl = from6to3x3(kl)
?>	partial_d_llll.<?=xij?>.<?=xkl?> = 0.
<?		if l <= solver.dim then
?>	+ .5 * (	// 1/2 d_kij,l
		U[solver->stepsize.<?=xl?>].d_lll.<?=xk?>.<?=xij?>
		- U[-solver->stepsize.<?=xl?>].d_lll.<?=xk?>.<?=xij?>
	) / (2. * solver->grid_dx.<?=xl?>);
<?
		end
		if k <= solver.dim then
?>	+ .5 * (
		U[solver->stepsize.<?=xk?>].d_lll.<?=xl?>.<?=xij?>
		- U[-solver->stepsize.<?=xk?>].d_lll.<?=xl?>.<?=xij?>
	) / (2. * solver->grid_dx.<?=xk?>);
<?		end
?>	;
<?	end
end
?>

	// R_ll.ij := R_ij
	//	= gamma^kl (-gamma_ij,kl - gamma_kl,ij + gamma_ik,jl + gamma_jl,ik)
	//		+ conn^k_ij (d_k - 2 e_k)
	//		- 2 d^l_ki d^k_lj 
	//		+ 2 d^l_ki d_lj^k 
	//		+ d_il^k d_jk^l
	sym3 R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * (U->V_l.<?=xk?> - e_l.<?=xk?>)
<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_ull.<?=xl?>.<?=sym(k,i)?> * (d_llu.<?=xl?>.<?=xj?>.<?=xk?> - d_ull.<?=xk?>.<?=sym(l,j)?>)
			+ d_llu.<?=xi?>.<?=xl?>.<?=xk?> * d_llu.<?=xj?>.<?=xk?>.<?=xl?>
			+ gamma_uu.<?=sym(k,l)?> * (
				- partial_d_llll.<?=xij?>.<?=sym(k,l)?>
				- partial_d_llll.<?=sym(k,l)?>.<?=xij?>
				+ partial_d_llll.<?=sym(i,k)?>.<?=sym(j,l)?>
				+ partial_d_llll.<?=sym(j,l)?>.<?=sym(i,k)?>
			)
<? 		end
	end
?>		,
<? end
?>	};

	//scaled down by 1/8 to match B&S BSSNOK equations ... maybe I'll scale theirs up by 8 ...
	//B&S eqn 2.125 ... divded by two
	//Alcubierre eqn 2.5.9
	//H = 1/2 (R + K^2 - K_ij K^ij) - 8 pi rho
	real R = sym3_dot(R_ll, gamma_uu);
	real tr_KSq = sym3_dot(KSq_ll, gamma_uu);
	U->H = .5 * (R + tr_K * tr_K - tr_KSq) <? 
if eqn.useStressEnergyTerms then ?>
	- 8. * M_PI * U->rho <? 
end ?>;

<?=eqn:makePartial1'K_ll'?>	

	/*
	momentum constraint
	Alcubierre eqn 2.4.11
	M^i = K^ij_;j 
		- gamma^ij K_,j 
		- 8 pi S^i
	M^i = gamma^im gamma^jn (
			K_mn,j 
			- conn^k_mj K_kn 
			- conn^k_nj K_mk
		) 
		- gamma^ij (gamma^mn_,j K_mn + gamma^mn K_mn,j)
		- 8 pi S^i
	M^i = gamma^ij (
			gamma^mn (
				K_jn,m 
				- K_mn,j
				- conn^k_jm K_kn 
				- conn^k_nm K_jk
			) 
			+ 2 d_jmn K^mn 
		)
		- 8 pi S^i
	*/
<? for i,xi in ipairs(xNames) do 
?>	U->M_u.<?=xi?> = 
<?	for j,xj in ipairs(xNames) do
		for m,xm in ipairs(xNames) do
			for n,xn in ipairs(xNames) do
?>		+ gamma_uu.<?=sym(i,j)?> * (
			gamma_uu.<?=sym(m,n)?> * (0.
				+ partial_K_lll[<?=m-1?>].<?=sym(j,n)?>
				- partial_K_lll[<?=j-1?>].<?=sym(m,n)?>
<?				for k,xk in ipairs(xNames) do
?>				- conn_ull.<?=xk?>.<?=sym(j,m)?> * U->K_ll.<?=sym(k,n)?>
				- conn_ull.<?=xk?>.<?=sym(n,m)?> * U->K_ll.<?=sym(j,k)?>
<?				end			
?>			)
			+ 2. * U->d_lll.<?=xj?>.<?=sym(m,n)?> * K_uu.<?=sym(m,n)?>
		)
<?			end
		end
	end
	if eqn.useStressEnergyTerms then
?>		- 8. * M_PI * U->S_u.<?=xi?>
<?	end
?>	;
<? end
?>

}

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
