//// MODULE_NAME: <?=eqn_common?>

//pressure function for ideal gas
real calc_P(
	constant <?=solver_t?> const * const solver,
	real const rho,
	real const eInt
) {
	return (solver->heatCapacityRatio - 1.) * rho * eInt;
}	

//chi in most papers
real calc_dP_drho(
	constant <?=solver_t?> const * const solver,
	real const rho,
	real const eInt
) {
	return (solver->heatCapacityRatio - 1.) * eInt;
}

//kappa in most papers
real calc_dP_deInt(
	constant <?=solver_t?> const * const solver,
	real const rho,
	real const eInt
) {
	return (solver->heatCapacityRatio - 1.) * rho;
}

real calc_eInt_from_P(
	constant <?=solver_t?> const * const solver,
	real const rho,
	real const P
) {
	return P / ((solver->heatCapacityRatio - 1.) * rho);
}

real calc_h(real rho, real P, real eInt) {
	return 1. + eInt + P / rho;
}

<?=eqn.cons_only_t?> consOnlyFromPrim(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const prim,
	real const alpha,
	real3 const beta,
	real3s3 const gamma
) {
	//2008 Font eqn 31 etc 
	real det_gamma = real3s3_det(gamma);
	real3s3 gammaU = real3s3_inv(gamma, det_gamma);
	real3 vU = real3s3_real3_mul(gammaU, prim.v);
	real vSq = real3_dot(prim.v, vU);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(solver, prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 28-30:
	
	real D = prim.rho * W;
	real3 S = real3_real_mul(prim.v, prim.rho * h * WSq);
	real tau = prim.rho * h * WSq - P - D;

	return (<?=eqn.cons_only_t?>){.D=D, .S=S, .tau=tau};
}

//// MODULE_NAME: <?=applyInitCondCell?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> * const cell<?=
	solver:getADMArgs()?>
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool const lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	//ignored:
	real3 B = real3_zero;

	<?=solver:getADMVarCode()?>

	<?=code?>
	
	real eInt = calc_eInt_from_P(solver, rho, P);

	<?=prim_only_t?> prim = {.rho=rho, .v=v, .eInt=eInt};
	U->prim = prim;
	U->cons = consOnlyFromPrim(solver, prim, alpha, beta, gamma);
}


//// MODULE_NAME: <?=calcDTCell?>

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
void <?=calcDTCell?>(
	global real * const dt,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell<?=
	solver:getADMArgs()?>
) {
	<?=prim_only_t?> prim = primOnlyFromPrim(U, solver, x);
	<?=solver:getADMVarCode()?>

	real const det_gamma = real3s3_det(gamma);
	real3s3 const gammaU = real3s3_inv(gamma, det_gamma);
	real3 const vU = real3s3_real3_mul(gammaU, prim.v);

	real const rho = prim.rho;
	real const eInt = prim.eInt;
	//2008 Font Eqn 31: v^2 = gamma_ij v^i v^j
	real const vSq = real3_dot(prim.v, vU);
	real const P = calc_P(solver, rho, eInt);
	real const h = calc_h(rho, P, eInt);
	real const csSq = solver->heatCapacityRatio * P / (rho * h);
	real const cs = sqrt(csSq);
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		//for the particular direction
		real vUi = vU.s<?=side?>;
		real vUiSq = vUi * vUi;
		
		//Font 2008 eqn 106
		real const betaUi = beta.s<?=side?>;
		real discr = sqrt((1. - vSq) * (gammaU.xx * (1. - vSq * csSq) - vUiSq * (1. - csSq)));
		real lambdaMin = (vUi * (1. - csSq) - cs * discr) * alpha / (1. - vSq * csSq) - betaUi;
		real lambdaMax = (vUi * (1. - csSq) + cs * discr) * alpha / (1. - vSq * csSq) - betaUi;
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		*(dt) = (real)min(*(dt), solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
}

//// MODULE_NAME: <?=fluxFromCons?>

<? if false then ?>
<?=cons_t?> <?=fluxFromCons?>(
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const U,
	<?=cell_t?> const * const cell,
	<?=normal_t?> const n,<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	real const det_gamma = real3s3_det(gamma);
	real3s3 const gammaU = real3s3_inv(gamma, det_gamma);
	
	real const vUi = gammaU.<?=sym(side+1,1)?> * U.prim.v.x
			+ gammaU.<?=sym(side+1,2)?> * U.prim.v.y
			+ gammaU.<?=sym(side+1,3)?> * U.prim.v.z;
	real const vUi_shift = vUi - beta.s<?=side?> / alpha;

	<?=cons_t?> F;
	F.cons.D = U.cons.D * vUi_shift;
	F.S = real3_real_mul(U.cons.S, vUi_shift);
	F.S.s<?=side?> += U.prim.p;
	F.tau = U.cons.tau * vUi_shift + p * vUi;
	F.prim = (<?=prim_t?>){
		.rho = 0,
		.v = real3_zero,
		.eInt = 0,
	};

	return F;
}
<? end ?>

//// MODULE_NAME: <?=eigen_forCell?>

//used by PLM
//TODO SRHD PLM needs to do this:
//1) calcLR for the <?=prim_t?> (that means put calcLR in its own file, and a new primLR buf)
//2) have a new kernel for calc consLR from primLR, since calcDeltaUEig and calcFlux both need this
//or does the eigenbasis need to be derived from the variables being transformed?
//shoud I PLM the U's then converge the prims ... and therefore track the prims on edges as well?
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*real3 const */n\
)

//// MODULE_NAME: <?=calcEigenBasis?>

#error <?=calcEigenBasis?> has been removed, and eigen_t structs are now calculated inline ... soooo ... convert this to something compatible
kernel void <?=calcEigenBasis?>(
	constant <?=solver_t?> const * const solver,
	global <?=eigen_t?> * const eigenBuf,
	
	//TODO 
	//turn this into a LR extrapolation
	//actually make use of PLM somehow 
	//right now only primBuf is being used for getting neighbor values
	//so SRHD should perform the PLM stuff on the primBuf instead of the UBUf?
	// or do the PLM on the UBuf and do the cons->prim on the ULR edge values
	global <?=cons_t?> const * const UBuf<?=
	solver:getADMArgs()?>
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost - 1);
	
	int const indexR = index;
	<?=prim_t?> primR = UBuf[indexR].prim;
	
	<?=solver:getADMVarCode{suffix='R'} --[[ produce alphaR, betaR, gammaR at indexR ]] ?>
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		
		int const indexL = index - solver->stepsize.s<?=side?>;
		<?=prim_t?> primL = UBuf[indexL].prim;
	
		<?=solver:getADMVarCode{suffix='L'} --[[ produce alphaL, betaL, gammaL at indexL ]] ?>

<? if true then -- arithmetic averaging ?>
		<?=prim_t?> avg = (<?=prim_t?>){
			.rho = .5 * (primL.rho + primR.rho),
			.v = real3_real_mul(real3_add(primL.v, primR.v), .5),
			.eInt = .5 * (primL.eInt + primR.eInt),
		};
		real alpha = .5 * (alphaL + alphaR);
		real3 beta = real3_real_mul(real3_add(betaL, betaR), .5);
		real3s3 gamma = real3s3_real_mul(real3s3_add(gammaL, gammaR), .5);
<? -- else -- Roe-averaging, Font 2008 eqn 38 ?>
<? end ?>
		
		real rho = avg.rho;
		real3 vL = avg.v;
		real eInt = avg.eInt;
			
		//these match <?=eigen_leftTransform?>
		<? if side == 1 then ?>
		//v' = P [vx,vy,vz] = [vy,-vx,vz]
		/*
		should gamma' = P gamma P' or P' gamma P ?
		I'm going to use P gamma P' because the signs match how 'v' changes as well
		however, because gamma the metric of v, maybe it should be opposite? 
		Also, the norms work out, so (P gamma P') (P v) has the same norm as (gamma v)
		Same with (P v)' (P gamma P') (P v) works vs (P v)' (P' gamma P) (P v) which doesn't
		*/
		vL = _real3(vL.y, -vL.x, vL.z);	// -90' rotation to put the y axis contents into the x axis
		beta = _real3(beta.y, -beta.x, beta.z);
		gamma = (real3s3){
			.xx = gamma.yy,
			.xy = -gamma.xy,
			.xz = gamma.yz,
			.yy = gamma.xx,
			.yz = -gamma.xz,
			.zz = gamma.zz,
		};
		<? elseif side == 2 then ?>
		//x,z -> z,-x
		vL = _real3(vL.z, vL.y, -vL.x);	//-90' rotation to put the z axis in the x axis
		beta = _real3(beta.z, beta.y, -beta.x);
		gamma = (real3s3){
			.xx = gamma.zz,
			.xy = gamma.yz,
			.xz = -gamma.xz,
			.yy = gamma.yy,
			.yz = -gamma.xy,
			.zz = gamma.xx,
		};
		<? end ?>
		
		real det_gamma = real3s3_det(gamma);
		real3s3 gammaU = real3s3_inv(gamma, det_gamma);

		real3 vU = real3s3_real3_mul(gammaU, vL);
		real vSq = real3_dot(vL, vU);
		real oneOverW2 = 1. - vSq;
		real oneOverW = sqrt(oneOverW2);
		real W = 1. / oneOverW;	//alpha?
		real W2 = 1. / oneOverW2;
		real P = (solver->heatCapacityRatio - 1.) * rho * eInt;
		real h = 1. + eInt + P / rho;

		real hW = h * W;

		//just after 2008 Font eqn 107:
		//h cs^2 = chi + P / rho^2 kappa = dp/drho + p / rho^2 dp/deInt
		// = (gamma-1) eInt + P/rho^2 (gamma-1) rho  for an ideal gas
		// = (gamma-1) (eInt + P/rho)
		// = 1/rho ( (gamma-1) rho eInt + (gamma-1) P )
		// = 1/rho ( P + (gamma-1) P)
		// = gamma P / rho
		real vUxSq = vU.x * vU.x;
		real csSq = solver->heatCapacityRatio * P / (rho * h);
		real cs = sqrt(csSq);

		//Font 2008 eqn 106 -- matches calcDTCell
		real const betaUi = beta.s<?=side?>;
		real discr = sqrt((1. - vSq) * (gammaU.xx * (1. - vSq * csSq) - vUxSq * (1. - csSq)));
		real lambdaMin = (vU.x * (1. - csSq) - cs * discr) * alpha / (1. - vSq * csSq) - betaUi;
		real lambdaMax = (vU.x * (1. - csSq) + cs * discr) * alpha / (1. - vSq * csSq) - betaUi;

		real LambdaMin = (lambdaMin + betaUi) / alpha;	//2008 Font eqn 114
		real LambdaMax = (lambdaMax + betaUi) / alpha;	//2008 Font eqn 114
		
		//used by evL and evR
		real ATildeMinus = (gammaU.xx - vUxSq) / (gammaU.xx - vU.x * LambdaMin);	//2008 Font eqn 113
		real ATildePlus  = (gammaU.xx - vUxSq) / (gammaU.xx - vU.x * LambdaMax);	//2008 Font eqn 113
		
		//used by evL
		real VMinus = (vU.x - LambdaMin) / (gammaU.xx - vU.x * LambdaMin);	//2008 Font eqn 113
		real VPlus = (vU.x - LambdaMax) / (gammaU.xx - vU.x * LambdaMax);		//2008 Font eqn 113
	
		//used by evL and evR
		//hmm, should these be lower?  Time to derive the equations?
		real CMinus = vL.x - VMinus;	//2008 Font eqn 112
		real CPlus = vL.x - VPlus;		//2008 Font eqn 112

		real kappa = calc_dP_deInt(solver, rho, eInt);	//2008 Font note just after eqn 107
		real kappaTilde = kappa / rho;	//2008 Font eqn 112.  
		//used by evL and evR
		real Kappa = kappaTilde / (kappaTilde - csSq);	//2008 Font eqn 112.  
		//Kappa = h;	//approx for ideal gas
	
		int indexInt = side + dim * index;
		global <?=eigen_t?>* eig = eigenBuf + indexInt;
<?
for _,var in ipairs(eqn.eigenVars) do
?>	eig-><?=var.name?> = <?=var.name?>;
<? end
?>
	}<? end ?>
}


<? 
local prefix = require 'ext.table'.map(eqn.eigenVars, function(var)
	return '\t'..var.type..' '..var.name..' = eig.'..var.name..';\n'
end):concat()
?>

//// MODULE_NAME: <?=eigen_leftTransform?>

void <?=eigen_leftTransform?>(
	<?=waves_t?> * const Y,
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const eig,
	<?=cons_t?> const X_,
	real3 const x
) { 
	//rotate incoming v's in X
	//this should match calcEigenBasis
	//eig.beta and eig.gamma should already be rotated
	<? if side==0 then ?>
	<?=cons_t?> X = X_;
	<? elseif side == 1 then ?>
	<?=cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[2], -X_.ptr[1], X_.ptr[3], X_.ptr[4]}};
	<? elseif side == 2 then ?>
	<?=cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[3], X_.ptr[2], -X_.ptr[1], X_.ptr[4]}};
	<? end ?>
	
	<?=prefix?>
	
	real const det_gamma = real3s3_det(gamma);
	real3s3 const gammaU = real3s3_inv(gamma, det_gamma);

	real const vUxSq = vU.x * vU.x;
	real const hSq = h * h;
	real const hW = h * W;
	real const W2 = W * W;

	real const gamma_gammaUxx = det_gamma * gammaU.xx;
	real const gamma_gammaUxy = det_gamma * gammaU.xy;
	real const gamma_gammaUxz = det_gamma * gammaU.xz;
	real const xi = det_gamma * (gammaU.xx - vUxSq);//2008 Font eqn 121
	real const Delta = hSq * hW * (Kappa - 1.) * (CPlus - CMinus) * xi;	//2008 Font eqn 121
	
	//min row	2008 Font eqn 118
	real scale;
	scale = hSq / Delta;
	real l5minus = (1 - Kappa) * (-det_gamma * vU.x + VPlus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VPlus * xi;
	Y->ptr[0] = (
		X.ptr[0] * (hW * VPlus * xi + l5minus)
		+ X.ptr[1] * (gamma_gammaUxx * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * vU.x * xi - gamma_gammaUxx * vU.x))
		+ X.ptr[2] * (gamma_gammaUxy * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * vU.y * xi - gamma_gammaUxy * vU.x))
		+ X.ptr[3] * (gamma_gammaUxz * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * vU.z * xi - gamma_gammaUxz * vU.x))
		+ X.ptr[4] * l5minus
	) * scale;
	//mid normal row	2008 Font eqn 115
	scale = W / (Kappa - 1.);
	Y->ptr[1] = (
		X.ptr[0] * (h - W) 
		+ X.ptr[1] * (W * vU.x) 
		+ X.ptr[2] * (W * vU.y) 
		+ X.ptr[3] * (W * vU.z) 
		+ X.ptr[4] * (-W)
	) * scale;
	//mid tangent A row	2008 Font eqn 116
	scale = 1. / (h * xi);
	Y->ptr[2] = (
		X.ptr[0] * (-gamma.zz * vL.y + gamma.yz * vL.z) 
		+ X.ptr[1] * vU.x * (gamma.zz * vL.y - gamma.yz * vL.z)
		+ X.ptr[2] * (gamma.zz * (1. - vL.x * vU.x) + gamma.xz * vL.z * vU.x)
		+ X.ptr[3] * (-gamma.yz * (1. - vL.x * vU.x) - gamma.xz * vL.y * vU.x)
		+ X.ptr[4] * (-gamma.zz * vL.y + gamma.yz * vL.z)
	) * scale;
	//mid tangent B row	2008 Font eqn 117
	Y->ptr[3] = (
		X.ptr[0] * (-gamma.yy * vL.z + gamma.yz * vL.y)
		+ X.ptr[1] * vU.x * (gamma.yy * vL.z - gamma.yz * vL.y)
		+ X.ptr[2] * (-gamma.yz * (1. - vL.x * vU.x) - gamma.xy * vL.z * vU.x)
		+ X.ptr[3] * (gamma.yy * (1. - vL.x * vU.x) + gamma.xy * vL.y * vU.x)
		+ X.ptr[4] * (-gamma.yy * vL.z + gamma.yz * vL.y)
	) * scale;
	//max row	2008 Font eqn 118
	scale = -hSq / Delta;
	real const l5plus = (1 - Kappa) * (-det_gamma * vU.x + VMinus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VMinus * xi;
	Y->ptr[4] = (
		X.ptr[0] * (h * W * VMinus * xi + l5plus)
		+ X.ptr[1] * (gamma_gammaUxx * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * vU.x * xi - gamma_gammaUxx * vU.x))
		+ X.ptr[2] * (gamma_gammaUxy * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * vU.y * xi - gamma_gammaUxy * vU.x))
		+ X.ptr[3] * (gamma_gammaUxz * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * vU.z * xi - gamma_gammaUxz * vU.x))
		+ X.ptr[4] * l5plus
	) * scale;
}

//// MODULE_NAME: <?=eigen_rightTransform?>

void <?=eigen_rightTransform?>(
	<?=cons_t?> * const Y,
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const eig,
	<?=waves_t?> const X,
	real3 const x,
	normal_t const n
) {
	<?=prefix?>
	
	real hW = h * W;
	real W2 = W * W;

	//2008 Font eqns 108-111
	Y->ptr[0] = X.ptr[0]
		+ X.ptr[1] * (Kappa / hW)
		+ X.ptr[2] * (W * vL.y)
		+ X.ptr[3] * (W * vL.z)
		+ X.ptr[4];
	Y->ptr[1] = X.ptr[0] * (hW * CMinus)
		+ X.ptr[1] * (vL.x)
		+ X.ptr[2] * (h * (gamma.xy + 2. * W2 * vL.y * vL.x))
		+ X.ptr[3] * (h * (gamma.xz + 2. * W2 * vL.x * vL.z))
		+ X.ptr[4] * (hW * CPlus);
	Y->ptr[2] = X.ptr[0] * (hW * vL.y)
		+ X.ptr[1] * (vL.y)
		+ X.ptr[2] * (h * (gamma.yy + 2. * W2 * vL.y * vL.y))
		+ X.ptr[3] * (h * (gamma.yz + 2. * W2 * vL.y * vL.z))
		+ X.ptr[4] * (hW * vL.y);
	Y->ptr[3] = X.ptr[0] * (hW * vL.z)
		+ X.ptr[1] * (vL.z)
		+ X.ptr[2] * (h * (gamma.yz + 2. * W2 * vL.y * vL.z))
		+ X.ptr[3] * (h * (gamma.zz + 2. * W2 * vL.z * vL.z))
		+ X.ptr[4] * (hW * vL.z);
	Y->ptr[4] =X.ptr[0] * (hW * ATildeMinus - 1.)
		+ X.ptr[1] * (1. - Kappa / hW)
		+ X.ptr[2] * (W * vL.y * (2. * hW - 1.))
		+ X.ptr[3] * (W * vL.z * (2. * hW - 1.))
		+ X.ptr[4] * (hW * ATildePlus - 1.);
	
	//rotate outgoing y's x's into side
	real const tmp = Y->ptr[1];
	Y->ptr[1] = -Y->ptr[1+n->side];
	Y->ptr[1+n->side] = tmp;
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

void <?=eigen_fluxTransform?>(
	<?=cons_t?> * const result,
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const eig,
	<?=cons_t?> const X_,
	<?=cell_t?> const * const cell,
	<?=normal_t?> const n
) {
#if 0
	//rotate incoming v's in x
	<? if side==0 then ?>
	<?=cons_t?> X = X_;
	<? elseif side == 1 then ?>
	<?=cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[2], -X_.ptr[1], X_.ptr[3], X_.ptr[4]}};
	<? elseif side == 2 then ?>
	<?=cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[3], X_.ptr[2], -X_.ptr[1], X_.ptr[4]}};
	<? end ?>

	//TODO do the matrix multiply here

	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real tmp = Y.ptr[1];
	Y.ptr[1] = Y[1+<?=side?>];
	Y.ptr[1+<?=side?>] = tmp;
	<? end ?>
#else
	//default
	<?=waves_t?> waves;
	<?=eigen_leftTransform?>(&waves, solver, eig, X_, (cell)->pos, n);
	<?=eqn:eigenWaveCodePrefix("n", "eig", "(cell)->pos")?>
<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode("n", "eig", "(cell)->pos", j)?>;
<? end 
?>	eigen_rightTransform(result, solver, eig, waves, (cell)->pos, n);
#endif
}

//// MODULE_NAME: <?=addSource?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf<?=
	solver:getADMArgs()?>
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	<?=solver:getADMVarCode()?>
}

//// MODULE_NAME: <?=constrainU?>

/*
This is from 2008 Alcubierre eqn 7.3.11
Also in 2013 Rezzolla, Zanotti eqn 7.17-7.22
TODO verify this matches with the rest of the Font 2008 stuff
because I remember seeing multiple definitions of u^i from v^i in GRHD, 
between Marti & Muller, Font, Alcubierre, and Baumgarte & Shapiro

W = alpha u^0
v^i = (u^i / u^0 + beta^i) / alpha

therefore 
u^0 = W / alpha
u^i = W (v^i - beta^i / alpha)
W = sqrt(1 - v^i v^j gamma_ij)
*/
kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf<?=
	solver:getADMArgs()?>
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost - 1);
	
	<?=solver:getADMVarCode()?>
	
	real det_gamma = real3s3_det(gamma);
	real3s3 gammaU = real3s3_inv(gamma, det_gamma);

	global <?=eqn.cons_only_t?> const * const U = &UBuf[index].cons;
	
	U->D = max(U->D, (real)solver->DMin);
	U->tau = max(U->tau, (real)solver->tauMin);

	U->D = min(U->D, (real)solver->DMax);
	U->tau = min(U->tau, (real)solver->tauMax);


	real D = U->D;
	real3 S = U->S;
	real tau = U->tau;

	global <?=prim_t?>* prim = &UBuf[index].prim;
	real3 v = prim->v;

	real SLen = real3_weightedLen(S, gammaU);
	//TODO update this to work with GRHD
	real PMin = max(SLen - tau - D + SLen * solver->solvePrimVelEpsilon, solver->solvePrimPMinEpsilon);
	real PMax = (solver->heatCapacityRatio - 1.) * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real vLen = SLen / (tau + D + P);	//tau + D + P = rho h W^2
		real vSq = vLen * vLen;
		real W = 1. / sqrt(1. - vSq);
		real eInt = (tau + D * (1. - W) + P * (1. - W*W)) / (D * W);
		real rho = D / W;
		real f = (solver->heatCapacityRatio - 1.) * rho * eInt - P;
		real csSq = (solver->heatCapacityRatio - 1.) * (tau + D * (1. - W) + P) / (tau + D + P);
		real df_dP = vSq * csSq - 1.;
		real newP = P - f / df_dP;
		newP = max(newP, PMin);
		real PError = fabs(1. - newP / P);
		P = newP;
		if (PError < solver->solvePrimStopEpsilon) {
			v = real3_real_mul(S, 1. / (tau + D + P));
			vSq = real3_weightedLenSq(v, gamma);
			W = 1. / sqrt(1. - vSq);
			rho = D / W;
			rho = max(rho, (real)solver->rhoMin);
			rho = min(rho, (real)solver->rhoMax);
			eInt = P / (rho * (solver->heatCapacityRatio - 1.));
			eInt = min(eInt, (real)solver->eIntMax);
			*prim = (<?=prim_t?>){
				.rho = rho,
				.v = v,
				.eInt = eInt,
			};
//printf("cell %d finished with prims = %f %f %f\n", index, rho, v.x, eInt);
			return;
		}
	}
}
