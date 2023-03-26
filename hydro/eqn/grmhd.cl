//// MODULE_NAME: <?=eqn_common?>

//I'm going to fix metric coordinates at first
//then later the transition to the evolved metric will be easier
constant real const alpha = 1;
constant real3 const betaU = real3_zero;

//pressure function for ideal gas
real calc_P(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * rho * eInt;
}	

//chi in most papers
real calc_dP_drho(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * eInt;
}

//kappa in most papers
real calc_dP_deInt(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * rho;
}

real calc_eInt_from_P(real rho, real P) {
	return P / ((heatCapacityRatio - 1.) * rho);
}

real calc_h(real rho, real P, real eInt) {
	return 1. + eInt + P / rho;
}

<?=cons_t?> consFromPrim(<?=prim_t?> prim, real3 x) {
	real vSq = coordLenSq(prim.v, x);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 40-42:
	
	//rest-mass density = J^0 = rho u^0
	real D = prim.rho * W;	
	
	//momentum = T^0i = rho h u^0 u^i + P g^0i
	real3 S = real3_add(
		real3_real_mul(prim.v, prim.rho * h * WSq),
		real3_real_mul(betaU, P / (alpha * alpha)));
	
	//energy = T^00 = rho h u^0 u^0 + P g^00
	real tau = prim.rho * h * WSq - D - P / (alpha * alpha);
	
	return (<?=cons_t?>){.D=D, .S=S, .tau=tau};
}

//// MODULE_NAME: <?=applyInitCondCell?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const cons,
	global <?=cell_t?> const * const cell,
	global <?=prim_t?> * const prim
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(mins, maxs), .5);
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

	<?=code?>
	
	real eInt = calc_eInt_from_P(rho, P);
	real vSq = coordLenSq(v, x);
	real W = 1. / sqrt(1. - vSq);
	real h = calc_h(rho, P, eInt);

	*prim = (<?=prim_t?>){.rho=rho, .v=v, .eInt=eInt};
	*cons = consFromPrim(*prim, x);
}


//// MODULE_NAME: <?=fluxFromCons?>

void <?=fluxFromCons?>(
	<?=cons_t?> * const resultF,
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const U,
	<?=cell_t?> const * const cell,
	<?=normal_t?> const n
) {
	real vi = W->v.s<?=side?>;
	real vi_shift = vi - betaU.s<?=side?> / alpha;

	/* 2008 Font eqn 34 */
	(resultF)->D = U->D * vi_shift;
	(resultF)->S = real3_real_mul(U->S, vi_shift);
	(resultF)->S.s<?=side?> += W->p;
	(resultF)->tau = U->tau * vi_shift + p * vi;
}

//// MODULE_NAME: calcDTCell

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
#error replace this with cons_t using the prim_only_t 
void calcDTCell(
	global real * const dt,
	constant <?=solver_t?> const * const solver,
	global <?=prim_t?> const * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3s3 const gammaU = coord_g_uu(x);

	<?=prim_only_t?> primOnly;
	primOnlyFromCons(&primOnly, solver, U, x);
	real const rho = prim.rho;
	real const eInt = prim.eInt;
	real const vSq = coordLenSq(prim.v, x);
	real const P = calc_P(rho, eInt);
	real const h = calc_h(rho, P, eInt);
	real const csSq = heatCapacityRatio * P / (rho * h);
	real const cs = sqrt(csSq);
	
	real dt = INFINITY;
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		//for the particular direction
		real const vi = prim.v.s<?=side?>;
		real const viSq = vi * vi;
		
		// Marti 1998 eqn 19
		// also Marti & Muller 2008 eqn 68
		// also Font 2008 eqn 106
		real const betaUi = betaU.s<?=side?>;
		real const discr = sqrt((1. - vSq) * (gammaU.xx * (1. - vSq * csSq) - viSq * (1. - csSq)));
		real const lambdaMin = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq) * alpha - betaUi;
		real const lambdaMax = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq) * alpha - betaUi;
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
	
	dtBuf[index] = dt; 
}

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

//// MODULE_NAME: calcEigenBasis

#error calcEigenBasis has been removed, and eigen_t structs are now calculated inline ... soooo ... convert this to something compatible
kernel void calcEigenBasis(
	constant <?=solver_t?> const * const solver,
	global real * const waveBuf,
	global <?=eigen_t?> * const eigenBuf,
	
	//TODO 
	//turn this into a LR extrapolation
	//actually make use of PLM somehow 
	//right now only primBuf is being used for getting neighbor values
	//so SRHD should perform the PLM stuff on the primBuf instead of the UBUf?
	// or do the PLM on the UBuf and do the cons->prim on the ULR edge values
	global <?=prim_t?> const * const primBuf	
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost - 1);
	real3 const x = cellBuf[index].pos;
	
	int const indexR = index;
	<?=prim_t?> primR = primBuf[indexR];
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		
		int const indexL = index - solver->stepsize.s<?=side?>;
		<?=prim_t?> const primL = primBuf[indexL];
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real3s3 gammaU = coord_g_uu(xInt);

<? if true then -- arithmetic averaging ?>
		<?=prim_t?> avg = (<?=prim_t?>){
			.rho = .5 * (primL.rho + primR.rho),
			.v = real3_real_mul(real3_add(primL.v, primR.v), .5),
			.eInt = .5 * (primL.eInt + primR.eInt),
		};
<? -- else -- Roe-averaging, Font 2008 eqn 38 ?>
<? end ?>

		real const rho = avg.rho;
		real3 const v = avg.v;
		real const eInt = avg.eInt;
	
		<? if side == 1 then ?>
		v = real3(v.y, -v.x, v.z);	// -90' rotation to put the y axis contents into the x axis
		<? elseif side == 2 then ?>
		v = real3(v.z, v.y, -v.x);	//-90' rotation to put the z axis in the x axis
		<? end ?>

//TODO NOTE if you're swapping vector components, you have to swap metric components too 

		real3 const vL = coord_lower(v, xInt);
		real const vSq = real3_dot(v, vL);
		real const oneOverW2 = 1. - vSq;
		real const oneOverW = sqrt(oneOverW2);
		real const W = 1. / oneOverW;
		real const W2 = 1. / oneOverW2;
		real const P = (heatCapacityRatio - 1.) * rho * eInt;
		real const h = 1. + eInt + P / rho;

		real const hW = h * W;

		//just after 2008 Font eqn 107:
		//h cs^2 = chi + P / rho^2 kappa = dp/drho + p / rho^2 dp/deInt
		// = (gamma-1) eInt + P/rho^2 (gamma-1) rho  for an ideal gas
		// = (gamma-1) (eInt + P/rho)
		// = 1/rho ( (gamma-1) rho eInt + (gamma-1) P )
		// = 1/rho ( P + (gamma-1) P)
		// = gamma P / rho
		real const vxSq = v.x * v.x;
		real const csSq = heatCapacityRatio * P / (rho * h);
		real const cs = sqrt(csSq);

		real const betaUi = betaU.s<?=side?>;
		real const discr = sqrt((1. - vSq) * ((1. - vSq * csSq) - vxSq * (1. - csSq)));
		real const lambdaMin = (v.x * (1. - csSq) - cs * discr) / (1. - vSq * csSq) * alpha * alpha - betaUi;
		real const lambdaMax = (v.x * (1. - csSq) + cs * discr) / (1. - vSq * csSq) * alpha * alpha - betaUi;

		int indexInt = side + dim * index;	
		global real* wave = waveBuf + numWaves * indexInt;
		wave[0] = lambdaMin;
		wave[1] = v.x * alpha - betaUi;
		wave[2] = v.x * alpha - betaUi;
		wave[3] = v.x * alpha - betaUi;
		wave[4] = lambdaMax;

		real const LambdaMin = (lambdaMin + betaUi) / alpha;	//2008 Font eqn 114
		real const LambdaMax = (lambdaMax + betaUi) / alpha;	//2008 Font eqn 114
		
		//used by evL and evR
		real const ATildeMinus = (gammaU.xx - vxSq) / (gammaU.xx - v.x * LambdaMin);	//2008 Font eqn 113
		real const ATildePlus  = (gammaU.xx - vxSq) / (gammaU.xx - v.x * LambdaMax);	//2008 Font eqn 113
		
		//used by evL
		real const VMinus = (v.x - LambdaMin) / (gammaU.xx - v.x * LambdaMin);	//2008 Font eqn 113
		real const VPlus = (v.x - LambdaMax) / (gammaU.xx - v.x * LambdaMax);	//2008 Font eqn 113
	
		//used by evL and evR
		real const CMinus = vL.x - VMinus;	//2008 Font eqn 112
		real const CPlus = vL.x - VPlus;	//2008 Font eqn 112

		real const kappa = calc_dP_deInt(rho, eInt);	//2008 Font note just after eqn 107
		real const kappaTilde = kappa / rho;	//2008 Font eqn 112.  
		//used by evL and evR
		real const Kappa = kappaTilde / (kappaTilde - csSq);	//2008 Font eqn 112.  
		//Kappa = h;	//approx for ideal gas
		
		global <?=eigen_t?> * const eig = eigenBuf + indexInt;	

<?
for _,field in ipairs(eqn.eigenVars) do
	local name,ctype = next(field)
?>
		eig-><?=name?> = <?=name?>;
<? end ?>

	}<? end ?>
}

//// MODULE_NAME: eigen_leftTransform
//// MODULE_DEPENDS: <?=coord_sqrt_det_g?>

void eigen_leftTransform(
	constant <?=solver_t?> const * const solver,
	real * const Y,
	<?=eigen_t?> const * const eig,
	real const * const X_,
	real3 const x
) { 
	//rotate incoming v's in X
	//TODO do the same for gamma_ij
	<? if side==0 then ?>
	<?=addr2?> const real* X = X_;
	<? elseif side == 1 then ?>
	real X[5] = {X_[0], X_[2], -X_[1], X_[3], X_[4]};
	<? elseif side == 2 then ?>
	real X[5] = {X_[0], X_[3], X_[2], -X_[1], X_[4]};
	<? end ?>

	<?=prefix?>
	real const gammaDet = coord_sqrt_det_g(solver, x);
	real3s3 const gammaL = coord_g_ll(x);
	real3s3 const gammaU = coord_g_uu(x);

	real3 const vL = coord_lower(v, x);
	real const vxSq = v.x * v.x;
	real const hSq = h * h;
	real const hW = h * W;
	real const W2 = W * W;

	real const gamma_gammaUxx = gammaDet * gammaU.xx;
	real const gamma_gammaUxy = gammaDet * gammaU.xy;
	real const gamma_gammaUxz = gammaDet * gammaU.xz;
	real const xi = gammaDet * (gammaU.xx - vxSq);//2008 Font eqn 121
	real const Delta = hSq * hW * (Kappa - 1.) * (CPlus - CMinus) * xi;	//2008 Font eqn 121
	
	//min row	2008 Font eqn 118
	real scale;
	scale = hSq / Delta;
	real l5minus = (1 - Kappa) * (-gammaDet * v.x + VPlus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VPlus * xi;
	Y[0] = (
		X[0] * (hW * VPlus * xi + l5minus)
		+ X[1] * (gamma_gammaUxx * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * v.x * xi - gamma_gammaUxx * v.x))
		+ X[2] * (gamma_gammaUxy * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * v.y * xi - gamma_gammaUxy * v.x))
		+ X[3] * (gamma_gammaUxz * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * v.z * xi - gamma_gammaUxz * v.x))
		+ X[4] * l5minus
	) * scale;
	//mid normal row	2008 Font eqn 115
	scale = W / (Kappa - 1.);
	Y[1] = (
		X[0] * (h - W) 
		+ X[1] * (W * v.x) 
		+ X[2] * (W * v.y) 
		+ X[3] * (W * v.z) 
		+ X[4] * (-W)
	) * scale;
	//mid tangent A row	2008 Font eqn 116
	scale = 1. / (h * xi);
	Y[2] = (
		X[0] * (-gammaL.zz * vL.y + gammaL.yz * vL.z) 
		+ X[1] * v.x * (gammaL.zz * vL.y - gammaL.yz * vL.z)
		+ X[2] * (gammaL.zz * (1. - v.x * vL.x) + gammaL.xz * vL.z * v.x)
		+ X[3] * (-gammaL.yz * (1. - vL.x * v.x) - gammaL.xz * vL.y * v.x)
		+ X[4] * (-gammaL.zz * vL.y + gammaL.yz * vL.z)
	) * scale;
	//mid tangent B row	2008 Font eqn 117
	Y[3] = (
		X[0] * (-gammaL.yy * vL.z + gammaL.yz * vL.y)
		+ X[1] * v.x * (gammaL.yy * vL.z - gammaL.yz * vL.y)
		+ X[2] * (-gammaL.yz * (1. - vL.x * v.x) - gammaL.xy * vL.z * v.x)
		+ X[3] * (gammaL.yy * (1. - vL.x * v.x) + gammaL.xy * vL.y * v.x)
		+ X[4] * (-gammaL.yy * vL.z + gammaL.yz * vL.y)
	) * scale;
	//max row	2008 Font eqn 118
	scale = -hSq / Delta;
	real l5plus = (1 - Kappa) * (-gammaDet * v.x + VMinus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VMinus * xi;
	Y[4] = (
		X[0] * (h * W * VMinus * xi + l5plus)
		+ X[1] * (gamma_gammaUxx * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * v.x * xi - gamma_gammaUxx * v.x))
		+ X[2] * (gamma_gammaUxy * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * v.y * xi - gamma_gammaUxy * v.x))
		+ X[3] * (gamma_gammaUxz * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * v.z * xi - gamma_gammaUxz * v.x))
		+ X[4] * l5plus
	) * scale;
}

//// MODULE_NAME: eigen_rightTransform
//// MODULE_DEPENDS: <?=coord_sqrt_det_g?>

void eigen_rightTransform(
	constant <?=solver_t?> const * const solver,
	real * const Y,
	<?=eigen_t?> const * const eig,
	real const * const X,
	real3 const x
) {
	<?=prefix?>
	real3s3 const gammaL = coord_g_ll(x);
	
	real3 const vL = coord_lower(v, x);
	real const hW = h * W;
	real const W2 = W * W;

	//2008 Font eqns 108-111
	Y[0] = X[0]
		+ X[1] * (Kappa / hW)
		+ X[2] * (W * vL.y)
		+ X[3] * (W * vL.z)
		+ X[4];
	Y[1] = X[0] * (hW * CMinus)
		+ X[1] * (vL.x)
		+ X[2] * (h * (gammaL.xy + 2. * W2 * vL.y * vL.x))
		+ X[3] * (h * (gammaL.xz + 2. * W2 * vL.x * vL.z))
		+ X[4] * (hW * CPlus);
	Y[2] = X[0] * (hW * vL.y)
		+ X[1] * (vL.y)
		+ X[2] * (h * (gammaL.yy + 2. * W2 * vL.y * vL.y))
		+ X[3] * (h * (gammaL.yz + 2. * W2 * vL.y * vL.z))
		+ X[4] * (hW * vL.y);
	Y[3] = X[0] * (hW * vL.z)
		+ X[1] * (vL.z)
		+ X[2] * (h * (gammaL.yz + 2. * W2 * vL.y * vL.z))
		+ X[3] * (h * (gammaL.zz + 2. * W2 * vL.z * vL.z))
		+ X[4] * (hW * vL.z);
	Y[4] =X[0] * (hW * ATildeMinus - 1.)
		+ X[1] * (1. - Kappa / hW)
		+ X[2] * (W * vL.y * (2. * hW - 1.))
		+ X[3] * (W * vL.z * (2. * hW - 1.))
		+ X[4] * (hW * ATildePlus - 1.);
	
	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real const tmp = Y[1];
	Y[1] = -Y[1+<?=side?>];
	Y[1+<?=side?>] = tmp;
	<? end ?>
}

//// MODULE_NAME: eigen_fluxTransform

void eigen_fluxTransform(
	constant <?=solver_t?> const * const solver,
	real * const Y,
	<?=eigen_t?> const * const eig,
	real const * const X_,
	real3 const x
) {
	//rotate incoming v's in x
	<? if side==0 then ?>
	<?=addr2?> const real* X = X_;
	<? elseif side == 1 then ?>
	real X[5] = {X_[0], X_[2], -X_[1], X_[3], X_[4]};
	<? elseif side == 2 then ?>
	real X[5] = {X_[0], X_[3], X_[2], -X_[1], X_[4]};
	<? end ?>

	//TODO do the matrix multiply here

	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real const tmp = Y[1];
	Y[1] = Y[1+<?=side?>];
	Y[1+<?=side?>] = tmp;
	<? end ?>
}

//// MODULE_NAME: constrainU

kernel void constrainU(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf
) {
	<?=SETBOUNDS?>(0,0);

	global <?=cons_t?> * const U = UBuf + index;
	
	U->D = max(U->D, (real)DMin);
	U->tau = max(U->tau, (real)tauMin);

	U->D = min(U->D, (real)DMax);
	U->tau = min(U->tau, (real)tauMax);
}

//// MODULE_NAME: updatePrims

//TODO update to include alphas, betas, and gammas
kernel void updatePrims(
	constant <?=solver_t?> const * const solver,
	global <?=prim_t?> * const primBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost - 1);
	real3 const x = cellBuf[index].pos;

	global <?=cons_t?> const * const U = UBuf + index;
	real const D = U->D;
	real3 const S = U->S;
	real const tau = U->tau;

	global <?=prim_t?> * const prim = primBuf + index;
	real3 v = prim->v;

	real const SLen = coordLen(S, x);
	real const PMin = max(SLen - tau - D + SLen * solvePrimVelEpsilon, solvePrimPMinEpsilon);
	real PMax = (heatCapacityRatio - 1.) * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real const vLen = SLen / (tau + D + P);
		real vSq = vLen * vLen;
		real W = 1. / sqrt(1. - vSq);
		real eInt = (tau + D * (1. - W) + P * (1. - W*W)) / (D * W);
		real rho = D / W;
		real const f = (heatCapacityRatio - 1.) * rho * eInt - P;
		real const csSq = (heatCapacityRatio - 1.) * (tau + D * (1. - W) + P) / (tau + D + P);
		real const df_dP = vSq * csSq - 1.;
		real newP = P - f / df_dP;
		newP = max(newP, PMin);
		real const PError = fabs(1. - newP / P);
		P = newP;
		if (PError < solvePrimStopEpsilon) {
			v = real3_real_mul(S, 1. / (tau + D + P));
			vSq = coordLenSq(v, x);
			W = 1. / sqrt(1. - vSq);
			rho = D / W;
			rho = max(rho, (real)rhoMin);
			rho = min(rho, (real)rhoMax);
			eInt = P / (rho * (heatCapacityRatio - 1.));
			eInt = min(eInt, (real)eIntMax);
			prim->rho = rho;
			prim->v = v;
			prim->eInt = eInt;
//printf("cell %d finished with prims = %f %f %f\n", index, rho, v.x, eInt);
			return;
		}
	}
}
