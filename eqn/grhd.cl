/*
2008 Font "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity"
*/

<? 
local table = require 'ext.table'

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym
?>

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
kernel void calcDT(
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf<?=
	solver:getADMArgs()?>
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	<?=eqn.prim_t?> prim = UBuf[index].prim;
	<?=solver:getADMVarCode()?>

	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	real3 vU = sym3_real3_mul(gammaU, prim.v);

	real rho = prim.rho;
	real eInt = prim.eInt;
	//2008 Font Eqn 31: v^2 = gamma_ij v^i v^j
	real vSq = real3_dot(prim.v, vU);
	real P = calc_P(rho, eInt);
	real h = calc_h(rho, P, eInt);
	real csSq = heatCapacityRatio * P / (rho * h);
	real cs = sqrt(csSq);
	
	real dt = INFINITY;
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		//for the particular direction
		real vUi = vU.s<?=side?>;
		real vUiSq = vUi * vUi;
		
		//Font 2008 eqn 106
		const real betaUi = beta.s<?=side?>;
		real discr = sqrt((1. - vSq) * (gammaU.xx * (1. - vSq * csSq) - vUiSq * (1. - csSq)));
		real lambdaMin = (vUi * (1. - csSq) - cs * discr) * alpha / (1. - vSq * csSq) - betaUi;
		real lambdaMax = (vUi * (1. - csSq) + cs * discr) * alpha / (1. - vSq * csSq) - betaUi;
		lambdaMin = min((real)0., lambdaMin);
		lambdaMax = max((real)0., lambdaMax);
		dt = min(dt, (real)dx<?=side?>_at(i) / (fabs(lambdaMax - lambdaMin) + (real)1e-9));
	}<? end ?>
	
	dtBuf[index] = dt; 
}

<? if false then ?>
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	
	real vUi = gammaU.<?=sym(side+1,1)?> * U.prim.v.x
			+ gammaU.<?=sym(side+1,2)?> * U.prim.v.y
			+ gammaU.<?=sym(side+1,3)?> * U.prim.v.z;
	real vUi_shift = vUi - beta.s<?=side?> / alpha;

	<?=eqn.cons_t?> F;
	F.cons.D = U.cons.D * vUi_shift;
	F.S = real3_scale(U.cons.S, vUi_shift);
	F.S.s<?=side?> += U.prim.p;
	F.tau = U.cons.tau * vUi_shift + p * vUi;
	F.prim = (<?=eqn.prim_t?>){
		.rho = 0,
		.v = {.s={0,0,0}},
		.eInt = 0,
	};
	F.ePot = 0;

	return F;
}
<? end ?>
<? end ?>

//used by PLM
//TODO SRHD PLM needs to do this:
//1) calcLR for the <?=eqn.prim_t?> (that means put calcLR in its own file, and a new primLR buf)
//2) have a new kernel for calc consLR from primLR, since calcDeltaUEig and calcFlux both need this
//or does the eigenbasis need to be derived from the variables being transformed?
//shoud I PLM the U's then converge the prims ... and therefore track the prims on edges as well?
<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){};
}
<? end ?>


kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,
	
	//TODO 
	//turn this into a LR extrapolation
	//actually make use of PLM somehow 
	//right now only primBuf is being used for getting neighbor values
	//so SRHD should perform the PLM stuff on the primBuf instead of the UBUf?
	// or do the PLM on the UBuf and do the cons->prim on the ULR edge values
	const global <?=eqn.cons_t?>* UBuf<?=
	solver:getADMArgs()?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	int indexR = index;
	<?=eqn.prim_t?> primR = UBuf[indexR].prim;
	
	<?=solver:getADMVarCode{suffix='R'} --[[ produce alphaR, betaR, gammaR at indexR ]] ?>
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
		<?=eqn.prim_t?> primL = UBuf[indexL].prim;
	
		<?=solver:getADMVarCode{suffix='L'} --[[ produce alphaL, betaL, gammaL at indexL ]] ?>

<? if true then -- arithmetic averaging ?>
		<?=eqn.prim_t?> avg = (<?=eqn.prim_t?>){
			.rho = .5 * (primL.rho + primR.rho),
			.v = real3_scale(real3_add(primL.v, primR.v), .5),
			.eInt = .5 * (primL.eInt + primR.eInt),
		};
		real alpha = .5 * (alphaL + alphaR);
		real3 beta = real3_scale(real3_add(betaL, betaR), .5);
		sym3 gamma = sym3_scale(sym3_add(gammaL, gammaR), .5);
<? -- else -- Roe-averaging, Font 2008 eqn 38 ?>
<? end ?>
		
		real rho = avg.rho;
		real3 vL = avg.v;
		real eInt = avg.eInt;
			
		//these match eigen_leftTransform
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
		gamma = (sym3){
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
		gamma = (sym3){
			.xx = gamma.zz,
			.xy = gamma.yz,
			.xz = -gamma.xz,
			.yy = gamma.yy,
			.yz = -gamma.xy,
			.zz = gamma.xx,
		};
		<? end ?>
		
		real det_gamma = sym3_det(gamma);
		sym3 gammaU = sym3_inv(gamma, det_gamma);

		real3 vU = sym3_real3_mul(gammaU, vL);
		real vSq = real3_dot(vL, vU);
		real oneOverW2 = 1. - vSq;
		real oneOverW = sqrt(oneOverW2);
		real W = 1. / oneOverW;	//alpha?
		real W2 = 1. / oneOverW2;
		real P = (heatCapacityRatio - 1.) * rho * eInt;
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
		real csSq = heatCapacityRatio * P / (rho * h);
		real cs = sqrt(csSq);

		//Font 2008 eqn 106 -- matches calcDT
		const real betaUi = beta.s<?=side?>;
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

		real kappa = calc_dP_deInt(rho, eInt);	//2008 Font note just after eqn 107
		real kappaTilde = kappa / rho;	//2008 Font eqn 112.  
		//used by evL and evR
		real Kappa = kappaTilde / (kappaTilde - csSq);	//2008 Font eqn 112.  
		//Kappa = h;	//approx for ideal gas
	
		int indexInt = side + dim * index;
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
<?
for _,field in ipairs(eqn.eigenVars) do
	local name,ctype = next(field)
?>	eig-><?=name?> = <?=name?>;
<? end
?>
	}<? end ?>
}

<? for side=0,solver.dim-1 do 
	local prefix = table.map(eqn.eigenVars, function(field)
		local name,ctype = next(field)
		return '\t'..ctype..' '..name..' = eig.'..name..';\n'
	end):concat()
?>
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X_,
	real3 x
) { 
	<?=eqn.waves_t?> Y;

	//rotate incoming v's in X
	//this should match calcEigenBasis
	//eig.beta and eig.gamma should already be rotated
	<? if side==0 then ?>
	<?=eqn.cons_t?> X = X_;
	<? elseif side == 1 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[2], -X_.ptr[1], X_.ptr[3], X_.ptr[4]}};
	<? elseif side == 2 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[3], X_.ptr[2], -X_.ptr[1], X_.ptr[4]}};
	<? end ?>
	
	<?=prefix?>
	
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);

	real vUxSq = vU.x * vU.x;
	real hSq = h * h;
	real hW = h * W;
	real W2 = W * W;

	real gamma_gammaUxx = det_gamma * gammaU.xx;
	real gamma_gammaUxy = det_gamma * gammaU.xy;
	real gamma_gammaUxz = det_gamma * gammaU.xz;
	real xi = det_gamma * (gammaU.xx - vUxSq);//2008 Font eqn 121
	real Delta = hSq * hW * (Kappa - 1.) * (CPlus - CMinus) * xi;	//2008 Font eqn 121
	
	//min row	2008 Font eqn 118
	real scale;
	scale = hSq / Delta;
	real l5minus = (1 - Kappa) * (-det_gamma * vU.x + VPlus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VPlus * xi;
	Y.ptr[0] = (
		X.ptr[0] * (hW * VPlus * xi + l5minus)
		+ X.ptr[1] * (gamma_gammaUxx * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * vU.x * xi - gamma_gammaUxx * vU.x))
		+ X.ptr[2] * (gamma_gammaUxy * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * vU.y * xi - gamma_gammaUxy * vU.x))
		+ X.ptr[3] * (gamma_gammaUxz * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * vU.z * xi - gamma_gammaUxz * vU.x))
		+ X.ptr[4] * l5minus
	) * scale;
	//mid normal row	2008 Font eqn 115
	scale = W / (Kappa - 1.);
	Y.ptr[1] = (
		X.ptr[0] * (h - W) 
		+ X.ptr[1] * (W * vU.x) 
		+ X.ptr[2] * (W * vU.y) 
		+ X.ptr[3] * (W * vU.z) 
		+ X.ptr[4] * (-W)
	) * scale;
	//mid tangent A row	2008 Font eqn 116
	scale = 1. / (h * xi);
	Y.ptr[2] = (
		X.ptr[0] * (-gamma.zz * vL.y + gamma.yz * vL.z) 
		+ X.ptr[1] * vU.x * (gamma.zz * vL.y - gamma.yz * vL.z)
		+ X.ptr[2] * (gamma.zz * (1. - vL.x * vU.x) + gamma.xz * vL.z * vU.x)
		+ X.ptr[3] * (-gamma.yz * (1. - vL.x * vU.x) - gamma.xz * vL.y * vU.x)
		+ X.ptr[4] * (-gamma.zz * vL.y + gamma.yz * vL.z)
	) * scale;
	//mid tangent B row	2008 Font eqn 117
	Y.ptr[3] = (
		X.ptr[0] * (-gamma.yy * vL.z + gamma.yz * vL.y)
		+ X.ptr[1] * vU.x * (gamma.yy * vL.z - gamma.yz * vL.y)
		+ X.ptr[2] * (-gamma.yz * (1. - vL.x * vU.x) - gamma.xy * vL.z * vU.x)
		+ X.ptr[3] * (gamma.yy * (1. - vL.x * vU.x) + gamma.xy * vL.y * vU.x)
		+ X.ptr[4] * (-gamma.yy * vL.z + gamma.yz * vL.y)
	) * scale;
	//max row	2008 Font eqn 118
	scale = -hSq / Delta;
	real l5plus = (1 - Kappa) * (-det_gamma * vU.x + VMinus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VMinus * xi;
	Y.ptr[4] = (
		X.ptr[0] * (h * W * VMinus * xi + l5plus)
		+ X.ptr[1] * (gamma_gammaUxx * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * vU.x * xi - gamma_gammaUxx * vU.x))
		+ X.ptr[2] * (gamma_gammaUxy * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * vU.y * xi - gamma_gammaUxy * vU.x))
		+ X.ptr[3] * (gamma_gammaUxz * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * vU.z * xi - gamma_gammaUxz * vU.x))
		+ X.ptr[4] * l5plus
	) * scale;
	
	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=prefix?>
	
	real hW = h * W;
	real W2 = W * W;

	<?=eqn.cons_t?> Y;

	//2008 Font eqns 108-111
	Y.ptr[0] = X.ptr[0]
		+ X.ptr[1] * (Kappa / hW)
		+ X.ptr[2] * (W * vL.y)
		+ X.ptr[3] * (W * vL.z)
		+ X.ptr[4];
	Y.ptr[1] = X.ptr[0] * (hW * CMinus)
		+ X.ptr[1] * (vL.x)
		+ X.ptr[2] * (h * (gamma.xy + 2. * W2 * vL.y * vL.x))
		+ X.ptr[3] * (h * (gamma.xz + 2. * W2 * vL.x * vL.z))
		+ X.ptr[4] * (hW * CPlus);
	Y.ptr[2] = X.ptr[0] * (hW * vL.y)
		+ X.ptr[1] * (vL.y)
		+ X.ptr[2] * (h * (gamma.yy + 2. * W2 * vL.y * vL.y))
		+ X.ptr[3] * (h * (gamma.yz + 2. * W2 * vL.y * vL.z))
		+ X.ptr[4] * (hW * vL.y);
	Y.ptr[3] = X.ptr[0] * (hW * vL.z)
		+ X.ptr[1] * (vL.z)
		+ X.ptr[2] * (h * (gamma.yz + 2. * W2 * vL.y * vL.z))
		+ X.ptr[3] * (h * (gamma.zz + 2. * W2 * vL.z * vL.z))
		+ X.ptr[4] * (hW * vL.z);
	Y.ptr[4] =X.ptr[0] * (hW * ATildeMinus - 1.)
		+ X.ptr[1] * (1. - Kappa / hW)
		+ X.ptr[2] * (W * vL.y * (2. * hW - 1.))
		+ X.ptr[3] * (W * vL.z * (2. * hW - 1.))
		+ X.ptr[4] * (hW * ATildePlus - 1.);
	
	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real tmp = Y.ptr[1];
	Y.ptr[1] = -Y.ptr[1+<?=side?>];
	Y.ptr[1+<?=side?>] = tmp;
	<? end ?>

	return Y;
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X_,
	real3 x
) {
#if 0
	//rotate incoming v's in x
	<? if side==0 then ?>
	<?=eqn.cons_t?> X = X_;
	<? elseif side == 1 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[2], -X_.ptr[1], X_.ptr[3], X_.ptr[4]}};
	<? elseif side == 2 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[3], X_.ptr[2], -X_.ptr[1], X_.ptr[4]}};
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
	<?=eqn.waves_t?> waves = eigen_leftTransform_<?=side?>(eig, X_, x);
	<?=eqn:eigenWaveCodePrefix(side, 'eig', 'x')?>
<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode(side, 'eig', 'x', j)?>;
<? end 
?>	return eigen_rightTransform_<?=side?>(eig, waves, x);
#endif
}
<? end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf<?=
	solver:getADMArgs()?>
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=solver:getADMVarCode()?>
}

kernel void constrainU(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);

	global <?=eqn.cons_only_t?>* U = &UBuf[index].cons;
	
	U->D = max(U->D, (real)DMin);
	U->tau = max(U->tau, (real)tauMin);

	U->D = min(U->D, (real)DMax);
	U->tau = min(U->tau, (real)tauMax);
}

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
kernel void updatePrims(
	global <?=eqn.cons_t?>* UBuf<?=
	solver:getADMArgs()?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	<?=solver:getADMVarCode()?>
	
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);

	const global <?=eqn.cons_only_t?>* U = &UBuf[index].cons;
	real D = U->D;
	real3 S = U->S;
	real tau = U->tau;

	global <?=eqn.prim_t?>* prim = &UBuf[index].prim;
	real3 v = prim->v;

	real SLen = real3_weightedLen(S, gammaU);
	//TODO update this to work with GRHD
	real PMin = max(SLen - tau - D + SLen * solvePrimVelEpsilon, solvePrimPMinEpsilon);
	real PMax = (heatCapacityRatio - 1.) * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real vLen = SLen / (tau + D + P);	//tau + D + P = rho h W^2
		real vSq = vLen * vLen;
		real W = 1. / sqrt(1. - vSq);
		real eInt = (tau + D * (1. - W) + P * (1. - W*W)) / (D * W);
		real rho = D / W;
		real f = (heatCapacityRatio - 1.) * rho * eInt - P;
		real csSq = (heatCapacityRatio - 1.) * (tau + D * (1. - W) + P) / (tau + D + P);
		real df_dP = vSq * csSq - 1.;
		real newP = P - f / df_dP;
		newP = max(newP, PMin);
		real PError = fabs(1. - newP / P);
		P = newP;
		if (PError < solvePrimStopEpsilon) {
			v = real3_scale(S, 1. / (tau + D + P));
			vSq = real3_weightedLenSq(v, gamma);
			W = 1. / sqrt(1. - vSq);
			rho = D / W;
			rho = max(rho, (real)rhoMin);
			rho = min(rho, (real)rhoMax);
			eInt = P / (rho * (heatCapacityRatio - 1.));
			eInt = min(eInt, (real)eIntMax);
			*prim = (<?=eqn.prim_t?>){
				.rho = rho,
				.v = v,
				.eInt = eInt,
			};
//printf("cell %d finished with prims = %f %f %f\n", index, rho, v.x, eInt);
			return;
		}
	}
}
