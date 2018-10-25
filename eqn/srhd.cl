/*
Marti & Muller 2008
Marti 1998
Font "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity" 2008 
*/

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);

	<?=eqn.prim_t?> prim = UBuf[index].prim;
	real rho = prim.rho;
	real eInt = prim.eInt;
	real vSq = coordLenSq(prim.v, x);
	real P = calc_P(solver, rho, eInt);
	real h = calc_h(rho, P, eInt);
	real csSq = solver->heatCapacityRatio * P / (rho * h);
	real cs = sqrt(csSq);
	
	real dt = INFINITY;
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		//for the particular direction
		real vi = prim.v.s<?=side?>;
		real viSq = vi * vi;
		
		// Marti 1998 eqn 19
		// also Marti & Muller 2008 eqn 68
		// also Font 2008 eqn 106
		real discr = sqrt((1. - vSq) * (1. - vSq * csSq - viSq * (1. - csSq)));
		real lambdaMin = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq);
		real lambdaMax = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq);
		lambdaMin = min((real)0., lambdaMin);
		lambdaMax = max((real)0., lambdaMax);
		dt = min(dt, (real)dx<?=side?>_at(i) / (fabs(lambdaMax - lambdaMin) + (real)1e-9));
	}<? end ?>
	
	dtBuf[index] = dt; 
}

<? if false then ?>
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real vi = U.prim.v.s<?=side?>;
	real P = calc_P(solver, U.prim.rho, U.prim.eInt);

	<?=eqn.cons_t?> F;
	F.cons.D = U.cons.D * vi;
	F.cons.S = real3_real_mul(U.cons.S, vi);
	F.cons.S.s<?=side?> += P;
	F.cons.tau = U.cons.tau * vi + P * vi;
	
	//make sure the rest is zero ...
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

/*
used by PLM
TODO SRHD PLM needs to do this:
1) calcLR for the <?=eqn.prim_t?> (that means put calcLR in its own file, and a new primLR buf)
2) have a new kernel for calc consLR from primLR, since calcDeltaUEig and calcFlux both need this
or does the eigenbasis need to be derived from the variables being transformed?
shoud I PLM the U's then converge the prims ... and therefore track the prims on edges as well?

NOTICE this is only going to use U->prim
... but that won't help the PLM, since it operates based on numIntStates
... and numIntStates only covers the cons_t vars ...
which means this function won't work with the PLM code
*/
<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){};
}
<? end ?>

kernel void calcEigenBasis(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.eigen_t?>* eigenBuf,
	
	//TODO 
	//turn this into a LR extrapolation
	//actually make use of PLM somehow 
	//right now only primBuf is being used for getting neighbor values
	//so SRHD should perform the PLM stuff on the primBuf instead of the UBUf?
	// or do the PLM on the UBuf and do the cons->prim on the ULR edge values
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;
	<?=eqn.prim_t?> primR = UBuf[indexR].prim;
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
		<?=eqn.prim_t?> primL = UBuf[indexL].prim;
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

<? if true then -- arithmetic averaging ?>
		<?=eqn.prim_t?> avg = (<?=eqn.prim_t?>){
			.rho = .5 * (primL.rho + primR.rho),
			.v = real3_real_mul(real3_add(primL.v, primR.v), .5),
			.eInt = .5 * (primL.eInt + primR.eInt),
		};
<? -- else -- Roe-averaging, Font 2008 eqn 38 ?>
<? end ?>

		avg.v = real3_swap<?=side?>(avg.v);

		real rho = avg.rho;
		real3 v = avg.v;
		real eInt = avg.eInt;

//TODO NOTE if you're swapping vector components, you have to swap metric components too 
		real3 vL = coord_lower(v, xInt);
		real vSq = real3_dot(v, vL);
		real oneOverW2 = 1. - vSq;
		real oneOverW = sqrt(oneOverW2);
		real W = 1. / oneOverW;
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
		real vxSq = v.x * v.x;
		real csSq = solver->heatCapacityRatio * P / (rho * h);
		real cs = sqrt(csSq);

		real discr = sqrt((1. - vSq) * ((1. - vSq * csSq) - vxSq * (1. - csSq)));
		real lambdaMin = (v.x * (1. - csSq) - cs * discr) / (1. - vSq * csSq);
		real lambdaMax = (v.x * (1. - csSq) + cs * discr) / (1. - vSq * csSq);

		//used by evL and evR
		real ATildeMinus = (1. - vxSq) / (1. - v.x * lambdaMin);	//2008 Font eqn 113
		real ATildePlus  = (1. - vxSq) / (1. - v.x * lambdaMax);	//2008 Font eqn 113
		
		//used by evL
		real VMinus = (v.x - lambdaMin) / (1. - v.x * lambdaMin);	//2008 Font eqn 113
		real VPlus = (v.x - lambdaMax) / (1. - v.x * lambdaMax);	//2008 Font eqn 113
	
		//used by evL and evR
		real CMinus = vL.x - VMinus;	//2008 Font eqn 112
		real CPlus = vL.x - VPlus;	//2008 Font eqn 112

		real kappa = calc_dP_deInt(solver, rho, eInt);	//2008 Font note just after eqn 107
		real kappaTilde = kappa / rho;	//2008 Font eqn 112.  
		//used by evL and evR
		real Kappa = kappaTilde / (kappaTilde - csSq);	//2008 Font eqn 112.  
		//Kappa = h;	//approx for ideal gas
		
		int indexInt = side + dim * index;	
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;	

<?
for _,field in ipairs(eqn.eigenVars) do
	local name,ctype = next(field)
?>
		eig-><?=name?> = <?=name?>;
<? end ?>

	}<? end ?>
}

<? -- create code to initialize local vars of all the eig vars
local eigVarCode = require 'ext.table'.map(eqn.eigenVars, function(field)
	local name,ctype = next(field)
	return '\t'..ctype..' '..name..' = eig.'..name..';\n'
end):concat()
?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X_,
	real3 x
) { 
	<?=eqn.waves_t?> Y;

	//rotate incoming v's in X
	<? if side==0 then ?>
	<?=eqn.cons_t?> X = X_;
	<? elseif side == 1 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[2], X_.ptr[1], X_.ptr[3], X_.ptr[4]}};
	<? elseif side == 2 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[3], X_.ptr[2], X_.ptr[1], X_.ptr[4]}};
	<? end ?>

	<?=eigVarCode?>

	real3 vL = coord_lower(v, x);
	real vxSq = v.x * v.x;
	real hSq = h * h;
	real hW = h * W;
	real W2 = W * W;

	real xi = 1. - vxSq;	//2008 Font eqn 121
	real Delta = hSq * hW * (Kappa - 1.) * (CPlus - CMinus) * xi;	//2008 Font eqn 121
	
	//min row	2008 Font eqn 118
	real scale;
	scale = hSq / Delta;
	real l5minus = (1 - Kappa) * (-v.x + VPlus * (W2 * xi - 1.)) - Kappa * W2 * VPlus * xi;
	Y.ptr[0] = (
		X.ptr[0] * (hW * VPlus * xi + l5minus)
		+ X.ptr[1] * (1 - Kappa * ATildePlus + (2. * Kappa - 1.) * VPlus * (W2 * v.x * xi - v.x))
		+ X.ptr[2] * ((2. * Kappa - 1.) * VPlus * W2 * v.y * xi)
		+ X.ptr[3] * ((2. * Kappa - 1.) * VPlus * W2 * v.z * xi)
		+ X.ptr[4] * l5minus
	) * scale;
	//mid normal row	2008 Font eqn 115
	scale = W / (Kappa - 1.);
	Y.ptr[1] = (
		X.ptr[0] * (h - W) 
		+ X.ptr[1] * (W * v.x) 
		+ X.ptr[2] * (W * v.y) 
		+ X.ptr[3] * (W * v.z) 
		+ X.ptr[4] * (-W)
	) * scale;
	//mid tangent A row	2008 Font eqn 116
	scale = 1. / (h * xi);
	Y.ptr[2] = (
		X.ptr[0] * -vL.y
		+ X.ptr[1] * v.x * vL.y
		+ X.ptr[2] * ((1. - v.x * vL.x))
		+ X.ptr[4] * -vL.y
	) * scale;
	//mid tangent B row	2008 Font eqn 117
	Y.ptr[3] = (
		X.ptr[0] * -vL.z
		+ X.ptr[1] * v.x * vL.z
		+ X.ptr[3] * (1. - vL.x * v.x)
		+ X.ptr[4] * -vL.z
	) * scale;
	//max row	2008 Font eqn 118
	scale = -hSq / Delta;
	real l5plus = (1 - Kappa) * (-v.x + VMinus * (W2 * xi - 1.)) - Kappa * W2 * VMinus * xi;
	Y.ptr[4] = (
		X.ptr[0] * (h * W * VMinus * xi + l5plus)
		+ X.ptr[1] * (1 - Kappa * ATildeMinus + (2. * Kappa - 1.) * VMinus * (W2 * v.x * xi - v.x))
		+ X.ptr[2] * ((2. * Kappa - 1.) * VMinus * W2 * v.y * xi)
		+ X.ptr[3] * ((2. * Kappa - 1.) * VMinus * W2 * v.z * xi)
		+ X.ptr[4] * l5plus
	) * scale;
	
	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=eigVarCode?>
	
	real3 vL = coord_lower(v, x);
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
		+ X.ptr[2] * (h * 2. * W2 * vL.y * vL.x)
		+ X.ptr[3] * (h * 2. * W2 * vL.x * vL.z)
		+ X.ptr[4] * (hW * CPlus);
	Y.ptr[2] = X.ptr[0] * (hW * vL.y)
		+ X.ptr[1] * (vL.y)
		+ X.ptr[2] * (h * (1. + 2. * W2 * vL.y * vL.y))
		+ X.ptr[3] * (h * (2. * W2 * vL.y * vL.z))
		+ X.ptr[4] * (hW * vL.y);
	Y.ptr[3] = X.ptr[0] * (hW * vL.z)
		+ X.ptr[1] * (vL.z)
		+ X.ptr[2] * (h * (2. * W2 * vL.y * vL.z))
		+ X.ptr[3] * (h * (1. + 2. * W2 * vL.z * vL.z))
		+ X.ptr[4] * (hW * vL.z);
	Y.ptr[4] =X.ptr[0] * (hW * ATildeMinus - 1.)
		+ X.ptr[1] * (1. - Kappa / hW)
		+ X.ptr[2] * (W * vL.y * (2. * hW - 1.))
		+ X.ptr[3] * (W * vL.z * (2. * hW - 1.))
		+ X.ptr[4] * (hW * ATildePlus - 1.);
	
	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real tmp = Y.ptr[1];
	Y.ptr[1] = Y.ptr[1+<?=side?>];
	Y.ptr[1+<?=side?>] = tmp;
	<? end ?>

	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X_,
	real3 x
) {
#if 0
	//rotate incoming v's in x
	<? if side==0 then ?>
	<?=eqn.cons_t?> X = X_;
	<? elseif side == 1 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[2], X_.ptr[1], X_.ptr[3], X_.ptr[4]}};
	<? elseif side == 2 then ?>
	<?=eqn.cons_t?> X = {.ptr={X_.ptr[0], X_.ptr[3], X_.ptr[2], X_.ptr[1], X_.ptr[4]}};
	<? end ?>

	//TODO do the matrix multiply here

	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real tmp = Y.ptr[1];
	Y.ptr[1] = Y.ptr[1+<?=side?>];
	Y.ptr[1+<?=side?>] = tmp;
	<? end ?>
#else
	//default
	<?=eqn.waves_t?> waves = eigen_leftTransform_<?=side?>(solver, eig, X_, x);
	<?=eqn:eigenWaveCodePrefix(side, 'eig', 'x')?>
<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode(side, 'eig', 'x', j)?>;
<? end 
?>	return eigen_rightTransform_<?=side?>(solver, eig, waves, x);
#endif
}
<? end ?>

kernel void constrainU(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);

	global <?=eqn.cons_only_t?>* U = &UBuf[index].cons;
	
	U->D = max(U->D, (real)solver->DMin);
	U->tau = max(U->tau, (real)solver->tauMin);

	U->D = min(U->D, (real)solver->DMax);
	U->tau = min(U->tau, (real)solver->tauMax);
}

kernel void updatePrims(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);

	const global <?=eqn.cons_only_t?>* U = &UBuf[index].cons;
	real D = U->D;
	real3 S = U->S;
	real tau = U->tau;

	global <?=eqn.prim_t?>* prim = &UBuf[index].prim;
	real3 v = prim->v;

	real SLen = coordLen(S, x);
	real PMin = max(SLen - tau - D + SLen * solver->solvePrimVelEpsilon, solver->solvePrimPMinEpsilon);
	real PMax = (solver->heatCapacityRatio - 1.) * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real vLen = SLen / (tau + D + P);
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
			vSq = coordLenSq(v, x);
			W = 1. / sqrt(1. - vSq);
			rho = D / W;
			rho = max(rho, (real)solver->rhoMin);
			rho = min(rho, (real)solver->rhoMax);
			eInt = P / (rho * (solver->heatCapacityRatio - 1.));
			eInt = min(eInt, (real)solver->eIntMax);
			prim->rho = rho;
			prim->v = v;
			prim->eInt = eInt;
//printf("cell %d finished with prims = %f %f %f\n", index, rho, v.x, eInt);
			return;
		}
	}
}
