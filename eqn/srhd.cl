/*
Font "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity" 2008 
*/

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
kernel void calcDT(
	global real* dtBuf,
	
	//TODO actually make use of PLM somehow 
	//right now only primBuf is being used for getting neighbor values
	//so SRHD should perform the PLM stuff on the primBuf instead of the UBUf?
	// or do the PLM on the UBuf and do the cons->prim on the ULR edge values
	<?= solver.getULRArg ?>,	
	
	const global prim_t* primBuf
) {
	SETBOUNDS(0,0);
	if (OOB(2,2)) {
		dtBuf[index] = INFINITY;
		return;
	}
	
	prim_t prim = primBuf[index];
	real rho = prim.rho;
	real eInt = prim.eInt;
	real vSq = metricLenSq(prim.v);
	real P = calc_P(rho, eInt);
	real h = calc_h(rho, P, eInt);
	real csSq = heatCapacityRatio * P / (rho * h);
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
		const real betaUi = betaU.s<?=side?>;
		real discr = sqrt((1. - vSq) * (gammaU.xx * (1. - vSq * csSq) - viSq * (1. - csSq)));
		real lambdaMin = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq) * alpha - betaUi;
		real lambdaMax = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq) * alpha - betaUi;
		lambdaMin = min((real)0., lambdaMin);
		lambdaMax = max((real)0., lambdaMax);
		dt = min(dt, dx<?=side?>_at(i) / (fabs(lambdaMax - lambdaMin) + (real)1e-9));
	}<? end ?>
	
	dtBuf[index] = dt; 
}

kernel void calcEigenBasis(
	global real* waveBuf,
	global eigen_t* eigenBuf,
	//TODO turn this into a LR extrapolation
	const global prim_t* primBuf	
) {
	SETBOUNDS(2,1);
	
	int indexR = index;
	prim_t primR = primBuf[indexR];
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize[side];
		prim_t primL = primBuf[indexL];

<? if true then -- arithmetic averaging ?>
		prim_t avg = (prim_t){
			.rho = .5 * (primL.rho + primR.rho),
			.v = real3_scale(real3_add(primL.v, primR.v), .5),
			.eInt = .5 * (primL.eInt + primR.eInt),
		};
<? -- else -- Roe-averaging, Font 2008 eqn 38 ?>
<? end ?>

		real rho = avg.rho;
		real3 v = avg.v;
		real eInt = avg.eInt;
	
		<? if side == 1 then ?>
		v = _real3(v.y, v.x, v.z);	// -90' rotation to put the y axis contents into the x axis
		<? elseif side == 2 then ?>
		v = _real3(v.z, v.y, v.x);	//-90' rotation to put the z axis in the x axis
		<? end ?>

//TODO NOTE if you're swapping vector components, you have to swap metric components too 

		real3 vL = lower(v);
		real vSq = real3_dot(v, vL);
		real oneOverW2 = 1. - vSq;
		real oneOverW = sqrt(oneOverW2);
		real W = 1. / oneOverW;
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
		real vxSq = v.x * v.x;
		real csSq = heatCapacityRatio * P / (rho * h);
		real cs = sqrt(csSq);

		const real betaUi = betaU.s<?=side?>;
		real discr = sqrt((1. - vSq) * ((1. - vSq * csSq) - vxSq * (1. - csSq)));
		real lambdaMin = (v.x * (1. - csSq) - cs * discr) / (1. - vSq * csSq) * alpha * alpha - betaUi;
		real lambdaMax = (v.x * (1. - csSq) + cs * discr) / (1. - vSq * csSq) * alpha * alpha - betaUi;

		int intindex = side + dim * index;	
		global real* wave = waveBuf + numWaves * intindex;
		wave[0] = lambdaMin;
		wave[1] = v.x * alpha - betaUi;
		wave[2] = v.x * alpha - betaUi;
		wave[3] = v.x * alpha - betaUi;
		wave[4] = lambdaMax;

		real LambdaMin = (lambdaMin + betaUi) / alpha;	//2008 Font eqn 114
		real LambdaMax = (lambdaMax + betaUi) / alpha;	//2008 Font eqn 114
		
		//used by evL and evR
		real ATildeMinus = (gammaU.xx - vxSq) / (gammaU.xx - v.x * LambdaMin);	//2008 Font eqn 113
		real ATildePlus  = (gammaU.xx - vxSq) / (gammaU.xx - v.x * LambdaMax);	//2008 Font eqn 113
		
		//used by evL
		real VMinus = (v.x - LambdaMin) / (gammaU.xx - v.x * LambdaMin);	//2008 Font eqn 113
		real VPlus = (v.x - LambdaMax) / (gammaU.xx - v.x * LambdaMax);	//2008 Font eqn 113
	
		//used by evL and evR
		real CMinus = vL.x - VMinus;	//2008 Font eqn 112
		real CPlus = vL.x - VPlus;	//2008 Font eqn 112

		real kappa = calc_dP_deInt(rho, eInt);	//2008 Font note just after eqn 107
		real kappaTilde = kappa / rho;	//2008 Font eqn 112.  
		//used by evL and evR
		real Kappa = kappaTilde / (kappaTilde - csSq);	//2008 Font eqn 112.  
		//Kappa = h;	//approx for ideal gas
		
		eigenBuf[intindex] = (eigen_t){
<?
for _,field in ipairs(eqn.eigenStructFields) do
	local name,ctype = next(field)
?>
			.<?=name?> = <?=name?>,
<? end ?>
		};
	}<? end ?>
}

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,solver.dim-1 do 
				local prefix = require 'ext.table'.map(eqn.eigenStructFields, function(field)
					local name,ctype = next(field)
					return '\t'..ctype..' '..name..' = eig->'..name..';\n'
				end):concat()
?>
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x_
) { 
	//rotate incoming v's in x
	//TODO do the same for gamma_ij
	<? if side==0 then ?>
	<?=addr2?> const real* x = x_;
	<? elseif side == 1 then ?>
	real x[5] = {x_[0], x_[2], x_[1], x_[3], x_[4]};
	<? elseif side == 2 then ?>
	real x[5] = {x_[0], x_[3], x_[2], x_[1], x_[4]};
	<? end ?>

	<?=prefix?>
	
	real3 vL = lower(v);
	real vxSq = v.x * v.x;
	real hSq = h * h;
	real hW = h * W;
	real W2 = W * W;

	real gamma_gammaUxx = gammaDet * gammaU.xx;
	real gamma_gammaUxy = gammaDet * gammaU.xy;
	real gamma_gammaUxz = gammaDet * gammaU.xz;
	real xi = gammaDet * (gammaU.xx - vxSq);//2008 Font eqn 121
	real Delta = hSq * hW * (Kappa - 1.) * (CPlus - CMinus) * xi;	//2008 Font eqn 121
	
	//min row	2008 Font eqn 118
	real scale;
	scale = hSq / Delta;
	real l5minus = (1 - Kappa) * (-gammaDet * v.x + VPlus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VPlus * xi;
	y[0] = (
		x[0] * (hW * VPlus * xi + l5minus)
		+ x[1] * (gamma_gammaUxx * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * v.x * xi - gamma_gammaUxx * v.x))
		+ x[2] * (gamma_gammaUxy * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * v.y * xi - gamma_gammaUxy * v.x))
		+ x[3] * (gamma_gammaUxz * (1 - Kappa * ATildePlus) + (2. * Kappa - 1.) * VPlus * (W2 * v.z * xi - gamma_gammaUxz * v.x))
		+ x[4] * l5minus
	) * scale;
	//mid normal row	2008 Font eqn 115
	scale = W / (Kappa - 1.);
	y[1] = (
		x[0] * (h - W) 
		+ x[1] * (W * v.x) 
		+ x[2] * (W * v.y) 
		+ x[3] * (W * v.z) 
		+ x[4] * (-W)
	) * scale;
	//mid tangent A row	2008 Font eqn 116
	scale = 1. / (h * xi);
	y[2] = (
		x[0] * (-gammaL.zz * vL.y + gammaL.yz * vL.z) 
		+ x[1] * v.x * (gammaL.zz * vL.y - gammaL.yz * vL.z)
		+ x[2] * (gammaL.zz * (1. - v.x * vL.x) + gammaL.xz * vL.z * v.x)
		+ x[3] * (-gammaL.yz * (1. - vL.x * v.x) - gammaL.xz * vL.y * v.x)
		+ x[4] * (-gammaL.zz * vL.y + gammaL.yz * vL.z)
	) * scale;
	//mid tangent B row	2008 Font eqn 117
	y[3] = (
		x[0] * (-gammaL.yy * vL.z + gammaL.yz * vL.y)
		+ x[1] * v.x * (gammaL.yy * vL.z - gammaL.yz * vL.y)
		+ x[2] * (-gammaL.yz * (1. - vL.x * v.x) - gammaL.xy * vL.z * v.x)
		+ x[3] * (gammaL.yy * (1. - vL.x * v.x) + gammaL.xy * vL.y * v.x)
		+ x[4] * (-gammaL.yy * vL.z + gammaL.yz * vL.y)
	) * scale;
	//max row	2008 Font eqn 118
	scale = -hSq / Delta;
	real l5plus = (1 - Kappa) * (-gammaDet * v.x + VMinus * (W2 * xi - gamma_gammaUxx)) - Kappa * W2 * VMinus * xi;
	y[4] = (
		x[0] * (h * W * VMinus * xi + l5plus)
		+ x[1] * (gamma_gammaUxx * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * v.x * xi - gamma_gammaUxx * v.x))
		+ x[2] * (gamma_gammaUxy * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * v.y * xi - gamma_gammaUxy * v.x))
		+ x[3] * (gamma_gammaUxz * (1 - Kappa * ATildeMinus) + (2. * Kappa - 1.) * VMinus * (W2 * v.z * xi - gamma_gammaUxz * v.x))
		+ x[4] * l5plus
	) * scale;
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	<?=prefix?>
	
	real3 vL = lower(v);
	real hW = h * W;
	real W2 = W * W;

	//2008 Font eqns 108-111
	y[0] = x[0]
		+ x[1] * (Kappa / hW)
		+ x[2] * (W * vL.y)
		+ x[3] * (W * vL.z)
		+ x[4];
	y[1] = x[0] * (hW * CMinus)
		+ x[1] * (vL.x)
		+ x[2] * (h * (gammaL.xy + 2. * W2 * vL.y * vL.x))
		+ x[3] * (h * (gammaL.xz + 2. * W2 * vL.x * vL.z))
		+ x[4] * (hW * CPlus);
	y[2] = x[0] * (hW * vL.y)
		+ x[1] * (vL.y)
		+ x[2] * (h * (gammaL.yy + 2. * W2 * vL.y * vL.y))
		+ x[3] * (h * (gammaL.yz + 2. * W2 * vL.y * vL.z))
		+ x[4] * (hW * vL.y);
	y[3] = x[0] * (hW * vL.z)
		+ x[1] * (vL.z)
		+ x[2] * (h * (gammaL.yz + 2. * W2 * vL.y * vL.z))
		+ x[3] * (h * (gammaL.zz + 2. * W2 * vL.z * vL.z))
		+ x[4] * (hW * vL.z);
	y[4] =x[0] * (hW * ATildeMinus - 1.)
		+ x[1] * (1. - Kappa / hW)
		+ x[2] * (W * vL.y * (2. * hW - 1.))
		+ x[3] * (W * vL.z * (2. * hW - 1.))
		+ x[4] * (hW * ATildePlus - 1.);
	
	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real tmp = y[1];
	y[1] = y[1+<?=side?>];
	y[1+<?=side?>] = tmp;
	<? end ?>
}

<?	if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x_
) {
	//rotate incoming v's in x
	<? if side==0 then ?>
	<?=addr2?> const real* x = x_;
	<? elseif side == 1 then ?>
	real x[5] = {x_[0], x_[2], x_[1], x_[3], x_[4]};
	<? elseif side == 2 then ?>
	real x[5] = {x_[0], x_[3], x_[2], x_[1], x_[4]};
	<? end ?>

	//rotate outgoing y's x's into side
	<? if side ~= 0 then ?>
	real tmp = y[1];
	y[1] = y[1+<?=side?>];
	y[1+<?=side?>] = tmp;
	<? end ?>
}
<?				end
			end
		end
	end
end ?>

kernel void constrainU(
	global cons_t* UBuf
) {
	SETBOUNDS(0,0);

	global cons_t* U = UBuf + index;
	
	U->D = max(U->D, (real)DMin);
	U->tau = max(U->tau, (real)tauMin);

	U->D = min(U->D, (real)DMax);
	U->tau = min(U->tau, (real)tauMax);
}

//TODO update to include alphas, betas, and gammas
kernel void updatePrims(
	global prim_t* primBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(2,1);
	
	const global cons_t* U = UBuf + index;
	real D = U->D;
	real3 S = U->S;
	real tau = U->tau;

	global prim_t* prim = primBuf + index;
	real3 v = prim->v;

	real SLen = metricLen(S);
	real PMin = max(SLen - tau - D + SLen * solvePrimVelEpsilon, solvePrimPMinEpsilon);
	real PMax = (heatCapacityRatio - 1.) * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real vLen = SLen / (tau + D + P);
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
			vSq = metricLenSq(v);
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
