range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	const __global prim_t* prim,
	int side
) {
	real rho = prim->rho;
	real eInt = prim->eInt;
	real vSq = coordLenSq(prim->v);
	real P = calc_P(rho, eInt);
	real h = calc_h(rho, P, eInt);
	real csSq = gamma * P / (rho * h);
	real cs = sqrt(csSq);
	//for the particular direction
	real vi = prim->v.s[side];
	real viSq = vi * vi;

	// Marti 1998 eqn 19
	// also Marti & Muller 2008 eqn 68
	// also Font 2008 eqn 106
	real discr = sqrt((1. - vSq) * ((1. - vSq * csSq) - viSq * (1. - csSq)));
	return (range_t){
		.min = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq),
		.max = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq)
	};
}

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
__kernel void calcDT(
	__global real* dtBuf,
	const __global cons_t* UBuf,
	const __global prim_t* primBuf
) {
	SETBOUNDS(0,0);
	if (OOB(2,2)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const __global cons_t* U = UBuf + index;
	const __global prim_t* prim = primBuf + index;

	real dt = INFINITY;
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		range_t lambda = calcCellMinMaxEigenvalues(U, prim, <?=side?>); 
		lambda.min = min((real)0., lambda.min);
		lambda.max = max((real)0., lambda.max);
		dt = min(dt, dx_at<?=side?>(i) / (fabs(lambda.max - lambda.min) + (real)1e-9));
	}<?end?>
	dtBuf[index] = dt; 
}

prim_t calcEigenBasisSide(prim_t primL, prim_t primR) {
	return (prim_t){
		.rho = .5 * (primL.rho + primR.rho),
		.v = real3_scale(real3_add(primL.v, primR.v), .5),
		.eInt = .5 * (primL.eInt + primR.eInt),
	};
}

__kernel void calcEigenBasis(
	__global real* waveBuf,
	__global eigen_t* eigenBuf,
	const __global prim_t* primBuf
#if defined(checkFluxError)
	, __global fluxXform_t* fluxXformBuf
#endif
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
	
		prim_t primL = primBuf[indexL];
		prim_t primR = primBuf[indexR];

		prim_t avg = calcEigenBasisSide(primL, primR);

		real rho = avg.rho;
		real3 v = avg.v;
		real eInt = avg.eInt;
		
#if dim > 1
		if (side == 1) {
			v = _real3(v.y, -v.x, v.z);	// -90' rotation to put the y axis contents into the x axis
		} 
#endif
#if dim > 2
		else if (side == 2) {
			v = _real3(v.z, v.y, -v.x);	//-90' rotation to put the z axis in the x axis
		}
#endif
		
		real vSq = coordLenSq(v);
		real oneOverW2 = 1. - vSq;
		real oneOverW = sqrt(oneOverW2);
		real W = 1. / oneOverW;
		real W2 = 1. / oneOverW2;
		real P = gamma_1 * rho * eInt;
		real h = 1. + eInt + P / rho;
		real P_over_rho_h = P / (rho * h);

		real hW = h * W;
		real hSq = h * h;

		//TODO check out Font's paper on where the metric coefficients go ...
		real vxSq = v.x * v.x;
		real csSq = gamma * P_over_rho_h;
		real cs = sqrt(csSq);

		real discr = sqrt((1. - vSq) * ((1. - vSq * csSq) - vxSq * (1. - csSq)));
		real lambdaMin = (v.x * (1. - csSq) - cs * discr) / (1. - vSq * csSq);
		real lambdaMax = (v.x * (1. - csSq) + cs * discr) / (1. - vSq * csSq);

		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		wave[0] = lambdaMin;
		wave[1] = v.x;
		wave[2] = v.x;
		wave[3] = v.x;
		wave[4] = lambdaMax;

		__global eigen_t* eigen = eigenBuf + intindex;
		
		real Kappa = h;	//true for ideal gas. otherwise the general equation for Kappa gets instable at high Lorentz factors 
		real AMinus = (1. - vxSq) / (1. - v.x * lambdaMin);
		real APlus  = (1. - vxSq) / (1. - v.x * lambdaMax);
		__global real* evR = eigen->evR;
		//min col
		evR[0 + numStates * 0] = 1.;
		evR[1 + numStates * 0] = hW * AMinus * lambdaMin;	//inf
		evR[2 + numStates * 0] = hW * v.y;
		evR[3 + numStates * 0] = hW * v.z;
		evR[4 + numStates * 0] = hW * AMinus - 1.;
		//mid col (normal)
		evR[0 + numStates * 1] = oneOverW;	// = Kappa / hW
		evR[1 + numStates * 1] = v.x;
		evR[2 + numStates * 1] = v.y;
		evR[3 + numStates * 1] = v.z;
		evR[4 + numStates * 1] = 1. - oneOverW;	// = 1. - Kappa / hW;
		//mid col (tangent A)
		evR[0 + numStates * 2] = W * v.y;
		evR[1 + numStates * 2] = 2. * h * W2 * v.x * v.y;
		evR[2 + numStates * 2] = h * (1. + 2. * W2 * v.y * v.y);
		evR[3 + numStates * 2] = 2. * h * W2 * v.y * v.z;
		evR[4 + numStates * 2] = (2. * hW - 1.) * W * v.y;
		//mid col (tangent B)
		evR[0 + numStates * 3] = W * v.z;
		evR[1 + numStates * 3] = 2. * h * W2 * v.x * v.z;
		evR[2 + numStates * 3] = 2. * h * W2 * v.y * v.z;
		evR[3 + numStates * 3] = h * (1. + 2. * W2 * v.z * v.z);
		evR[4 + numStates * 3] = (2. * hW - 1.) * W * v.z;
		//max col 
		evR[0 + numStates * 4] = 1.;
		evR[1 + numStates * 4] = hW * APlus * lambdaMax;	//inf
		evR[2 + numStates * 4] = hW * v.y;
		evR[3 + numStates * 4] = hW * v.z;
		evR[4 + numStates * 4] = hW * APlus - 1.;

		__global real* evL = eigen->evL; 
		real Delta = hSq * hW * (Kappa - 1.) * (1. - vxSq) * (APlus * lambdaMax - AMinus * lambdaMin);
		//min row
		real scale;
		scale = hSq / Delta;
		evL[0 + numStates * 0] = scale * (hW * APlus * (v.x - lambdaMax) - v.x - W2 * (vSq - vxSq) * (2. * Kappa - 1.) * (v.x - APlus * lambdaMax) + Kappa * APlus * lambdaMax);
		evL[0 + numStates * 1] = scale * (1. + W2 * (vSq - vxSq) * (2. * Kappa - 1.) * (1. - APlus) - Kappa * APlus);
		evL[0 + numStates * 2] = scale * (W2 * v.y * (2. * Kappa - 1.) * APlus * (v.x - lambdaMax));
		evL[0 + numStates * 3] = scale * (W2 * v.z * (2. * Kappa - 1.) * APlus * (v.x - lambdaMax));
		evL[0 + numStates * 4] = scale * (-v.x - W2 * (vSq - vxSq) * (2. * Kappa - 1.) * (v.x - APlus * lambdaMax) + Kappa * APlus * lambdaMax);
		//mid normal row
		scale = W / (Kappa - 1.);
		evL[1 + numStates * 0] = scale * (h - W);
		evL[1 + numStates * 1] = scale * (W * v.x);
		evL[1 + numStates * 2] = scale * (W * v.y);
		evL[1 + numStates * 3] = scale * (W * v.z);
		evL[1 + numStates * 4] = scale * (-W);
		//mid tangent A row
		scale = 1. / (h * (1. - vxSq));
		evL[2 + numStates * 0] = scale * (-v.y);
		evL[2 + numStates * 1] = scale * (v.x * v.y);
		evL[2 + numStates * 2] = scale * (1. - vxSq);
		evL[2 + numStates * 3] = 0.;
		evL[2 + numStates * 4] = scale * (-v.y);
		//mid tangent B row
		evL[3 + numStates * 0] = scale * (-v.z);
		evL[3 + numStates * 1] = scale * (v.x * v.z);
		evL[3 + numStates * 2] = 0.;
		evL[3 + numStates * 3] = scale * (1. - vxSq);
		evL[3 + numStates * 4] = scale * (-v.z);
		//max row
		scale = -hSq / Delta;
		evL[4 + numStates * 0] = scale * (hW * AMinus * (v.x - lambdaMin) - v.x - W2 * (vSq - vxSq) * (2. * Kappa - 1.) * (v.x - AMinus * lambdaMin) + Kappa * AMinus * lambdaMin);
		evL[4 + numStates * 1] = scale * (1. + W2 * (vSq - vxSq) * (2. * Kappa - 1.) * (1. - AMinus) - Kappa * AMinus);
		evL[4 + numStates * 2] = scale * (W2 * v.y * (2. * Kappa - 1.) * AMinus * (v.x - lambdaMin));
		evL[4 + numStates * 3] = scale * (W2 * v.z * (2. * Kappa - 1.) * AMinus * (v.x - lambdaMin));
		evL[4 + numStates * 4] = scale * (-v.x - W2 * (vSq - vxSq) * (2. * Kappa - 1.) * (v.x - AMinus * lambdaMin) + Kappa * AMinus * lambdaMin);
	
#if dim > 1
		if (side == 1) {
			for (int i = 0; i < numStates; ++i) {
				real tmp;
				//-90' rotation applied to the LHS of incoming velocity vectors, to move their y axis into the x axis
				// is equivalent of a -90' rotation applied to the RHS of the flux jacobian A
				// and A = Q V Q-1 for Q = the right eigenvectors and Q-1 the left eigenvectors
				// so a -90' rotation applied to the RHS of A is a +90' rotation applied to the RHS of Q-1 the left eigenvectors
				//and while a rotation applied to the LHS of a vector rotates the elements of its column vectors, a rotation applied to the RHS rotates the elements of its row vectors 
				//each row's y <- x, x <- -y
				tmp = evL[i + numStates * cons_Sx];
				evL[i + numStates * cons_Sx] = -evL[i + numStates * cons_Sy];
				evL[i + numStates * cons_Sy] = tmp;
				//a -90' rotation applied to the RHS of A must be corrected with a 90' rotation on the LHS of A
				//this rotates the elements of the column vectors by 90'
				//each column's x <- y, y <- -x
				tmp = evR[cons_Sx + numStates * i];
				evR[cons_Sx + numStates * i] = -evR[cons_Sy + numStates * i];
				evR[cons_Sy + numStates * i] = tmp;
			}
		}
#endif
#if dim > 2
		else if (side == 2) {
			for (int i = 0; i < numStates; ++i) {
				real tmp;
				tmp = evL[i + numStates * cons_Sx];
				evL[i + numStates * cons_Sx] = -evL[i + numStates * cons_Sz];
				evL[i + numStates * cons_Sz] = tmp;
				tmp = evR[cons_Sx + numStates * i];
				evR[cons_Sx + numStates * i] = -evR[cons_Sz + numStates * i];
				evR[cons_Sz + numStates * i] = tmp;
			}
		}
#endif
	}
}

__kernel void constrainU(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);

	__global cons_t* U = UBuf + index;
	
	U->D = max(U->D, (real)DMin);
	U->tau = max(U->tau, (real)tauMin);

	U->D = min(U->D, (real)DMax);
	U->tau = min(U->tau, (real)tauMax);
}

__kernel void updatePrims(
	__global prim_t* primBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(2,1);
	
	const __global cons_t* U = UBuf + index;
	real D = U->D;
	real3 S = U->S;
	real tau = U->tau;

	__global prim_t* prim = primBuf + index;
	real3 v = prim->v;

	real SLen = coordLen(S);
	real PMin = max(SLen - tau - D + SLen * solvePrimVelEpsilon, solvePrimPMinEpsilon);
	real PMax = gamma_1 * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real vLen = SLen / (tau + D + P);
		real vSq = vLen * vLen;
		real W = 1. / sqrt(1. - vSq);
		real eInt = (tau + D * (1. - W) + P * (1. - W*W)) / (D * W);
		real rho = D / W;
		real f = gamma_1 * rho * eInt - P;
		real csSq = gamma_1 * (tau + D * (1. - W) + P) / (tau + D + P);
		real df_dP = vSq * csSq - 1.;
		real newP = P - f / df_dP;
		newP = max(newP, PMin);
		real PError = fabs(1. - newP / P);
		P = newP;
		if (PError < solvePrimStopEpsilon) {
			v = real3_scale(S, 1. / (tau + D + P));
			vSq = coordLenSq(v);
			W = 1. / sqrt(1. - vSq);
			rho = D / W;
			rho = max(rho, (real)rhoMin);
			rho = min(rho, (real)rhoMax);
			eInt = P / (rho * gamma_1);
			eInt = min(eInt, (real)eIntMax);
			prim->rho = rho;
			prim->v = v;
			prim->eInt = eInt;
//printf("cell %d finished with prims = %f %f %f\n", index, rho, v.x, eInt);
			return;
		}
	}
}
