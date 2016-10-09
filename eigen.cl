// the default eigen transforms, using eigen struct as a dense matrix:

void fluxTransform(
	real* y,
	const __global fluxXform_t* f,
	const real* x,
	int side
) {
	const __global real* A = f->A;
	for (int i = 0; i < numStates; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numStates*j] * x[j];
		}
		y[i] = sum;
	}
}


void eigen_leftTransform(
	real* y,
	const __global eigen_t* eigen,
	const real* x,
	int side
) {
	const __global real* A = eigen->evL;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

void eigen_rightTransform(
	real* y,
	const __global eigen_t* eigen,
	const real* x,
	int side
) {
	const __global real* A = eigen->evR;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

real eigen_calcDisplayVar(
	int displayVar,
	const __global eigen_t* eigen
) {
	displayVar -= displayFirst_eigen;
	if (displayVar < numStates * numWaves) return eigen->evL[displayVar];
	displayVar -= numStates * numWaves;
	if (displayVar < numStates * numWaves) return eigen->evR[displayVar];
	return 0.;
}
