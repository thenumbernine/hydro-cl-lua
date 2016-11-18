// the default eig transforms, using eig struct as a dense matrix:

<? for side=0,2 do ?>
void eigen_leftTransform_<?=side?>(
	real* y,
	eigen_t eig,
	const real* x
) {
	const real* A = eig.evL;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

void eigen_rightTransform_<?=side?>(
	real* y,
	eigen_t eig,
	const real* x
) {
	const real* A = eig.evR;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

<? if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>(
	real* y,
	eigen_t eig,
	const real* x
) {
	const real* A = eig.A;
	for (int i = 0; i < numStates; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numStates*j] * x[j];
		}
		y[i] = sum;
	}
}
<?	end ?>
<? end ?>
