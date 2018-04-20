// the default eig transforms, using eig struct as a dense matrix:

<?
for side=0,2 do ?>
/*
eig provides the variables used to compute the left eigenvector linear transformation.
	Notice I'm not holding the entire matrix, and I'm not multiplying the entire matrix -- only the variables used to reconstruct it.  Why waste memory and instructions?
*/
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> x
) {
	<?=eqn.waves_t?> y;
	<?=addr1?> const real* A = eig.evL;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numIntStates; ++j) {
			sum += A[i+numWaves*j] * x.ptr[j];
		}
		y.ptr[i] = sum;
	}
	return y;
}

/*
eig provides the variables used to compute the right eigenvector linear transformation.
*/
<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> x
) {
	<?=eqn.cons_t?> y;
	<?=addr1?> const real* A = eig.evR;
	for (int i = 0; i < numIntStates; ++i) {
		real sum = 0;
		for (int j = 0; j < numWaves; ++j) {
			sum += A[i+numIntStates*j] * x.ptr[j];
		}
		y.ptr[i] = sum;
	}
	for (int i = numIntStates; i < numStates; ++i) {
		y.ptr[i] = 0;
	}
	return y;
}

<? 	if solver.checkFluxError then ?>
/*
This is the dA/dU matrix linear transformation.

Y stores the result in real[numStates]
eig provides the variables used to compute the flux Jacobian linear transformation.
X accepts real[numWaves]

How about also implementing dA/dW and then implementing this as dA/dW * dW/dU * x = dA_dW(dW_dU(x))?  would it be less instructions?.
*/
<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> x
) {
	<?=eqn.cons_t?> y;
	<?=addr1?> const real* A = eig.A;
	for (int i = 0; i < numIntStates; ++i) {
		real sum = 0;
		for (int j = 0; j < numIntStates; ++j) {
			sum += A[i+numIntStates*j] * x.ptr[j];
		}
		y.ptr[i] = sum;
	}
	for (int i = numIntStates; i < numStates; ++i) {
		y.ptr[i] = 0;
	}
	return y;
}
<?	end
end ?>
