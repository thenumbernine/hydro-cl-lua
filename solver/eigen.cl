// the default eig transforms, using eig struct as a dense matrix:

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,2 do ?>
/*
Y stores the results in real[numWaves]
eig provides the variables used to compute the left eigenvector linear transformation.
	Notice I'm not holding the entire matrix, and I'm not multiplying the entire matrix -- only the variables used to reconstruct it.  Why waste memory and instructions?
X points to the cons_t structure, stored in real[numStates] 
*/
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y/*[numWaves]*/,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x/*[numStates]*/
) {
	<?=addr1?> const real* A = eig->evL;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numIntStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

/*
Y stores the result in real[numIntStates], meaning if your output is a cons_t then you will have to write Y[numIntStates..numStates-1] yourself
eig provides the variables used to compute the right eigenvector linear transformation.
X points to the waves, stored in real[numWaves] 
*/
void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y/*[numIntStates]*/,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x/*[numWaves]*/
) {
	<?=addr1?> const real* A = eig->evR;
	for (int i = 0; i < numIntStates; ++i) {
		real sum = 0;
		for (int j = 0; j < numWaves; ++j) {
			sum += A[i+numIntStates*j] * x[j];
		}
		y[i] = sum;
	}
}

<? 				if solver.checkFluxError then ?>
/*
This is the dA/dU matrix linear transformation.

Y stores the result in real[numIntStates] following the eigen_rightTransform convention
eig provides the variables used to compute the flux Jacobian linear transformation.
X accepts real[numWaves]

How about also implementing dA/dW and then implementing this as dA/dW * dW/dU * x = dA_dW(dW_dU(x))?  would it be less instructions?.
*/
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y/*[numIntStates]*/,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x/*[numStates]*/
) {
	<?=addr1?> const real* A = eig->A;
	for (int i = 0; i < numIntStates; ++i) {
		real sum = 0;
		for (int j = 0; j < numIntStates; ++j) {
			sum += A[i+numIntStates*j] * x[j];
		}
		y[i] = sum;
	}
}
<?				end
			end
		end
	end
end ?>
