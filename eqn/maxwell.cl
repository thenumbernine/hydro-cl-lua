#define sqrt_1_2 0.70710678118654757273731092936941422522068023681641

real calcEigenvalue() { 
	return 1./(sqrt_eps0 * sqrt_mu0);
}

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global cons_t* U
) {
	real lambda = calcEigenvalue();
	return (range_t){-lambda, lambda};
}
<? end ?>

kernel void calcEigenBasis(
	global real* waveBuf,
	global eigen_t* eigenBuf,
	const global consLR_t* ULRBuf
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
		real lambda = calcEigenvalue();
		
		int intindex = side + dim * index;	
		global real* wave = waveBuf + numWaves * intindex;
		wave[0] = -lambda;
		wave[1] = -lambda;
		wave[2] = 0;
		wave[3] = 0;
		wave[4] = lambda;
		wave[5] = lambda;
	
		//no eigenbuf info since waves are unrelated to state
	}
}

<? for side=0,2 do ?>

void eigen_leftTransform_<?=side?>(
	real* y,
	const global eigen_t* eigen,
	const real* x
) {
	const real ise = sqrt_1_2 / sqrt_eps0;
	const real isu = sqrt_1_2 / sqrt_mu0;

	<? if side==0 then ?>
	
	y[0] = x[2] *  ise + x[4] * isu;
	y[1] = x[1] * -ise + x[5] * isu;
	y[2] = x[0] * -ise + x[3] * isu;
	y[3] = x[0] *  ise + x[3] * isu;
	y[4] = x[1] *  ise + x[5] * isu;
	y[5] = x[2] * -ise + x[4] * isu;
	
	<? elseif side==1 then ?>
	
	y[0] = x[0] *  ise + x[5] * isu;
	y[1] = x[2] * -ise + x[3] * isu;
	y[2] = x[1] * -ise + x[4] * isu;
	y[3] = x[1] *  ise + x[4] * isu;
	y[4] = x[2] *  ise + x[3] * isu;
	y[5] = x[0] * -ise + x[5] * isu;
	
	<? elseif side==2 then ?>
	
	y[0] = x[1] *  ise + x[3] * isu;
	y[1] = x[0] * -ise + x[4] * isu;
	y[2] = x[2] * -ise + x[5] * isu;
	y[3] = x[2] *  ise + x[5] * isu;
	y[4] = x[0] *  ise + x[4] * isu;
	y[5] = x[1] * -ise + x[3] * isu;
	
	<? end ?>
}

void eigen_rightTransform_<?=side?>(
	real* y,
	const global eigen_t* eigen,
	const real* x
) {
	const real se = sqrt_1_2 * sqrt_eps0;
	const real su = sqrt_1_2 * sqrt_mu0;

	<? if side==0 then ?>
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
	y[0] = se * (-x[2] + x[3]);
	y[1] = se * (-x[1] + x[4]);
	y[2] = se * (x[0] + -x[5]);
	y[3] = su * (x[2] + x[3]);
	y[4] = su * (x[0] + x[5]);
	y[5] = su * (x[1] + x[4]);
	
	<? elseif side==1 then ?>

/*
x, -z, -y, y, z, -x
z,  x,  y, y, x,  z

1  0  0 0 0 -1
0  0 -1 1 0  0
0 -1  0 0 1  0
0  1  0 0 1  0
0  0  1 1 0  0
1  0  0 0 0  1
*/
	y[0] = se * (x[0] - x[5]);
	y[1] = se * (-x[2] + x[3]);
	y[2] = se * (-x[1] + x[4]);
	y[3] = su * (x[1] + x[4]);
	y[4] = su * (x[2] + x[3]);
	y[5] = su * (x[0] + x[5]);
	
	<? elseif side==2 then ?>

/*
y, -x, -z, z, x, -y
x,  y,  z, z,  y,  x

0 -1  0 0 1  0
1  0  0 0 0 -1
0  0 -1 1 0  0
1  0  0 0 0  1
0  1  0 0 1  0
0  0  1 1 0  0
*/
	y[0] = se * (-x[1] + x[4]);
	y[1] = se * (x[0] - x[5]);
	y[2] = se * (-x[2] + x[3]);
	y[3] = su * (x[0] + x[5]);
	y[4] = su * (x[1] + x[4]);
	y[5] = su * (x[2] + x[3]);
	
	<? end ?>
}

<? if solver.checkFluxError then ?>
void fluxTransform_<?=side?>(
	real* y,
	const global eigen_t* eigen,
	const real* x_
) {
	//swap input dim x<->side
	const cons_t* x = (const cons_t*)x_;
	real3 epsE = x->epsE;
	real3 B = x->B;

	<? if side==0 then ?>
	
	y[0] = 0;
	y[1] = B.z / mu0;
	y[2] = -B.y / mu0;
	y[3] = 0;
	y[4] = -epsE.z / eps0;
	y[5] = epsE.y / eps0;

	<? elseif side==1 then ?>
		
	y[0] = -B.z / mu0;
	y[1] = 0;
	y[2] = B.x / mu0;
	y[3] = epsE.z / eps0;
	y[4] = 0;
	y[5] = -epsE.x / eps0;
		
	<? elseif side==2 then ?>
		
	y[0] = B.y / mu0;
	y[1] = -B.x / mu0;
	y[2] = 0;
	y[3] = -epsE.y / eps0;
	y[4] = epsE.x / eps0;
	y[5] = 0;
		
	<? end ?>
}
<?	end ?>
<? end ?>

kernel void addSource(
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(2,2);
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(U->epsE, 1. / eps0 * sigma));
}
