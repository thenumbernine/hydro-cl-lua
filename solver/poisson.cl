/*
rho is a Dirac delta function

-4 pi rho delta3(r) = del^2 rho/|r| = del^2 phi
del rho/|r| = del phi = -rho r/|r|  points inwards towards the rho


let E = del rho/|r| = del phi
then -4 pi rho delta3(r) = del^2 phi = del^2 rho/|r| = del . E
then phi = -4 pi del^-2 (rho delta3(r))

TODO curved space divergence? g^ab phi,ab - phi,c Gamma^c_ab g^ab
or is this already solved in the discrete case? 

1/sqrt(g) (sqrt(g) f^,i)_,i

discrete evaluation:
1/sqrt(g) (sqrt(g(x+dxi/2)) g^ij(x+dxj/2) (f(x+dxj) - f(x)) / dx(x+dxj/2))_,i
*/

kernel void initPotential(
	global <?=poisson:getPotBufType()?>* UBuf
) {
	SETBOUNDS(2,2);
	real rho = 0;
	<?=calcRho?>
	UBuf[index].<?=poisson.potentialField?> = -rho;
}

kernel void solvePoisson(
	global <?=poisson:getPotBufType()?>* UBuf
) {
	SETBOUNDS(2,2);

#if 1
<? for j=0,solver.dim-1 do ?>
	real dx<?=j?> = dx<?=j?>_at(i);
<? end ?>

	real skewSum = 0.
<? for j=0,solver.dim-1 do ?>
		+ UBuf[index + stepsize.s<?=j?>].<?=poisson.potentialField?> / (dx<?=j?> * dx<?=j?>)
		+ UBuf[index - stepsize.s<?=j?>].<?=poisson.potentialField?> / (dx<?=j?> * dx<?=j?>)
<? end 
?>	;

	const real diag = -2. * (0.
<? for j=0,solver.dim-1 do ?>
		+ 1. / (dx<?=j?> * dx<?=j?>)
<? end ?>
	);
#else
	//f_;a^a = g^ab (f_,ab - Gamma^c_ab f,c)
	// = 1/sqrt|g| (sqrt|g| g^ab f_,a)_,b
	//I think I'm gonna use finite-differencing with the second one

#endif

	real rho = 0;
	<?=calcRho?>

	UBuf[index].<?=poisson.potentialField?> = (rho - skewSum) / diag;
}
