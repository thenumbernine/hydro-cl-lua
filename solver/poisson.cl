/*
used with Euler
but works with anything that has a cons_t::rho, m[], ETotal

rho is a delta function

-4 pi rho delta3(r) = del^2 rho/|r| = del^2 phi
del rho/|r| = del phi = -rho r/|r|  points inwards towards the rho


let E = del rho/|r| = del phi
then -4 pi rho delta3(r) = del^2 phi = del^2 rho/|r| = del . E
then phi = -4 pi del^-2 (rho delta3(r))

*/

kernel void initPotential(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	real rho = 0;
	<?=calcRho?>
	UBuf[index].<?=potentialField?> = -rho;
}

kernel void solvePoisson(
	<?=args?>
) {
	SETBOUNDS(2,2);

<? for j=0,solver.dim-1 do ?>
	real dx<?=j?> = dx<?=j?>_at(i);
<? end ?>

	real skewSum = 0;
<? for j=0,solver.dim-1 do ?>
	if (i.s<?=j?> < gridSize.s<?=j?>-3) {
		skewSum += UBuf[index + stepsize.s<?=j?>].<?=potentialField?> / (dx<?=j?> * dx<?=j?>);
	}
	if (i.s<?=j?> > 2) {
		skewSum += UBuf[index - stepsize.s<?=j?>].<?=potentialField?> / (dx<?=j?> * dx<?=j?>);
	}
<? end ?>

	const real diag = -2. * (0
<? for j=0,solver.dim-1 do ?>
		+ 1. / (dx<?=j?> * dx<?=j?>)
<? end ?>
	);

	real rho = 0;
	<?=calcRho?>

	UBuf[index].<?=potentialField?> = (rho - skewSum) / diag;
}
