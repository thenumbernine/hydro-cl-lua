<?
local makePartials = require 'eqn.makepartial'
local derivOrder = 2 * solver.numGhost
local makePartial = function(...) return makePartials.makePartial(derivOrder, solver, ...) end
local makePartial2 = function(...) return makePartials.makePartial2(derivOrder, solver, ...) end
?>

kernel void solveMinimalDistortionEllipticShift<?=op.suffix?>(
	global <?=op:getPotBufType()?>* UBuf<?
if op.stopOnEpsilon then ?>,
	global real* reduceBuf<?
end ?>
) {
<? if not op.stopOnEpsilon then ?>
	SETBOUNDS(numGhost,numGhost);
<? else ?>
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		reduceBuf[index] = 0.;
		return;
	}
<? end ?>

	global <?=op:getPotBufType()?>* U = UBuf + index;

<? for j=0,solver.dim-1 do ?>
	real dx<?=j?> = dx<?=j?>_at(i);
<? end ?>

	real3 intIndex = _real3(i.x, i.y, i.z);
	real3 volL, volR;
<? for j=0,solver.dim-1 do ?>
	intIndex.s<?=j?> = i.s<?=j?> - .5;
	volL.s<?=j?> = volume_at(cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?> + .5;
	volR.s<?=j?> = volume_at(cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?>;
<? end ?>
	real volAtX = volume_at(cell_x(i));

<?
local scalar = op.scalar
local zero = scalar..'_zero'
local add3 = scalar..'_add3'
local sub = scalar..'_sub'
local mul = scalar..'_mul'
local lenSq = scalar..'_lenSq'
local real_mul = scalar..'_real_mul'
?>


	<?=scalar?> skewSum = <?=scalar?>_zero;

<? for j=0,solver.dim-1 do ?>
	skewSum = <?=add3?>(skewSum,
		<?=real_mul?>(U[stepsize.s<?=j?>].<?=op.potentialField?>, volR.s<?=j?> / (dx<?=j?> * dx<?=j?>)),
		<?=real_mul?>(U[-stepsize.s<?=j?>].<?=op.potentialField?>, volL.s<?=j?> / (dx<?=j?> * dx<?=j?>)));
<? end ?>
	skewSum = <?=real_mul?>(skewSum, 1. / volAtX);

	const real diag = (0.
<? for j=0,solver.dim-1 do ?>
		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end ?>
	) / volAtX;


	<?=scalar?> rho = <?=zero?>;
	<?=op:getCalcRhoCode() or ''?>

	<?=scalar?> oldU = U-><?=op.potentialField?>;
	
	//Gauss-Seidel iteration: x_i = (b_i - A_ij x_j) / A_ii
	<?=scalar?> newU = <?=real_mul?>(<?=sub?>(rho, skewSum), 1. / diag);

<?
if op.stopOnEpsilon then
?>	<?=scalar?> dU = <?=sub?>(newU, oldU);
	reduceBuf[index] = <?=lenSq?>(dU);
<?
end
?>
	U-><?=op.potentialField?> = newU;
}
