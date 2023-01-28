//// MODULE_NAME: <?=op.symbolPrefix?>comments

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
<?
local scalar = op.scalar 
?>


//// MODULE_NAME: <?=initPotential?>
//// MODULE_DEPENDS: <?=table.concat(op.codeDepends or {}, ' ')?>
/*
used by hydro/op/poisson_krylov.lua and hydro/op/relaxation.lua
initialize the relaxation solver field 
this is only called upon solver reset
each iteration uses the previous iteration's results as the starting point
*/
kernel void <?=initPotential?>(
	constant <?=solver_t?> const * const psolver,
	global <?=op:getPotBufType()?> * const UBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(solver.numGhost, solver.numGhost);
	global <?=op:getPotBufType()?> * const U = UBuf + index;
	<?=scalar?> source = {};
<?=op:getPoissonDivCode() or ""?>
	
<? if cmdline.selfGravInitPotential == "+" then
?>	UBuf[index].<?=op.potentialField?> = source;
<? else
?>	UBuf[index].<?=op.potentialField?> = -source;
<? end
?>
}

//// MODULE_NAME: <?=copyWriteToPotentialNoGhost?>
//// MODULE_DEPENDS: <?=SETBOUNDS_NOGHOST?>
//used by hydro/op/relaxation.lua
kernel void <?=copyWriteToPotentialNoGhost?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const UBuf,
	global real const * const writeBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS_NOGHOST?>();
	UBuf[index].<?=op.potentialField?> = writeBuf[index];
}

//// MODULE_NAME: <?=setReduceToPotentialSquared?>
//// MODULE_DEPENDS: <?=SETBOUNDS_NOGHOST?>
//used by hydro/op/relaxation.lua
kernel void <?=setReduceToPotentialSquared?>(
	constant <?=solver_t?> const * const psolver,
	global real * const reduceBuf,
	global <?=cons_t?> const * const UBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS_NOGHOST?>();
	reduceBuf[index] = lenSq(UBuf[index].<?=op.potentialField?>);
}
