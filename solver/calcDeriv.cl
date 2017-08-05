/*
used by all the finite volume solvers
*/
kernel void calcDerivFromFlux(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	
	real3 x = cell_x(i);
	real volume = volume_at(x);

	<? for side=0,solver.dim-1 do ?>{
		int indexIntL = <?=side?> + dim * index;
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * stepsize[<?=side?>]; 
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + indexIntR;

<? -- trying to apply Trangenstein to my holonomic coordinate code:
if false
and require 'geom.cylinder'.is(solver.geometry) 
and solver.dim == 2
then
?>
		real rL = x.x - .5 * grid_dx0;
		real rR = x.x + .5 * grid_dx0;
		real thetaL = x.y - .5 * grid_dx1;
		real thetaR = x.y + .5 * grid_dx1;
		real dr = rR - rL;
		real dtheta = thetaR - thetaL;
		
		//the volume of a discrete cell
		// is the integral of the jacobian determinant across the coord range of the cell
		//so integral from r1 to r2 of r = r^2/2 at r2 - r1 = (r2^2-r1^2)/2
<?	if side==0 then	-- r+ and r-
?>		real areaL = rL * dtheta;
		real areaR = rR * dtheta;
<?	elseif side==1 then	-- θ+ and θ-	
?>		real areaL = dr;
		real areaR = dr;
<?	end
?>
		//so the next question is -- how to incorporate the d/dxi ln(sqrt(g)) ...
		// is it magically already considered?
		//or do I have to also scale by that here (in addition to scaling by the sides) ?
		real volume = .5 * (rR*rR - rL*rL) * (thetaR - thetaL);

<? else
?>		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real volumeIntL = volume_at(xIntL);
		real areaL = volumeIntL / grid_dx<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		real volumeIntR = volume_at(xIntR);
		real areaR = volumeIntR / grid_dx<?=side?>;
<? end
?>
		for (int j = 0; j < numStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
	}<? end ?>
}
