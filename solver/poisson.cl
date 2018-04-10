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

//initialize the poisson solver field 
//this is only called upon solver reset
//each iteration uses the previous iteration's results as the starting point
kernel void initPoissonPotential<?=poisson.suffix?>(
	global <?=poisson:getPotBufType()?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=poisson:getPotBufType()?>* U = UBuf + index;
	real rho = 0;
	<?=poisson:getCalcRhoCode() or ''?>
	UBuf[index].<?=poisson.potentialField?> = -rho;
}

/*
called every Jacobi method iteration

del phi = f
(d/dx^2 + d/dy^2 + ...) phi = f
sum_i ((phi[x+e[i] ] - 2 phi[x] + phi[x-e[i] ]) / dx[i]^2) = f
sum_i (1 / dx[i]^2) phi[x+e[i] ] 
	+ sum_i (-2 / dx[i]^2) phi[x] 
	+ sum_i (1 / dx[i]^2) phi[x-e[i] ]
	= f

a_kk = sum_i (-2 / dx[i]^2)
a_jk = sum_i (1 / dx[i]^2) for j != k

jacobi update:
phi[x,k+1] = (f[x] - sum_i,j!=k (phi[x+e[i],k] / dx[i]^2))
	/ sum_i (-2 / dx[i]^2)
*/
kernel void solvePoissonJacobi<?=poisson.suffix?>(
	global <?=poisson:getPotBufType()?>* UBuf<?
if poisson.stopOnEpsilon then ?>,
	global real* reduceBuf<?
end ?>
) {
<? if not poisson.stopOnEpsilon then ?>
	SETBOUNDS(numGhost,numGhost);
<? else ?>
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		reduceBuf[index] = 0.;
		return;
	}
<? end ?>

	global <?=poisson:getPotBufType()?>* U = UBuf + index;

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

#if 1	//defined(geometry_cartesian)
	real skewSum = (0.
<? for j=0,solver.dim-1 do ?>
		+ volR.s<?=j?> * U[stepsize.s<?=j?>].<?=poisson.potentialField?> / (dx<?=j?> * dx<?=j?>)
		+ volL.s<?=j?> * U[-stepsize.s<?=j?>].<?=poisson.potentialField?> / (dx<?=j?> * dx<?=j?>)
<? end 
?>	) / volAtX;

	const real diag = (0.
<? for j=0,solver.dim-1 do ?>
		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end ?>
	) / volAtX;
#else
	//f_;a^a = g^ab (f_,ab - Gamma^c_ab f,c)
	// = 1/sqrt|g| (sqrt|g| g^ab f_,a)_,b
	//I think I'm gonna use finite-differencing with the second one
	// = 1/sqrt|g|(x) ((sqrt|g| g^ab)(x+h/2) (f(x+h)-f(x))/h - (sqrt|g| g^ab)(x-h/2) (f(x)-f(x-h))/h)/h
	// = 1/sqrt|g|(x) ((sqrt|g| g^ab)(x+h/2) f(x+h)/h - ((sqrt|g| g^ab)(x+h/2) + (sqrt|g| g^ab)(x-h/2)) f(x)/h + (sqrt|g| g^ab)(x-h/2) f(x-h)/h)/h
	//that's a lot of metric matrix multiplies ...

	//sqrt|g| = r, g_ab = diag(1, r^2)
	// = 1/sqrt|g|(x) ((sqrt|g| g^ab)(x+h/2) (f(x+h)-f(x))/h - (sqrt|g| g^ab)(x-h/2) (f(x)-f(x-h))/h)/h
	// = 1/r (
	//		(diag(r+dr/2, 1/(r+dr/2)) (f(r+dr,th)-f(r,th))/dr - diag(r-dr/2, 1/(r-dr/2)) (f(r,th)-f(r-dr,th))/dr)/dr
	//		+ (diag(r, 1/r) (f(r,th+dth)-f(r,th))/dth - diag(r, 1/r) (f(r,th)-f(r,th-dth))/dth)/dth			)
	// = 1/r (
	//		(diag(r+dr/2, 1/(r+dr/2)) (f(r+dr,th)-f(r,th))/dr - diag(r-dr/2, 1/(r-dr/2)) (f(r,th)-f(r-dr,th))/dr)/dr
	//		+ (diag(r, 1/r) (f(r,th+dth)-f(r,th))/dth - diag(r, 1/r) (f(r,th)-f(r,th-dth))/dth)/dth			)

	//ok now it's looking like i'll have to apply one component to the the half-step inverse metric evaluated at one point and the other component at another point  ... bad idea
	// just do the easy way
	//f_;a^b = g^ab (f_,ab - Gamma^c_ab f,c)
	// = diag(1, 1/r^2) * (diag(f,rr; f,theta_theta) - [[0; 1/r f,theta], [1/r f,theta; -r f,r]])
	// f,rr + 1/r^2 f,theta_theta - 1/r f,r
	// (f,r r),r / r + (f,theta r),theta / r
	//so where do we put the r's?
#endif

	real rho = 0;
	<?=poisson:getCalcRhoCode() or ''?>

	real oldU = U-><?=poisson.potentialField?>;
	
	//Gauss-Seidel iteration: x_i = (b_i - A_ij x_j) / A_ii
	real newU = (rho - skewSum) / diag;

<?
if poisson.stopOnEpsilon then
?>	real dU = newU - oldU;
	reduceBuf[index] = dU * dU;
<?
end
?>
	U-><?=poisson.potentialField?> = newU;
}
