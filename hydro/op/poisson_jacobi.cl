<?
local scalar = op.scalar 
?>

//// MODULE_NAME: <?=solveJacobi?>
//// MODULE_DEPENDS: <?=cell_dx_i?> <?=table.concat(op.codeDepends or {}, ' ')?>
/*
called every Jacobi method iteration

reads from UBuf
writes to writeBuf
optionally if stopOnEpsilon is enabled, writes residual component-wise squared to reduceBuf

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

input is poisson source divergence, in 1/s^2
output is potentialField, in m^2/s^2
*/
kernel void <?=solveJacobi?>(
	constant <?=solver_t?> const * const psolver,
	global real * const writeBuf,
	global <?=op:getPotBufType()?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf<?
if op.stopOnEpsilon then ?>,
	global real * const reduceBuf<? 
end ?>
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(0,0);
	if (<?=OOB?>(solver.numGhost, solver.numGhost)) {
		writeBuf[index] = UBuf[index].<?=op.potentialField?>;
<? if op.stopOnEpsilon then	
?>		reduceBuf[index] = 0;
<? end
?>		return;
	}
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;

	global <?=op:getPotBufType()?> const * const U = UBuf + index;

<? for j=0,solver.dim-1 do
?>	real const dx<?=j?> = cell_dx<?=j?>(solver, x);
<? end
?>
	
	real3 xInt = x;
	real3 volL, volR;
<? for j=0,solver.dim-1 do 
?>	xInt.s<?=j?> = x.s<?=j?> - .5 * solver.grid_dx.s<?=j?>;
	// TODO instead of volume_intL as the avg between two cell volumes, and then divide by dx to get the face, instead, just store the face.
	real const volume_intL<?=j?> = .5 * (cell->volume + cell[-solver.stepsize.s<?=j?>].volume);
	volL.s<?=j?> = volume_intL<?=j?>;
	xInt.s<?=j?> = x.s<?=j?> + .5 * solver.grid_dx.s<?=j?>;
	// TODO instead of volume_intL as the avg between two cell volumes, and then divide by dx to get the face, instead, just store the face.
	real const volume_intR<?=j?> = .5 * (cell->volume + cell[solver.stepsize.s<?=j?>].volume);
	volR.s<?=j?> = volume_intR<?=j?>;
	xInt.s<?=j?> = x.s<?=j?>;
<? end 
?>	real const volAtX = cell->volume;

<? 
--[=[
volume-weighted ... however volume-weighted laplace beltrami looks like this:
lap phi = 1/sqrt|g| ( sqrt|g| g^ij phi_,j )_,i
...so I should be sampling sqrt|g| g^ij phi_,j at the + and - on each dimension
= 1/sqrt|g| (
	[( sqrt|g| g^ij phi_,j )|(x+dx_i) 
	- ( sqrt|g| g^ij phi_,j )|(x-dx_i)] / (2*dx_i)
)
= 1/sqrt|g|(x) (
	[
		( sqrt|g|(x+dx_i) g^ij(x+dx_i) 
			* (phi(x+dx_i+dx_j) - phi(x+dx_i-dx_j)) / (2*dx_j) 
		)|(x+dx_i) 
		- ( sqrt|g|(x+dx_i) g^ij(x+dx_i) 
			* (phi(x-dx_i+dx_j) - phi(x-dx_i-dx_j)) / (2*dx_j) 
		)|(x-dx_i)
	] / (2*dx_i)
)
--]=]
if true 
then 
?>	
	<?=scalar?> skewSum = {};

<? for j=0,solver.dim-1 do 
?>	skewSum += U[solver.stepsize.s<?=j?>].<?=op.potentialField?> * (volR.s<?=j?> / (dx<?=j?> * dx<?=j?>))
			+ U[-solver.stepsize.s<?=j?>].<?=op.potentialField?> * (volL.s<?=j?> / (dx<?=j?> * dx<?=j?>));
<? end 
?>	skewSum /= volAtX;

	real const diag = (0.
<? for j=0,solver.dim-1 do 
?>		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end 
?>	) / volAtX;

<? 
else 	-- not cartesian
?>
/*
for scalars:
f_;a^a = g^ab (f_,ab - Gamma^c_ab f,c)
 = 1/sqrt|g| (sqrt|g| g^ab f_,a)_,b
I think I'm gonna use finite-differencing with the second one
 = 1/sqrt|g|(x) ((sqrt|g| g^ab)(x+h/2) (f(x+h)-f(x))/h - (sqrt|g| g^ab)(x-h/2) (f(x)-f(x-h))/h)/h
 = 1/sqrt|g|(x) ((sqrt|g| g^ab)(x+h/2) f(x+h)/h - ((sqrt|g| g^ab)(x+h/2) + (sqrt|g| g^ab)(x-h/2)) f(x)/h + (sqrt|g| g^ab)(x-h/2) f(x-h)/h)/h

or for arbitrary tensors:
(wiki says (T_;ab - T_;ba) g^ab)
t^i1..ip_j1..jq^;a_;a
= (t^i1..ip_j1..jq_,a 
	+ Sum_I=1..p Conn^iI_k_a t^i1..k..ip_j1..jq 
	- Sum_J=1..q Conn^k_jJ_a t^i1..ip_j1..k..jq )_;b g^ab
= (t^i1..ip_j1..jq_,ab 
	+ (Sum_I=1..p Conn^iI_k_a t^i1..k..ip_j1..jq)_;b
	- (Sum_J=1..q Conn^k_jJ_a t^i1..ip_j1..k..jq)_;b ) g^ab
= (t^i1..ip_j1..jq_,ab 
	+ Sum_I=1..p Conn^iI_k_a,b t^i1..k..ip_j1..jq 
	- Sum_J=1..q Conn^k_jJ_a,b t^i1..ip_j1..k..jq 
	+ Sum_I=1..p Conn^Ii_k_a t^i1..k..ip_j1..jq,b
	- Sum_J=1..q Conn^k_jJ_a t^i1..ip_j1..k..jq,b 
) g^ab

*/
<?
end
?>

	//source is 4 pi G rho delta(x) is the laplacian of the gravitational potential field, which is integrated across discretely here
	//in units of 1/s^2
	<?=scalar?> source = {};
<?=op:getPoissonDivCode() or ""?>

	<?=scalar?> oldU = U-><?=op.potentialField?>;
	
	//Jacobi iteration: x_i = sum i!=j of (b_i - A_ij x_j) / A_ii
	<?=scalar?> newU = (source - skewSum) * (1. / diag);

	writeBuf[index] = newU;	
<? if op.stopOnEpsilon then
?>	//residual = b_i - A_ij x_j
	<?=scalar?> residual = source - (skewSum + diag * U-><?=op.potentialField?>);
	reduceBuf[index] = lenSq(residual);
<? end
?>
}
