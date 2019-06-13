local table = require 'ext.table'
local class = require 'ext.class'
local template = require 'template'
local Relaxation = require 'op.relaxation'

local PoissonJacobi = class(Relaxation)

PoissonJacobi.name = 'PoissonJacobi'

PoissonJacobi.solverCodeFile = 'op/poisson.cl'

local poissonJacobiCode = [[
<?
local solver = op.solver
local eqn = solver.eqn

local scalar = op.scalar 
local neg = scalar..'_neg'
local zero = scalar..'_zero'
local add3 = scalar..'_add3'
local sub = scalar..'_sub'
local mul = scalar..'_mul'
local lenSq = scalar..'_lenSq'
local real_mul = scalar..'_real_mul'
?>

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
kernel void solveJacobi<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global real* writeBuf,
	const global <?=op:getPotBufType()?>* UBuf<?
if op.stopOnEpsilon then ?>,
	global real* reduceBuf<? 
end ?>
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost, numGhost)) {
		writeBuf[index] = UBuf[index].<?=op.potentialField?>;
<? if op.stopOnEpsilon then	
?>		reduceBuf[index] = 0;
<? end
?>		return;
	}
	
	real3 x = cell_x(i);

	global <?=op:getPotBufType()?>* U = UBuf + index;

<? for j=0,solver.dim-1 do
?>	real dx<?=j?> = cell_dx<?=j?>(x);
<? end
?>
	
	real3 xInt = x;
	real3 volL, volR;
<? for j=0,solver.dim-1 do 
?>	xInt.s<?=j?> = x.s<?=j?> - .5 * solver->grid_dx.s<?=j?>;
	volL.s<?=j?> = cell_volume(solver, xInt);
	xInt.s<?=j?> = x.s<?=j?> + .5 * solver->grid_dx.s<?=j?>;
	volR.s<?=j?> = cell_volume(solver, xInt);
	xInt.s<?=j?> = x.s<?=j?>;
<? end 
?>	real volAtX = cell_volume(solver, x);

<? if true -- require 'coord.cartesian'.is(solver.coord) 
then 
?>	<?=scalar?> skewSum = <?=scalar?>_zero;

<? for j=0,solver.dim-1 do 
?>	skewSum = <?=add3?>(skewSum,
		<?=real_mul?>(U[solver->stepsize.s<?=j?>].<?=op.potentialField?>, volR.s<?=j?> / (dx<?=j?> * dx<?=j?>)),
		<?=real_mul?>(U[-solver->stepsize.s<?=j?>].<?=op.potentialField?>, volL.s<?=j?> / (dx<?=j?> * dx<?=j?>)));
<? end 
?>	skewSum = <?=real_mul?>(skewSum, 1. / volAtX);

	const real diag = (0.
<? for j=0,solver.dim-1 do 
?>		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end 
?>	) / volAtX;
<? 
else 
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
	<?=scalar?> source = <?=zero?>;
<?=op:getPoissonDivCode() or ''?>

	<?=scalar?> oldU = U-><?=op.potentialField?>;
	
	//Jacobi iteration: x_i = sum i!=j of (b_i - A_ij x_j) / A_ii
	<?=scalar?> newU = <?=real_mul?>(<?=sub?>(source, skewSum), 1. / diag);

	writeBuf[index] = newU;	
<? if op.stopOnEpsilon then
?>	//residual = b_i - A_ij x_j
	<?=scalar?> residual = <?=sub?>(<?=sub?>(source, skewSum), <?=mul?>(diag, U-><?=op.potentialField?>));
	reduceBuf[index] = <?=lenSq?>(residual);
<? end
?>
}
]]

function PoissonJacobi:getCode()
	return table{
		PoissonJacobi.super.getCode(self),
		template(poissonJacobiCode, {op = self}),
		self:getPoissonCode() or '',
	}:concat'\n'
end

return PoissonJacobi
