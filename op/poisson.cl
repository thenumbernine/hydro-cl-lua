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
local solver = op.solver
local eqn = solver.eqn
?>

//initialize the relaxation solver field 
//this is only called upon solver reset
//each iteration uses the previous iteration's results as the starting point
kernel void initPotential<?=op.name?><?=op.suffix?>(
	global <?=op:getPotBufType()?>* UBuf
) {
<?
local scalar = op.scalar 
local neg = scalar..'_neg'
?>
	SETBOUNDS(numGhost,numGhost);
	global <?=op:getPotBufType()?>* U = UBuf + index;
	<?=scalar?> rho = <?=scalar?>_zero;
	<?=op:getCalcRhoCode() or ''?>
	UBuf[index].<?=op.potentialField?> = <?=neg?>(rho);
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
kernel void solveJacobi<?=op.name?><?=op.suffix?>(
	constant <?=solver.solver_t?>* solver,
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
	volL.s<?=j?> = volume_at(solver, cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?> + .5;
	volR.s<?=j?> = volume_at(solver, cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?>;
<? end ?>
	real volAtX = volume_at(solver, cell_x(i));

<?
local scalar = op.scalar
local zero = scalar..'_zero'
local add3 = scalar..'_add3'
local sub = scalar..'_sub'
local mul = scalar..'_mul'
local lenSq = scalar..'_lenSq'
local real_mul = scalar..'_real_mul'
?>

<? 
if true -- require 'coord.cartesian'.is(solver.coord) 
then 
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
