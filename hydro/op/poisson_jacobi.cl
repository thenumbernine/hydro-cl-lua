<?
local scalar = op.scalar 
local neg = scalar.."_neg"
local zero = scalar.."_zero"
local add3 = scalar.."_add3"
local sub = scalar.."_sub"
local mul = scalar.."_mul"
local lenSq = scalar.."_lenSq"
local real_mul = scalar.."_real_mul"
?>

//// MODULE_NAME: <?=solveJacobi?>
//// MODULE_DEPENDS: <?=cell_dx_i?> <?=table.concat(op.codeDepends or {}, ' ')?>
/*
called every Jacobi method iteration

reads from UBuf
writes to writeBuf
optionally if stopOnEpsilon is enabled, writes residual component-wise squared to reduceBuf

Δφ = f
(∂/∂x^2 + ∂/∂y^2 + ... (plus connections...)) φ = f
Σ_i ((φ[x+e[i] ] - 2 φ[x] + φ[x-e[i] ]) / dx[i]^2) = f
Σ_i (1 / dx[i]^2) φ[x+e[i] ] 
	+ Σ_i (-2 / dx[i]^2) φ[x] 
	+ Σ_i (1 / dx[i]^2) φ[x-e[i] ]
	= f

a_kk = Σ_i (-2 / dx[i]^2)
a_jk = Σ_i (1 / dx[i]^2) for j≠k

jacobi update:
φ[x,k+1] = (f[x] - Σ_i,j≠k (φ[x+e[i],k] / dx[i]^2))
	/ Σ_i (-2 / dx[i]^2)

input is poisson source divergence, in 1/s^2
output is potentialField, in m^2/s^2
*/
kernel void <?=solveJacobi?>(
	constant <?=solver_t?> const * const solver,
	global real * const writeBuf,
	global <?=op:getPotBufType()?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf<?
if op.stopOnEpsilon then ?>,
	global real * const reduceBuf<? 
end ?>
) {
	<?=SETBOUNDS?>(0,0);
	if (<?=OOB?>(solver->numGhost, solver->numGhost)) {
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
?>	real const dx<?=j?> = cell_dx<?=j?>(x);
<? end
?>
	
	real3 xInt = x;
	real3 volL, volR;
<? for j=0,solver.dim-1 do 
?>	xInt.s<?=j?> = x.s<?=j?> - .5 * solver->grid_dx.s<?=j?>;
	// TODO instead of volume_intL as the avg between two cell volumes, and then divide by dx to get the face, instead, just store the face.
	real const volume_intL<?=j?> = .5 * (cell->volume + cell[-solver->stepsize.s<?=j?>].volume);
	volL.s<?=j?> = volume_intL<?=j?>;
	xInt.s<?=j?> = x.s<?=j?> + .5 * solver->grid_dx.s<?=j?>;
	// TODO instead of volume_intL as the avg between two cell volumes, and then divide by dx to get the face, instead, just store the face.
	real const volume_intR<?=j?> = .5 * (cell->volume + cell[solver->stepsize.s<?=j?>].volume);
	volR.s<?=j?> = volume_intR<?=j?>;
	xInt.s<?=j?> = x.s<?=j?>;
<? end 
?>	real const volAtX = cell->volume;

<? 
--[=[
volume-weighted ... however volume-weighted laplace beltrami looks like this:
Δφ = 1/√|g| ( √|g| g^ij φ_,j )_,i
...so I should be sampling √|g| g^ij φ_,j at the + and - on each dimension
= 1/√|g| (
	[( √|g| g^ij φ_,j )|(x+dx_i) 
	- ( √|g| g^ij φ_,j )|(x-dx_i)] / (2*dx_i)
)
= 1/√|g|(x) (
	[
		( √|g|(x+dx_i) g^ij(x+dx_i) 
			* (φ(x+dx_i+dx_j) - φ(x+dx_i-dx_j)) / (2*dx_j) 
		)|(x+dx_i) 
		- ( √|g|(x+dx_i) g^ij(x+dx_i) 
			* (φ(x-dx_i+dx_j) - φ(x-dx_i-dx_j)) / (2*dx_j) 
		)|(x-dx_i)
	] / (2*dx_i)
)
--]=]
if true 
then 
?>	
	<?=scalar?> skewSum = <?=scalar?>_zero;

<? for j=0,solver.dim-1 do 
?>	skewSum = <?=add3?>(skewSum,
		<?=real_mul?>(U[solver->stepsize.s<?=j?>].<?=op.potentialField?>, volR.s<?=j?> / (dx<?=j?> * dx<?=j?>)),
		<?=real_mul?>(U[-solver->stepsize.s<?=j?>].<?=op.potentialField?>, volL.s<?=j?> / (dx<?=j?> * dx<?=j?>)));
<? end 
?>	skewSum = <?=real_mul?>(skewSum, 1. / volAtX);

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
f_;a^a = g^ab (f_,ab - Γ^c_ab f,c)
 = 1/√|g| (√|g| g^ab f_,a)_,b
I think I'm gonna use finite-differencing with the second one
 = 1/√|g|(x) ((√|g| g^ab)(x+h/2) (f(x+h)-f(x))/h - (√|g| g^ab)(x-h/2) (f(x)-f(x-h))/h)/h
 = 1/√|g|(x) ((√|g| g^ab)(x+h/2) f(x+h)/h - ((√|g| g^ab)(x+h/2) + (√|g| g^ab)(x-h/2)) f(x)/h + (√|g| g^ab)(x-h/2) f(x-h)/h)/h

or for arbitrary tensors:
(wiki says (T_;ab - T_;ba) g^ab)
t^i1..ip_j1..jq^;a_;a
= (t^i1..ip_j1..jq_,a 
	+ Σ_I=1..p Γ^iI_k_a t^i1..k..ip_j1..jq 
	- Σ_J=1..q Γ^k_jJ_a t^i1..ip_j1..k..jq )_;b g^ab
= (t^i1..ip_j1..jq_,ab 
	+ (Σ_I=1..p Γ^iI_k_a t^i1..k..ip_j1..jq)_;b
	- (Σ_J=1..q Γ^k_jJ_a t^i1..ip_j1..k..jq)_;b ) g^ab
= (t^i1..ip_j1..jq_,ab 
	+ Σ_I=1..p Γ^iI_k_a,b t^i1..k..ip_j1..jq 
	- Σ_J=1..q Γ^k_jJ_a,b t^i1..ip_j1..k..jq 
	+ Σ_I=1..p Γ^Ii_k_a t^i1..k..ip_j1..jq,b
	- Σ_J=1..q Γ^k_jJ_a t^i1..ip_j1..k..jq,b 
) g^ab

*/
<?
end
?>

	//source is 4 π G ρ δ(x) is the laplacian of the gravitational potential field, which is integrated across discretely here
	//in units of 1/s^2
	<?=scalar?> source = <?=zero?>;
<?=op:getPoissonDivCode() or ""?>

	<?=scalar?> oldU = U-><?=op.potentialField?>;
	
	//Jacobi iteration: x_i = Σ i≠j of (b_i - A_ij x_j) / A_ii
	<?=scalar?> newU = <?=real_mul?>(<?=sub?>(source, skewSum), 1. / diag);

	writeBuf[index] = newU;	
<? if op.stopOnEpsilon then
?>	//residual = b_i - A_ij x_j
	<?=scalar?> residual = <?=sub?>(<?=sub?>(source, skewSum), <?=mul?>(diag, U-><?=op.potentialField?>));
	reduceBuf[index] = <?=lenSq?>(residual);
<? end
?>
}
