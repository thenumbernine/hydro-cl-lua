<?
local solver = op.solver
?>
/*
this will use an inverse Laplace equation solver to converge the shift 

using eqns Alubierre's book 4.3.14, 4.3.15 ...

D^2 beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j = 2 D_j (alpha A^ij)
A^ij = K^ij - 1/3 K gamma^ij
solve for beta^i

here's my thoughts ...
1) calculate Gamma^i_jk
2) calculate R^i_j
3) solve inverse linear system

this means the equation will look like ...
gamma^jk D_j D_k beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j = 2 D_j (alpha A^ij)
gamma^jk D_j (beta^i_,k + Gamma^i_lk beta^l) + 1/3 gamma^ik D_k (beta^j_,j + Gamma^j_lj beta^l) + R^i_j beta^j = 2 (D_j alpha A^ij + alpha D_j A^ij)
gamma^jk (D_j beta^i_,k + D_j Gamma^i_lk beta^l + Gamma^i_lk D_j beta^l) + 1/3 gamma^ik (D_k beta^j_,j + D_k Gamma^j_lj beta^l + Gamma^j_lj D_k beta^l) + R^i_j beta^j = 2 (alpha_,j A^ij + alpha (A^ij_,j + Gamma^i_kj A^kj + Gamma^j_kj A^ik))

(
	+ gamma^jk beta^i_,jk 
	+ 1/3 gamma^ik beta^j_,jk 
	+ gamma^jk Gamma^i_lk beta^l_,j 
	- gamma^jk Gamma^l_kj beta^i_,l 
	+ 1/3 gamma^ik Gamma^j_lj beta^l_,k 
	+ 1/3 gamma^ik beta^l Gamma^j_lj,k 
	+ gamma^jk beta^l Gamma^i_lk,j 
	+ gamma^jk Gamma^i_lj beta^l_k 
	+ gamma^jk beta^l Gamma^i_mj Gamma^m_lk 
	- gamma^jk beta^l Gamma^m_lj Gamma^i_mk 
	- gamma^jk beta^l Gamma^m_kj Gamma^i_lm 
	+ gamma^jk Gamma^i_lk Gamma^l_mj beta^m
	+ 1/3 gamma^ik beta^l Gamma^j_mk Gamma^m_lj 
	- 1/3 gamma^ik beta^l Gamma^m_lk Gamma^j_mj 
	- 1/3 gamma^ik beta^l Gamma^m_jk Gamma^j_lm
	+ 1/3 gamma^ik Gamma^j_lj Gamma^l_mk beta^m
	+ R^i_j beta^j 
) = 2 alpha_,j (K^ij - 1/3 gamma^ij K) 
	+ 2 alpha (K^ij - 1/3 gamma^ij K)_,j 
	+ 2 alpha Gamma^i_kj (K^kj - 1/3 gamma^kj K) 
	+ 2 alpha Gamma^j_kj (K^ik - 1/3 gamma^ik K)

(
	+ beta^i_,jk gamma^jk
	+ 1/3 beta^j_,jk gamma^ik
	- 4 d^ij_l beta^l_,j
	+ 4 d_l^ji beta^l_,j 
	+ 4 d^ji_l beta^l_,j 
	- 4 d_j^jl beta^i_,l
	+ 2 d^lj_j beta^i_,l
	+ 2/3 gamma^ik d_lj^j beta^l_,k 
	+ 4/3 gamma^ik beta^l gamma^mj d_ljm,k 
	- 2/3 gamma^ik beta^l gamma^jm d_mjl,k
	+ 2 gamma^jk beta^l gamma^im d_klm,j 
	+ 2 gamma^jk beta^l gamma^im d_lkm,j 
	- 2 gamma^jk beta^l gamma^im d_mlk,j
	
	- 8/3 d^ijk d_ljk beta^l
	+ 4/3 d^ijk d_kjl beta^l 
	
	- 8 d^jmi d_lmj beta^l
	+ 4 d^jmi d_ljm beta^l
	
	- 4 d^jmi d_jml beta^l
	+ 4 d^jmi d_mjl beta^l
	
	- 4 gamma^jk (d_kj^m + d_jk^m - d^m_kj) (d_lm^i + d_ml^i - d^i_lm) beta^l
	+ 4 gamma^jk (d_lk^i + d_kl^i - d^i_lk) (d_mj^l + d_jm^l - d^l_mj) beta^m
	+ 4/3 gamma^ik (d_mk^j + d_km^j - d^j_mk) (d_lj^m + d_jl^m - d^m_lj) beta^l
	- 4/3 gamma^ik (d_lk^m + d_kl^m - d^m_lk) (d_mj^j + d_jm^j - d^j_mj) beta^l
	- 4/3 gamma^ik (d_jk^m + d_kj^m - d^m_jk) (d_lm^j + d_ml^j - d^j_lm) beta^l
	+ 4/3 gamma^ik (d_lj^j + d_jl^j - d^j_lj) (d_mk^l + d_km^l - d^l_mk) beta^m
	
	+ R^i_j beta^j 
) = 2 alpha_,j K^ij 
	- 2/3 alpha_,j gamma^ij K
	+ 2 alpha gamma^im K_mn,j gamma^nj 
	- 2/3 alpha gamma^ij K_mn,j gamma^mn 
	
	- 4 alpha d_j^im K_m^j
	+ 4 alpha d_kj^i K^kj 
	+ 4 alpha d_jk^i K^kj 
	- 4 alpha d^i_kj K^kj
	+ 4/3 alpha d^imn K_mn
	
	- 4 alpha d_j^nj K^i_n
	+ 4 alpha d_kj^j K^ik 
	+ 4 alpha d_jk^j K^ik 
	- 4 alpha d^j_kj K^ik
	
	+ 4/3 alpha d_j^ij gamma^kl K_kl 
	- 4/3 alpha d_kj^i gamma^kj K_mn gamma^mn
	- 4/3 alpha d_jk^i gamma^kj K_mn gamma^mn
	+ 4/3 alpha d^i_kj gamma^kj K_mn gamma^mn
	- 4/3 alpha d_kj^j gamma^ik K_mn gamma^mn 
	- 4/3 alpha d_jk^j gamma^ik K_mn gamma^mn 
	+ 4/3 alpha d^j_kj gamma^ik K_mn gamma^mn


Gamma_ijk = d_jki + d_kji - d_ijk
gamma^ij_,k = -gamma^im gamma_mn,k gamma^nj = -2 d_k^ij

*/

<?
local makePartials = require 'eqn.makepartial'
local derivOrder = 2 * solver.numGhost
local makePartial = function(...) return makePartials.makePartial(derivOrder, solver, ...) end
local makePartial2 = function(...) return makePartials.makePartial2(derivOrder, solver, ...) end
?>

kernel void solveMinimalDistortionEllipticShift<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global real* writeBuf,
	global <?=op:getPotBufType()?>* UBuf<?
if op.stopOnEpsilon then ?>,
	global real* reduceBuf<?
end ?>
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		writeBuf[index] = 0.;
<? if op.stopOnEpsilon then 
?>		reduceBuf[index] = 0.;
<? end 
?>		return;
	}

	global <?=op:getPotBufType()?>* U = UBuf + index;

<? for j=0,solver.dim-1 do ?>
	real dx<?=j?> = cell_dx<?=j?>(x);
<? end ?>

	real3 intIndex = _real3(i.x, i.y, i.z);
	real3 volL, volR;
<? for j=0,solver.dim-1 do ?>
	intIndex.s<?=j?> = i.s<?=j?> - .5;
	volL.s<?=j?> = coord_sqrt_det_g(solver, cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?> + .5;
	volR.s<?=j?> = coord_sqrt_det_g(solver, ell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?>;
<? end ?>
	real volAtX = coord_sqrt_det_g(solver, cell_x(i));

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
		<?=real_mul?>(U[solver->stepsize.s<?=j?>].<?=op.potentialField?>, volR.s<?=j?> / (dx<?=j?> * dx<?=j?>)),
		<?=real_mul?>(U[-solver->stepsize.s<?=j?>].<?=op.potentialField?>, volL.s<?=j?> / (dx<?=j?> * dx<?=j?>)));
<? end ?>
	skewSum = <?=real_mul?>(skewSum, 1. / volAtX);

	const real diag = (0.
<? for j=0,solver.dim-1 do ?>
		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end ?>
	) / volAtX;


	<?=scalar?> source = <?=zero?>;
<?=op:getPoissonDivCode() or ''?>

	<?=scalar?> oldU = U-><?=op.potentialField?>;
	
	//Jacobi iteration: x_i = (b_i - A_ij x_j) / A_ii
	<?=scalar?> newU = <?=real_mul?>(<?=sub?>(source, skewSum), 1. / diag);

	writeBuf[index] = newU;
<? if op.stopOnEpsilon then
?>	reduceBuf[index] = <?=lenSq?>(newU);
<? end 
?>
}
