<?
local table = require 'ext.table'
local from3x3to6_table = {
	{1, 2, 3},
	{2, 4, 5},
	{3, 5, 6},
}
local function from3x3to6(i,j)
	return from3x3to6_table[i][j]
end
local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
local function from6to3x3(i)
	return table.unpack(from6to3x3_table[i])
end

local function sym(a,b)
	if a > b then a,b = b,a end
	return xNames[a]..xNames[b]
end
?>
kernel void calcDeriv(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
		
	const global <?=eqn.cons_t?>* Up[dim];
	const global <?=eqn.cons_t?>* Um[dim];
	for (int i = 0; i < dim; ++i) {
		Up[i] = U + index + stepsize[i];
		Um[i] = U + index - stepsize[i];
	}

	real3 dx_alpha;			//alpha,i = dx_alpha[i]
	real3 dx_phi;			//phi,i = dx_phi[i]
	real3 dx_beta[3];		//beta^i_,j = dx_beta[j][i]
	real tr_dx_beta = 0.;	//beta^i_,i
	sym3 dx_gammaTilde[3];
<? for i=0,solver.dim-1 do
?>	dx_alpha.s<?=i?> = (Up[<?=i?>]->alpha - Um[<?=i?>]->alpha) / (2. * dx<?=i?>_at(i));
	dx_phi.s<?=i?> = (Up[<?=i?>]->phi - Um[<?=i?>]->phi) / (2. * dx<?=i?>_at(i));
	dx_beta[<?=i?>] = real3_scale( real3_sub(Up[<?=i?>]->beta, Um[<?=i?>]->beta), 1. / (2. * dx<?=i?>_at(i) ));
	dx_gammaTilde[<?=i?>] = sym3_scale( sym3_sub(Up[<?=i?>]->gammaTilde, Um[<?=i?>]->gammaTilde), 1. / (2. * dx<?=i?>_at(i) ));
	tr_dx_beta += dx_beta[<?=i?>].s<?=i?>;
<? end
for i=solver.dim,2 do
?>	dx_alpha.s<?=i?> = 0;
	dx_phi.s<?=i?> = 0;
	dx_beta[<?=i?>] = _real3(0,0,0);
	dx_gammaTilde[<?=i?>] = (sym3){.s={0,0,0,0,0,0}};
<? end
?>

	real exp_4phi = exp(4. * U->phi);
	real exp_neg4phi = 1. / exp_4phi;

	//gamma_ij = exp(4 phi) gammaTilde_ij
	sym3 gamma = sym3_scale(U->gammaTilde, exp_4phi);

	//gammaTilde^ij = inv gammaTilde_ij
	real det_gammaTilde = sym3_det(U->gammaTilde);
	sym3 gammaTildeInv = sym3_inv(det_gammaTilde, U->gammaTilde);

	//gamma^ij = inv gamma_ij	
	//sym3 gammaU = sym3_inv(gamma);
	//gamma^ij = exp(-4 phi) gammaTilde^ij
	sym3 gammaInv = sym3_scale(gammaTildeInv, exp_neg4phi);

	//Q = f(alpha) tr K
	real Q = calc_f(U->alpha) * U->tr_K;
	
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q + real3_dot(dx_alpha, U->beta);

	//manuall update elsewhere?
	//deriv->beta += _real3(0,0,0);

	//B&S 11.50
	//phi,t = -1/6 alpha tr K + beta^i phi,i + 1/6 beta^i_,i
	deriv->phi += -U->alpha * U->tr_K / 6. + real3_dot(U->beta, dx_phi) + tr_dx_beta / 6.;

	//B&S 11.51
	//gammaTilde_ij,t = -2 alpha ATilde_ij + beta^k gammaTilde_ij,k + gammaTilde_ik beta^k_,j + gammaTilde_kj beta^k_,i - 2/3 gammaTilde_ij beta^k_,k
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->gammaTilde.<?=xij?> += -2 * U->alpha * U->ATilde.<?=xij?>	//-2 alpha ATilde_ij 
<? 	for k,xk in ipairs(xNames) do
		local ik = from3x3to6(i,k)
		local jk = from3x3to6(j,k)
?>		+ dx_gammaTilde[<?=k?>].<?=xij?> * U->beta.<?=xk?>		//+ beta^k gammaTilde_ij,k 
		+ U->gammaTilde.<?=sym(i,k)?> * dx_beta[<?=j+1?>].<?=xk?> 	//+ gammaTilde_ik beta^k_,j
		+ U->gammaTilde.<?=sym(j,k)?> * dx_beta[<?=j+1?>].<?=xi?>	//+ gammaTilde_kj beta^k_,i 
<? 	end
?>		- 2./3. * U->gammaTilde.<?=sym(i,j)?> * tr_dx_beta;		//- 2/3 gammaTilde_ij beta^k_,k
<?
end
?>
}
