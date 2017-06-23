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
	assert(a >= 1 and a <= 3, "tried to index sym with "..a..", "..b)
	assert(b >= 1 and b <= 3, "tried to index sym with "..a..", "..b)
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

	real3 partial_alpha_l;		//alpha_,i = partial_alpha_l[i]
	real3 partial_phi_l;		//phi_,i = partial_phi_l[i]
	real3 partial_K_l;			//K,i = partial_K_l[i]
	real3 partial_beta_ul[3];		//beta^i_,j = partial_beta_ul[j][i]
	real tr_partial_beta = 0.;	//beta^i_,i
	sym3 partial_gammaBar_ll[3];	//gammaBar_ij,k = gammaBar[k].ij
<? for i=0,solver.dim-1 do
?>	partial_alpha_l.s<?=i?> = (Up[<?=i?>]->alpha - Um[<?=i?>]->alpha) / (2. * grid_dx<?=i?>);
	partial_phi_l.s<?=i?> = (Up[<?=i?>]->phi - Um[<?=i?>]->phi) / (2. * grid_dx<?=i?>);
	partial_beta_ul[<?=i?>] = real3_scale( real3_sub(Up[<?=i?>]->beta_u, Um[<?=i?>]->beta_u), 1. / (2. * grid_dx<?=i?>));
	partial_gammaBar_ll[<?=i?>] = sym3_scale( sym3_sub(Up[<?=i?>]->gammaBar_ll, Um[<?=i?>]->gammaBar_ll), 1. / (2. * grid_dx<?=i?>));
	partial_K_l.s<?=i?> = (Up[<?=i?>]->K - Um[<?=i?>]->K) / (2. * grid_dx<?=i?>);
	tr_partial_beta += partial_beta_ul[<?=i?>].s<?=i?>;
<? end
for i=solver.dim,2 do
?>	partial_alpha_l.s<?=i?> = 0;
	partial_phi_l.s<?=i?> = 0;
	partial_beta_ul[<?=i?>] = _real3(0,0,0);
	partial_gammaBar_ll[<?=i?>] = (sym3){.s={0,0,0,0,0,0}};
	partial_K_l.s<?=i?> = 0;
<? end
?>

	//alpha_,ij = partial2_alpha_ll.ij
	sym3 partial2_alpha_ll = (sym3){.s={0,0,0,0,0,0}};
<? 
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	if i==j then
?>	partial2_alpha_ll.<?=xij?> = (Up[<?=i-1?>]->alpha - 2. * U->alpha + Um[<?=i-1?>]->alpha) / (grid_dx<?=i-1?> * grid_dx<?=i-1?>);
<?	end
?>	partial2_alpha_ll.<?=xij?> = (
		U[index + stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].alpha 
		- U[index - stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].alpha 
		- U[index + stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].alpha
		+ U[index - stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].alpha 
	) / (4 * grid_dx<?=i-1?> * grid_dx<?=j-1?>);
<? 
end
?>


	real exp_4phi = exp(4. * U->phi);
	real exp_neg4phi = 1. / exp_4phi;

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma = sym3_scale(U->gammaBar_ll, exp_4phi);

	//gammaBar^ij = inv gammaBar_ij
	//shouldn't this be 1 anyways?
	real det_gammaBar = sym3_det(U->gammaBar_ll);
	sym3 gammaBar_uu = sym3_inv(det_gammaBar, U->gammaBar_ll);

	//gamma^ij = inv gamma_ij	
	//sym3 gammaU = sym3_inv(gamma);
	//gamma^ij = exp(-4 phi) gammaBar^ij
	sym3 gamma_uu = sym3_scale(gammaBar_uu, exp_neg4phi);

	sym3 connBar_lll[3];	//connBar_ijk = connBar_lll[i].jk
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll[<?=i-1?>].<?=xjk?> = .5 * (partial_gammaBar_ll[<?=k-1?>].<?=sym(i,j)?> + partial_gammaBar_ll[<?=j-1?>].<?=sym(i,k)?> - partial_gammaBar_ll[<?=i-1?>].<?=xjk?>);
<?	end
end
?>	

	sym3 connBar_ull[3];	//connBar^i_jk = connBar_ull[i].jk
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
?>	connBar_ull[<?=i-1?>].<?=xjk?> = 0. <?
		for l,xl in ipairs(xNames) do
?> + gammaBar_uu.<?=sym(i,l)?> * connBar_lll[<?=l-1?>].<?=xjk?>
<?		end
?>;
<?	end
end
?>

	sym3 conn_ull[3];	//conn^i_jk = conn_ull[i].jk
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	conn_ull[<?=i-1?>].<?=xjk?> = connBar_ull[<?=i-1?>].<?=xjk?> <?
		if i==j then
?> + 2 * partial_phi_l.<?=xNames[k]?><?
		end
		if i==k then
?> + 2 * partial_phi_l.<?=xNames[j]?><?
		end
		for l,xl in ipairs(xNames) do
?> - 2 * U->gammaBar_ll.<?=xjk?> * gammaBar_uu.<?=sym(i,l)?> * partial_phi_l.<?=xl?><?
		end
?>;
<?	end
end
?>

	//D2_alpha.ij = D_i D_j alpha = partial_i partial_j alpha - conn^k_ij partial_k alpha
	sym3 D2_alpha = (sym3){.s={0,0,0,0,0,0}};
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	D2_alpha.<?=xij?> = partial2_alpha_ll.<?=xij?> <?
	for k,xk in ipairs(xNames) do 
?> - conn_ull[<?=k-1?>].<?=xij?> * partial_alpha_l.<?=xk?><?
	end ?>;<?
end
?>

	//Q = f(alpha) K
	real Q = calc_f(U->alpha) * U->K;
	
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q + real3_dot(partial_alpha_l, U->beta_u);

	//manuall update elsewhere?
	//deriv->beta += _real3(0,0,0);

	//B&S 11.50
	//phi,t = -1/6 alpha K + beta^i phi,i + 1/6 beta^i_,i
	deriv->phi += -U->alpha * U->K / 6. + real3_dot(U->beta_u, partial_phi_l) + tr_partial_beta / 6.;

	//B&S 11.51
	//gammaBar_ij,t = -2 alpha ATilde_ij + beta^k gammaBar_ij,k + gammaBar_ik beta^k_,j + gammaBar_kj beta^k_,i - 2/3 gammaBar_ij beta^k_,k
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->gammaBar_ll.<?=xij?> += -2 * U->alpha * U->ATilde_ll.<?=xij?>	//-2 alpha ATilde_ij 
<? 	for k,xk in ipairs(xNames) do
?>		+ partial_gammaBar_ll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>		//+ beta^k gammaBar_ij,k 
		+ U->gammaBar_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?> 	//+ gammaBar_ik beta^k_,j
		+ U->gammaBar_ll.<?=sym(j,k)?> * partial_beta_ul[<?=j-1?>].<?=xi?>	//+ gammaBar_kj beta^k_,i 
<? 	end
?>		- 2./3. * U->gammaBar_ll.<?=sym(i,j)?> * tr_partial_beta;		//- 2/3 gammaBar_ij beta^k_,k
<?
end
?>

	mat3 ATilde_lu = sym3_sym3_mul(U->ATilde_ll, gammaBar_uu);
	mat3 ATildeSq_lu = mat3_mat3_mul(ATilde_lu, ATilde_lu); 
	real ATildeSq = mat3_trace(ATildeSq_lu);

	real S = sym3_dot(U->S_ll, gamma_uu);

	//B&S 11.52
	//K_,t = -gamma^ij D_ij alpha + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	deriv->K += -sym3_dot(gammaBar_uu, D2_alpha) + U->alpha * (ATildeSq + U->K * U->K / 3.) + 4 * M_PI * U->alpha * (U->rho + S) + real3_dot(U->beta_u, partial_K_l);

}
