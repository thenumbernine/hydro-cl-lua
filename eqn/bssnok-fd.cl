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
	assert(a >= 1 and a <= 3, "tried to index sym with "..tostring(a)..", "..tostring(b))
	assert(b >= 1 and b <= 3, "tried to index sym with "..tostring(a)..", "..tostring(b))
	if a > b then a,b = b,a end
	return xNames[a]..xNames[b]
end
?>

kernel void constrainU(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	/*
	det(gammaBar_ij) 
	= det(gamma^-1/3 gamma_ij)
	= gamma^-1 gamma
	= 1
	*/
	real det_gammaBar = sym3_det(U->gammaBar_ll);
	real cbrt_det_gammaBar = 1. / cbrt(det_gammaBar);
	U->gammaBar_ll = sym3_scale(U->gammaBar_ll, cbrt_det_gammaBar);

	sym3 gammaBar_uu = sym3_inv(1., U->gammaBar_ll);

	/*
	tr(A_ij)
	= tr(K_ij - 1/3 gamma_ij K)
	= gamma^ij K_ij - 1/3 gamma^ij gamma_ij K
	= K - 1/3 3 K
	= 0

	tr(ATilde_ij) = 3 psi^-4 tr(A_ij) = 3 psi^-4 * 0 
	= 0
	*/
	real tr_ATilde = sym3_dot(gammaBar_uu, U->ATilde_ll);
	U->ATilde_ll = sym3_sub(
		U->ATilde_ll, 
		sym3_scale(U->gammaBar_ll, -1./3. * tr_ATilde));
}


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

	real3 partial_alpha_l;			//alpha_,i = partial_alpha_l.i
	real3 partial_phi_l;			//phi_,i = partial_phi_l.i
	real3 partial_K_l;				//K,i = partial_K_l.i
	real3 partial_beta_ul[3];		//beta^i_,j = partial_beta_ul[j].i
	real3 partial_connBar_ul[3];	//connBar^i_,j = partial_connBar_ul[j].i
	real tr_partial_beta = 0.;		//beta^i_,i
	sym3 partial_gammaBar_lll[3];	//gammaBar_ij,k = partial_gammaBar[k].ij
	sym3 partial_ATilde_lll[3];		//ATilde_ij,k = partial_ATilde_lll[k].ij
<? for i=1,solver.dim do
	local xi = xNames[i]
?>	partial_alpha_l.<?=xi?> = (Up[<?=i-1?>]->alpha - Um[<?=i-1?>]->alpha) / (2. * grid_dx<?=i-1?>);
	partial_phi_l.<?=xi?> = (Up[<?=i-1?>]->phi - Um[<?=i-1?>]->phi) / (2. * grid_dx<?=i-1?>);
	partial_beta_ul[<?=i-1?>] = real3_scale( real3_sub(Up[<?=i-1?>]->beta_u, Um[<?=i-1?>]->beta_u), 1. / (2. * grid_dx<?=i-1?>));
	partial_connBar_ul[<?=i-1?>] = real3_scale( real3_sub(Up[<?=i-1?>]->connBar_u, Um[<?=i-1?>]->connBar_u), 1. / (2. * grid_dx<?=i-1?>));
	partial_gammaBar_lll[<?=i-1?>] = sym3_scale( sym3_sub(Up[<?=i-1?>]->gammaBar_ll, Um[<?=i-1?>]->gammaBar_ll), 1. / (2. * grid_dx<?=i-1?>));
	partial_ATilde_lll[<?=i-1?>] = sym3_scale( sym3_sub(Up[<?=i-1?>]->ATilde_ll, Um[<?=i-1?>]->ATilde_ll), 1. / (2. * grid_dx<?=i-1?>));
	partial_K_l.<?=xi?> = (Up[<?=i-1?>]->K - Um[<?=i-1?>]->K) / (2. * grid_dx<?=i-1?>);
	tr_partial_beta += partial_beta_ul[<?=i-1?>].<?=xi?>;
<? end
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>	partial_alpha_l.<?=xi?> = 0;
	partial_phi_l.<?=xi?> = 0;
	partial_beta_ul[<?=i-1?>] = _real3(0,0,0);
	partial_connBar_ul[<?=i-1?>] = _real3(0,0,0);
	partial_gammaBar_lll[<?=i-1?>] = (sym3){.s={0,0,0,0,0,0}};
	partial_ATilde_lll[<?=i-1?>] = (sym3){.s={0,0,0,0,0,0}};
	partial_K_l.<?=xi?> = 0;
<? end
?>

	//alpha_,ij = partial2_alpha_ll.ij
	sym3 partial2_alpha_ll = (sym3){.s={0,0,0,0,0,0}};
<? 
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	if i==j then
?>	partial2_alpha_ll.<?=xij?> = (Up[<?=i-1?>]->alpha - 2. * U->alpha + Um[<?=i-1?>]->alpha) / (grid_dx<?=i-1?> * grid_dx<?=i-1?>);
<?	else
?>	partial2_alpha_ll.<?=xij?> = (
		U[index + stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].alpha 
		- U[index - stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].alpha 
		- U[index + stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].alpha
		+ U[index - stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].alpha 
	) / (4 * grid_dx<?=i-1?> * grid_dx<?=j-1?>);
<? 	end
end
?>
	//phi_,ij = partial2_phi_ll.ij
	sym3 partial2_phi_ll = (sym3){.s={0,0,0,0,0,0}};
<? 
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	if i==j then
?>	partial2_phi_ll.<?=xij?> = (Up[<?=i-1?>]->phi - 2. * U->phi + Um[<?=i-1?>]->phi) / (grid_dx<?=i-1?> * grid_dx<?=i-1?>);
<?	else
?>	partial2_phi_ll.<?=xij?> = (
		U[index + stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].phi 
		- U[index - stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].phi 
		- U[index + stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].phi
		+ U[index - stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].phi 
	) / (4 * grid_dx<?=i-1?> * grid_dx<?=j-1?>);
<? 	end
end
?>

	real exp_4phi = exp(4. * U->phi);
	real exp_neg4phi = 1. / exp_4phi;

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma_ll = sym3_scale(U->gammaBar_ll, exp_4phi);

	//gammaBar^ij = inv gammaBar_ij
	sym3 gammaBar_uu = sym3_inv(1., U->gammaBar_ll);

	//gamma^ij = inv gamma_ij	
	//sym3 gamma_uu = sym3_inv(gamma_ll);
	//gamma^ij = exp(-4 phi) gammaBar^ij
	sym3 gamma_uu = sym3_scale(gammaBar_uu, exp_neg4phi);

	sym3 connBar_lll[3];	//connBar_ijk = connBar_lll[i].jk
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll[<?=i-1?>].<?=xjk?> = .5 * (partial_gammaBar_lll[<?=k-1?>].<?=sym(i,j)?> + partial_gammaBar_lll[<?=j-1?>].<?=sym(i,k)?> - partial_gammaBar_lll[<?=i-1?>].<?=xjk?>);
<?	end
end
?>	

	sym3 connBar_ull[3];	//connBar^i_jk = connBar_ull[i].jk
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
?>	connBar_ull[<?=i-1?>].<?=xjk?> = 0. <?
		for l,xl in ipairs(xNames) do
?> + gammaBar_uu.<?=sym(i,l)?> * connBar_lll[<?=l-1?>].<?=xjk?><?
		end
?>;
<?	end
end
?>

	sym3 conn_ull[3];	//conn^i_jk = conn_ull[i].jk
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	conn_ull[<?=i-1?>].<?=xjk?> = connBar_ull[<?=i-1?>].<?=xjk?><?
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

	//D2_alpha_ll.ij = D_i D_j alpha = partial_i partial_j alpha - conn^k_ij partial_k alpha
	sym3 D2_alpha_ll = (sym3){.s={0,0,0,0,0,0}};
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	D2_alpha_ll.<?=xij?> = partial2_alpha_ll.<?=xij?> <?
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
?>		+ partial_gammaBar_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>		//+ beta^k gammaBar_ij,k 
		+ U->gammaBar_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?> 	//+ gammaBar_ik beta^k_,j
		+ U->gammaBar_ll.<?=sym(j,k)?> * partial_beta_ul[<?=j-1?>].<?=xi?>	//+ gammaBar_kj beta^k_,i 
<? 	end
?>		- 2./3. * U->gammaBar_ll.<?=sym(i,j)?> * tr_partial_beta;		//- 2/3 gammaBar_ij beta^k_,k
<?
end
?>

	mat3 ATilde_lu = sym3_sym3_mul(U->ATilde_ll, gammaBar_uu);	//ATilde_i^j
	mat3 ATildeSq_lu = mat3_mat3_mul(ATilde_lu, ATilde_lu); 	//ATilde_i^k ATilde_k^j
	real tr_ATildeSq = mat3_trace(ATildeSq_lu);					//ATilde_ij ATilde^ij
	mat3 ATildeSq_ll = mat3_sym3_mul(ATilde_lu, U->ATilde_ll);	//ATilde_i^k ATilde_kj

	real S = sym3_dot(U->S_ll, gamma_uu);

	//B&S 11.52
	//K_,t = -gamma^ij D_ij alpha + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	deriv->K += -sym3_dot(gammaBar_uu, D2_alpha_ll) + U->alpha * (tr_ATildeSq + U->K * U->K / 3.) + 4 * M_PI * U->alpha * (U->rho + S) + real3_dot(U->beta_u, partial_K_l);

	sym3 partial2_gammaBar_llll[6];	//partial2_gammaBar_llll[kl].ij = gammaBar_ij,kl
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	if i==j then
?>	partial2_gammaBar_llll[<?=ij-1?>] = sym3_scale(
		sym3_add(
			sym3_scale(U->gammaBar_ll, -2.),
			sym3_add(
				Up[<?=i-1?>]->gammaBar_ll,
				Um[<?=i-1?>]->gammaBar_ll)),
			1. / (grid_dx<?=i-1?> * grid_dx<?=i-1?>));
<?	else
?>	partial2_gammaBar_llll[<?=ij-1?>] = sym3_scale(
		sym3_sub(
			sym3_add(
				U[index + stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].gammaBar_ll,
				U[index - stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].gammaBar_ll),
			sym3_add(
				U[index - stepsize[<?=i-1?>] + stepsize[<?=j-1?>]].gammaBar_ll,
				U[index + stepsize[<?=i-1?>] - stepsize[<?=j-1?>]].gammaBar_ll)
		), 1. / (4. * grid_dx<?=i-1?> * grid_dx<?=j-1?>));
<?	end
end
?>
	sym3 partial2_gammaBar_ll;	//partial2_gammaBar_ll.ij = gammaBar^kl gammaBar_ij,kl
<? for ij,xij in ipairs(symNames) do
?>	partial2_gammaBar_ll.<?=xij?> = 0. <?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?> + gammaBar_uu.<?=sym(k,l)?> * partial2_gammaBar_llll[<?=from3x3to6(k,l)-1?>].<?=xij?><?
		end
	end
?>;
<? end
?>

	//B&S 11.54
	sym3 RBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	RBar_ll.<?=xij?> = -1./2. * partial2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		+ .5 * U->gammaBar_ll.<?=sym(k,i)?> * partial_connBar_ul[<?=j-1?>].<?=xk?>
		+ .5 * U->gammaBar_ll.<?=sym(k,j)?> * partial_connBar_ul[<?=i-1?>].<?=xk?>
		+ .5 * U->connBar_u.<?=xk?> * connBar_lll[<?=i-1?>].<?=sym(j,k)?>
		+ .5 * U->connBar_u.<?=xk?> * connBar_lll[<?=j-1?>].<?=sym(i,k)?>
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>		+ connBar_ull[<?=k-1?>].<?=sym(l,i)?> * connBar_ull[<?=j-1?>].<?=sym(k,m)?>
		+ connBar_ull[<?=k-1?>].<?=sym(l,j)?> * connBar_ull[<?=i-1?>].<?=sym(k,m)?>
		+ connBar_ull[<?=k-1?>].<?=sym(i,m)?> * connBar_ull[<?=k-1?>].<?=sym(l,j)?>
<?			end
		end
	end
?>;
<? end
?>

	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll.<?=xij?> <?
	for k,xk in ipairs(xNames) do
?> - connBar_ull[<?=k?>].<?=xij?> * partial_phi_l.<?=xk?><?
	end
?>;
<? end
?>

	real tr_DBar2_phi = sym3_dot(gammaBar_uu, DBar2_phi_ll);

	real DBar_phi_norm = real3_weighted_norm(partial_phi_l, gammaBar_uu);

	//Baumgarte & Shapiro p.57 eqn 3.10
	//R_ll(i,j) := R_ij = RBar_ij - 2 (DBar_i DBar_j ln(psi) + gammaBar_ij gammaBar^lm DBar_l DBar_m ln(psi)) + 4((DBar_i ln(psi)) (DBar_j ln(psi)) - gammaBar_ij gammaBar^lm (DBar_l ln(psi)) (DBar_m ln(psi)))
	//Then Baumgarte & Shapiro on p.390 say RPhi_ij is the same as p.57 substituting phi for ln(psi)
	// ... but I thought phi was ln(psi)?  Then why would you need to separate R_ij = RBar_ij + RPhi_ij ?  I thought the substitution showed that R_ij was RPhi_ij?
	//phi = ln(psi), so DBar_i ln psi = partial_phi_i
	// so DBar_i DBar_j = 
	sym3 RPhi_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	RPhi_ll.<?=xij?> + RBar_ll.<?=xij?> - 2. * (DBar2_phi_ll.<?=xij?> + gammaBar_uu.<?=xij?> * tr_DBar2_phi) + 4. * (partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?> - U->gammaBar_ll.<?=xij?> * DBar_phi_norm);
<? end 
?>
	sym3 R_ll = sym3_add(RPhi_ll, RBar_ll);

	//traceless portion of -D^2 alpha + alpha (R_ij - 8 pi S_ij)
	sym3 tracelessPart_ll = sym3_sub(
		sym3_scale(
			sym3_add(R_ll, sym3_scale(U->S_ll, -8. * M_PI)), 
			U->alpha),
		D2_alpha_ll);
	real tr_tracelessPart = sym3_dot(gamma_uu, tracelessPart_ll);
	tracelessPart_ll = sym3_sub(tracelessPart_ll, sym3_scale(gamma_ll, -1./3. * tr_tracelessPart));

	//B&S 11.53
	//ATilde_ij,t = exp(-4phi) (-(D_ij alpha)^TF + alpha (R_ij^TF - 8 pi S_ij^TF) ) 
	//	+ alpha (K ATilde_ij - 2 ATilde_il ATilde^l_j)
	//	+ beta^k ATilde_ij,k + ATilde_ik beta^k_,j + ATidle_kj beta^k_,i - 2/3 ATilde_ij beta^k_,k
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->ATilde_ll.<?=xij?> = exp_neg4phi * tracelessPart_ll.<?=xij?>
		+ U->alpha * (U->K * U->ATilde_ll.<?=xij?> - 2. * ATildeSq_ll.<?=xi?>.<?=xj?>)
<?	for k,xk in ipairs(xNames) do
?>		+ partial_ATilde_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>
		+ U->ATilde_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
		+ U->ATilde_ll.<?=sym(j,k)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
<?	end
?>		- 2./3. * U->ATilde_ll.<?=xij?> * tr_partial_beta;
<? end
?>
}
