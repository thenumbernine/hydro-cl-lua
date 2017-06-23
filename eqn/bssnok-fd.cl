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
	
	const real rho = 0;
	const real S = 0;
		
	const global <?=eqn.cons_t?>* Up[dim];
	const global <?=eqn.cons_t?>* Um[dim];
	for (int i = 0; i < dim; ++i) {
		Up[i] = U + index + stepsize[i];
		Um[i] = U + index - stepsize[i];
	}

	real3 partial_alpha;		//alpha_,i = partial_alpha[i]
	real3 partial_phi;			//phi_,i = partial_phi[i]
	real3 partial_tr_K;			//tr K,i = partial_tr_K[i]
	real3 partial_beta[3];		//beta^i_,j = partial_beta[j][i]
	real tr_partial_beta = 0.;	//beta^i_,i
	sym3 partial_gammaTilde[3];	//gammaTilde_ij,k = gammaTilde[k].ij
<? for i=0,solver.dim-1 do
?>	partial_alpha.s<?=i?> = (Up[<?=i?>]->alpha - Um[<?=i?>]->alpha) / (2. * grid_dx<?=i?>);
	partial_phi.s<?=i?> = (Up[<?=i?>]->phi - Um[<?=i?>]->phi) / (2. * grid_dx<?=i?>);
	partial_beta[<?=i?>] = real3_scale( real3_sub(Up[<?=i?>]->beta, Um[<?=i?>]->beta), 1. / (2. * grid_dx<?=i?>));
	partial_gammaTilde[<?=i?>] = sym3_scale( sym3_sub(Up[<?=i?>]->gammaTilde, Um[<?=i?>]->gammaTilde), 1. / (2. * grid_dx<?=i?>));
	partial_tr_K.s<?=i?> = (Up[<?=i?>]->tr_K - Um[<?=i?>]->tr_K) / (2. * grid_dx<?=i?>);
	tr_partial_beta += partial_beta[<?=i?>].s<?=i?>;
<? end
for i=solver.dim,2 do
?>	partial_alpha.s<?=i?> = 0;
	partial_phi.s<?=i?> = 0;
	partial_beta[<?=i?>] = _real3(0,0,0);
	partial_gammaTilde[<?=i?>] = (sym3){.s={0,0,0,0,0,0}};
	partial_tr_K.s<?=i?> = 0;
<? end
?>

	//alpha_,ij = partial2_alpha.ij
	sym3 partial2_alpha = (sym3){.s={0,0,0,0,0,0}};
<? for i=0,solver.dim-1 do
?>	partial2_alpha.<?=sym(i+1,i+1)?> = (Up[<?=i?>]->alpha - 2. * U->alpha + Um[<?=i?>]->alpha) / (grid_dx<?=i?> * grid_dx<?=i?>);
<? 	for j=i+1,solver.dim-1 do
?>	partial2_alpha.<?=sym(i+1,j+1)?> = (
		U[index + stepsize[<?=i?>] + stepsize[<?=j?>]]->alpha 
		- Um[index - stepsize[<?=i?>] + stepsize[<?=j?>]]->alpha 
		- Um[index + stepsize[<?=i?>] - stepsize[<?=j?>]]->alpha
		+ Up[index - stepsize[<?=i?>] - stepsize[<?=j?>]]->alpha 
	) / (4 * grid_dx<?=i?> * grid_dx<?=j?>);
<? 	end
end
?>


	real exp_4phi = exp(4. * U->phi);
	real exp_neg4phi = 1. / exp_4phi;

	//gamma_ij = exp(4 phi) gammaTilde_ij
	sym3 gamma = sym3_scale(U->gammaTilde, exp_4phi);

	//gammaTilde^ij = inv gammaTilde_ij
	//shouldn't this be 1 anyways?
	real det_gammaTilde = sym3_det(U->gammaTilde);
	sym3 gammaTildeInv = sym3_inv(det_gammaTilde, U->gammaTilde);

	//gamma^ij = inv gamma_ij	
	//sym3 gammaU = sym3_inv(gamma);
	//gamma^ij = exp(-4 phi) gammaTilde^ij
	sym3 gammaInv = sym3_scale(gammaTildeInv, exp_neg4phi);


	sym3 connTildeL[3];	//connTilde_ijk = connTildeL[i].jk
<? 
for i=0,2 do
	for j=0,2 do
		for k=j,2 do
?>	connTildeL[<?=i?>].<?=sym(j+1,k+1)?> = .5 * (
		partial_gammaTilde[<?=i?>].<?=sym(j+1,k+1)?>
		+ partial_gammaTilde[<?=i?>].<?=sym(k+1,j+1)?>
		- partial_gammaTilde[<?=j?>].<?=sym(k+1,i+1)?>);
<?		end
	end
end
?>	

	sym3 connTilde[3];	//connTilde^i_jk = connTilde[i].jk
<? 
for i=0,2 do
	for j=0,2 do
		for k=j,2 do
?>	connTilde[<?=i?>].<?=sym(j+1,k+1)?> = 0.
<?			for l=0,2 do
?>		+ gammaTildeInv.<?=sym(i+1,l+1)?> * connTildeL[<?=l?>].<?=sym(j+1,k+1)?>
<?			end
?>	;
<?		end
	end
end
?>

	sym3 conn[3];	//conn^i_jk = conn[i].jk
<? for i=0,2 do
	for j=0,2 do
		for k=0,2 do
?>	conn[<?=i?>].<?=sym(j+1,k+1)?> = connTilde[<?=i?>].<?=sym(j+1,k+1)?>
<?			if i==j then
?>		+ 2 * partial_phi.s<?=k?>
<?			end
			if i==k then
?>		+ 2 * partial_phi.s<?=j?>
<?			end
			for l=0,2 do
?>		- 2 * U->gammaTilde.<?=sym(j+1,k+1)?> * gammaTildeInv.<?=sym(i+1,l+1)?> * partial_phi.s<?=l?>
<?			end
?>	;
<?		end
	end
end
?>

	//D2_alpha.ij = D_i D_j alpha = partial_i partial_j alpha - conn^k_ij partial_k alpha
	sym3 D2_alpha = (sym3){.s={0,0,0,0,0,0}};
<? for i=0,solver.dim-1 do
	for j=i,solver.dim-1 do
?>	D2_alpha.<?=sym(i+1,j+1)?> = partial2_alpha.<?=sym(i+1,j+1)?>
<?		for k=0,2 do
?>		- conn[<?=k?>].<?=sym(i+1,j+1)?> * partial_alpha.s<?=k?>
<?		end ?>;
<?	end
end
?>

	//Q = f(alpha) tr K
	real Q = calc_f(U->alpha) * U->tr_K;
	
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q + real3_dot(partial_alpha, U->beta);

	//manuall update elsewhere?
	//deriv->beta += _real3(0,0,0);

	//B&S 11.50
	//phi,t = -1/6 alpha tr K + beta^i phi,i + 1/6 beta^i_,i
	deriv->phi += -U->alpha * U->tr_K / 6. + real3_dot(U->beta, partial_phi) + tr_partial_beta / 6.;

	//B&S 11.51
	//gammaTilde_ij,t = -2 alpha ATilde_ij + beta^k gammaTilde_ij,k + gammaTilde_ik beta^k_,j + gammaTilde_kj beta^k_,i - 2/3 gammaTilde_ij beta^k_,k
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->gammaTilde.<?=xij?> += -2 * U->alpha * U->ATilde.<?=xij?>	//-2 alpha ATilde_ij 
<? 	for k,xk in ipairs(xNames) do
?>		+ partial_gammaTilde[<?=k-1?>].<?=xij?> * U->beta.<?=xk?>		//+ beta^k gammaTilde_ij,k 
		+ U->gammaTilde.<?=sym(i,k)?> * partial_beta[<?=j-1?>].<?=xk?> 	//+ gammaTilde_ik beta^k_,j
		+ U->gammaTilde.<?=sym(j,k)?> * partial_beta[<?=j-1?>].<?=xi?>	//+ gammaTilde_kj beta^k_,i 
<? 	end
?>		- 2./3. * U->gammaTilde.<?=sym(i,j)?> * tr_partial_beta;		//- 2/3 gammaTilde_ij beta^k_,k
<?
end
?>

	mat3 ATilde_lu = sym3_sym3_mul(U->ATilde, gammaTildeInv);
	mat3 ATildeSq_lu = mat3_mat3_mul(ATilde_lu, ATilde_lu); 
	real ATildeSq = mat3_trace(ATildeSq_lu);

	//B&S 11.52
	//K_,t = -gamma^ij D_ij alpha + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	deriv->tr_K += -sym3_dot(gammaTildeInv, D2_alpha) + U->alpha * (ATildeSq + U->tr_K * U->tr_K / 3.) + 4 * M_PI * U->alpha * (rho + S) + real3_dot(U->beta, partial_tr_K);

}
