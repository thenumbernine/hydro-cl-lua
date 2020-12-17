//// MODULE_NAME: calc_gammaHat_ll
//// MODULE_DEPENDS: <?=coord_g_ll?>

#define calc_gammaHat_ll	coord_g_ll

//// MODULE_NAME: calc_gammaHat_uu
//// MODULE_DEPENDS: <?=coord_g_uu?>

#define calc_gammaHat_uu 	coord_g_uu

//// MODULE_NAME: calc_det_gammaHat
//// MODULE_DEPENDS: coord_det_g

#define calc_det_gammaHat 	coord_det_g

//// MODULE_NAME: calc_gammaHat_LL

#define calc_gammaHat_LL(x) (sym3_ident)

//// MODULE_NAME: calc_gammaHat_UU

#define calc_gammaHat_UU(x) (sym3_ident)

//// MODULE_NAME: calc_gammaBar_LL
//// MODULE_DEPENDS: <?=cons_t?> calc_gammaHat_LL

sym3 calc_gammaBar_LL(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaHat_LL = calc_gammaHat_LL(x);
	sym3 gammaBar_LL = sym3_add(gammaHat_LL, U->epsilon_LL);
	return gammaBar_LL;
}

//// MODULE_NAME: calc_det_gammaBarLL
/*
det(epsilon_IJ + gammaHat_IJ) 
= det(epsilon_IJ + delta_IJ) 
= det(e^i_I e^j_J (gammaHat_ij + epsilon_ij))
= det(e^i_I) det(e^j_J) det(gammaHat_ij + epsilon_ij)
= det(gammaHat^ij) det(gammaBar_ij)
= det(gammaBar_ij) / det(gammaHat_ij)
= 1
TODO detg ... unless we want to change the constraint
*/
#if 0	//use the value
real calc_det_gammaBarLL(global const <?=cons_t?>* U, ral3 x) {
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = sym3_det(gammaBar_LL);
	return det_gammaBarLL;
}
#else	//use the constraint
#define calc_det_gammaBarLL(x) 1.
#endif

//// MODULE_NAME: calc_gammaBar_UU
//// MODULE_DEPENDS: <?=cons_t?> calc_gammaBar_LL calc_det_gammaBarLL

sym3 calc_gammaBar_UU(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	return gammaBar_UU;
}

//// MODULE_NAME: calc_gammaBar_ll
//// MODULE_DEPENDS: <?=cons_t?> calc_gammaHat_ll rescaleFromCoord/rescaleToCoord

//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaHat_ll = calc_gammaHat_ll(x);
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);
	return gammaBar_ll;
}

//// MODULE_NAME: calc_det_gammaBar
//// MODULE_DEPENDS: calc_det_gammaHat

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//...except sometimes, according to 2012 Baumgarte et al, last paragraph of II B
real calc_det_gammaBar(real3 const x) {
	//TODO detg ...
	real det_gammaHat = calc_det_gammaHat(x);
	real detg = 1.;
	real det_gammaBar = det_gammaHat * detg;
	return det_gammaBar;
}

//// MODULE_NAME: calc_exp_neg4phi

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

//// MODULE_NAME: calc_gammaBar_uu
//// MODULE_DEPENDS: <?=cons_t?> calc_gammaBar_ll calc_det_gammaBar

sym3 calc_gammaBar_uu(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}

//// MODULE_NAME: <?=calc_gamma_ll?>
//// MODULE_DEPENDS: <?=cons_t?> calc_gammaBar_ll calc_exp_neg4phi

sym3 <?=calc_gamma_ll?>(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	return gamma_ll;
}

//// MODULE_NAME: <?=calc_gamma_uu?>
//// MODULE_DEPENDS: <?=cons_t?> calc_gammaBar_ll calc_exp_neg4phi calc_det_gammaBar

sym3 <?=calc_gamma_uu?>(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	real det_gamma = calc_det_gammaBar(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma); 
	return gamma_uu;
}

//// MODULE_NAME: mystery_C_U

//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero

//// MODULE_NAME: <?=setFlatSpace?>
//// MODULE_DEPENDS: mystery_C_U

void <?=setFlatSpace?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	U->alpha = 1.;
	U->beta_U = real3_zero;
	U->epsilon_LL = sym3_zero;
	U->W = 1;
	U->K = 0;
	U->ABar_LL = sym3_zero;

	//LambdaBar^i = Delta^i + C^i = Delta^i_jk gammaBar^jk = (connBar^i_jk - connHat^i_jk) gammaBar^jk + C^i
	//but when space is flat we have connBar^i_jk = connHat^i_jk and therefore Delta^i_jk = 0, Delta^i = 0, and LambdaBar^i = 0
	U->LambdaBar_U = mystery_C_U;

<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_U = real3_zero;
<? end
?>

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_U = real3_zero;
}
