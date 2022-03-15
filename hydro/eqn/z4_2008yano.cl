//// MODULE_NAME: <?=calc_gamma_ll?>

#define <?=calc_gamma_ll?>(U, x)	((U)->gamma_ll)

//// MODULE_NAME: <?=calc_gamma_uu?>
//// MODULE_DEPENDS: <?=cons_t?>

static inline sym3 <?=calc_gamma_uu?>(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	return gamma_uu;
}

<?
if eqn.initCond.initAnalytical then
	error("TODO - can't handle analytical initial conditions yet")
end
?>

<? if eqn.initCond.useBSSNVars then ?>

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=coordMap?> <?=coord_gHol_ll?> <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real W = 1.;
	real K = 0.;
	real3 LambdaBar_U = real3_zero;
	real3 beta_U = real3_zero;
	real3 B_U = real3_zero;
	sym3 epsilon_LL = sym3_zero;
	sym3 ABar_LL = sym3_zero;

	real rho = 0.;

	<?=initCode()?>

	//for (int i = 0; i < numStates; ++i) {
	//	U->ptr[i] = 0.;
	//}
	*U = (<?=cons_t?>){.ptr={ 0. / 0. }};

	U->alpha = alpha;

//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?>
	// gammaHat_IJ = delta_IJ
	// gamma_ij = e_i^I e_j^J (epsilon_IJ + gammaHat_IJ) / W^2
	sym3 const gammaBar_LL = sym3_add(epsilon_LL, sym3_ident);
	sym3 const gamma_LL = sym3_real_mul(gammaBar_LL, 1. / (W*W));
	U->gamma_ll = sym3_rescaleToCoord_LL(gamma_LL, x);
	
	// K_ij = e_i^I e_j^J (ABar_IJ + gammaBar_IJ K/3) / W^2
	U->K_ll = sym3_rescaleToCoord_LL(
		sym3_add(
			sym3_real_mul(ABar_LL, 1. / (W*W)),
			sym3_real_mul(gamma_LL, K / 3.)
		), x);

	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if eqn.useShift ~= "none" then
?>	U->beta_u = real3_rescaleFromCoord_U(beta_U);
<? end -- TODO support for hyperbolic gamma driver, so we can read B_U
?>

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	
	U->H = 0;
	U->M_u = real3_zero;
}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf 
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (
		log(U[solver->stepsize.<?=xi?>].alpha) 
		- log(U[-solver->stepsize.<?=xi?>].alpha)
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (
		U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> 
		- U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>
}

<? else	-- not eqn.initCond.useBSSNVars ?>

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=coordMap?> <?=coord_gHol_ll?> <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real3 beta_u = real3_zero;

	sym3 gamma_ll = coord_gHol_ll(x);

	sym3 K_ll = sym3_zero;

	//TODO more stress-energy vars 
	real rho = 0.;

	<?=initCode()?>

	*U = (<?=cons_t?>){.ptr={ 0. / 0. }};
	
	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;

	//Z_u n^u = 0
	//Theta = alpha n_u Z^u = alpha Z^u
	//for n_a = (-alpha, 0)
	//n^a_l = (1/alpha, -beta^i/alpha)
	//(Z_t - Z_i beta^i) / alpha = Theta ... = ?
	//Z^t n_t + Z^i n_i = -alpha Z^t = Theta
	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if eqn.useShift ~= "none" then
?>	U->beta_u = beta_u;
<? end
?>

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	
	U->H = 0;
	U->M_u = real3_zero;
}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	
	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (
		log(U[solver->stepsize.<?=xi?>].alpha) 
		- log(U[-solver->stepsize.<?=xi?>].alpha)
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (
		U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> 
		- U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? end ?>
<?
end
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>
}

<? end	-- eqn.initCond.useBSSNVars ?>

//// MODULE_NAME: <?=setFlatSpace?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?>

static inline void <?=setFlatSpace?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	(U)->alpha = 1.;
	(U)->gamma_ll = sym3_ident;
	(U)->a_l = real3_zero;
	(U)->d_lll.x = sym3_zero;
	(U)->d_lll.y = sym3_zero;
	(U)->d_lll.z = sym3_zero;
	(U)->K_ll = sym3_zero;
	(U)->Theta = 0.;
	(U)->Z_l = real3_zero;
<? if eqn.useShift ~= "none" then 
?>	(U)->beta_u = real3_zero;
<? end 
?>	
<? if eqn.useStressEnergyTerms then ?>
	//what to do with the constraint vars and the source vars?
	(U)->rho = 0;
	(U)->S_u = real3_zero;
	(U)->S_ll = sym3_zero;
<? end ?>
	
	(U)->H = 0;
	(U)->M_u = real3_zero;
}

//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=SETBOUNDS?> <?=cons_t?> <?=initCond_codeprefix?>

#define <?=calcDTCell?>(\
	/*global real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell\
) {\
	/* the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma_ij) is only called once */\
	real const f_alphaSq = calc_f_alphaSq(U->alpha);\
	real const det_gamma = sym3_det(U->gamma_ll);\
	real const alpha_sqrt_f = sqrt(f_alphaSq);\
\
	<? for side=0,solver.dim-1 do ?>{\
\
		<? if side == 0 then ?>\
		real const gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;\
		<? elseif side == 1 then ?>\
		real const gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;\
		<? elseif side == 2 then ?>\
		real const gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;\
		<? end ?>	\
		real const sqrt_gammaUjj = sqrt(gammaUjj);\
		real const lambdaLight = sqrt_gammaUjj * U->alpha;\
		real const lambdaGauge = sqrt_gammaUjj * alpha_sqrt_f;\
		real const lambda = (real)max(lambdaGauge, lambdaLight);\
\
		<? if eqn.useShift ~= "none" then ?>\
		real const betaUi = U->beta_u.s<?=side?>;\
		<? else ?>\
		real const betaUi = 0.;\
		<? end ?>\
\
		real const lambdaMin = (real)min((real)0., -betaUi - lambda);\
		real const lambdaMax = (real)max((real)0., -betaUi + lambda);\
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
		absLambdaMax = max((real)1e-9, absLambdaMax);\
		*(dt) = (real)min(*(dt), solver->grid_dx.s<?=side?> / absLambdaMax);\
	}<? end ?>\
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=cons_t?> <?=solver_t?> <?=normal_t?> rotate sym3_rotate _3sym3_rotate

/*
Going by the 2008 Yano et al paper flux term.
Actually I'm going by my attempts to reproduce it, in my symmath/tests/output/Z4.html file, since I found math errors in the 2008 Yano paper.
But my own worksheet jumps right into the flux jacobian eigensystem, doesn't stop at isolating the flux, so I'm sort of winging this part ...
TODO in that worksheet separate the flux out and codegen for it, and put that code HERE.

Notice in my Z4.html worksheet, ^γ_ij is the background metric, which is static, 
and Δγ_ij is the deviation from the background metric, which is evolved,
so γ_ij = Δγ_ij + ^γ_ij
same with d_ijk

a_k,t 
	+ ((-2 Θ α f + α f K) δ^r_k)_,r -- flux without shift
	+ (-β^r a_k)_,r 				-- flux with shift
	= 
	+ b^l_k a_l - b^l_l a_k			-- source with shift
... compared to 2008 yano seems I lost an a_m β^m term in the flux 
... unless maybe 2008 yano didn't want to start with the assumption that α_,t = α_,i β^i - α^2 f (K - 2 Θ)
... maybe instead 2008 yano drops the shift term and starts with α_,t = α^2 f (K - 2 Θ), which is what it says in eqn 14
... also in its eqn 15 it seems to have lumped the α_,i β^i into the Q term of the lapse gauge.
... so I'll go ahead and consider this to be done ...
*** CHECK ***

Δd_kij,t - ((+ ^d_tij - 1/2 ^γ_il b^l_j - 1/2 ^γ_jl b^l_i - 1/2 Δγ_il b^l_j - 1/2 Δγ_jl b^l_i + α K_ij - β^l ^d_lij) δ^r_k - β^r Δd_kij)_,r = b^l_k Δd_lij - b^l_l Δd_kij
... using simplified terms ...
Δd_kij,t 
	+ (α K_ij δ^r_k)_,r 								-- flux without shift
	+ (-(b_(ij) + β^l ^d_lij) δ^r_k - β^r Δd_kij)_,r 	-- flux with shift
	= 
	- ^d_kij,t											-- source
	+ b^l_k Δd_lij - b^l_l Δd_kij						-- source with shift
... this could be true or not.  the 2008 Yano paper uses a mystery Q_ij gauage var.
... it cites its def in 2005 Bona et al "Geometrically..." which uses the same Q_ij var with no def
... so I don't have any way to know wtf this var is.
... but other Z4 defs without shift give d_kij,t + α K_ij,k = ... as a def, so I'm close 
*** CLOSE ENOUGH ***

K_ij,t's definition using d_ijk's:
K_ij,t = -1/2 (α a_i δ^k_j + α a_j δ^k_i)_,k
	+ α a_k γ^km (- d_mij + d_jmi + d_imj)
	+ α (
		( -1/2 γ^pq (d_ipq δ^k_j + d_jpq δ^k_i) + γ^kl (1/2 d_ilj + 1/2 d_jli - d_lij) )_,k
		+ γ^mn (d_inj + d_jni - d_nij) γ^kl d_mlk
		- γ^kl (d_mlj + d_jlm - d_lmj) γ^mn (d_ink + d_kni - d_nik)
		- γ^ak γ^bl d_kab (d_ijl + d_jil)
		+ 1/2 γ^mp (d_mpj δ^k_i + d_mpi δ^k_j)_,k
	)
	+ α (
		-2 Z_k γ^km (d_imj + d_jmi - d_mij)
		+ K_ij (K - 2 Θ)
		- 2 K_ik γ^kl K_lj
	)
	+ 4 π α (γ_ij (S - ρ) - 2 S_ij)
	+ K_ij,k β^k
	+ K_ki b^k_j
	+ K_kj b^k_i
	+ (α (Z_i δ^k_j + Z_j δ^k_i))_,k
	- α a_k (Z_i δ^k_j + Z_j δ^k_i)
	+ (β^k K_ij)_,k
	- β^k K_ij,k
	- K_ij b^k_k

... separating out terms ...
K_ij,t 
	+ (
		α (												-- flux without shift ... CHECK
			- 1/2 γ^kl (d_ilj + d_jli - 2 d_lij)
			- 1/2 γ^mp (d_mpj δ^k_i + d_mpi δ^k_j)
			+ 1/2 γ^pq (d_ipq δ^k_j + d_jpq δ^k_i)
			+ 1/2 (a_i δ^k_j + a_j δ^k_i)
			- Z_i δ^k_j - Z_j δ^k_i
		)
		- β^k K_ij										-- flux with shift ... CHECK
	)_,k
	= 
	+ α (
		+ a_k (											-- source without shift ... bleh ... have to check in 2005 Bona et al ... 
			+ 1/2 γ^mp (d_ipj δ^k_m + d_jpi δ^k_m)
			- 1/2 γ^mp (d_mpj δ^k_i + d_mpi δ^k_j)
			+ 1/2 γ^pq (d_ipq δ^k_j + d_jpq δ^k_i)
		)
		- γ^ak γ^bl d_kab (d_ijl + d_jil)
		+ 1/2 γ^ma γ_ab,k γ^pb (d_mpj δ^k_i + d_mpi δ^k_j)
		+ γ^mn (d_inj + d_jni - d_nij) γ^kl d_mlk
		- 2 Z_k γ^km (d_imj + d_jmi - d_mij)
		- γ^kl (d_mlj + d_jlm - d_lmj) γ^mn (d_ink + d_kni - d_nik)
		- a_j Z_i - a_i Z_j
		+ K_ij (K - 2 Θ)
		- 2 K_ik γ^kl K_lj
	)
	+ K_ki b^k_j										-- source with shift ... CHECK
	+ K_kj b^k_i
	- K_ij b^k_k
	+ 4 π α (γ_ij (S - ρ) - 2 S_ij)						-- source with stress-energy ... CHECK
... ok from 2008 Yano, flux matches.
... but source term S(K_ij) aren't specified in the 2008 Yano paper
... so how about the source ...

Θ_,t = (Θ β^k)_,k
	+ (α Z_l γ^kl)_,k
	+ 1/2 (
		+ (2 α γ^ij γ^kl (d_ijl - d_lij))_,k
		- α_,k γ^ij γ^kl (d_ijl - d_lij)
		+ α (
			+ 2 γ^ci γ^dj d_kcd (-1/2 γ^pq (d_ipq δ^k_j + d_jpq δ^k_i) + 1/2 γkl (d_ilj + d_jli - 2 d_lij))
			+ γ^ij γ^mn (d_inj + d_jni - d_nij) γ^kl d_mlk
			- γ^ij γ^kl (d_mlj + d_jlm - d_lmj) γ^mn (d_ink + d_kni - d_nik)
			- γ^ij γ^ak γ^bl d_kab (d_ijl + d_jil)
			+ (γ^mp γ^ci γ^dj d_kcd + γ^ij γ^cm γ^dp d_kcd) (d_mpj δ^k_i + d_mpi δ^k_j)
		)
	)
	+ 2 α γ^ak γ^bl d_kab Z_l 
	- 2 α γ^kl Z_k a_l
	- α γ^kl Γ^m_kl Z_m
	+ 1/2 α (
		K (K - 2 Θ)
		- K^kl K_kl
		- 16 π ρ
	)
	- Θ b^k_k

Θ_,t 
	+ (
		- α Z_l γ^kl								-- flux without shift ... CHECK
		- α γ^ij γ^kl (d_ijl - d_lij)
		- Θ β^k										-- flux with shift ... CHECK
	)_,k
	=
	+ 1/2 (											-- source without shift
		- α a_k γ^ij γ^kl (d_ijl - d_lij)
		+ α (
			+ 2 γ^ci γ^dj d_kcd (-1/2 γ^pq (d_ipq δ^k_j + d_jpq δ^k_i) + 1/2 γkl (d_ilj + d_jli - 2 d_lij))
			+ γ^ij γ^mn (d_inj + d_jni - d_nij) γ^kl d_mlk
			- γ^ij γ^kl (d_mlj + d_jlm - d_lmj) γ^mn (d_ink + d_kni - d_nik)
			- γ^ij γ^ak γ^bl d_kab (d_ijl + d_jil)
			+ (γ^mp γ^ci γ^dj d_kcd + γ^ij γ^cm γ^dp d_kcd) (d_mpj δ^k_i + d_mpi δ^k_j)
		)
	)
	+ 2 α γ^ak γ^bl d_kab Z_l 
	- 2 α γ^kl Z_k a_l
	- α γ^kl Γ^m_kl Z_m
	+ 1/2 α (
		K (K - 2 Θ)
		- K^kl K_kl
	)
	- Θ b^k_k										-- source with shift ... CHECK
	+ 8 α π ρ										-- source with stress-energy ... CHECK

Z_k,t 
	+ (
		+ α (										-- flux without shift
			+ γ^lr K_kl
			+ (K - Θ) δ^r_k
		)
		- β^r Z_k									-- flux with shift
	)_,r
	= 
	+ α γ^mn_,k K_mn
	+ α a_k γ^mn K_mn
	+ α a_m γ^lm K_kl
	+ α γ^lm_,m K_kl
	+ 2 α K_mn γ^lm γ^pn d_klp
	- α K_kn Γ^n_lm γ^lm
	- α K_nl Γ^n_km γ^lm
	- 2 α Z_m γ^lm K_kl
	- 2 Θ α a_k 
	+ b^l_k Z_l										-- source with shift ... CHECK
	- b^l_l Z_k
	- 8 α π S_k										-- source with stress-energy ... CHECK

MDE:
β^l_,t = 
	+ α^2 γ^kl (2 γ^jm d_mjk - γ^jm d_kjm - a_k)	-- source without shift
	+ β^k b^l_k 									-- source with shift
*** CHECK *** i think

MDE:
b^l_k,t 
	+ (-α^2 γ^il (2 γ^jm d_mji - γ^jm d_ijm - a_i))_,k	-- flux without shift
	+ (-β^i b^l_i)_,k									-- flux with shift
	= 0
*/
#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const f_alpha = calc_f_alpha((U)->alpha);\
\
	real const det_gamma = sym3_det((U)->gamma_ll);\
	sym3 const gamma_uu = sym3_inv((U)->gamma_ll, det_gamma);\
	real const K = sym3_dot((U)->K_ll, gamma_uu);\
	real3 const d_l = _3sym3_sym3_dot23((U)->d_lll, gamma_uu);\
	real3 const e_l = sym3_3sym3_dot12(gamma_uu, (U)->d_lll);\
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, (U)->d_lll);\
	real3x3x3 const d_llu = _3sym3_sym3_mul((U)->d_lll, gamma_uu);\
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);\
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);\
	real3 const Z_u = sym3_real3_mul(gamma_uu, (U)->Z_l);\
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, (U)->K_ll);\
\
	(resultFlux)->H = 0.;\
	(resultFlux)->M_u = real3_zero;\
	(resultFlux)->alpha = 0.;\
	(resultFlux)->gamma_ll = sym3_zero;\
\
/* a_k,t		+ (α f (K - 2 Θ) δ^r_k)_,r */\
	(resultFlux)->a_l.x = f_alpha * (K - 2 * (U)->Theta);\
\
/* Δd_kij,t		+ (α K_ij δ^r_k)_,r */\
	(resultFlux)->d_lll.x = sym3_real_mul((U)->K_ll, (U)->alpha);\
\
/* K_ij,t 		+ (1/2 α (									\
					+ δ^r_i (d_j - e_j + a_j - 2 Z_j)		\
					+ δ^r_j (d_i - e_i + a_i - 2 Z_i)		\
					+ (2 d^r_ij - d_ij^r - d_ji^r)			\
				))_,r 										*/	\
<? --\
for ij,xij in ipairs(symNames) do --\
	local i,j,xi,xj = from6to3x3(ij) --\
?>	(resultFlux)->K_ll.<?=xij?> = .5 * (U)->alpha * (\
<?	if i == 1 then --\
?>		+ d_l.<?=xj?> - e_l.<?=xj?> + (U)->a_l.<?=xj?> - 2. * (U)->Z_l.<?=xj?>\
<? 	end --\
	if j == 1 then --\
?>		+ d_l.<?=xi?> - e_l.<?=xi?> + (U)->a_l.<?=xi?> - 2. * (U)->Z_l.<?=xi?>\
<? 	end --\
?>		+ 2. * d_ull.x.<?=xij?> - d_llu.<?=xi?>.<?=xj?>.x - d_llu.<?=xj?>.<?=xi?>.x\
	);\
<? --\
end ?>\
\
/* Θ_,t			+ (α (d^r - e^r - Z^r))_,r */	\
	(resultFlux)->Theta = (U)->alpha * (d_u.x - e_u.x - Z_u.x);\
\
/* Z_k,t 		+ (α (K_k^r + (K - Θ) δ^r_k))_,r */	\
	(resultFlux)->Z_l = real3_real_mul(K_ul.x, (U)->alpha);\
	(resultFlux)->Z_l.x += (U)->alpha * (K - (U)->Theta);\
\
<? if eqn.useShift ~= "none" then ?>/*									\
flux without shift:												\
b^l_k,t 	+ (-α^2 γ^il (2 γ^jm d_mji - γ^jm d_ijm - a_i))_,k	\
																\
flux with shift:												\
a_k,t		+ (-β^r a_k)_,r 									\
Δd_kij,t 	+ (-(b_(ij) + β^l ^d_lij) δ^r_k - β^r Δd_kij)_,r 	\
K_ij,t		+ (- β^k K_ij)_,k									\
Θ_,t 		+ (- Θ β^k)_,k										\
Z_k,t 		+ (													\
				- β^r Z_k										\
			)_,r												\
b^l_k,t 	+ (-β^i b^l_i)_,k									*/	\
	(resultFlux)->beta_l = real3_zero;\
	(resultFlux)->b_ul = real3x3_zero; /* TODO FINISHME */\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=solver_t?> <?=eigen_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>
// used by hll, roe, weno, plm ... anything that uses eigenvalues or eigenvector transforms

//used for interface eigen basis
#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	(resultEig)->alpha = .5 * ((UL)->alpha + (UR)->alpha);\
	(resultEig)->alpha_sqrt_f = sqrt(calc_f_alphaSq((resultEig)->alpha));\
\
	sym3 const avg_gamma = sym3_real_mul(sym3_add((UL)->gamma_ll, (UR)->gamma_ll), .5);\
	real const det_avg_gamma = sym3_det(avg_gamma);\
	(resultEig)->gamma_ll = avg_gamma;\
	(resultEig)->gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
/*  I'm using .side for holonomic(coordinate) and anholonomic(orthonormal) */\
/* but for cartesian vector componets there is no .side, just .n, which is covariant iirc */\
/* and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes) */\
/* so here I'm going to just wing it */\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, (resultEig)->gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = (resultEig)->gamma_uu.xx;\
	} else if (n.side == 1) {\
		gammaUnn = (resultEig)->gamma_uu.yy;\
	} else if (n.side == 2) {\
		gammaUnn = (resultEig)->gamma_uu.zz;\
	}\
<? end ?>\
\
	(resultEig)->sqrt_gammaUnn = sqrt(gammaUnn);\
\
<? if eqn.useShift ~= "none" then ?>\
	(resultEig)->beta_u = real3_real_mul(real3_add((UL)->beta_u, (UR)->beta_u), .5);\
<? end ?>\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: <?=range_t?> <?=normal_t?> cons_t <?=initCond_codeprefix?>
// not used anymore, replaced in calcDT by eqn:consMinWaveCode/eqn:consMaxWaveCode eigenvalue inlining

#define <?=calcCellMinMaxEigenvalues?>(\
	/*<?=range_t?> * const */result,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const det_gamma = sym3_det((U)->gamma_ll);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
	sym3 const gamma_uu = sym3_inv((U)->gamma_ll, det_gamma);\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = ((U)->gamma_ll.yy * (U)->gamma_ll.zz - (U)->gamma_ll.yz * (U)->gamma_ll.yz) / det_gamma;\
	} else if (n.side == 1) {\
		gammaUnn = ((U)->gamma_ll.xx * (U)->gamma_ll.zz - (U)->gamma_ll.xz * (U)->gamma_ll.xz) / det_gamma;\
	} else if (n.side == 2) {\
		gammaUnn = ((U)->gamma_ll.xx * (U)->gamma_ll.yy - (U)->gamma_ll.xy * (U)->gamma_ll.xy) / det_gamma;\
	}\
<? end ?>\
	real const sqrt_gammaUnn = sqrt(gammaUnn);\
	real const lambdaLight = (U)->alpha * sqrt_gammaUnn;\
\
	real const f_alphaSq = calc_f_alphaSq((U)->alpha);\
	real const lambdaGauge = sqrt(f_alphaSq) * sqrt_gammaUnn;\
\
	real lambdaMax = max(lambdaGauge, lambdaLight);\
	real lambdaMin = -lambdaMin;\
\
	<? if eqn.useShift ~= "none" then ?>\
	lambdaMin -= normal_vecDotN1(n, (U)->beta_u);\
	lambdaMax -= normal_vecDotN1(n, (U)->beta_u);\
	<? end ?>\
\
	(result)->min = lambdaMin;\
	(result)->max = lambdaMax;\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>

//used by PLM, and by the default <?=fluxFromCons?> (used by hll, or roe when roeUseFluxFromCons is set)
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	(resultEig)->alpha = (U)->alpha;\
	(resultEig)->alpha_sqrt_f = sqrt(calc_f_alphaSq((U)->alpha));\
	(resultEig)->gamma_ll = (U)->gamma_ll;\
	real const det_gamma = sym3_det((U)->gamma_ll);\
	(resultEig)->gamma_uu = sym3_inv((U)->gamma_ll, det_gamma);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
/*  I'm using .side for holonomic(coordinate) and anholonomic(orthonormal) */\
/* but for cartesian vector componets there is no .side, just .n, which is covariant iirc */\
/* and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes) */\
/* so here I'm going to just wing it */\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, (resultEig)->gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = (resultEig)->gamma_uu.xx;\
	} else if (n.side == 1) {\
		gammaUnn = (resultEig)->gamma_uu.yy;\
	} else if (n.side == 2) {\
		gammaUnn = (resultEig)->gamma_uu.zz;\
	}\
<? end ?>\
\
	(resultEig)->sqrt_gammaUnn = sqrt(gammaUnn);\
\
	<? if eqn.useShift ~= "none" then ?>\
	(resultEig)->beta_u = (U)->beta_u;\
	<? end ?>\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> rotate <?=initCond_codeprefix?>
// used by roe, weno, some plm

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */results,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */inputU,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	/* input */\
	real3 const a_l = real3_swap((inputU)->a_l, n.side);							/* 0-2 */\
	_3sym3 const d_lll = _3sym3_swap((inputU)->d_lll, n.side);						/* 3-20 */\
	sym3 const K_ll = sym3_swap((inputU)->K_ll, n.side);							/* 21-26 */\
	real const Theta = (inputU)->Theta;												/* 27 */\
	real3 const Z_l = real3_swap((inputU)->Z_l, n.side);							/* 28-30 */\
\
	/* eig */\
	real const sqrt_f = (eig)->alpha_sqrt_f / (eig)->alpha;\
	real const f = sqrt_f * sqrt_f;	\
	real const _1_sqrt_f = 1. / sqrt_f;\
	real const _1_f = _1_sqrt_f * _1_sqrt_f;\
	real const fSq = f * f;\
	real const f_toThe_3_2 = f * sqrt_f;\
	real const f_toThe_5_2 = f * f * sqrt_f;\
\
	real const f_minus_1 = f - 1.;\
	real const f_minus_1_sq = f_minus_1 * f_minus_1;\
\
	sym3 const gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 const gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
\
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);\
\
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, K_ll);\
\
	real const tr_K = real3x3_trace(K_ul);\
\
	real sqrt_gUxx = (eig)->sqrt_gammaUnn;\
\
	real const _1_sqrt_gUxx = 1. / sqrt_gUxx;\
	real const _1_gammaUxx = _1_sqrt_gUxx * _1_sqrt_gUxx;\
\
	real const a_u_x = a_l.x * gamma_uu.xx + a_l.y * gamma_uu.xy + a_l.z * gamma_uu.xz;\
	real const Z_u_x = Z_l.x * gamma_uu.xx + Z_l.y * gamma_uu.xy + Z_l.z * gamma_uu.xz;\
\
	/* d_i = d_ijk gamma^jk */\
	real3 const d_l = _real3(\
		sym3_dot(d_lll.x, gamma_uu),\
		sym3_dot(d_lll.y, gamma_uu),\
		sym3_dot(d_lll.z, gamma_uu));\
\
	real const d_u_x = d_l.x * gamma_uu.xx + d_l.y * gamma_uu.xy + d_l.z * gamma_uu.xz;\
\
	/* e_i = d_jki gamma^jk */\
	real3 const e_l = (real3){\
<? for i,xi in ipairs(xNames) do --\
?>		.<?=xi?> = 0. <? --\
	for j,xj in ipairs(xNames) do --\
		for k,xk in ipairs(xNames) do --\
?> + d_lll.<?=xj?>.<?=sym(k,i)?> * gamma_uu.<?=sym(j,k)?><? --\
		end --\
	end --\
?>,\
<? end --\
?>	};\
\
	real const e_u_x = e_l.x * gamma_uu.xx + e_l.y * gamma_uu.xy + e_l.z * gamma_uu.xz;\
\
	real const m = solver->m;\
\
	(results)->ptr[0] = (\
		-(2. * gamma_uu.xy * a_l.y * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * a_l.y * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xy * a_l.y * f * _1_sqrt_gUxx \
		+ 8. * gamma_uu.xy * K_ll.xy * f_toThe_3_2 \
		- 4. * gamma_uu.xy * K_ll.xy * sqrt_f \
		- 4. * gamma_uu.xy * K_ll.xy * f_toThe_5_2 \
		+ 8. * gamma_uu.xz * K_ll.xz * f_toThe_3_2 \
		- 4. * gamma_uu.xz * K_ll.xz * f_toThe_5_2 \
		- 4. * gamma_uu.xz * K_ll.xz * sqrt_f \
		+ 2. * gamma_uu.xz * a_l.z * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xz * a_l.z * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xz * a_l.z * f * _1_sqrt_gUxx \
		+ 4. * gamma_uu.yy * K_ll.yy * f_toThe_3_2 \
		- 2. * gamma_uu.yy * K_ll.yy * sqrt_f \
		- 2. * gamma_uu.yy * K_ll.yy * f_toThe_5_2 \
		- 4. * gamma_uu.yz * K_ll.yz * f_toThe_5_2 \
		- 4. * gamma_uu.yz * K_ll.yz * sqrt_f \
		+ 8. * gamma_uu.yz * K_ll.yz * f_toThe_3_2 \
		+ 4. * gamma_uu.zz * K_ll.zz * f_toThe_3_2 \
		- 2. * gamma_uu.zz * K_ll.zz * sqrt_f \
		- 2. * gamma_uu.zz * K_ll.zz * f_toThe_5_2 \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		- 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 8. * f * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 8. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		+ 8. * f * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 8. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * Z_l.x * sqrt_gUxx \
		- 4. * fSq * Z_l.x * sqrt_gUxx \
		+ 2. * f_toThe_5_2 * m * Theta \
		- 2. * f_toThe_3_2 * m * Theta \
		- 2. * f * m * Z_l.x * sqrt_gUxx \
		+ 2. * fSq * m * Z_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx * fSq \
		- 4. * a_l.x * sqrt_gUxx * f \
		+ 4. * K_ll.xx * gamma_uu.xx * f_toThe_3_2 \
		- 2. * K_ll.xx * gamma_uu.xx * sqrt_f \
		- 2. * K_ll.xx * f_toThe_5_2 * gamma_uu.xx \
		+ 4. * Theta * sqrt_f \
		- 4. * Theta * f_toThe_3_2)) / (gamma_uu.xx * sqrt_f * (4. \
		- 8. * f \
		+ 4. * fSq));\
	(results)->ptr[1] = (2. * gamma_uu.xx * sqrt_gUxx * K_ll.xz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.zz \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.z.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yz * sqrt_gUxx \
		+ 2. * gamma_uu.xz * K_ll.zz * sqrt_gUxx \
		+ gamma_uu.yy * gamma_uu.xx * d_lll.y.yz \
		- gamma_uu.yy * gamma_uu.xx * d_lll.z.yy \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.y.zz \
		- gamma_uu.yz * gamma_uu.xx * d_lll.z.yz \
		- a_l.z * gamma_uu.xx \
		+ 2. * Z_l.z * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[2] = (2. * gamma_uu.xx * sqrt_gUxx * K_ll.xy \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yy * sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.yz \
		+ 2. * gamma_uu.xz * K_ll.yz * sqrt_gUxx \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.zz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.zz * gamma_uu.xx * d_lll.z.yz \
		- a_l.y * gamma_uu.xx \
		+ 2. * Z_l.y * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[3] = (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		+ gamma_uu.xy * a_l.y \
		- 2. * gamma_uu.xy * K_ll.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		- gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		- 2. * gamma_uu.xz * K_ll.xz * sqrt_gUxx \
		+ gamma_uu.xz * a_l.z \
		- 2. * gamma_uu.yy * K_ll.yy * sqrt_gUxx \
		- 4. * gamma_uu.yz * K_ll.yz * sqrt_gUxx \
		- 2. * gamma_uu.zz * K_ll.zz * sqrt_gUxx \
		+ 2. * Theta * sqrt_gUxx \
		+ 2. * Z_l.x * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[4] = (gamma_uu.xx * gamma_uu.yy * d_lll.y.xy \
		- gamma_uu.xx * gamma_uu.yy * d_lll.x.yy \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.y.xz \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.z.xy \
		- 2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz \
		+ gamma_uu.xx * gamma_uu.zz * d_lll.z.xz \
		- gamma_uu.xx * gamma_uu.zz * d_lll.x.zz \
		+ gamma_uu.xx * Z_l.x \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.xz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.xy \
		+ 2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		- gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.y.xy \
		+ gamma_uu.xy * Z_l.y \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.x.yy \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		- gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.xz \
		+ gamma_uu.xz * Z_l.z \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.x.zz \
		+ Theta * sqrt_gUxx) / (2. * sqrt_gUxx);\
	(results)->ptr[5] = (gamma_uu.xx * d_lll.z.xz \
		- gamma_uu.xx * d_lll.x.zz \
		- gamma_uu.xy * d_lll.y.zz \
		+ gamma_uu.xy * d_lll.z.yz \
		+ K_ll.zz * sqrt_gUxx) / (2. * sqrt_gUxx);\
	(results)->ptr[6] = (gamma_uu.xx * d_lll.y.xz \
		+ gamma_uu.xx * d_lll.z.xy \
		- 2. * gamma_uu.xx * d_lll.x.yz \
		- gamma_uu.xy * d_lll.y.yz \
		+ gamma_uu.xy * d_lll.z.yy \
		+ gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * d_lll.z.yz \
		+ 2. * K_ll.yz * sqrt_gUxx) / (4. * sqrt_gUxx);\
	(results)->ptr[7] = (\
		-(gamma_uu.xx * gamma_uu.yy * d_lll.y.xy * f \
		- gamma_uu.xx * gamma_uu.yy * d_lll.x.yy * f \
		- gamma_uu.xx * gamma_uu.yy * m * d_lll.y.xy * f \
		+ gamma_uu.xx * gamma_uu.yy * m * d_lll.x.yy * f \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.y.xz * f \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.z.xy * f \
		- 2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz * f \
		- gamma_uu.xx * gamma_uu.yz * m * d_lll.y.xz * f \
		- gamma_uu.xx * gamma_uu.yz * m * d_lll.z.xy * f \
		+ 2. * gamma_uu.xx * gamma_uu.yz * m * d_lll.x.yz * f \
		+ gamma_uu.xx * gamma_uu.zz * d_lll.z.xz * f \
		- gamma_uu.xx * gamma_uu.zz * d_lll.x.zz * f \
		- gamma_uu.xx * gamma_uu.zz * m * d_lll.z.xz * f \
		+ gamma_uu.xx * gamma_uu.zz * m * d_lll.x.zz * f \
		- 2. * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * f \
		- 2. * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * f \
		+ 4. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * f \
		+ gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * f \
		+ gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * f \
		- 2. * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * f \
		- 2. * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * f \
		+ 2. * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * f \
		+ gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * f \
		- gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * f \
		- 2. * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * f \
		+ 2. * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * f \
		+ gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * f \
		- gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * f \
		- 2. * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * f \
		+ gamma_uu.xy * a_l.y \
		+ 2. * gamma_uu.xy * Z_l.y * f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * f \
		+ gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * f \
		- gamma_uu.xy * m * Z_l.y * f \
		- gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * f \
		+ 2. * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * f \
		- 2. * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * f \
		- gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * f \
		+ gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * f \
		+ 2. * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * f \
		- 2. * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * f \
		- gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * f \
		+ gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * f \
		- 2. * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * f \
		+ gamma_uu.xz * a_l.z \
		+ 2. * gamma_uu.xz * Z_l.z * f \
		+ 2. * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * f \
		+ gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * f \
		- gamma_uu.xz * m * Z_l.z * f \
		- gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * f \
		- f * gamma_uu.xy * a_l.y \
		- f * gamma_uu.xz * a_l.z \
		+ a_l.x * gamma_uu.xx \
		- d_lll.x.xx * gamma_uu.xx * gamma_uu.xx * f \
		- m * Z_l.x * gamma_uu.xx * f)) / (gamma_uu.xx * gamma_uu.xx * f);\
	(results)->ptr[8] = (\
		-(2. * gamma_uu.xy * d_lll.y.xy \
		- 2. * gamma_uu.xy * d_lll.x.yy \
		+ gamma_uu.xz * d_lll.y.xz \
		+ gamma_uu.xz * d_lll.z.xy \
		- 2. * gamma_uu.xz * d_lll.x.yz \
		+ gamma_uu.yz * d_lll.y.yz \
		- gamma_uu.yz * d_lll.z.yy \
		+ gamma_uu.zz * d_lll.y.zz \
		- gamma_uu.zz * d_lll.z.yz \
		+ a_l.y \
		- 2. * Z_l.y \
		- 2. * d_lll.x.xy * gamma_uu.xx)) / (2. * gamma_uu.xx);\
	(results)->ptr[9] = (\
		-(gamma_uu.xy * d_lll.y.xz \
		+ gamma_uu.xy * d_lll.z.xy \
		- 2. * gamma_uu.xy * d_lll.x.yz \
		+ 2. * gamma_uu.xz * d_lll.z.xz \
		- 2. * gamma_uu.xz * d_lll.x.zz \
		- gamma_uu.yy * d_lll.y.yz \
		+ gamma_uu.yy * d_lll.z.yy \
		- gamma_uu.yz * d_lll.y.zz \
		+ gamma_uu.yz * d_lll.z.yz \
		+ a_l.z \
		- 2. * Z_l.z \
		- 2. * d_lll.x.xz * gamma_uu.xx)) / (2. * gamma_uu.xx);\
	(results)->ptr[10] = d_lll.y.xx;\
	(results)->ptr[11] = d_lll.y.xy;\
	(results)->ptr[12] = d_lll.y.xz;\
	(results)->ptr[13] = d_lll.y.yy;\
	(results)->ptr[14] = d_lll.y.yz;\
	(results)->ptr[15] = d_lll.y.zz;\
	(results)->ptr[16] = d_lll.z.xx;\
	(results)->ptr[17] = d_lll.z.xy;\
	(results)->ptr[18] = d_lll.z.xz;\
	(results)->ptr[19] = d_lll.z.yy;\
	(results)->ptr[20] = d_lll.z.yz;\
	(results)->ptr[21] = d_lll.z.zz;\
	(results)->ptr[22] = (\
		-(gamma_uu.xy * gamma_uu.xz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yy \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.yz \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.zz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.zz * gamma_uu.xx * d_lll.z.yz \
		- a_l.y * gamma_uu.xx)) / (2. * gamma_uu.xx);\
	(results)->ptr[23] = (gamma_uu.xy * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yz \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.z.yy \
		- gamma_uu.yy * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yy * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yz \
		+ a_l.z * gamma_uu.xx) / (2. * gamma_uu.xx);\
	(results)->ptr[24] = (\
		-(2. * gamma_uu.xx * sqrt_gUxx * K_ll.xz \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.yz \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yz * sqrt_gUxx \
		+ 2. * gamma_uu.xz * K_ll.zz * sqrt_gUxx \
		- gamma_uu.yy * gamma_uu.xx * d_lll.y.yz \
		+ gamma_uu.yy * gamma_uu.xx * d_lll.z.yy \
		- gamma_uu.yz * gamma_uu.xx * d_lll.y.zz \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.z.yz \
		+ a_l.z * gamma_uu.xx \
		- 2. * Z_l.z * gamma_uu.xx)) / (4. * gamma_uu.xx);\
	(results)->ptr[25] = (\
		-(2. * gamma_uu.xx * sqrt_gUxx * K_ll.xy \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.xz * d_lll.z.yy \
		+ 2. * gamma_uu.xy * K_ll.yy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.xz * d_lll.y.zz \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.z.yz \
		+ 2. * gamma_uu.xz * K_ll.yz * sqrt_gUxx \
		+ gamma_uu.yz * gamma_uu.xx * d_lll.y.yz \
		- gamma_uu.yz * gamma_uu.xx * d_lll.z.yy \
		+ gamma_uu.zz * gamma_uu.xx * d_lll.y.zz \
		- gamma_uu.zz * gamma_uu.xx * d_lll.z.yz \
		+ a_l.y * gamma_uu.xx \
		- 2. * Z_l.y * gamma_uu.xx)) / (4. * gamma_uu.xx);\
	(results)->ptr[26] = (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		- gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		+ gamma_uu.xy * a_l.y \
		+ 2. * gamma_uu.xy * K_ll.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		- gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		+ 2. * gamma_uu.xz * K_ll.xz * sqrt_gUxx \
		+ gamma_uu.xz * a_l.z \
		+ 2. * gamma_uu.yy * K_ll.yy * sqrt_gUxx \
		+ 4. * gamma_uu.yz * K_ll.yz * sqrt_gUxx \
		+ 2. * gamma_uu.zz * K_ll.zz * sqrt_gUxx \
		- 2. * Theta * sqrt_gUxx \
		+ 2. * Z_l.x * gamma_uu.xx) / (4. * gamma_uu.xx);\
	(results)->ptr[27] = (\
		-(gamma_uu.xx * gamma_uu.yy * d_lll.y.xy \
		- gamma_uu.xx * gamma_uu.yy * d_lll.x.yy \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.y.xz \
		+ gamma_uu.xx * gamma_uu.yz * d_lll.z.xy \
		- 2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz \
		+ gamma_uu.xx * gamma_uu.zz * d_lll.z.xz \
		- gamma_uu.xx * gamma_uu.zz * d_lll.x.zz \
		+ gamma_uu.xx * Z_l.x \
		- gamma_uu.xy * gamma_uu.xz * d_lll.y.xz \
		- gamma_uu.xy * gamma_uu.xz * d_lll.z.xy \
		+ 2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz \
		- gamma_uu.xy * gamma_uu.yz * d_lll.y.yz \
		+ gamma_uu.xy * gamma_uu.yz * d_lll.z.yy \
		- gamma_uu.xy * gamma_uu.zz * d_lll.y.zz \
		+ gamma_uu.xy * gamma_uu.zz * d_lll.z.yz \
		- gamma_uu.xy * gamma_uu.xy * d_lll.y.xy \
		+ gamma_uu.xy * Z_l.y \
		+ gamma_uu.xy * gamma_uu.xy * d_lll.x.yy \
		+ gamma_uu.xz * gamma_uu.yy * d_lll.y.yz \
		- gamma_uu.xz * gamma_uu.yy * d_lll.z.yy \
		+ gamma_uu.xz * gamma_uu.yz * d_lll.y.zz \
		- gamma_uu.xz * gamma_uu.yz * d_lll.z.yz \
		- gamma_uu.xz * gamma_uu.xz * d_lll.z.xz \
		+ gamma_uu.xz * Z_l.z \
		+ gamma_uu.xz * gamma_uu.xz * d_lll.x.zz \
		- Theta * sqrt_gUxx)) / (2. * sqrt_gUxx);\
	(results)->ptr[28] = (\
		-(gamma_uu.xx * d_lll.z.xz \
		- gamma_uu.xx * d_lll.x.zz \
		- gamma_uu.xy * d_lll.y.zz \
		+ gamma_uu.xy * d_lll.z.yz \
		- K_ll.zz * sqrt_gUxx)) / (2. * sqrt_gUxx);\
	(results)->ptr[29] = (\
		-(gamma_uu.xx * d_lll.y.xz \
		+ gamma_uu.xx * d_lll.z.xy \
		- 2. * gamma_uu.xx * d_lll.x.yz \
		- gamma_uu.xy * d_lll.y.yz \
		+ gamma_uu.xy * d_lll.z.yy \
		+ gamma_uu.xz * d_lll.y.zz \
		- gamma_uu.xz * d_lll.z.yz \
		- 2. * K_ll.yz * sqrt_gUxx)) / (4. * sqrt_gUxx);\
	(results)->ptr[30] = (2. * gamma_uu.xy * a_l.y * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * a_l.y * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xy * a_l.y * f * _1_sqrt_gUxx \
		- 8. * gamma_uu.xy * K_ll.xy * f_toThe_3_2 \
		+ 4. * gamma_uu.xy * K_ll.xy * sqrt_f \
		+ 4. * gamma_uu.xy * K_ll.xy * f_toThe_5_2 \
		- 8. * gamma_uu.xz * K_ll.xz * f_toThe_3_2 \
		+ 4. * gamma_uu.xz * K_ll.xz * f_toThe_5_2 \
		+ 4. * gamma_uu.xz * K_ll.xz * sqrt_f \
		+ 2. * gamma_uu.xz * a_l.z * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xz * a_l.z * fSq * _1_sqrt_gUxx \
		- 4. * gamma_uu.xz * a_l.z * f * _1_sqrt_gUxx \
		- 4. * gamma_uu.yy * K_ll.yy * f_toThe_3_2 \
		+ 2. * gamma_uu.yy * K_ll.yy * sqrt_f \
		+ 2. * gamma_uu.yy * K_ll.yy * f_toThe_5_2 \
		+ 4. * gamma_uu.yz * K_ll.yz * f_toThe_5_2 \
		+ 4. * gamma_uu.yz * K_ll.yz * sqrt_f \
		- 8. * gamma_uu.yz * K_ll.yz * f_toThe_3_2 \
		- 4. * gamma_uu.zz * K_ll.zz * f_toThe_3_2 \
		+ 2. * gamma_uu.zz * K_ll.zz * sqrt_f \
		+ 2. * gamma_uu.zz * K_ll.zz * f_toThe_5_2 \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.y.xy \
		- 4. * f * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.yy * d_lll.x.yy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.y.xy \
		- 2. * fSq * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 2. * f * sqrt_gUxx * gamma_uu.yy * m * d_lll.x.yy \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.y.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.z.xy \
		- 8. * f * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 8. * fSq * sqrt_gUxx * gamma_uu.yz * d_lll.x.yz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.y.xz \
		- 2. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.z.xy \
		- 4. * fSq * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		+ 4. * f * sqrt_gUxx * gamma_uu.yz * m * d_lll.x.yz \
		- 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.z.xz \
		+ 4. * fSq * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 4. * f * sqrt_gUxx * gamma_uu.zz * d_lll.x.zz \
		- 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.z.xz \
		+ 2. * f * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 2. * fSq * sqrt_gUxx * gamma_uu.zz * m * d_lll.x.zz \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy * _1_sqrt_gUxx \
		+ 8. * f * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 8. * fSq * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.y.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.z.xy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xz * m * d_lll.x.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.yz * m * d_lll.z.yy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.zz * m * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * Z_l.y * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.y.xy * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * m * Z_l.y * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xy * gamma_uu.xy * m * d_lll.x.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.y.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yy * m * d_lll.z.yy * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.y.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.yz * m * d_lll.z.yz * _1_sqrt_gUxx \
		- 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * Z_l.z * _1_sqrt_gUxx \
		- 4. * fSq * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.z.xz * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * m * Z_l.z * _1_sqrt_gUxx \
		- 2. * f * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 2. * fSq * gamma_uu.xz * gamma_uu.xz * m * d_lll.x.zz * _1_sqrt_gUxx \
		+ 4. * f * Z_l.x * sqrt_gUxx \
		- 4. * fSq * Z_l.x * sqrt_gUxx \
		- 2. * f_toThe_5_2 * m * Theta \
		+ 2. * f_toThe_3_2 * m * Theta \
		- 2. * f * m * Z_l.x * sqrt_gUxx \
		+ 2. * fSq * m * Z_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx \
		+ 2. * a_l.x * sqrt_gUxx * fSq \
		- 4. * a_l.x * sqrt_gUxx * f \
		- 4. * K_ll.xx * gamma_uu.xx * f_toThe_3_2 \
		+ 2. * K_ll.xx * gamma_uu.xx * sqrt_f \
		+ 2. * K_ll.xx * f_toThe_5_2 * gamma_uu.xx \
		- 4. * Theta * sqrt_f \
		+ 4. * Theta * f_toThe_3_2) / (gamma_uu.xx * sqrt_f * (4. \
		- 8. * f \
		+ 4. * fSq));\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> rotate <?=initCond_codeprefix?> _3sym3_rotate sym3_rotate

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */resultU,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */input,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	for (int j = 0; j < numStates; ++j) {\
		(resultU)->ptr[j] = 0;\
	}\
	\
	sym3 gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
\
	real sqrt_f = (eig)->alpha_sqrt_f / (eig)->alpha;\
	real f = sqrt_f * sqrt_f;\
	real fSq = f * f;\
	real f_toThe_3_2 = f * sqrt_f;\
\
	real f_minus_1 = f - 1.;\
	real sqrt_gUxx = sqrt(gamma_uu.xx);\
	real _1_sqrt_gUxx = 1. / sqrt_gUxx;\
	real _1_gammaUxx = _1_sqrt_gUxx * _1_sqrt_gUxx;\
	\
	real m = solver->m;\
\
	(resultU)->a_l.x = (gamma_uu.xy * gamma_uu.yz * (input)->ptr[14] * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yz * (input)->ptr[14] * f * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yz * (input)->ptr[19] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.yz * (input)->ptr[19] * f * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.zz * (input)->ptr[15] * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.zz * (input)->ptr[15] * f * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.zz * (input)->ptr[20] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.zz * (input)->ptr[20] * f * _1_sqrt_gUxx \
		- 2. * gamma_uu.xy * (input)->ptr[22] * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * (input)->ptr[22] * f * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[14] * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[14] * f * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[19] * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[19] * f * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yz * (input)->ptr[15] * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yz * (input)->ptr[15] * f * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yz * (input)->ptr[20] * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yz * (input)->ptr[20] * f * _1_sqrt_gUxx \
		- 2. * gamma_uu.xz * (input)->ptr[23] * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xz * (input)->ptr[23] * f * _1_sqrt_gUxx \
		+ f_toThe_3_2 * gamma_uu.xx * (input)->ptr[0] \
		- sqrt_f * gamma_uu.xx * (input)->ptr[0] \
		+ sqrt_f * gamma_uu.xx * (input)->ptr[30] \
		- f_toThe_3_2 * gamma_uu.xx * (input)->ptr[30] \
		+ 2. * f * (input)->ptr[27] \
		- 2. * f * (input)->ptr[4] \
		- f * m * (input)->ptr[27] \
		+ f * m * (input)->ptr[4]) / (sqrt_gUxx * (1. \
		- f));\
	(resultU)->a_l.y = (gamma_uu.xy * gamma_uu.xz * (input)->ptr[14] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[19] \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[15] \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[20] \
		- gamma_uu.yz * gamma_uu.xx * (input)->ptr[14] \
		+ gamma_uu.yz * gamma_uu.xx * (input)->ptr[19] \
		- gamma_uu.zz * gamma_uu.xx * (input)->ptr[15] \
		+ gamma_uu.zz * gamma_uu.xx * (input)->ptr[20] \
		+ 2. * (input)->ptr[22] * gamma_uu.xx) / gamma_uu.xx;\
	(resultU)->a_l.z = (\
		-(gamma_uu.xy * gamma_uu.xz * (input)->ptr[15] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[20] \
		+ gamma_uu.xy * gamma_uu.xy * (input)->ptr[14] \
		- gamma_uu.xy * gamma_uu.xy * (input)->ptr[19] \
		- gamma_uu.yy * gamma_uu.xx * (input)->ptr[14] \
		+ gamma_uu.yy * gamma_uu.xx * (input)->ptr[19] \
		- gamma_uu.yz * gamma_uu.xx * (input)->ptr[15] \
		+ gamma_uu.yz * gamma_uu.xx * (input)->ptr[20] \
		- 2. * (input)->ptr[23] * gamma_uu.xx)) / gamma_uu.xx;\
	(resultU)->d_lll.x.xx = (gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f_toThe_3_2 * (input)->ptr[4] \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f_toThe_3_2 * m * (input)->ptr[4] \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * sqrt_f \
		- sqrt_gUxx * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] * _1_sqrt_gUxx * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] * f_toThe_3_2 * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * f_toThe_3_2 * sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] * _1_sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] * f_toThe_3_2 * _1_sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * f_toThe_3_2 * sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * f_toThe_3_2 \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * sqrt_f * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * f_toThe_3_2 \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * sqrt_f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * _1_sqrt_gUxx * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f_toThe_3_2 * _1_sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * _1_sqrt_gUxx * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f_toThe_3_2 * _1_sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * _1_sqrt_gUxx * sqrt_f \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f_toThe_3_2 * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * _1_sqrt_gUxx * sqrt_f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f_toThe_3_2 * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f_toThe_3_2 * sqrt_gUxx * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f_toThe_3_2 * (input)->ptr[27] * _1_sqrt_gUxx \
		+ 3. * gamma_uu.xy * gamma_uu.xy * f_toThe_3_2 * (input)->ptr[27] * sqrt_gUxx * gamma_uu.yy \
		- 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.yy * sqrt_f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * f_toThe_3_2 \
		+ 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * f_toThe_3_2 * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx * gamma_uu.yy * sqrt_f \
		+ 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.yy * f_toThe_3_2 * gamma_uu.xx \
		- 3. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.yy * sqrt_f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * sqrt_f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f_toThe_3_2 \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * sqrt_f * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * sqrt_f * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * f_toThe_3_2 * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * sqrt_f * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * sqrt_gUxx * sqrt_f \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * f_toThe_3_2 * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * f_toThe_3_2 * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * sqrt_f \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * sqrt_gUxx * sqrt_f \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * f_toThe_3_2 * sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * f_toThe_3_2 * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] * sqrt_f * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] * f_toThe_3_2 * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] * f_toThe_3_2 * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] * sqrt_f * gamma_uu.xx * gamma_uu.xx \
		- f_toThe_3_2 * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] \
		+ 2. * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * _1_sqrt_gUxx \
		- 3. * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * sqrt_gUxx * gamma_uu.yy \
		- f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * m * (input)->ptr[4] * _1_sqrt_gUxx \
		+ 2. * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * m * (input)->ptr[4] * sqrt_gUxx * gamma_uu.yy \
		+ f_toThe_3_2 * m * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] \
		+ f_toThe_3_2 * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * _1_sqrt_gUxx \
		- 2. * f_toThe_3_2 * m * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx * gamma_uu.yy \
		+ (input)->ptr[0] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * (input)->ptr[0] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ 2. * (input)->ptr[0] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * f \
		- (input)->ptr[0] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f \
		+ (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		- (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f \
		- (input)->ptr[30] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * (input)->ptr[30] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		- (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ (input)->ptr[30] * sqrt_gUxx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * (input)->ptr[30] * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.yy * gamma_uu.yy \
		- (input)->ptr[7] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * sqrt_f \
		+ (input)->ptr[7] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * f_toThe_3_2 * gamma_uu.yy * gamma_uu.yy \
		+ 2. * (input)->ptr[7] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * sqrt_f * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[7] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_f * gamma_uu.xx \
		+ (input)->ptr[7] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		- 2. * (input)->ptr[7] * f_toThe_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx) / (\
		-gamma_uu.xx * sqrt_f * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ 2. * gamma_uu.xx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.xx * f * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->d_lll.x.xy = (2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] \
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		- 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] \
		+ 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] \
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] \
		+ gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] \
		+ gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] \
		+ gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx \
		+ gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[27] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx \
		+ gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[4] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] \
		- gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] \
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * (input)->ptr[25] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy \
		- (input)->ptr[8] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy \
		- (input)->ptr[8] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ 2. * (input)->ptr[8] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy) / (\
		-gamma_uu.xx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->d_lll.x.xz = (\
		-(gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		- gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		+ (input)->ptr[1] \
		+ (input)->ptr[24] \
		- (input)->ptr[9] * gamma_uu.xx)) / gamma_uu.xx;\
	(resultU)->d_lll.x.yy = (2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		- gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.yy \
		- gamma_uu.xx * (input)->ptr[27] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * (input)->ptr[27] * gamma_uu.yy \
		- gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.yy \
		+ gamma_uu.xx * (input)->ptr[4] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * (input)->ptr[4] * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] \
		+ 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] * gamma_uu.xx * gamma_uu.yy \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] \
		- 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[25] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[25] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[2] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[2] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * (input)->ptr[14] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * (input)->ptr[14] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xz * (input)->ptr[14] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xz * (input)->ptr[19] * _1_sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xz * (input)->ptr[19] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		- gamma_uu.xz * (input)->ptr[19] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ gamma_uu.xz * (input)->ptr[1] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xz * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx \
		- gamma_uu.xz * (input)->ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xz * (input)->ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy \
		+ (input)->ptr[11] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * (input)->ptr[11] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ (input)->ptr[11] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy) / (sqrt_gUxx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->d_lll.x.yz = (\
		-(gamma_uu.xy * (input)->ptr[14] * _1_sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[19] * _1_sqrt_gUxx \
		- gamma_uu.xz * (input)->ptr[15] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[20] * _1_sqrt_gUxx \
		- (input)->ptr[12] * sqrt_gUxx \
		- (input)->ptr[17] * sqrt_gUxx \
		- 2. * (input)->ptr[29] \
		+ 2. * (input)->ptr[6])) / (2. * sqrt_gUxx);\
	(resultU)->d_lll.x.zz = (\
		-(gamma_uu.xy * (input)->ptr[15] * _1_sqrt_gUxx \
		- gamma_uu.xy * (input)->ptr[20] * _1_sqrt_gUxx \
		- (input)->ptr[18] * sqrt_gUxx \
		- (input)->ptr[28] \
		+ (input)->ptr[5])) / sqrt_gUxx;\
	(resultU)->d_lll.y.xx = (input)->ptr[10];\
	(resultU)->d_lll.y.xy = (input)->ptr[11];\
	(resultU)->d_lll.y.xz = (input)->ptr[12];\
	(resultU)->d_lll.y.yy = (input)->ptr[13];\
	(resultU)->d_lll.y.yz = (input)->ptr[14];\
	(resultU)->d_lll.y.zz = (input)->ptr[15];\
	(resultU)->d_lll.z.xx = (input)->ptr[16];\
	(resultU)->d_lll.z.xy = (input)->ptr[17];\
	(resultU)->d_lll.z.xz = (input)->ptr[18];\
	(resultU)->d_lll.z.yy = (input)->ptr[19];\
	(resultU)->d_lll.z.yz = (input)->ptr[20];\
	(resultU)->d_lll.z.zz = (input)->ptr[21];\
	(resultU)->K_ll.xx = (gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[26] * f \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[3] * f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[29] * f \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * f * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] * f * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * (input)->ptr[6] * f \
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * f * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] * f * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[25] * sqrt_gUxx * f \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * sqrt_gUxx * f \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f \
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xx * sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * (input)->ptr[2] * f * sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] * f * gamma_uu.xx * sqrt_gUxx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy * gamma_uu.yy * f * gamma_uu.xx * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * f * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy * gamma_uu.yy * f * gamma_uu.xx * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f \
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * f * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[28] * f * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy * gamma_uu.xx \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * (input)->ptr[5] * f * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[27] \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[27] * gamma_uu.yy * gamma_uu.xx \
		+ 4. * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[4] \
		+ 4. * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * (input)->ptr[4] * gamma_uu.yy * gamma_uu.xx \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * sqrt_gUxx \
		- 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * sqrt_gUxx * f \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * sqrt_gUxx \
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f * sqrt_gUxx \
		- 4. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy \
		+ 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] * f * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.xx \
		+ 2. * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		- 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * sqrt_gUxx \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		+ 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[1] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[1] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx \
		+ 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * f \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f \
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[24] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * f * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[28] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[28] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * f * gamma_uu.xx * gamma_uu.xx \
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * (input)->ptr[5] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * (input)->ptr[5] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- f * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- f * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] \
		- 3. * f * m * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 3. * f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * gamma_uu.yy * gamma_uu.xx \
		- f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] \
		- 3. * f * m * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		+ 3. * f * m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * gamma_uu.yy * gamma_uu.xx \
		+ f * m * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		+ f * m * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * (input)->ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy \
		+ (input)->ptr[0] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f \
		+ (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.xx \
		+ 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- 3. * (input)->ptr[0] * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy \
		+ (input)->ptr[30] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f \
		+ (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx \
		+ 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		- 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx \
		- (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.xx \
		- 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx \
		+ 3. * (input)->ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx) \
		/ (gamma_uu.xx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 3. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ 3. * gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f \
		+ 3. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy \
		- 3. * gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f));\
	(resultU)->K_ll.xy = (2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] \
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		+ 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] \
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[26] \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * (input)->ptr[26] * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[3] \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * (input)->ptr[3] * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[25] \
		- gamma_uu.xx * gamma_uu.yy * (input)->ptr[25] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[2] \
		+ gamma_uu.xx * gamma_uu.yy * (input)->ptr[2] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] \
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] \
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[1] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] \
		+ gamma_uu.xy * gamma_uu.xz * (input)->ptr[24] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[27] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[27] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[4] * sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[4] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[29] \
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * (input)->ptr[6] \
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx) / (\
		-sqrt_gUxx * (gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy));\
	(resultU)->K_ll.xz = (\
		-(gamma_uu.xy * (input)->ptr[29] * _1_sqrt_gUxx \
		+ gamma_uu.xy * (input)->ptr[6] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[28] * _1_sqrt_gUxx \
		+ gamma_uu.xz * (input)->ptr[5] * _1_sqrt_gUxx \
		- (input)->ptr[1] \
		+ (input)->ptr[24])) / sqrt_gUxx;\
	(resultU)->K_ll.yy = (sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[25] \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * (input)->ptr[25] * gamma_uu.yy \
		- sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * (input)->ptr[2] \
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * (input)->ptr[2] * gamma_uu.yy \
		- sqrt_gUxx * gamma_uu.xz * (input)->ptr[1] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.xz * (input)->ptr[1] * gamma_uu.yy \
		+ sqrt_gUxx * gamma_uu.xz * (input)->ptr[24] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.xz * (input)->ptr[24] * gamma_uu.yy \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[29] * gamma_uu.yy \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.xy * gamma_uu.xy \
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * (input)->ptr[6] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[28] * gamma_uu.yy \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * (input)->ptr[5] * gamma_uu.yy \
		+ gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[26] * gamma_uu.yy \
		+ gamma_uu.xx * (input)->ptr[27] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * (input)->ptr[27] * gamma_uu.yy \
		- gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.xy * gamma_uu.xy \
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * (input)->ptr[3] * gamma_uu.yy \
		+ gamma_uu.xx * (input)->ptr[4] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xx * gamma_uu.xx * (input)->ptr[4] * gamma_uu.yy \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] \
		- 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[29] * gamma_uu.xx * gamma_uu.yy \
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] \
		- 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[6] * gamma_uu.xx * gamma_uu.yy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[28] * gamma_uu.xy * gamma_uu.xy \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xx * gamma_uu.yy \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[5] * gamma_uu.xy * gamma_uu.xy) / (\
		-(gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy \
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy \
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.yy));\
	(resultU)->K_ll.yz = (input)->ptr[29] \
		+ (input)->ptr[6];\
	(resultU)->K_ll.zz = (input)->ptr[28] \
		+ (input)->ptr[5];\
	(resultU)->Theta = (input)->ptr[27] \
		+ (input)->ptr[4];\
	(resultU)->Z_l.x = (\
		-(gamma_uu.xy * (input)->ptr[22] \
		+ gamma_uu.xz * (input)->ptr[23] \
		- (input)->ptr[26] * gamma_uu.xx \
		- (input)->ptr[3] * gamma_uu.xx)) / gamma_uu.xx;\
	(resultU)->Z_l.y = (input)->ptr[22] \
		+ (input)->ptr[25] \
		+ (input)->ptr[2];\
	(resultU)->Z_l.z = (input)->ptr[1] \
		+ (input)->ptr[23] \
		+ (input)->ptr[24];\
\
	(resultU)->a_l = real3_swap((resultU)->a_l, n.side);							/* 0-2 */\
	(resultU)->d_lll = _3sym3_swap((resultU)->d_lll, n.side);						/* 3-20 */\
	(resultU)->K_ll = sym3_swap((resultU)->K_ll, n.side);							/* 21-26 */\
	(resultU)->Z_l = real3_swap((resultU)->Z_l, n.side);							/* 28-30 */\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> _3sym3_rotate sym3_rotate
// used by roe, some plm
//so long as roeUseFluxFromCons isn't set for the roe solver, 
// and fluxFromCons is provided/unused,
// eigen_fluxTransform isn't needed.
// but some solvers do use a boilerplate right(lambda(left(U)))
//however if you want to use the HLL solver then fluxFromCons is needed
//...however fluxFromCons is not provided by this eqn.	


//notice this paper uses the decomposition alpha A = R Lambda L
// so this computation is for alpha A
#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */input,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	for (int i = 0; i < numStates; ++i) {\
		(resultFlux)->ptr[i] = 0;\
	}\
\
	/* input */\
	real3 a_l = real3_swap((input)->a_l, n.side);			/* 0-2 */\
	_3sym3 d_lll = _3sym3_swap((input)->d_lll, n.side);		/* 3-20 */\
	sym3 K_ll = sym3_swap((input)->K_ll, n.side);			/* 21-26 */\
	real Theta = (input)->Theta;							/* 27 */\
	real3 Z_l = real3_swap((input)->Z_l, n.side);			/* 28-30 */\
\
	sym3 gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
\
	real sqrt_f = (eig)->alpha_sqrt_f / (eig)->alpha;\
	real f = sqrt_f * sqrt_f;\
\
	real m = solver->m;\
\
	/* TODO this still needs to know what vars are going into it */\
	/* so it can use the correctly swapped vars */\
\
	(resultFlux)->ptr[0] = f * (gamma_uu.xx * (input)->ptr[21] \
		+ 2. * gamma_uu.xy * (input)->ptr[22] \
		+ 2. * gamma_uu.xz * (input)->ptr[23] \
		+ gamma_uu.yy * (input)->ptr[24] \
		+ 2. * gamma_uu.yz * (input)->ptr[25] \
		+ gamma_uu.zz * (input)->ptr[26] \
		- m * (input)->ptr[27]);\
	(resultFlux)->ptr[1] = 0.;\
	(resultFlux)->ptr[2] = 0.;\
	(resultFlux)->ptr[3] = (input)->ptr[21];\
	(resultFlux)->ptr[4] = (input)->ptr[22];\
	(resultFlux)->ptr[5] = (input)->ptr[23];\
	(resultFlux)->ptr[6] = (input)->ptr[24];\
	(resultFlux)->ptr[7] = (input)->ptr[25];\
	(resultFlux)->ptr[8] = (input)->ptr[26];\
	(resultFlux)->ptr[9] = 0.;\
	(resultFlux)->ptr[10] = 0.;\
	(resultFlux)->ptr[11] = 0.;\
	(resultFlux)->ptr[12] = 0.;\
	(resultFlux)->ptr[13] = 0.;\
	(resultFlux)->ptr[14] = 0.;\
	(resultFlux)->ptr[15] = 0.;\
	(resultFlux)->ptr[16] = 0.;\
	(resultFlux)->ptr[17] = 0.;\
	(resultFlux)->ptr[18] = 0.;\
	(resultFlux)->ptr[19] = 0.;\
	(resultFlux)->ptr[20] = 0.;\
	(resultFlux)->ptr[21] = \
		-(gamma_uu.yy * (input)->ptr[10] \
		- gamma_uu.yy * (input)->ptr[6] \
		+ gamma_uu.yz * (input)->ptr[11] \
		+ gamma_uu.yz * (input)->ptr[16] \
		- 2. * gamma_uu.yz * (input)->ptr[7] \
		+ gamma_uu.zz * (input)->ptr[17] \
		- gamma_uu.zz * (input)->ptr[8] \
		- (input)->ptr[0] \
		+ 2. * (input)->ptr[28]);\
	(resultFlux)->ptr[22] = (2. * gamma_uu.xy * (input)->ptr[10] \
		- 2. * gamma_uu.xy * (input)->ptr[6] \
		+ gamma_uu.xz * (input)->ptr[11] \
		+ gamma_uu.xz * (input)->ptr[16] \
		- 2. * gamma_uu.xz * (input)->ptr[7] \
		+ gamma_uu.yz * (input)->ptr[13] \
		- gamma_uu.yz * (input)->ptr[18] \
		+ gamma_uu.zz * (input)->ptr[14] \
		- gamma_uu.zz * (input)->ptr[19] \
		+ (input)->ptr[1] \
		- 2. * (input)->ptr[29]) / 2.;\
	(resultFlux)->ptr[23] = (gamma_uu.xy * (input)->ptr[11] \
		+ gamma_uu.xy * (input)->ptr[16] \
		- 2. * gamma_uu.xy * (input)->ptr[7] \
		+ 2. * gamma_uu.xz * (input)->ptr[17] \
		- 2. * gamma_uu.xz * (input)->ptr[8] \
		- gamma_uu.yy * (input)->ptr[13] \
		+ gamma_uu.yy * (input)->ptr[18] \
		- gamma_uu.yz * (input)->ptr[14] \
		+ gamma_uu.yz * (input)->ptr[19] \
		+ (input)->ptr[2] \
		- 2. * (input)->ptr[30]) / 2.;\
	(resultFlux)->ptr[24] = \
		-(gamma_uu.xx * (input)->ptr[10] \
		- gamma_uu.xx * (input)->ptr[6] \
		+ gamma_uu.xz * (input)->ptr[13] \
		- gamma_uu.xz * (input)->ptr[18]);\
	(resultFlux)->ptr[25] = (\
		-(gamma_uu.xx * (input)->ptr[11] \
		+ gamma_uu.xx * (input)->ptr[16] \
		- 2. * gamma_uu.xx * (input)->ptr[7] \
		- gamma_uu.xy * (input)->ptr[13] \
		+ gamma_uu.xy * (input)->ptr[18] \
		+ gamma_uu.xz * (input)->ptr[14] \
		- gamma_uu.xz * (input)->ptr[19])) / 2.;\
	(resultFlux)->ptr[26] = \
		-(gamma_uu.xx * (input)->ptr[17] \
		- gamma_uu.xx * (input)->ptr[8] \
		- gamma_uu.xy * (input)->ptr[14] \
		+ gamma_uu.xy * (input)->ptr[19]);\
	(resultFlux)->ptr[27] = \
		-(gamma_uu.xx * gamma_uu.yy * (input)->ptr[10] \
		- gamma_uu.xx * gamma_uu.yy * (input)->ptr[6] \
		+ gamma_uu.xx * gamma_uu.yz * (input)->ptr[11] \
		+ gamma_uu.xx * gamma_uu.yz * (input)->ptr[16] \
		- 2. * gamma_uu.xx * gamma_uu.yz * (input)->ptr[7] \
		+ gamma_uu.xx * gamma_uu.zz * (input)->ptr[17] \
		- gamma_uu.xx * gamma_uu.zz * (input)->ptr[8] \
		+ gamma_uu.xx * (input)->ptr[28] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[11] \
		- gamma_uu.xy * gamma_uu.xz * (input)->ptr[16] \
		+ 2. * gamma_uu.xy * gamma_uu.xz * (input)->ptr[7] \
		- gamma_uu.xy * gamma_uu.yz * (input)->ptr[13] \
		+ gamma_uu.xy * gamma_uu.yz * (input)->ptr[18] \
		- gamma_uu.xy * gamma_uu.zz * (input)->ptr[14] \
		+ gamma_uu.xy * gamma_uu.zz * (input)->ptr[19] \
		- gamma_uu.xy * gamma_uu.xy * (input)->ptr[10] \
		+ gamma_uu.xy * (input)->ptr[29] \
		+ gamma_uu.xy * gamma_uu.xy * (input)->ptr[6] \
		+ gamma_uu.xz * gamma_uu.yy * (input)->ptr[13] \
		- gamma_uu.xz * gamma_uu.yy * (input)->ptr[18] \
		+ gamma_uu.xz * gamma_uu.yz * (input)->ptr[14] \
		- gamma_uu.xz * gamma_uu.yz * (input)->ptr[19] \
		- gamma_uu.xz * gamma_uu.xz * (input)->ptr[17] \
		+ gamma_uu.xz * (input)->ptr[30] \
		+ gamma_uu.xz * gamma_uu.xz * (input)->ptr[8]);\
	(resultFlux)->ptr[28] = gamma_uu.xy * (input)->ptr[22] \
		+ gamma_uu.xz * (input)->ptr[23] \
		+ gamma_uu.yy * (input)->ptr[24] \
		+ 2. * gamma_uu.yz * (input)->ptr[25] \
		+ gamma_uu.zz * (input)->ptr[26] \
		- (input)->ptr[27];\
	(resultFlux)->ptr[29] = \
		-(gamma_uu.xx * (input)->ptr[22] \
		+ gamma_uu.xy * (input)->ptr[24] \
		+ gamma_uu.xz * (input)->ptr[25]);\
	(resultFlux)->ptr[30] = \
		-(gamma_uu.xx * (input)->ptr[23] \
		+ gamma_uu.xy * (input)->ptr[25] \
		+ gamma_uu.xz * (input)->ptr[26]);\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=SETBOUNDS_NOGHOST?> <?=initCond_codeprefix?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const deriv = derivBuf + index;

	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real const f = calc_f(U->alpha);	/* could be based on alpha... */

	real3 const S_l = real3_zero;
	sym3 const S_ll = sym3_zero;
	real const S = 0.;
	real const rho = 0.;

	/*  source terms */
	
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			/* K^i_j */
	real const trK = real3x3_trace(K_ul);								/* K^k_k */
	sym3 const KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		/* KSq_ij = K_ik K^k_j */

	/* d_llu = d_ij^k = d_ijl * gamma^lk */
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};

	/* d_ull = d^i_jk = gamma^il d_ljk */
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	/* e_l = d^j_ji */
	real3 const e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end	?>,
<? end
?>	};

	/* conn^k_ij = d_ij^k + d_ji^k - d^k_ij */
	_3sym3 const conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d_lll.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};


	/* d_l = d_i = d_ij^j */
	real3 const d_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_trace(d_llu[<?=i-1?>]),
<? end
?>	};
	
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 const Z_u = sym3_real3_mul(gamma_uu, U->Z_l);

	/* d_luu = d_i^jk = gamma^jl d_il^k */
	_3sym3 const d_luu = (_3sym3){
<? for i,xi in ipairs(xNames) do		
?>		.<?=xi?> = sym3_real3x3_to_sym3_mul(gamma_uu, d_llu[<?=i-1?>]),
<? end
?>	};

	/* alpha_,t = shift terms - alpha^2 f (gamma^ij K_ij - m Theta) */
	real const f_alphaSq = calc_f_alphaSq(U->alpha);
	deriv->alpha += -f_alphaSq * (trK - solver->m * U->Theta);
	
	/* gamma_ij,t = shift terms - 2 alpha K_ij */
	deriv->gamma_ll = sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));

	/* 2005 Bona et al A.1 */
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->K_ll.<?=xij?> += U->alpha * (
		.5 * (U->a_l.<?=xi?> * d_l.<?=xj?> + U->a_l.<?=xj?> * d_l.<?=xi?>)
		- .25 * (U->a_l.<?=xj?> * (2. * e_l.<?=xi?> - d_l.<?=xi?>) + U->a_l.<?=xi?> * (2. * e_l.<?=xj?> - d_l.<?=xj?>))
		- U->a_l.<?=xi?> * U->Z_l.<?=xj?> - U->a_l.<?=xj?> * U->Z_l.<?=xi?>
		+ (trK - 2. * U->Theta) * U->K_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- .5 * U->a_l.<?=xk?> * conn_ull.<?=xk?>.<?=xij?>
		+ .5 * U->a_l.<?=xk?> * d_ull.<?=xk?>.<?=xij?>
		- 2. * e_l.<?=xk?> * (d_llu[<?=i-1?>].<?=xj?>.<?=xk?> + d_llu[<?=j-1?>].<?=xi?>.<?=xk?>)
		+ (d_l.<?=xk?> + U->a_l.<?=xk?> - 2. * U->Z_l.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		- 2. * K_ul.<?=xk?>.<?=xi?> * U->K_ll.<?=sym(k,j)?>
<?		for l,xl in ipairs(xNames) do
?>		+ 2. * (d_llu[<?=i-1?>].<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,j)?> + d_llu[<?=j-1?>].<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,i)?>)
		- conn_ull.<?=xk?>.<?=sym(l,j)?> * conn_ull.<?=xl?>.<?=sym(k,i)?>
<?		end
	end
?> - 8. * M_PI * (S_ll.<?=xij?> - .5 * U->gamma_ll.<?=xij?> * (S - rho)));
<? end
?>

	/* 2005 Bona et al A.2 */
<? for i,xi in ipairs(xNames) do
?>	deriv->Z_l.<?=xi?> += U->alpha * (
		U->a_l.<?=xi?> * (trK - 2. * U->Theta)
<?	for k,xk in ipairs(xNames) do
?>		- U->a_l.<?=xk?> * K_ul.<?=xk?>.<?=xi?>
		+ K_ul.<?=xk?>.<?=xi?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?		for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * conn_ull.<?=xr?>.<?=sym(k,i)?>
<?		end
	end
?>
	) - 8 * M_PI * U->alpha * S_l.<?=xi?>;
<? end 
?>
	/* 2005 Bona et al A.3 */
	deriv->Theta += U->alpha * .5 * ( trK * (trK - 2. * U->Theta)
<? 
for k,xk in ipairs(xNames) do 
?>		+ 2. * U->a_l.<?=xk?> * (d_u.<?=xk?> - e_u.<?=xk?> - 2. * Z_u.<?=xk?>) 
		- d_u.<?=xk?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?	for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * K_ul.<?=xr?>.<?=xk?>
<?		for s,xs in ipairs(xNames) do
?>		+ d_luu.<?=xk?>.<?=sym(r,k)?> * conn_ull.<?=xk?>.<?=sym(r,s)?>
<?		end
	end
end?>
	) - 8. * M_PI * U->alpha * rho;

}

//// MODULE_NAME: <?=constrainU?>

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;

	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real const rho = 0.;

	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			/* K^i_j */
	real const trK = real3x3_trace(K_ul);								/* K^k_k */
	sym3 const KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		/* KSq_ij = K_ik K^k_j */

	/* d_ull = d^i_jk = gamma^il d_ljk */
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	/* e_l = d^j_ji */
	real3 const e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end	?>,
<? end
?>	};

	/* d_llu = d_ij^k = d_ijl * gamma^lk */
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};
	
	/* d_l = d_ij^j */
	real3 const d_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_trace(d_llu[<?=i-1?>]),
<? end
?>	};

	/* conn^k_ij = d_ij^k + d_ji^k - d^k_ij */
	_3sym3 const conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d_lll.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	real3 const V_l = real3_sub(d_l, e_l);

	sym3 const R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>
			+ conn_ull.<?=xk?>.<?=xij?> * (V_l.<?=xk?> - e_l.<?=xk?>)

<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(j,l)?>
			- 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=l-1?>].<?=xj?>.<?=xk?>
			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=j-1?>].<?=xl?>.<?=xk?>
			+ 2. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=k-1?>].<?=xj?>.<?=xl?>
			- 3. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=j-1?>].<?=xk?>.<?=xl?>
<? 		end
	end
?>		,
<? end
?>	};


	/* calculate the Hamiltonian and momentum constraints */
	/* scaled down by 1/8 to match B&S BSSNOK equations ... maybe I'll scale theirs up by 8 ... */
	/* B&S eqn 2.125 ... divded by two */
	/* Alcubierre eqn 2.5.9 */
	/* H = 1/2 (R + K^2 - K_ij K^ij) - 8 pi rho */
	real const R = sym3_dot(R_ll, gamma_uu);
	real const tr_KSq = sym3_dot(KSq_ll, gamma_uu);
	U->H = .5 * (R + trK * trK - tr_KSq) - 8. * M_PI * rho;
}
