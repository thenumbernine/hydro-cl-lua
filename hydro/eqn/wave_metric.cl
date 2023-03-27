//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: real3x3

/*
background metric ADM decomposition
this assumes the ADM spatial metric gamma_ij is equal to the grid metric

for the wave equation d'Lambertian phi = 0
i.e. g^tt phi_,tt + 2 g^ti phi_;ti + g^ij phi_;ij = 0
ADM is defined such that 
 	g^tt = -alpha^-2
	g^ti = alpha^-2 beta^i 
	g^ij = gamma^ij - alpha^-2 beta^i beta^j
TODO make this configurable somehow
or make it modular enough to merge with BSSNOK
*/

#define /*real */metric_alpha(/*real3 const */pt) \
	(<?=eqn:compile(eqn.metric.alpha)?>)

#define /*real3 */metric_beta_u(/*real3 const */pt) \
	real3{\
<? for i,xi in ipairs(xNames) do --\
?>		<?=eqn:compile(eqn.metric.beta_u[i])?>,\
<? end --\
?>	}

#define /*real */metric_K(/*real3 const */pt) \
	(<?=eqn:compile(eqn.metric.K)?>)

#define /*real */metric_f(/*real3 const */pt) \
	(<?=eqn:compile(eqn.metric.f)?>)

#define /*real3 */metric_partial_alpha_l(/*real3 const */pt) \
	real3{\
<? for i,xi in ipairs(xNames) do --\
?>		<?=eqn:compile(eqn.metric.alpha:diff(eqn.metric.coords[i])())?>,\
<? end --\
?>	}

//partial_beta_ul[i][j] = beta^i_,j
#define /*real3x3 */metric_partial_beta_ul(/*real3 const */pt) \
	real3x3{\
<? for i,xi in ipairs(xNames) do --\
?>		{\
<?	for j,xj in ipairs(xNames) do --\
?>			<?=eqn:compile(eqn.metric.beta_u[i]:diff(eqn.metric.coords[j])())?>,\
<?	end --\
?>		},\
<?	end --\
?>	}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=cartesianToCoord?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool const lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;
	
	<?=initCode()?>

<? if eqn.usePressure then
?>	U->Pi = <?=scalar?>_from_real(P);
<? else		
?>	U->Pi = <?=scalar?>_from_real(rho);
<? end		
?>	U->Psi_l = <?=vec3?>_from_real3(cartesianToCoord(v, x));
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=normal_t?> <?=cons_t?> <?=eqn_common?>
// eqn_common has metric_alpha

// What's the difference between <?=eigen_fluxTransform?> and <?=fluxFromCons?>?
// The difference is that the flux matrix of this is based on 'eig', which is derived from U's ... especially UL & UR in the case of the Roe solver
// whereas that of <?=fluxFromCons?> is based purely on 'U'.
// Since hydro/eqn/wave has no <?=eigen_t?> info derived from U, the two functions are identical.
#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const alpha = metric_alpha((cell)->pos);\
	real const beta_n = normal_vecDotN1(n, metric_beta_u((cell)->pos));\
	\
	real3 const nL = normal_l1(n);\
	real3 const nU = normal_u1(n);\
	\
	/* F^Pi = -c (Pi beta_n + alpha Psi_i n^i) */\
	(resultFlux)->Pi = <?=scalar?>_real_mul(\
		/* Pi beta_n + alpha Psi_i n^i */\
		real_<?=scalar?>_add(\
		\
			/* Pi beta_n: */\
			<?=scalar?>_real_mul((U)->Pi, beta_n), \
			\
			/* alpha Psi_i n^i: */\
			<?=scalar?>_real_mul(\
				/* Psi_i n^i: */\
				<?=vec3?>_real3_dot(\
					(U)->Psi_l, \
					nU\
				), \
				alpha\
			)\
		), \
		-solver->wavespeed\
	);\
	\
	/* F^{Psi_j} = -c (alpha Pi n_j + Psi_j beta_n) */\
	(resultFlux)->Psi_l = <?=vec3?>_real_mul(\
		<?=vec3?>_add(\
			/* Psi_j beta_n: */\
			<?=vec3?>_real_mul(\
				(U)->Psi_l,\
				beta_n\
			),\
			\
			/* alpha Pi n_j: */\
			<?=vec3?>_<?=scalar?>_mul(\
				/* n: */\
				<?=vec3?>_from_real3(nL),\
				\
				/* alpha Pi: */\
				<?=scalar?>_real_mul(\
					(U)->Pi, \
					alpha\
				)\
			)\
		),\
		/* -c */\
		-solver->wavespeed\
	);\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>

// used by PLM
void <?=calcCellMinMaxEigenvalues?>(
	<?=range_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const x,
	<?=normal_t?> const n
) {
	real const alpha_nLen = metric_alpha(x) * normal_len(n);
	real const beta_n = normal_vecDotN1(n, metric_beta_u(x));
	result->min = solver->wavespeed * (-beta_n - alpha_nLen);
	result->max = solver->wavespeed * (-beta_n + alpha_nLen);
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=eigen_t?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) \
	*(resultEig) = (<?=eigen_t?>){}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=eigen_t?>

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*real3 const */n\
) \
	*(resultEig) = (<?=eigen_t?>){};

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=waves_t?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */x,\
	/*<?=normal_t?> const */n\
) { \
	real const nLen = normal_len(n);\
	real const invDenom = 1. / (nLen * nLen);\
	(result)->ptr[0] = .5 * invDenom * (\
		  (X)->ptr[0] * nLen\
		+ (X)->ptr[1] * normal_u1x(n)\
		+ (X)->ptr[2] * normal_u1y(n)\
		+ (X)->ptr[3] * normal_u1z(n)\
	);\
	(result)->ptr[1] = invDenom * (\
		  (X)->ptr[1] * normal_u2x(n)\
		+ (X)->ptr[2] * normal_u2y(n)\
		+ (X)->ptr[3] * normal_u2z(n)\
	);\
	(result)->ptr[2] = invDenom * (\
		  (X)->ptr[1] * normal_u3x(n)\
		+ (X)->ptr[2] * normal_u3y(n)\
		+ (X)->ptr[3] * normal_u3z(n)\
	);\
	(result)->ptr[3] = .5 * invDenom * (\
		- (X)->ptr[0] * nLen\
		+ (X)->ptr[1] * normal_u1x(n)\
		+ (X)->ptr[2] * normal_u1y(n)\
		+ (X)->ptr[3] * normal_u1z(n)\
	);\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=waves_t?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X_,\
	/*real3 const */x,\
	/*<?=normal_t?> const */n\
) {\
	<?=scalar?>* X = (<?=scalar?>*)(X_)->ptr;\
	real nLen = normal_len(n);\
	(result)->ptr[0] = \
		  (X[0] - X[3]) * nLen;\
	(result)->ptr[1] = \
		  X[0] * normal_l1x(n) \
		+ X[1] * normal_l2x(n) \
		+ X[2] * normal_l3x(n) \
		+ X[3] * normal_l1x(n);\
	(result)->ptr[2] = \
		  X[0] * normal_l1y(n) \
		+ X[1] * normal_l2y(n) \
		+ X[2] * normal_l3y(n) \
		+ X[3] * normal_l1y(n);\
	(result)->ptr[3] = \
		  X[0] * normal_l1z(n) \
		+ X[1] * normal_l2z(n) \
		+ X[2] * normal_l3z(n) \
		+ X[3] * normal_l1z(n);\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

// by default in hydro/eqn/eqn.lua, <?=fluxFromCons?> is defined by <?=eigen_fluxTransform?>
// but since eig is empty, we can define <?=eigen_fluxTransform?> with <?=fluxFromCons?>
#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X_,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
)\
	<?=fluxFromCons?>(resultFlux, solver, X_, cell, n)

// TODO only if metric ~= ident and coord ~= cartesian 

<? if true 
--if not (
--	solver.coord.vectorComponent == "cartesian" 
--	or require "hydro.coord.cartesian":isa(solver.coord) 
--) 
then ?>

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=eqn_common?> <?=coord_raise?> <?=coord_connHol_trace23?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	real3 const pt = cellBuf[index].pos;
	
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	//TODO make use of this
	//real const c = solver->wavespeed / unit_m_per_s;

	real const alpha = metric_alpha(pt);
	real const K = metric_K(pt);
	real3x3 const partial_beta_ul = metric_partial_beta_ul(pt);
	real3 const partial_alpha_l = metric_partial_alpha_l(pt);
	real3 const conn23 = coord_connHol_trace23(pt);
	real const f = metric_f(pt);

	real3 const Psi_u = coord_raise(U->Psi_l, pt);

	deriv->Pi += 
		partial_alpha_l.dot(Psi_u)
		+ alpha * K * U->Pi
		- alpha * U->Psi_l.dot(conn23)
		- alpha * f 						//... for □Φ=f
	;

	deriv->Psi_l = real3_add3(
		deriv->Psi_l,
		real3_real3x3_mul(
			deriv->Psi_l,
			partial_beta_ul
		),
		real3_real_mul(partial_alpha_l, U->Pi)
	);
}

<? end ?>
