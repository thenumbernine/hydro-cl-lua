<? 
local common = require 'common'
local sym = common.sym
local xNames = common.xNames

local solver = eqn.solver
local scalar = eqn.scalar
local vec3 = eqn.vec3
?>

//TODO I have the templated types so I can mix solver code.  In such a case, these typedefs shouldn't go in the getCommonCode.
typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;
//typedef <?=scalar?> scalar;
//typedef <?=vec3?> vec3;


<? if getCommonCode then ?>


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

real metric_alpha(real3 pt) { 
	return <?=eqn:compile(eqn.metric.alpha)?>; 
}

real3 metric_beta_u(real3 pt) { 
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = <?=eqn:compile(eqn.metric.beta_u[i])?>,
<? end
?>	};
}

real metric_K(real3 pt) { 
	return <?=eqn:compile(eqn.metric.K)?>;
}

real metric_f(real3 pt) { 
	return <?=eqn:compile(eqn.metric.f)?>;
}

real3 metric_partial_alpha_l(real3 pt) { 
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = <?=eqn:compile(eqn.metric.alpha:diff(eqn.metric.coords[i])())?>,
<? end
?>	};
}

//partial_beta_ul[i][j] = beta^i_,j
real3x3 metric_partial_beta_ul(real3 pt) { 
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = <?=eqn:compile(eqn.metric.beta_u[i]:diff(eqn.metric.coords[j])())?>,
<?	end
?>		},
<?	end
?>	};
}


<? else -- getCommonCode ?> 


// What's the difference between eigen_fluxTransform and fluxFromCons?
// The difference is that the flux matrix of this is based on 'eig', which is derived from U's ... especially UL & UR in the case of the Roe solver
// whereas that of fluxFromCons is based purely on 'U'.
// Since eqn/wave has no eigen_t info derived from U, the two functions are identical.
cons_t fluxFromCons(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	normalInfo_t n
) {
	real alpha = metric_alpha(x);
	real beta_n = normalInfo_vecDotN1(n, metric_beta_u(x));
	
	real3 nL = normalInfo_l1(n);
	real3 nU = normalInfo_u1(n);
	
	cons_t F;
	//F^Pi = -c (Pi beta_n + alpha Psi_i n^i)
	F.Pi = <?=scalar?>_real_mul(
		//Pi beta_n + alpha Psi_i n^i
		real_<?=scalar?>_add(
		
			//Pi beta_n:
			<?=scalar?>_real_mul(U.Pi, beta_n), 
			
			//alpha Psi_i n^i:
			<?=scalar?>_real_mul(
				//Psi_i n^i:
				<?=vec3?>_real3_dot(
					U.Psi_l, 
					nU
				), 
				alpha
			)
		), 
		-solver->wavespeed
	);
	
	//F^{Psi_j} = -c (alpha Pi n_j + Psi_j beta_n)
	F.Psi_l = <?=vec3?>_real_mul(
		<?=vec3?>_add(
			//Psi_j beta_n:
			<?=vec3?>_real_mul(
				U.Psi_l,
				beta_n
			),
			
			//alpha Pi n_j:
			<?=vec3?>_<?=scalar?>_mul(
				//n:
				<?=vec3?>_from_real3(nL),
				
				//alpha Pi:
				<?=scalar?>_real_mul(
					U.Pi, 
					alpha
				)
			)
		),
		//-c
		-solver->wavespeed
	);
	
	return F;
}


// used by PLM
range_t calcCellMinMaxEigenvalues(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x,
	normalInfo_t n
) {
	real alpha_nLen = metric_alpha(x) * normalInfo_len(n);
	real beta_n = normalInfo_vecDotN1(n, metric_beta_u(x));
	return (range_t){
		.min = solver->wavespeed * (-beta_n - alpha_nLen),
		.max = solver->wavespeed * (-beta_n + alpha_nLen),
	};
}


#define eigen_forInterface(solver, UL, UR, x, n) ((eigen_t){})
#define eigen_forCell(solver, UL, UR, x, n) ((eigen_t){})


waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x,
	normalInfo_t n
) { 
	waves_t Y;

	real nLen = normalInfo_len(n);
	real invDenom = 1. / (nLen * nLen);

	Y.ptr[0] = .5 * invDenom * (
		X.ptr[0] * nLen
		+ X.ptr[1] * normalInfo_u1x(n)
		+ X.ptr[2] * normalInfo_u1y(n)
		+ X.ptr[3] * normalInfo_u1z(n)
	);
	Y.ptr[1] = invDenom * (
		X.ptr[1] * normalInfo_u2x(n)
		+ X.ptr[2] * normalInfo_u2y(n)
		+ X.ptr[3] * normalInfo_u2z(n)
	);
	Y.ptr[2] = invDenom * (
		X.ptr[1] * normalInfo_u3x(n)
		+ X.ptr[2] * normalInfo_u3y(n)
		+ X.ptr[3] * normalInfo_u3z(n)
	);
	Y.ptr[3] = .5 * invDenom * (
		-X.ptr[0] * nLen
		+ X.ptr[1] * normalInfo_u1x(n)
		+ X.ptr[2] * normalInfo_u1y(n)
		+ X.ptr[3] * normalInfo_u1z(n)
	);

	return Y;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x,
	normalInfo_t n
) {
	<?=scalar?>* X = (<?=scalar?>*)X_.ptr;
	
	real nLen = normalInfo_len(n);

	cons_t Y;
	Y.ptr[0] = (X[0] - X[3]) * nLen;
	Y.ptr[1] = X[0] * normalInfo_l1x(n) 
			+ X[1] * normalInfo_l2x(n) 
			+ X[2] * normalInfo_l3x(n) 
			+ X[3] * normalInfo_l1x(n);
	Y.ptr[2] = X[0] * normalInfo_l1y(n) 
			+ X[1] * normalInfo_l2y(n) 
			+ X[2] * normalInfo_l3y(n) 
			+ X[3] * normalInfo_l1y(n);
	Y.ptr[3] = X[0] * normalInfo_l1z(n) 
			+ X[1] * normalInfo_l2z(n) 
			+ X[2] * normalInfo_l3z(n) 
			+ X[3] * normalInfo_l1z(n);
	
	return Y;
}

// by default in eqn/eqn.lua, fluxFromCons is defined by eigen_fluxTransform
// but since eig is empty, we can define eigen_fluxTransform with fluxFromCons
#define eigen_fluxTransform(solver, eig, X, x, n) fluxFromCons(solver, X, x, n)


kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
#if 0
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;
	//TODO make use of this
	//real c = solver->wavespeed / unit_m_per_s;

	real alpha = metric_alpha(x);
	real K = metric_K(x);
	real3x3 partial_beta_ul = metric_partial_beta_ul(x);
	real3 partial_alpha_l = metric_partial_alpha_l(x);
	real3 conn23 = coord_conn_trace23(x);
	real f = metric_f(x);

	real3 Psi_u = coord_raise(U->Psi_l, x);

	deriv->Pi += 
		real3_dot(partial_alpha_l, Psi_u)
		+ alpha * K * U->Pi
		- alpha * real3_dot(U->Psi_l, conn23)
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
#endif
}

<? end -- getCommonCode ?> 
