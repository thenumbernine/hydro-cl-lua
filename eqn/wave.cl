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


// TODO this can be automatically done based on flagging the state variables, whether their indexes are upper or lower
<? for side=0,solver.dim-1 do ?>
cons_t cons_parallelPropagate<?=side?>(cons_t U, real3 x, real dx) {
	U.Psi_l = coord_parallelPropagateL<?=side?>(U.Psi_l, x, dx);
	return U;
}
<? end ?>


<? else -- getCommonCode ?> 



cons_t fluxFromConsForNormal(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	real3 n	//covariant
) {
	real alpha = metric_alpha(x);
	real beta_n = real3_dot(metric_beta_u(x), n);
	
	real3 nU = coord_raise(n, x);
	
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
				<?=vec3?>_from_real3(n),
				
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
range_t calcCellMinMaxEigenvaluesForNormal(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x,
	real3 n
) {
	real3 nU = coord_raise(n, x);
	real nLen = sqrt(real3_dot(n, nU));
	real alpha_nLen = metric_alpha(x) * nLen;
	real beta_n = real3_dot(metric_beta_u(x), n);
	return (range_t){
		.min = solver->wavespeed * (-beta_n - alpha_nLen),
		.max = solver->wavespeed * (-beta_n + alpha_nLen),
	};
}


eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3 n
) {
	return (eigen_t){};
}

eigen_t eigen_forCellForNormal(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x
) {
	return (eigen_t){};
}




waves_t eigen_leftTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x,
	real3 n		// in covariant components
) { 
	waves_t Y;

	real3 nU = coord_raise(n, x);
	real nLenSq = real3_dot(n, nU);
	real nLen = sqrt(nLenSq);
	real invDenom = 1. / nLenSq;

	real3 n2U, n3U;
	getPerpendicularBasis(nU, &n2U, &n3U);

	Y.ptr[0] = .5 * invDenom * (
		X.ptr[0] * nLen
		+ X.ptr[1] * nU.x
		+ X.ptr[2] * nU.y
		+ X.ptr[3] * nU.z
	);
	Y.ptr[1] = invDenom * (
		X.ptr[1] * n2U.x
		+ X.ptr[2] * n2U.y
		+ X.ptr[3] * n2U.z
	);
	Y.ptr[2] = invDenom * (
		X.ptr[1] * n3U.x
		+ X.ptr[2] * n3U.y
		+ X.ptr[3] * n3U.z
	);
	Y.ptr[3] = .5 * invDenom * (
		-X.ptr[0] * nLen
		+ X.ptr[1] * nU.x
		+ X.ptr[2] * nU.y
		+ X.ptr[3] * nU.z
	);

	return Y;
}

cons_t eigen_rightTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x,
	real3 n
) {
	<?=scalar?>* X = (<?=scalar?>*)X_.ptr;
	cons_t Y;
	
	real3 nU = coord_raise(n, x);
	real nLenSq = real3_dot(n, nU);
	real nLen = sqrt(nLenSq);
	real invDenom = 1. / nLenSq;

	real3 n2, n3;
	getPerpendicularBasis(n, &n2, &n3);

	Y.ptr[0] = (X[0] - X[3]) * nLen;
	Y.ptr[1] = X[0] * n.x + X[1] * n2.x + X[2] * n3.x + X[3] * n.x;
	Y.ptr[2] = X[0] * n.y + X[1] * n2.y + X[2] * n3.y + X[3] * n.y;
	Y.ptr[3] = X[0] * n.z + X[1] * n2.z + X[2] * n3.z + X[3] * n.z;
	
	return Y;
}


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


<? for side=0,solver.dim-1 do ?>
#define fluxFromCons_<?=side?>(solver, U, x) 				fluxFromConsForNormal(solver, U, x, normalForSide<?=side?>())
#define calcCellMinMaxEigenvalues_<?=side?>(solver, U, x) 	calcCellMinMaxEigenvaluesForNormal(solver, U, x, normalForSide<?=side?>())
#define eigen_leftTransform_<?=side?>(solver, eig, X, x) 	eigen_leftTransformForNormal(solver, eig, X, x, normalForSide<?=side?>())
#define eigen_rightTransform_<?=side?>(solver, eig, X, x) 	eigen_rightTransformForNormal(solver, eig, X, x, normalForSide<?=side?>())

// What's the difference between eigen_fluxTransform and fluxFromCons?
// The difference is that the flux matrix of this is based on 'eig', which is derived from U's ... especially UL & UR in the case of the Roe solver
// whereas that of fluxFromCons is based purely on 'U'.
// Since eqn/wave has no eigen_t info derived from U, the two functions are identical.
#define eigen_fluxTransform_<?=side?>(solver, eig, X, x) 	fluxFromCons_<?=side?>(solver, X, x)

#define eigen_forCell_<?=side?>(solver, U, x) 				eigen_forCellForNormal(solver, U, x, normalForSide<?=side?>())
<? end ?>


<? end -- getCommonCode ?> 
