local class = require 'ext.class'
local template = require 'template'

-- TODO make this a ctor parameter 
local Poisson = require(
	cmdline.srhdSelfGravPoissonSolver 
	and 'op.poisson_'..cmdline.srhdSelfGravPoissonSolver
	or 'op.poisson_krylov'		-- Krylov
	--or 'op.poisson_jacobi'		-- Jacobi
)


local SRHDSelfGrav = class(Poisson)

SRHDSelfGrav.name = 'SRHDSelfGrav'
SRHDSelfGrav.enableField = 'useGravity'

SRHDSelfGrav.guiVars = {
	{name='gravitationalConstant', value=2},
}

function SRHDSelfGrav:init(args)
	args.verbose = cmdline.selfGravVerbose
	args.linearSolver = cmdline.selfGravLinearSolver or 'conjres'	-- so far works best for selfgrav
	SRHDSelfGrav.super.init(self, args)
	self.solver[self.enableField] = not not self.solver[self.enableField]
end

-- params for op/poisson.cl 
function SRHDSelfGrav:getPoissonDivCode()
	-- because op/poisson.cl assumes it's UBuf, 
	-- we gotta keep the name 'UBuf'
	-- even though it's the primBuf ...
	return template([[
<? local eqn = op.solver.eqn ?>
	global <?=eqn.prim_t?>* prim = &UBuf[index].prim;
	source = 4. * M_PI * (solver->gravitationalConstant / unit_m3_per_kg_s2) * prim->rho;
]], {
		op = self,
	})
end

function SRHDSelfGrav:getPoissonCode()
	return template([[
<?
local solver = op.solver
local eqn = solver.eqn
?>
kernel void calcGravityDeriv<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuffer,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_only_t?>* U = &UBuf[index].cons;
	const global <?=eqn.prim_t?>* prim = &UBuf[index].prim;

	real3 du_dt = real3_zero;
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - solver->stepsize.s<?=side?>;
		int indexR = index + solver->stepsize.s<?=side?>;

		du_dt.s<?=side?> = (
			UBuf[indexR].<?=op.potentialField?> 
			- UBuf[indexL].<?=op.potentialField?>
		) / (2. * cell_dx<?=side?>(x));
	}<? end ?>

	real Phi = UBuf[index].<?=op.potentialField?>;

	//u = W v
	real W = U->D / prim->rho;
	real3 u = real3_real_mul(prim->v, W);

	/*
	u,t = W,t v + W v,t
	W = 1/sqrt(1 - v^2)
	in general relativity, v^2 = v^i v^j gamma_ij
	ADM formalism: gamma_ij = g_ij
	W,t = -1/2 1/(1 - v^2)^(3/2) * (-2 v dot v,t)
	W,t = W^3 (v dot v,t)
	W,t v + W v,t = W^3 v (v dot v,t) + W v,t
	u,t = W (W^2 v (v dot v,t) + v,t)
	u,t = W (u outer u + I) * v,t

	(u outer u + I) (I - alpha u outer u)
	= u outer u - alpha u^2 u outer u + I - alpha u outer u
	= I + u outer u (1 - alpha (1 + u^2))
	is the inverse for 1 - alpha (1 + u^2) = 0
	 i.e. alpha = 1 / (1 + u^2)
	so the inverse of (I + u outer u) is (I - u outer u / (1 + u^2))

	v,t = 1/W (I - u outer u / (1 + u^2)) u,t
	v,t = (u,t - u (u dot u,t) / (1 + u^2)) / W
	
	
	
	in general relativity
	using a scalar metric g_ab = eta_ab - 2 Phi delta_ab so that the g_at = 0
	W_,t = -1/2 1/(1 - v^m v^n gamma_mn)^(3/2) * -(2 v^i_,t v_i + v^i v^j gamma_ij,t)
	W_,t = W^3 (v^i_,t v_i + v^i v^j gamma_ij,t / 2)
	W_,t v^j + W v^j_,t = v^j W^3 (v^i_,t v_i + v^i v^k gamma_ik,t / 2) + W v^j_,t
	u^j_,t = W (W^2 v^j v_i + delta^j_i) v^i_,t + v^j W^3 v^i v^k gamma_ik,t / 2
		v_i = gamma_ij v^j by v_a = g_ab v^b = (g_at v^t + g_aj v^j) 
	u^j_,t - u^i u^j u^k gamma_ik,t / 2 = (u_i u^j + delta_i^j) W v^i_,t 

	delta_i^k = (delta_i^j + u_i u^j) (delta_j^k - u_j u^k / (1 + u^m u_m))
	= delta_i^k + u_i u^k - u_i u^k / (1 + u^m u_m) - (u_i u^k (u^j u_j) / (1 + u^m u_m)
	= delta_i^k + u_i u^k - u_i u^k (1 + u_j u^j) / (1 + u_m u^m)
	= delta_i^k
	
	(u^j_,t - u^i u^j u^l gamma_il,t / 2) (delta_j^k - u_j u^k / (1 + u^m u_m))
		= (delta_j^k - u_j u^k / (1 + u^m u_m)) (u_i u^j + delta_i^j) W v^i_,t 
	v^k_,t = (delta^k_j - u^k u_j / (1 + u^m u_m)) (u^j_,t - u^i u^j u^l gamma_il,t / 2) / W
	v^i_,t = (delta^i_j - u^i u_j / (1 + u^m u_m)) (u^j_,t - u^j u^k u^l gamma_jk,t / 2) / W
		ignore gamma_ij,t for now
	v^i_,t = (delta^i_j - u^i u_j / (1 + u^m u_m)) u^j_,t / W
	v^i_,t = u^i_,t - u^i u_j / (1 + u^m u_m) u^j_,t / W
		using the scalar metric 
	v^i_,t = u^i_,t / W - u^i u^j u^j_,t / (W (1 - 2 Phi) (1 + u^m u_m))

	*/
	real uSq = real3_dot(u, u) / (1. - 2. * Phi);
	real3 dv_dt = real3_add(
		real3_real_mul(du_dt, 1. / W),
		real3_real_mul(u, real3_dot(u, du_dt) / (W * (1. - 2. * Phi) * (1. + uSq)))
	);

	//W_,t = W^3 (v^i_,t v_i + v^i v^j gamma_ij,t / 2)
	//W_,t = W^3 (v^i_,t v^i) / (1 - 2 Phi)
	real dW_dt = W * W * W * real3_dot(prim->v, dv_dt) / (1. - 2. * Phi);
	real h = 1. + solver->heatCapacityRatio * prim->eInt;

	//why am I integrating negative again?
	//why does "Hydrodynamics II" say to integrate negative for the Euler equations?	

	//D = W rho
	//D,t = W,t rho
	deriv->cons.D -= dW_dt * prim->rho;

	//S = rho h W^2 v = rho h W u
	//assuming rho and h are constant ... 
	//S,t = rho h (W,t u + W u,t)
	deriv->cons.S = real3_sub(deriv->cons.S,
		real3_add(
			real3_real_mul(u, prim->rho * h * dW_dt),
			real3_real_mul(du_dt, prim->rho * h * W)
		)
	);
	
	//tau = rho h W^2 - p - rho W
	//tau,t = rho h (2 W W,t) - rho W,t
	//tau,t = rho W,t (2 h W - 1)
	deriv->cons.tau -= prim->rho * dW_dt * (2. * h * W - 1.);
}

]],
	{
		op = self,
	})
end

function SRHDSelfGrav:refreshSolverProgram()
	SRHDSelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	self.calcGravityDerivKernelObj = solver.solverProgramObj:kernel('calcGravityDeriv'..self.name)
	self.calcGravityDerivKernelObj.obj:setArg(0, solver.solverBuf)
	self.calcGravityDerivKernelObj.obj:setArg(2, solver.UBuf)
end

function SRHDSelfGrav:updateGUI()
	SRHDSelfGrav.super.updateGUI(self)
	ig.igPushIDStr'SRHDSelfGrav behavior'
	tooltip.checkboxTable('use gravity', self.solver, self.enableField)
	ig.igPopID()
end

function SRHDSelfGrav:addSource(derivBufObj)
	local solver = self.solver
	if not solver[self.enableField] then return end
		
	self:relax()
	self.calcGravityDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcGravityDerivKernelObj()
end

return SRHDSelfGrav
