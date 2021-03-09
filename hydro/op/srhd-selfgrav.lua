local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local ig = require 'ffi.imgui'
local tooltip = require 'hydro.tooltip'

-- TODO make this a ctor parameter 
local Poisson = require(
	cmdline.srhdSelfGravPoissonSolver 
	and 'hydro.op.poisson_'..cmdline.srhdSelfGravPoissonSolver
	or 'hydro.op.poisson_krylov'		-- Krylov
	--or 'hydro.op.poisson_jacobi'		-- Jacobi
)


local SRHDSelfGrav = class(Poisson)

SRHDSelfGrav.name = 'srhd_selfgrav'
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

function SRHDSelfGrav:getSymbolFields()
	return SRHDSelfGrav.super.getSymbolFields(self):append{
		'copyPotentialToReduce',	-- TODO use this?
		'offsetPotential',			-- TODO this too?
		'calcGravityAccel',
		'calcGravityDeriv',
	}
end

-- params for hydro/op/poisson.cl 
function SRHDSelfGrav:getPoissonDivCode()
	return self.solver.eqn:template([[
	source = 4. * M_PI * U->rho 
		* solver->gravitationalConstant / unit_m3_per_kg_s2;
]], {
		op = self,
	})
end

--[[
This is all using a Cartesian background iirc
TODO make it work for any grid metric


u,t = W,t v + W v,t
W = 1/sqrt(1 - v^2)
in general relativity, v^2 = v^i v^j γ_ij
ADM formalism: γ_ij = g_ij
W,t = -1/2 1/(1 - v^2)^(3/2) * (-2 v dot v,t)
W,t = W^3 (v dot v,t)
W,t v + W v,t = W^3 v (v dot v,t) + W v,t
u,t = W (W^2 v (v dot v,t) + v,t)
u,t = W (u outer u + I) * v,t

(u outer u + I) (I - α u outer u)
= u outer u - α u^2 u outer u + I - α u outer u
= I + u outer u (1 - α (1 + u^2))
is the inverse for 1 - α (1 + u^2) = 0
 i.e. α = 1 / (1 + u^2)
so the inverse of (I + u outer u) is (I - u outer u / (1 + u^2))

v,t = 1/W (I - u outer u / (1 + u^2)) u,t
v,t = (u,t - u (u dot u,t) / (1 + u^2)) / W


in general relativity
using a scalar metric g_ab = η_ab - 2 Φ δ_ab so that the g_at = 0
NOTICE this assumes a Cartesian background.
TODO rewrite Newtonian limit of GR for arbitrary spatial background.

W_,t = -1/2 1/(1 - v^m v^n γ_mn)^(3/2) * -(2 v^i_,t v_i + v^i v^j γ_ij,t)
W_,t = W^3 (v^i_,t v_i + v^i v^j γ_ij,t / 2)
W_,t v^j + W v^j_,t = v^j W^3 (v^i_,t v_i + v^i v^k γ_ik,t / 2) + W v^j_,t
u^j_,t = W (W^2 v^j v_i + δ^j_i) v^i_,t + v^j W^3 v^i v^k γ_ik,t / 2
	v_i = γ_ij v^j by v_a = g_ab v^b = (g_at v^t + g_aj v^j) 
u^j_,t - u^i u^j u^k γ_ik,t / 2 = (u_i u^j + δ_i^j) W v^i_,t 

δ_i^k = (δ_i^j + u_i u^j) (δ_j^k - u_j u^k / (1 + u^m u_m))
= δ_i^k + u_i u^k - u_i u^k / (1 + u^m u_m) - (u_i u^k (u^j u_j) / (1 + u^m u_m)
= δ_i^k + u_i u^k - u_i u^k (1 + u_j u^j) / (1 + u_m u^m)
= δ_i^k

(u^j_,t - u^i u^j u^l γ_il,t / 2) (δ_j^k - u_j u^k / (1 + u^m u_m))
	= (δ_j^k - u_j u^k / (1 + u^m u_m)) (u_i u^j + δ_i^j) W v^i_,t 
v^k_,t = (δ^k_j - u^k u_j / (1 + u^m u_m)) (u^j_,t - u^i u^j u^l γ_il,t / 2) / W
v^i_,t = (δ^i_j - u^i u_j / (1 + u^m u_m)) (u^j_,t - u^j u^k u^l γ_jk,t / 2) / W
	ignore γ_ij,t for now <-> grid unchanging in time
v^i_,t = (δ^i_j - u^i u_j / (1 + u^m u_m)) u^j_,t / W
v^i_,t = u^i_,t - u^i u_j / (1 + u^m u_m) u^j_,t / W
	using the scalar metric 
v^i_,t = u^i_,t / W - u^i u^j u^j_,t / (W (1 - 2 Φ) (1 + u^m u_m))
--]]
function SRHDSelfGrav:getPoissonCode()
	return file['hydro/op/srhd-selfgrav.cl']
end

function SRHDSelfGrav:initCodeModules()
	SRHDSelfGrav.super.initCodeModules(self)

	local solver = self.solver
	solver.solverModulesEnabled[self.symbols.calcGravityDeriv] = true
	--solver.solverModulesEnabled[self.symbols.copyPotentialToReduce] = true
	--solver.solverModulesEnabled[self.symbols.offsetPotential] = true
end

function SRHDSelfGrav:refreshSolverProgram()
	SRHDSelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	self.calcGravityDerivKernelObj = solver.solverProgramObj:kernel(self.symbols.calcGravityDeriv)
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
