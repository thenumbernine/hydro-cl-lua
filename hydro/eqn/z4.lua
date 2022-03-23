--[[
based on whatever my numerical_relativity_codegen z4 is based on, which is probably a Bona-Masso paper,
probably 2004 Bona et al "A symmetry-breaking mechanism for the Z4 general-covariant evolution system"
and then 2005 Bona et al "Geometrically Motivated ..."
and then some 2008 Yano et al
and then my symmath/tests/Z4.lua ... which is trying to do the same thing but grid-background-indepdendent

scratch paper:

g_uv = spacetime 4-metric

-- state vars: gauge
α = lapse
γ_ij = spatial metric
β^i = shift

-- state vars
K_ij = extrinsic curvature

-- state vars - 1st deriv
a_k = log(α)_,k
d_kij = 1/2 γ_ij,k
b^i_j = β^i_,j

-- state vars - z4
Z_i = spatial component of killing vector
Θ = Z^u n_u

-- background metric
^γ_ij = background grid holonomic metric
^d_kij = partial of background grid holonomic metric

-- state vars - deviation from background metric
Δγ_ij = deviation from background grid metric = γ_ij - ^γ_ij
Δd_kij = deviation from partial of backgrond = d_kij - ^d_kij

-- hyperbolic shift vars
B^i = β^i_,t

A_ij = K_ij - 1/3 K γ_ij

d_i = d_ij^j
e_i = d^j_ji

Γ^i_jk = spatial connection = 1/2 γ^im (γ_mj,k + γ_mk,j - γ_jk,m)
Γ^i_jk = d_kj^i + d_jk^i - d^i_jk
Γ^i = Γ^ij_j = d_j^ji + d_j^ji - d^ij_j = 2 e^i - d^i
log(√γ)_,i = 1/2 γ_,i / γ = Γ^j_ji = d_ij^j + d^j_ji - d^j_ji = d_i


gauge ivp:
α_,t - α_,k β^k = = -α^2 Q = -α^2 f K Bona-Masso slicing lapse family
d/dt α = -α^2 f (K - m Θ) ... with Z4 vars ... I think I'm assuming m=2 in some calculations

d/dt T = (∂_t - L_β) T ... where L_β is the Lie-derivative in the shift direction
∂_0 T = T_,0 = (∂_t - ∂_k β^k) T ... which is only the scalar component of the Lie-derivative

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local sym = common.sym


local Z4_2004Bona = class(EinsteinEqn)
Z4_2004Bona.name = 'Z4_2004Bona'

--[[
args:
noZeroRowsInFlux = true by default.
	true = use 13 vars of a_x, d_xij, K_ij
	false = use all 30 or so hyperbolic conservation variables
useShift
	
	useShift = 'none'


	useShift = 'MinimalDistortionElliptic'
	minimal distortion, as an elliptic equation, via Poisson relaxation.
	2008 Alcubierre's book, eqn 4.3.14 and 4.3.15
	Δ_L(β^i) = β^i_;j^;j + 1/3 β^j_;j^;i + R^i_j β^j = 2 (α A^ij)_;j


	useShift = 'MinimalDistortionParabolic'
	minimal distortion, parabolic equation, via evolution.
	eqn 10 of 1996 Balakrishna et al "Coordinate Conditions and their Implementations in 3D Numerical Relativity"
		β_i,t = ε (β_i^;j_j + 1/3 β_j^;j_i + R^j_i β_j - 2 (α A_ij)^;j)
	2008 Alcubierre's book, "parabolic minimal distortion driver"
	eqn 4.3.26:
		β^i_,t = ε (β^i;j_j + 1/3 β_j^;ji + R^i_j β^j - 2 (α A^ij)_;j)
		Notice that Balakrishna's paper is in lower and Alcubierre's is in upper

	useShift = "MinimalDistortionHyperbolic"
	2008 Alcubierre's book, eqn 4.3.27
		β^i_,tt = α^2 ξ (β^i;j_j + 1/3 β_j^;ji + R^i_j β^j - 2 (α A^ij)_;j)
	so the ε is replaced with α^2 ξ , and the comment says something about ξ => ξ α^n, for integer n, so choose n=-2 and you recover the original
	then change to hyperbolic form:
		β^i_,t = B^i
		B^i_,t = α^2 ξ (β^i;j_j + 1/3 β_j^;ji + R^i_j β^j - 2 (α A^ij)_;j)
	then change to "shifting-shift" form:
	2008 Alcubierre's book, eqn 4.3.31 - 4.3.32
		β^i_,t - β^i_,j β^j = B^i
		B^i_,t - B^i_,j β^j = α^2 ξ (β^i;j_j + 1/3 β_j^;ji + R^i_j β^j - 2 (α A^ij)_;j)


	useShift = 'GammaDriver'
	2010 Baumgarte & Shapiro's book, eqn 4.82
		β^i_,t = k (_Γ^i_,t - η _Γ^i)
		k = 3/4
		η = 1 or 1/(2 M) for total mass M
	2008 Alcubierre's book ... which eqns?
	"parabolic" + "non-shifting-shift" form:
		β^i_,t = α^2 ξ ~Γ^i_,t - η ~Γ^i
	2008 Alcubierre later suggets ξ => ξ α^n, for integer n ...
	... so for n=-2 and ξ = k = 3/4 this brings us to the 2010 Baumgarte & Shapiro's book definition
	... with the exception of our exchange of B^i => ~Γ^i ... how about this ...

	useShift = 'HyperbolicGammaDriver'
	2017 Ruchlin, Etienne, Baumgarte - "SENR-NRPy- Numerical Relativity in Singular Curvilinear Coordinate Systems"
		eqn 14.a: β^i_,t = B^i
		eqn 14.b: B^i_,t - B^i_,j β^j = 3/4 (_Λ^i_,t - _Λ^i_,j β^j) - η B^i
		eqn 11.e: _Λ^i_,t - _Λ^i_,j β^j + _Λ^j β^i_,j = γ^jk ^D_j ^D_k β^i + 2/3 ΔΓ^i (_D_j β^j) + 1/3 _D^i _D_j β^j - 2 _A^ij (α_,j - 6 φ_,j) + 2 _A^jk ΔΓ^i_jk - 4/3 α _γ^ij K_,j
	2010 Baumgarte & Shapiro's book, eqn 4.83
		(non-shifting-shift)
		β^i_,t = k B^i ... k = 3/4
		B^i_,t = _Γ^i_,t - η B^i
		... so when you go from parabolic to hyperbolic, how come you replace the -η _Γ^i with a -η B^i
		... maybe because in the parabolic form, the 2nd term represented the coefficient converging β^i into _Γ^i
	2008 Alcubierre's book, eqn 4.3.33 4.3.34
	"hyperbolic" + "shifting-shift" form:
		β^i_,t - β^i_,j β^j = B^i
		B^i_,t - B^i_,j β^j = α^2 ξ (~Γ^i_,t - ~Γ^i_,j β^j) - η B^i

	useShift = 'HarmonicShiftParabolic'
	2008 Yano et al "Flux-Vector-Splitting..." eqn 16-17
	2005 Bona et al "Geometrically Motivated..." says "to convert the minimal distortion elliptic equations into time-dependent parabolic equations by means of the Hamilton-Jacobi method"
	2005 Bona mentions a few, but 2008 Yano picks the first one from the 2005 Bona paper.
	2005 Bona et al, section B.1, "to convert the minimal distortion elliptic equations into time-dependent parabolic equations by means of the Hamilton-Jacobi method"
	2008 Yano eqn 16-17:
		β^i_,t = -α Q^i
		β^i_,t = β^k β^i_,k + α^2 γ^ki (γ^jm γ_jk,m - 1/2 γ_,k / γ - α_,k / α))
		β^i_,t = β^k β^i_,k + α^2 (2 e^i - d^i - a^i)
		β^i_,t = β^i_,k β^k + α^2 (Γ^i - a^i)
		β^i_,0 = α^2 (Γ^i - a^i)
	2008 Alcubierre 4.3.37
		β^i_,t = β^j β^i_,j - α α^,i + α^2 Γ^i + β^i / α (α_,t - β^j α_,j + α^2 K)
		β^i_,0 = α^2 (Γ^i - a^i) + β^i / α (α_,t - β^j α_,j + α^2 K)
		β^i_,0 = α^2 (Γ^i - a^i) + β^i / α (-α^2 Q + α^2 K)
		β^i_,0 = α^2 (Γ^i - a^i) + α β^i (K - Q)
	... which for Q = f K and f=1 and Z4's m=0 or Θ=0 ...
		β^i_,0 = α^2 (Γ^i - a^i)
	... and this is identical to the HarmonicShiftParabolic shift
	... so back to the first def, but substituting the Z4 lapse ivp:
		β^i_,t = β^j β^i_,j + α^2 (Γ^i - a^i) + β^i / α (-α^2 f (K - m Θ) + α^2 K)
		β^i_,t = β^j β^i_,j + α^2 (Γ^i - a^i) + α β^i ((1 - f) K + f m Θ)
	... but wait, if the Alcubierre 4.3.37 eqns shouldve underwent K => K - m Θ for the Z4 impl (TODO determine if so)
	... then it looks like this:
		β^i_,t = β^j β^i_,j + α^2 (Γ^i - a^i) + β^i (K - m Θ) (α - α f)
	2008 Alcubierre's book, eqn 4.3.39
		σ^i = β^i / α
		h = h(α) = and to recover some slicing choose h = f
		σ^i_,t = α σ^i_,j σ^j - α^,i + α h (σ^i K + Γ^i)
		(β^i / α)_,t = α (β^i / α)_,j β^j / α - α_,j γ^ij + α f (β^i / α K + Γ^i)
		β^i_,t α^-1 - β^i α^-2 α_,t = α (β^i_,j α^-1 - β^i α^-2 α_,j) β^j / α - α_,j γ^ij + α f (β^i / α K + Γ^i)
		β^i_,t = β^i_,j β^j - α α_,j γ^ij + f (α β^i K + α^2 Γ^i) + β^i (α_,t - α_,j β^j) / α
		β^i_,t = β^i_,j β^j + α^2 (f Γ^i - a^i) + α f β^i K + β^i (α_,t - α_,j β^j) / α
	... then insert α_,t - α_,j β^j = -α^2 Q
		β^i_,t = β^i_,j β^j + α^2 (f Γ^i - a^i) + α β^i (f K - Q)
	... so this looks a lot like Alcubierre 4.3.37 except with some extra 'f's laying around
	... now insert Z4 Q = f (K - m Θ)
		β^i_,t = β^i_,j β^j + α^2 (f Γ^i - a^i) + α β^i (f K - f (K - m Θ))
		β^i_,t = β^i_,j β^j + α^2 (f Γ^i - a^i) + α β^i f m Θ
	... and if the first formulation involved a K => (K - m Θ) then it would look like ...
		β^i_,t = β^i_,j β^j + α^2 (f Γ^i - a^i)
	... which is just like the 2005 Bona et al / 2008 Yano et al papers' shift,
	    but with a 'f' inserted next to the Γ^i


so the breakdown of shift conditions ...
1) Elliptic, Laplacian(curved space) based for minimal-distortion
2) IVP ...
	2.a) minimal-distortion
	2.b) gamma-driver
	2.c) harmonic
	- option for non-shifting-shift vs shifting-shift (_,t => _,0 = _,t - _,k β^k)
	- option for parabolic (β^i_,t = ...) vs hyperbolic (β^i_,tt = ...)
	


	useShift = 'LagrangianCoordinates'
	--[=[
	Step backwards along shift vector and advect the state
	Idk how accurate this is ...
	Hmm, even if I implement the Lie derivative as Lagrangian coordinate advection
	I'll still be responsible for setting some β^i_,t gauge
	so for the L_β Lie derivative, we have some options:
	1) none (for no-shift)
	2) finite difference
	3) finite volume / absorb into the eigensystem
	4) Lagrangian coordinates
	and this should be a separate variable, separate of the shift gauge

	so
	one variable for what β^i_,t is
	another variable for how to
	--]=]
--]]
function Z4_2004Bona:init(args)
	local solver = args.solver

	local fluxVars = table{
		{name='a_l', type='real3'},
		{name='dDelta_lll', type='_3sym3'},
		{name='K_ll', type='sym3'},
		{name='Theta', type='real'},
		{name='Z_l', type='real3'},
	}

	self.consVars = table{
		{name='alpha', type='real'},
		{name='gammaDelta_ll', type='sym3'},
	}:append(fluxVars)


	--[[
	how are shift conditions impmlemented?
	options for determining β^i:
	1) by solving a constraint equation (minimal distortion elliptic solves it with a numerical Poisson solver)
	2) by solving an initial value problem

	Once β^i is determined, how are the variables iterated?
	This question is split into a) and b):
	a) How are the source-only variables iterated wrt β^i?
	options:
	1) put the Lie derivative terms into the source side of these variables
	2) give them -β^i eigenvalues?  this is the equivalent of rewriting the hyperbolic vars associated with these (a_i,d_kij) back into first-derivative 0th order vars (α, γ_ij)

	b) How are the flux variables iterated wrt β^i?
	options:
	1) this would be solved by offsetting the eigenvalues
		with the offset eigenvalues,
		the hyperbolic state vars' contributions get incorporated into the flux,
		but there are still some shift-based terms that end up in the source ...
		... so the shift is split between the flux and source ...
	2) ... and adding a few terms to the source
	--]]

	self.useShift = args.useShift or 'none'
	
	-- set to false to disable ρ, S_i, S^ij
	self.useStressEnergyTerms = true

	if self.useShift ~= 'none' then
		self.consVars:insert{name='beta_u', type='real3'}

		if self.useShift == 'MinimalDistortionElliptic' then
			self.consVars:insert{name='betaLap_u', type='real3'}
		elseif self.useShift == 'HarmonicShiftParabolic'
		or self.useShift == 'MinimalDistortionParabolic'
		then
			self.consVars:insert{name='b_ul', type='real3x3'}
		elseif self.useShift == 'HyperbolicGammaDriver' then
			self.consVars:insert{name='B_u', type='real3'}
			self.consVars:insert{name='b_ul', type='real3x3'}
		end
	end


	--[[
	solve a smaller eigendecomposition that doesn't include the rows of variables whose d/dt is zero.
	kind of like how ideal MHD is an 8-var system, but the flux jacobian solved is a 7x7 because Bx,t = 0.
	TODO make this a ctor arg - so solvers can run in parallel with and without this
	...or why don't I just scrap the old code, because this runs a lot faster.
	--]]
	self.noZeroRowsInFlux = true
	--self.noZeroRowsInFlux = false
	if args.noZeroRowsInFlux ~= nil then
		self.noZeroRowsInFlux = args.noZeroRowsInFlux
	end

	-- NOTE this doesn't work when using shift ... because then all the eigenvalues are -β^i, so none of them are zero (except the source-only α, β^i, γ_ij)
	-- with the exception of the lagrangian shift.  if we split that operator out then we can first solve the non-shifted system, then second advect it by the shift vector ...
	--if self.useShift ~= 'none' then
	--	self.noZeroRowsInFlux = false
	--end

	if not self.noZeroRowsInFlux then
		-- skip α and γ_ij
		self.numWaves = Struct.countScalars{vars=fluxVars}
		assert(self.numWaves == 31)
	else
		-- Z4 has a_x, d_xij, K_ij, Θ, Z_i ...
		-- which is 17 waves
		self.numWaves = 17
	end

	-- only count int vars after the shifts have been added
	self:cdefAllVarTypes(solver, self.consVars)	-- have to call before countScalars in eqn:init
	self.numIntStates = Struct.countScalars{vars=self.consVars}

	-- now add in the source terms (if you want them)
	if self.useStressEnergyTerms then
		self.consVars:append{
			--stress-energy variables:
			{name='rho', type='real'},			--1: n_a n_b T^ab
			{name='S_u', type='real3'},			--3: -γ^ij n_a T_aj
			{name='S_ll', type='sym3'},			--6: γ_i^c γ_j^d T_cd
		}
	end
	self.consVars:append{
		--constraints:
		{name='H', type='real'},				--1
		{name='M_u', type='real3'},				--3
	}

	self.eigenVars = table{
		{name='alpha', type='real'},
		{name='alpha_sqrt_f', type='real'},
		{name='gamma_ll', type='sym3'},
		{name='gamma_uu', type='sym3'},
		-- sqrt(n_i n_j γ^ij) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
		{name='sqrt_gammaUnn', type='real'},
	}

	-- hmm, only certain shift methods actually use β^i ...
	if self.useShift ~= 'none' then
		self.eigenVars:insert{name='beta_u', type='real3'}
	end


	-- build stuff around consVars
	Z4_2004Bona.super.init(self, args)


	if self.useShift == 'MinimalDistortionElliptic' then
		local MinimalDistortionEllipticShift = require 'hydro.op.gr-shift-mde'
		self.solver.ops:insert(MinimalDistortionEllipticShift{solver=self.solver})
	elseif self.useShift == 'LagrangianCoordinates' then
		local LagrangianCoordinateShift = require 'hydro.op.gr-shift-lc'
		self.solver.ops:insert(LagrangianCoordinateShift{solver=self.solver})
	end
end

function Z4_2004Bona:getSymbolFields()
	return Z4_2004Bona.super.getSymbolFields(self):append{
		'initDeriv_numeric_and_useBSSNVars',
		'calc_gammaHat_ll',		-- calc ^γ_ij = grid metric (hol.) derivative ... also in common with bssnok-fd-*
		'calc_dHat_lll',		-- calc ^d_kij = 1/2 ^γ_ij,k
		'calc_d_lll',			-- from U->dDelta_ijk and ^d(pt)_ijk
		'calcFromGrad_a_l',
		'calcFromGrad_d_lll',	-- finite difference from grid
		'calcFromGrad_b_ul',
	}
end

function Z4_2004Bona:createInitState()
	Z4_2004Bona.super.createInitState(self)
	self:addGuiVars{
		
		-- from 2004 Bona et al, "A symmetry breaking..." eqn A.20
		{name='m', value=2},
	
		-- convergence between finite-difference of log(α)_,i and a_i
		{name='a_convCoeff', value=0},
		
		-- convergence between finite-difference of 1/2 γ_ij,k and d_kij
		{name='d_convCoeff', value=0},

		{name='dissipationCoeff', value=cmdline.dissipationCoeff or 0},
	
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=-math.huge},


		-- TODO add shift option
		-- but that means moving the consVars construction to the :init()
		-- so until then, just add shift vars no matter what ...
		-- convergence between finite-difference of 1/2 γ_ij,k and d_kij
		{name='b_convCoeff', value=0},
	}
end

-- don't use default
function Z4_2004Bona:initCodeModule_calcDTCell() end
function Z4_2004Bona:initCodeModule_fluxFromCons() end

Z4_2004Bona.solverCodeFile = 'hydro/eqn/z4.cl'

Z4_2004Bona.predefinedDisplayVars = {
	'U alpha',
--[[ for x dir only
	'U gamma_ll x x',
	'U d_lll_x x x',
	'U K_ll x x',
	'U Theta',
	'U Z_l x',
--]]
-- [[ for all dirs
	'U gamma_ll norm',
	'U a_l mag',
	-- TODO a 3sym3 Frobenius norm - use for 'd'?
	'U d_lll x norm',
	'U d_lll y norm',
	'U d_lll z norm',
	'U K_ll norm',
	'U Theta',
	'U Z_l mag',
--]]
	'U beta_u mag',
	'U H',
	'U M_u mag',
	'U volume',
	'U f*alpha',
}

function Z4_2004Bona:getDisplayVars()
	local vars = Z4_2004Bona.super.getDisplayVars(self)

	vars
	:append(xNames:mapi(function(xi,i)
		return {
			name = 'd_lll '..xi,
			type = 'sym3',
			code = self:template([[
//// MODULE_DEPENDS: <?=calc_d_lll?>
value.vsym3 = <?=calc_d_lll?>(U, cell->pos).<?=xi?>;
]], {xi=xi})
		}
	end))
	:append{
		{
			name = 'gamma_ll',
			type = 'sym3',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gamma_ll?>
value.vsym3 = <?=calc_gamma_ll?>(U, cell->pos);
]],
		},
		{	-- this is spacetime volume
			-- spatial volume is just "U gamma_ll det"
			name = 'volume',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gamma_ll?>
sym3 const gamma_ll = <?=calc_gamma_ll?>(U, cell->pos);
value.vreal = U->alpha * sqrt(sym3_det(gamma_ll));
]],
		},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='f*alpha', code='value.vreal = calc_f_alpha(U->alpha);'},
		{name='f*alpha^2', code='value.vreal = calc_f_alphaSq(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		{name='alpha^2*df/dalpha', code='value.vreal = calc_alphaSq_dalpha_f(U->alpha);'},
		
		-- is expansion really just -K?  aren't there shift terms too?
		{
			name = 'expansion',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gamma_uu?>
sym3 const gamma_uu = <?=calc_gamma_uu?>(U, cell->pos);
value.vreal = -sym3_dot(U->K_ll, gamma_uu);
]],
		},
	}:append{
--[=[
--[[
Alcubierre 3.1.1
Baumgarte & Shapiro 2.90
H = R + K^2 - K^ij K_ij - 16 π ρ
for ρ = n_a n_b T^ab (B&S eqn 2.89)
and n_a = -α t_,a (B&S eqns 2.19, 2.22, 2.24)

momentum constraints
--]]
		{H = [[
	.5 *
]]		},
--]=]
	
		-- Γ^i = Γ^ij_j = d_j^ji + d_j^ji - d^ij_j = 2 e^i - d^i
		{
			name = 'Gamma^i',
			type = 'real3',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_d_lll?>
_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);
//// MODULE_DEPENDS: <?=calc_gamma_uu?>
sym3 const gamma_uu = <?=calc_gamma_uu?>(U, cell->pos);
value.vreal3 = real3_sub(
	real3_real_mul(
		sym3_3sym3_dot12(gamma_uu, d_lll),	//e_l
		2.
	),
	_3sym3_sym3_dot23(d_lll, gamma_uu)		//d_l
);
]],
		},
		-- log(√γ)_,i = 1/2 γ_,i / γ = Γ^j_ji = d_ij^j + d^j_ji - d^j_ji = d_i
		{
			name = 'log(sqrt(gamma))_,i',
			type = 'real3',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_d_lll?>
_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);
//// MODULE_DEPENDS: <?=calc_gamma_uu?>
sym3 const gamma_uu = <?=calc_gamma_uu?>(U, cell->pos);
value.vreal3 = _3sym3_sym3_dot23(d_lll, gamma_uu);		//d_l
]],
		},
	}

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
	vars:insert{
		name = 'gravity',
		type = 'real3',
		code = self:template[[
//// MODULE_DEPENDS: <?=calc_gamma_uu?>
sym3 const gamma_uu = <?=calc_gamma_uu?>(U, cell->pos);
value.vreal3 = real3_real_mul(sym3_real3_mul(gamma_uu, U->a_l), -U->alpha * U->alpha);
]],
	}

	-- a_i vs. log(α)_,i
	vars:insert{
		name = 'alpha vs a_i',
		type = 'real3',
		code = self:template[[
if (<?=OOB?>(1,1)) {
	value.vreal3 = real3_zero;
} else {
//// MODULE_DEPENDS: <?=calcFromGrad_a_l?>
	real3 const target_a_l = <?=calcFromGrad_a_l?>(solver, U);
	value.vreal3 = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = fabs(target_a_l.<?=xi?> - U->a_l.<?=xi?>),
<? end
?>	};
}
]]}

	-- d_kij vs. 1/2 γ_ij,k
	for i,xi in ipairs(xNames) do
		vars:insert{
			name = 'gamma_ij vs d_'..xi..'ij',
			type = 'sym3',
			code = self:template([[
if (<?=OOB?>(1,1)) {
	value.vsym3 = sym3_zero;
} else {
//// MODULE_DEPENDS: <?=calcFromGrad_d_lll?>
	_3sym3 const target_d_lll = <?=calcFromGrad_d_lll?>(solver, U, cell);
//// MODULE_DEPENDS: <?=calc_d_lll?>
	_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);
	value.vsym3 = (sym3){
<? for jk,xjk in ipairs(symNames) do
?>		.<?=xjk?> = fabs(target_d_lll.<?=xi?>.<?=xjk?> - d_lll.<?=xi?>.<?=xjk?>),
<? end
?>	};
}
]], 		{
				i = i,
				xi = xi,
			})
		}
	end

	-- b^i_j vs. β^i_,j
	if self.consStruct.vars:find(nil, function(var) return var.name == 'b_ul' end) then
		for i,xi in ipairs(xNames) do
			vars:insert{
				name = 'beta^j_,'..xi..' vs b^j_'..xi,
				type = 'real3',
				code = self:template([[
if (<?=OOB?>(1,1)) {
	value.vreal3 = real3_zero;
} else {
//// MODULE_DEPENDS: <?=calcFromGrad_b_ul?>
	real3x3 const target_b_ul = <?=calcFromGrad_b_ul?>(solver, U);
	value.vreal3 = (real3){
<? for j,xj in ipairs(xNames) do
?>		.<?=xj?> = fabs(target_b_ul.<?=xi?>.<?=xj?> - U->b_ul.<?=xi?>.<?=xj?>),
<? end
?>	};
}
]],				{
					i = i,
					xi = xi,
				}),
			}
		end
	end

	return vars
end

function Z4_2004Bona:eigenWaveCodePrefix(args)
	return self:template([[
real const eig_lambdaLight = (<?=eig?>)->sqrt_gammaUnn * (<?=eig?>)->alpha;
real const eig_lambdaGauge = (<?=eig?>)->sqrt_gammaUnn * (<?=eig?>)->alpha_sqrt_f;
<? if eqn.useShift ~= 'none' then ?>
real const betaUi = normal_vecDotN1(<?=n?>, (<?=eig?>)->beta_u);
<? end ?>
]], args)
end

function Z4_2004Bona:eigenWaveCode(args)
	-- TODO find out if -- if we use the lagrangian coordinate shift operation -- do we still need to offset the eigenvalues by -β^i?
	--local shiftingLambdas = self.useShift ~= 'none'
	--and self.useShift ~= 'LagrangianCoordinates'

	local betaUi = self.useShift ~= 'none' and 'betaUi' or '0'

	if not self.noZeroRowsInFlux then
		if args.waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif args.waveIndex >= 1 and args.waveIndex <= 6 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif args.waveIndex >= 7 and args.waveIndex <= 23 then
			return '-'..betaUi
		elseif args.waveIndex >= 24 and args.waveIndex <= 29 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif args.waveIndex == 30 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	else	-- noZeroRowsInFlux
		-- noZeroRowsInFlux implies useShift == 'none'
		if args.waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif args.waveIndex >= 1 and args.waveIndex <= 6 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif args.waveIndex >= 7 and args.waveIndex <= 9 then
			return '-'..betaUi
		elseif args.waveIndex >= 10 and args.waveIndex <= 15 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif args.waveIndex == 16 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	end
	error'got a bad waveIndex'
end

Z4_2004Bona.consWaveCode = Z4_2004Bona.eigenWaveCode
-- Z4_2004Bona.eigenWaveCodeMinMax -- use default

function Z4_2004Bona:consWaveCodePrefix(args)
	return self:template([[
//// MODULE_DEPENDS: <?=calc_gamma_uu?>
sym3 const gamma_uu = <?=calc_gamma_uu?>(<?=U?>, <?=pt?>);

<? if solver.coord.vectorComponent == 'cartesian' then ?>
real3 const n_l = normal_l1(<?=n?>);
real const gammaUnn = real3_weightedLenSq(n_l, gamma_uu);
<? else ?>
real gammaUnn = 0./0.;
if ((<?=n?>).side == 0) {
	gammaUnn = gamma_uu.xx;
} else if ((<?=n?>).side == 1) {
	gammaUnn = gamma_uu.yy;
} else if ((<?=n?>).side == 2) {
	gammaUnn = gamma_uu.zz;
}
<? end ?>

real const sqrt_gammaUnn = sqrt(gammaUnn);
real const eig_lambdaLight = sqrt_gammaUnn * (<?=U?>)->alpha;
real const alpha_sqrt_f = sqrt(calc_f_alphaSq((<?=U?>)->alpha));
real const eig_lambdaGauge = sqrt_gammaUnn * alpha_sqrt_f;
]], args)
end

-- Z4_2004Bona.consWaveCodeMinMax uses default
-- Z4_2004Bona.consWaveCodeMinMaxAllSidesPrefix uses default
-- Z4_2004Bona.consWaveCodeMinMaxAllSides uses default

return Z4_2004Bona
