local class = require 'ext.class'

local ViscousFlux = class()

function ViscousFlux:init(args)
	self.solver = assert(args.solver)
end

-- insert extra flux calculations into fvsolver's calcFlux() kernel
function ViscousFlux:addCalcFluxCode()
--[[
from "CFD - An Open Approach" section 14.2

dU^I/dt + dF^Ij/dx^j = 0
dU^I/dt = -dF^Ij/dx^j

i = flux direction
j = state variable vector component
Fv^ρ^i = 0
Fv^m^ij = -τ^ij
Fv^E^i = -τ^ij v_j - q^i

τ^ij = mu ((v^i;j + v^j;i) - 2/3 g^ij v^k_;k)
τ^ij_;j = mu ((v^i;j_j + v^j;i_j) - 2/3 g^ij v^k_;kj)

dU^I/dt + dFv^Ii/dx^i = 0
dU^I/dt = -dFv^Ii/dx^i

dρ/dt = 0
dm^j/dt = ∇_i τ^ij
dE/dt = ∇_i (τ^ij v_j + q^i)

dm^j/dt = ∇_i (mu ((∇^j v^i + ∇^i v^j) - 2/3 g^ij ∇_k v^k))
dm^j/dt = 
		∇_i mu ((∇^j v^i + ∇^i v^j) - 2/3 g^ij ∇_k v^k)
		+ mu (
			(
				∇^j (∇_i v^i)
				
				+ (∇_i ∇^i) v^j
			)
			- 2/3 g^ij ∇_i (∇_k v^k)
		)

dE/dt = ∇_i ((mu ((∇^j v^i + ∇^i v^j) - 2/3 g^ij ∇_k v^k)) v_j + q^i)

... do I really want to modify the flux vector?
or should I just put this in the source term?
flux vector here
source term somewhere else


per-component in cartesian metric
for rhs calcs:

Fv^m^ij = -τ^ij
-Fv^m^ij_,j = τ^ij_,j
-Fv^m^ij_,j = mu ((v_i,jj + v_j,ij) - 2/3 v_j,ji)
-Fv^m^ij_,j = mu (v_i,jj + 1/3 v_j,ij)
-Fv^m^ij_,j = mu (v_i,xx + v_i,yy + v_i,zz + 1/3 v_x,ix + 1/3 v_y,iy + 1/3 v_z,iz)

-Fv^m^xj_,j = mu (4/3 v_x,xx + v_x,yy + v_x,zz + 1/3 v_y,xy + 1/3 v_z,xz)
-Fv^m^yj_,j = mu (4/3 v_y,yy + v_y,xx + v_y,zz + 1/3 v_x,xy + 1/3 v_z,yz)
-Fv^m^zj_,j = mu (4/3 v_z,zz + v_z,xx + v_z,yy + 1/3 v_x,xz + 1/3 v_y,yz)

I could abstract my calculation of variables in my partial deriv code, and finite-difference 'v'
or I could factor our momentum and density ...

-Fv^m^xj_,j = mu (4/3 (m_x/rho),xx + (m_x/rho),yy + (m_x/rho),zz + 1/3 (m_y/rho),xy + 1/3 (m_z/rho),xz)
-Fv^m^yj_,j = mu (4/3 (m_y/rho),yy + (m_y/rho),xx + (m_y/rho),zz + 1/3 (m_x/rho),xy + 1/3 (m_z/rho),yz)
-Fv^m^zj_,j = mu (4/3 (m_z/rho),zz + (m_z/rho),xx + (m_z/rho),yy + 1/3 (m_x/rho),xz + 1/3 (m_y/rho),yz)

(m_i/rho),jk = (m_i,j/rho - m_i rho_,j/rho^2),k
(m_i/rho),jk = (
	m_i,jk / rho 
	- m_i,j rho_,k / rho^2
	- m_i,k rho_,j / rho^2
	- m_i rho_,jk / rho^2
	+ 2 m_i rho_,j rho_,k / rho^3
)

-Fv^m^ij_,j = mu (v_i,jj + 1/3 v_j,ij)
-Fv^m^ij_,j = mu ((m_i/rho),jj + 1/3 (m_j/rho),ij)
-Fv^m^ij_,j = mu (
	(
		m_i,jj / rho 
		- m_i,j rho_,j / rho^2
		- m_i,j rho_,j / rho^2
		- m_i rho_,jj / rho^2
		+ 2 m_i rho_,j rho_,j / rho^3
	)
	+ 1/3 (
		m_j,ij / rho 
		- m_j,i rho_,j / rho^2
		- m_j,j rho_,i / rho^2
		- m_j rho_,ij / rho^2
		+ 2 m_j rho_,i rho_,j / rho^3
	)
)
-Fv^m^ij_,j = mu (
	+ m_i,jj / rho + 1/3 m_j,ij / rho 
	- (
		+ 2 m_i,j rho_,j
		+ m_i rho_,jj
		+ 1/3 m_j,i rho_,j
		+ 1/3 m_j,j rho_,i
		+ 1/3 m_j rho_,ij
	) / rho^2
	+ (
		+ 2 m_i rho_,j rho_,j
		+ 2/3 m_j rho_,i rho_,j
	) / rho^3
)

... too many terms


-Fv^E^i = τ^ij v_j + q^i
-Fv^E^i_,i = (τ^ij v_j + q^i)_,i
-Fv^E^i_,i = τ^ij_,i v_j + τ^ij v_j_,i + q^i_,i

--]]
	return [[
	
#error TODO
	//how about the partials not across the flux axis?  time to construct a basis.
	//this works fine for grids, but what to do for unstructured?
	// even for grids, along the flux dir we will be using dx, but along the perp dirs we will use 2 dx
	// and even if I do 2dx along the non-flux axis then we have to recalc centers of faces in neighbors of our current face
	real3x3 dv;
	dv.x = real3_sub(
		real3_real_mul(ppUR->m, 1. / (ppUR->rho * dx)),
		real3_real_mul(ppUL->m, 1. / (ppUL->rho * dx))
	)

	sym3 tau_uu = 
]]
end

return ViscousFlux 
