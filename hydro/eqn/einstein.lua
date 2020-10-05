--[[
common functions for all Einstein field equation solvers
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'

local common = require 'hydro.common'
local xNames = common.xNames

local EinsteinEquation = class(Equation)

EinsteinEquation.initConds = table():append(
	require 'hydro.init.senr'
):append(
	require 'hydro.init.einstein'
)

function EinsteinEquation:createInitState()
	EinsteinEquation.super.createInitState(self)
	self:addGuiVars{
		{
			type = 'combo',
			name = 'f_eqn',
			options = {
				'2/alpha',							-- 1+log slicing
				'1 + 1/alpha^2', 					-- Alcubierre 10.2.24: "shock avoiding condition" for Toy 1+1 spacetimes 
				'1', 								-- Alcubierre 4.2.50 - harmonic slicing
				'0', '.49', '.5', '1.5', '1.69',
			}
		},
	}
end

function EinsteinEquation:initCodeModules()
	EinsteinEquation.super.initCodeModules(self)
	
	self:initCodeModule_setFlatSpace()

	local solver = self.solver

	-- rescaling from/to diagonalization from the grid metric
	-- used especially by the bssn solvers
	-- but anything that wants to use their initial conditions will also need this.
	solver.modules:add{
		name = 'rescaleFromCoord/rescaleToCoord',
		depends = {
			'coord_dx#',
		},
		code = self:template[[
#if 1

//rescaling, used for bssn finite-difference, but I am tempted to try it with other coordinate systems with singularities
//TODO for the initial conditions do this symbolically instead of numerically

//apply this to lower indexes to convert from coordinate metric to better metric
//apply this to upper indexes to convert from better metric to coordinate metric
real3 real3_rescaleFromCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x / coord_dx0(x),
		.y = v.y / coord_dx1(x),
		.z = v.z / coord_dx2(x),
	};
}
#define real3_rescaleToCoord_U real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_L(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y * coord_dx1(x),
		.z = v.z * coord_dx2(x),
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_L

sym3 sym3_rescaleFromCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleToCoord_UU sym3_rescaleFromCoord_ll

sym3 sym3_rescaleToCoord_LL(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleFromCoord_uu sym3_rescaleToCoord_LL

#else	//debugging -- turning it off

#define real3_rescaleFromCoord_l(a,x) a
#define real3_rescaleToCoord_U(a,x) a
#define real3_rescaleToCoord_L(a,x) a
#define real3_rescaleFromCoord_u(a,x) a
#define sym3_rescaleFromCoord_ll(a,x) a
#define sym3_rescaleToCoord_UU(a,x) a
#define sym3_rescaleToCoord_LL(a,x) a
#define sym3_rescaleFromCoord_uu(a,x) a

#endif
]],
	}

	-- only used by bssnok-fd-num with scalar field
	solver.modules:add{
		name = 'cplx3_rescaleFromCoord/cplx3_rescaleToCoord',
		depends = {
			'cplx3',
			'rescaleFromCoord/rescaleToCoord',
		},
		code = self:template[[
cplx3 cplx3_rescaleFromCoord_l(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		real3_rescaleFromCoord_l(cplx3_re(v), x),
		real3_rescaleFromCoord_l(cplx3_im(v), x));
}
]],
	}

	solver.modules:add{
		name = '_3sym3_rescaleFromCoord/_3sym3_rescaleToCoord',
		depends = {
			'_3sym3',
			--'rescaleFromCoord/rescaleToCoord',	-- I could use this for sub-member rescaling
			'coord_dx#',							-- ... but I just did it manually
		},
		code = self:template[[
#if 1

_3sym3 _3sym3_rescaleFromCoord_lll(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleToCoord_UUU _3sym3_rescaleFromCoord_lll

_3sym3 _3sym3_rescaleToCoord_LLL(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleFromCoord_uuu _3sym3_rescaleToCoord_LLL

_3sym3 _3sym3_rescaleFromCoord_ull(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * coord_dx<?=i-1?>(x) / (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}

_3sym3 _3sym3_rescaleToCoord_ULL(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)) / coord_dx<?=i-1?>(x),
<?	end
?>		},
<? end
?>	};
}

#else

#define _3sym3_rescaleFromCoord_lll(a,x) a
#define _3sym3_rescaleToCoord_UUU(a,x) a
#define _3sym3_rescaleToCoord_LLL(a,x) a
#define _3sym3_rescaleFromCoord_uuu(a,x) a

#endif
]]
	}

	-- only used by bssnok-fd-sym
	solver.modules:add{
		name = 'sym3sym3_rescaleFromCoord/sym3sym3_rescaleToCoord',
		depends = {
			'sym3sym3',
			--'rescaleFromCoord/rescaleToCoord',	-- I could use this for sub-member rescaling
			'coord_dx#',							-- ... but I just did it manually
		},
		code = self:template[[
#if 1

sym3sym3 sym3sym3_rescaleFromCoord_llll(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleToCoord_UUUU sym3sym3_rescaleFromCoord_llll

sym3sym3 sym3sym3_rescaleToCoord_LLLL(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleFromCoord_uuuu sym3sym3_rescaleToCoord_LLLL

#else

#define sym3sym3_rescaleFromCoord_llll(a,x) a
#define sym3sym3_rescaleToCoord_UUUU(a,x) a
#define sym3sym3_rescaleToCoord_LLLL(a,x) a
#define sym3sym3_rescaleFromCoord_uuuu (a,x) a

#endif
]],
	}
end

function EinsteinEquation:getModuleDependsSolver()
	return {
		'initCond.codeprefix',	-- eigen_forInterface & others uss calc_f
	}
end

-- add an option for fixed Minkowsky boundary spacetime
-- TODO now there is already a BoundaryFixed in hydro/solver/gridsolver, but no easy way to parameterize how to set what fixed values it is
function EinsteinEquation:createBoundaryOptions()
	local Boundary = self.solver.Boudary
	local BoundaryFixed = class(Boundary)
	BoundaryFixed.name = 'fixed'
	function BoundaryFixed:getCode(args)
		local lines = table()
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		for _,j in ipairs{'j', gridSizeSide..'-numGhost+j'} do
			local index = args.indexv(j)
			local U = 'buf[INDEX('..index..')]'
			lines:insert(self:template([[
	setFlatSpace(solver, &<?=U?>, cell_x((int4)(<?=index?>, 0)));
]], 		{
				U = U,
				index = index,
			}))
		end
		return lines:concat'\n'
	end
	
	self.solver:addBoundaryOption(BoundaryFixed)
end

function EinsteinEquation:getModuleDependsSolver() 
	return {
		-- for the addDisplayComponents 
		'calc_gamma_ll',
		'calc_gamma_uu',
	}
end

-- any modules this code needs, add to solver module dependencies
function EinsteinEquation:createDisplayComponents()
	local solver = self.solver
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma_ij',
		code = [[
const global <?=eqn.cons_t?>* U = buf + index;
sym3 gamma_ll = calc_gamma_ll(U, x);
value->vreal = real3_weightedLen(value->vreal3, gamma_ll);
]],
	})
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma^ij',
		code = [[
const global <?=eqn.cons_t?>* U = buf + index;
sym3 gamma_uu = calc_gamma_uu(U, x);
value->vreal = real3_weightedLen(value->vreal3, gamma_uu);
]],
	})
	solver:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^ij',
		code = [[
const global <?=eqn.cons_t?>* U = buf + index;
sym3 gamma_uu = calc_gamma_uu(U, x);
value->vreal = sym3_dot(value->vsym3, gamma_uu);]],
	})
end

return EinsteinEquation
