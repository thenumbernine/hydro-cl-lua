local class = require 'ext.class'
local table = require 'ext.table'

local ViscousExplicit = class()

ViscousExplicit.name = 'ViscousExplicit'

function ViscousExplicit:init(args)
	self.solver = assert(args.solver)
	
	require 'hydro.code.symbols'(self, self:getSymbolFields())
end

function ViscousExplicit:getSymbolFields()
	return table{
		'viscousExplicitUpdate',
	}
end

function ViscousExplicit:initCodeModules()
	local solver = self.solver
	solver.modules:addFromMarkup{
		code = solver.eqn:template[[
//// MODULE_NAME: <?=viscousExplicitUpdate?>

kernel void <?=viscousExplicitUpdate?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
// viscous terms as rhs.  work is in hydro/op/viscousflux.lua
// only for grid solvers at the moment
// notice that in the real world shear viscosity changes with temperature (which changes with internal energy, which changes with position and time)
// so add more derivative terms?

	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	real3 const x = cellBuf[index].pos;

	real3 v = <?=calc_v?>(U);

<? local function getV(offset) return eqn.symbols.calc_v.."(U + "..offset..")" end ?>
<?=eqn:makePartial2(getV, "real3", "d2v")?>

	real3 div_tau = real3_zero;		// tau^ij_,j
	<? for i,xi in ipairs(xNames) do ?>
		<? for j,xj in ipairs(xNames) do ?>
			div_tau.<?=xi?> += solver->shearViscosity * (
				d2v[<?=from3x3to6(j,j)-1?>].<?=xi?>
				+ (1./3.) * d2v[<?=from3x3to6(i,j)-1?>].<?=xj?>
			);
		<? end ?>
	<? end ?>

<?=eqn:makePartial1(getV, "real3", "dv")?>
	
	real div_v = dv.x.x + dv.y.y + dv.z.z;

	real tau_dot_dv = 0.;
	<? for ij,xij in ipairs(symNames) do 
		local i,j,xi,xj = from6to3x3(ij)
	?>{
		real tau_ij = dv.<?=xi?>.<?=xj?> + dv.<?=xj?>.<?=xj?>;
	<? if i==j then ?>
		tau_ij -= (2./3.) * div_v;
	<? end ?>
		tau_ij *= solver->shearViscosity;
		tau_dot_dv += tau_ij * dv.<?=xi?>.<?=xj?>;
	}<? end ?>

	//q_i = -k T_,i
	//q^i_;i = -(k T_,i g^ij)_;j = -g^ij (k T_;ij + k_,i T_,j) = -g^ij (k T_,ij - k T_,k Γ^k_ij + k_,i T_,j)
	//for an ideal gas: T = eInt / Cv = (ETotal - 1/2 g^ij m_i m_j / ρ) / ρ / Cv

// TODO calc_T's 'x' should vary with the offset
<? local function getT(offset) return eqn.symbols.calc_T.."(U + "..offset..", x)" end ?>
<?=eqn:makePartial2(getT, "real", "d2T")?>

	deriv->m = real3_add(deriv->m, div_tau);
	deriv->ETotal += real3_dot(div_tau, v) + tau_dot_dv - solver->heatConductivity * sym3_trace(d2T);
}
]],
	}
	solver.solverModulesEnabled[self.symbols.viscousExplicitUpdate] = true
end

function ViscousExplicit:refreshSolverProgram()
	local solver = self.solver
	self.viscousExplicitUpdateKernelObj = solver.solverProgramObj:kernel(self.symbols.viscousExplicitUpdate)
	self.viscousExplicitUpdateKernelObj.obj:setArg(0, solver.solverBuf)
	self.viscousExplicitUpdateKernelObj.obj:setArg(2, solver.UBuf)
	self.viscousExplicitUpdateKernelObj.obj:setArg(3, solver.cellBuf)
end

function ViscousExplicit:addSource(derivBufObj)
	self.viscousExplicitUpdateKernelObj.obj:setArg(1, derivBufObj.obj)
	self.viscousExplicitUpdateKernelObj()
end

return ViscousExplicit
