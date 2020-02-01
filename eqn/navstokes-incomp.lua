local class = require 'ext.class'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local NavierStokesDivFree = class(Equation)
NavierStokesDivFree.name = 'Navier-Stokes div-free' 
NavierStokesDivFree.hasCalcDTCode = true
NavierStokesDivFree.numStates = 4

NavierStokesDivFree.reflectVars = {
	mirror = {
		{'v.x'},
		{'v.y'},
		{'v.z'},
	},
}

NavierStokesDivFree.initStates = require 'init.euler'

function NavierStokesDivFree:getTypeCode()
	return template([[
typedef union {
	real ptr[4];
	struct {
		real rho;
		real3 v;
	};
} <?=eqn.prim_t?>;

typedef <?=eqn.prim_t?> <?=eqn.cons_t?>; 
]], {
		eqn = self,
	})
end

NavierStokesDivFree.initStateCode = [[
kernel void initState(
	global <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	real rho = 0;
	real3 v = real3_zero;
	
	// throw-away:
	real P = 0;
	real3 B = real3_zero;
	real ePot = 0;

	<?=code?>

	<?=eqn.cons_t?> U = {
		.rho = rho,
		.v = v,
	};
	UBuf[index] = U;
}
]]

NavierStokesDivFree.solverCodeFile = 'eqn/euler.cl'

function NavierStokesDivFree:getDisplayVars()
	-- k is 0,1,2
	local function vorticity(k,result)
		local xs = {'x','y','z'}
		local i = (k+1)%3
		local j = (i+1)%3
		return {['vorticity '..xs[k+1]] = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = buf + index - solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = buf + index + solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = buf + index - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = buf + index + solver->stepsize.s<?=j?>;
		
		real3 vim = Uim->v;
		real3 vip = Uip->v;
		real3 vjm = Ujm->v;
		real3 vjp = Ujp->v;
		
		<?=result?> = (vjp.s<?=i?> - vjm.s<?=i?>) / (2. * solver->grid_dx.s<?=i?>)
				- (vip.s<?=j?> - vim.s<?=j?>) / (2. * solver->grid_dx.s<?=j?>);
	}
]], 	{
			i = i,
			j = j,
			eqn = self,
		})}
	end
	local vars = table{
		{rho = 'value.vreal = U->rho;'},
		{v = 'value.vreal3 = U->v;', type='real3'},
		{m = 'value = real3_real_mul(U->v, U->rho);', type='real3'},
		--{P = 'value.vreal = W.P;'},
		--{eInt = 'value.vreal = calc_eInt(W);'},
		--{eKin = 'value.vreal = calc_eKin(W);'},
		--{ePot = 'value.vreal = U->ePot;'},
		--{eTotal = 'value.vreal = U->ETotal / W.rho;'},
		--{EInt = 'value.vreal = calc_EInt(W);'},
		{EKin = 'value.vreal = .5 * U->rho * coordLen(U->v, x);'},
		--{EPot = 'value.vreal = U->rho * U->ePot;'},
		--{ETotal = 'value.vreal = U->ETotal;'},
		--{S = 'value.vreal = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		--{H = 'value.vreal = calc_H(W.P);'},
		--{h = 'value.vreal = calc_h(W.rho, W.P);'},
		--{HTotal = 'value.vreal = calc_HTotal(W.P, U->ETotal);'},
		--{hTotal = 'value.vreal = calc_hTotal(W.rho, W.P, U->ETotal);'},
		--{['Speed of Sound'] = 'value.vreal = calc_Cs(&W);'},
		--{['Mach number'] = 'value.vreal = coordLen(W.v, x) / calc_Cs(&W);'},
	}
	
	if self.solver.dim == 2 then
	-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
		vars:insert(vorticity(2,'value.vreal'))
	elseif self.solver.dim == 3 then
		local v = range(0,2):map(function(i) return vorticity(self,i,'value['..i..']') end)
		vars:insert{vorticityVec = template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
		++value;
	}<? end ?>
	value -= 3;
]], {v=v}), type='real3'}
	end	
			
	return vars
end

return NavierStokesDivFree
