local class = require 'ext.class'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local NavierStokesDivFree = class(Equation)
NavierStokesDivFree.name = 'Navier-Stokes div-free' 
NavierStokesDivFree.numStates = 4

NavierStokesDivFree.mirrorVars = {{'v.x'}, {'v.y'}, {'v.z'}}

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

function NavierStokesDivFree:getInitStateCode()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)	
	local code = initState.init(self.solver)	
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	real rho = 0;
	real3 v = _real3(0,0,0);
	
	// throw-away:
	real P = 0;
	real3 B = _real3(0,0,0);
	real ePot = 0;

]]..code..[[

	<?=eqn.cons_t?> U = {
		.rho = rho,
		.v = v,
	};
	UBuf[index] = U;
}
]], {
		eqn = self,
	})
end

function NavierStokesDivFree:getSolverCode()
	return template(file['eqn/euler.cl'], {eqn=self, solver=self.solver})
end

function NavierStokesDivFree:getDisplayVarCodePrefix()
	return template([[
	<?=eqn.cons_t?> U = buf[index];
]], {
		eqn = self,
	})
end

function NavierStokesDivFree:getDisplayVars()
	-- k is 0,1,2
	local function vorticity(k)
		local xs = {'x','y','z'}
		local i = (k+1)%3
		local j = (i+1)%3
		return {['vorticity '..xs[k+1]] = template([[
	global const <?=eqn.cons_t?>* Uim = buf + index - stepsize.s<?=i?>;
	global const <?=eqn.cons_t?>* Uip = buf + index + stepsize.s<?=i?>;
	global const <?=eqn.cons_t?>* Ujm = buf + index - stepsize.s<?=j?>;
	global const <?=eqn.cons_t?>* Ujp = buf + index + stepsize.s<?=j?>;
	
	real3 vim = Uim->v;
	real3 vip = Uip->v;
	real3 vjm = Ujm->v;
	real3 vjp = Ujp->v;
	
	*value = (vjp.s<?=i?> - vjm.s<?=i?>) / (2. * grid_dx<?=i?>)
			- (vip.s<?=j?> - vim.s<?=j?>) / (2. * grid_dx<?=j?>);
]], 	{
			i = i,
			j = j,
			eqn = self,
		})}
	end
	return table{
		{rho = '*value = U.rho;'},
		{vx = '*value = U.v.x;'},
		{vy = '*value = U.v.y;'},
		{vz = '*value = U.v.z;'},
		{v = '*value = coordLen(U.v, x);'},
		{mx = '*value = U.rho * U.v.x;'},
		{my = '*value = U.rho * U.v.y;'},
		{mz = '*value = U.rho * U.v.z;'},
		{m = '*value = U.rho * coordLen(U.v, x);'},
		--{P = '*value = W.P;'},
		--{eInt = '*value = calc_eInt(W);'},
		--{eKin = '*value = calc_eKin(W);'},
		--{ePot = '*value = U.ePot;'},
		--{eTotal = '*value = U.ETotal / W.rho;'},
		--{EInt = '*value = calc_EInt(W);'},
		{EKin = '*value = .5 * U.rho * coordLen(U.v, x);'},
		--{EPot = '*value = U.rho * U.ePot;'},
		--{ETotal = '*value = U.ETotal;'},
		--{S = '*value = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		--{H = '*value = calc_H(W.P);'},
		--{h = '*value = calc_h(W.rho, W.P);'},
		--{HTotal = '*value = calc_HTotal(W.P, U.ETotal);'},
		--{hTotal = '*value = calc_hTotal(W.rho, W.P, U.ETotal);'},
		--{['Speed of Sound'] = '*value = calc_Cs(&W);'},
		--{['Mach number'] = '*value = coordLen(W.v, x) / calc_Cs(&W);'},
	}:append( ({
	-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
			[1] = {},
			[2] = {vorticity(2)},
			[3] = range(0,2):map(vorticity),
	})[self.solver.dim] )
end

return NavierStokesDivFree 
