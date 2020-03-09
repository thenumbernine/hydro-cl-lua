-- check out my 'wave equation hyperbolic form' worksheet
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local common = require 'common'
local xNames = common.xNames

local Wave = class(Equation)
Wave.name = 'wave'

-- depending on which you use, the source terms change
Wave.weightFluxByGridVolume = true
--Wave.weightFluxByGridVolume = false

Wave.useSourceTerm = true

Wave.hasEigenCode = true
Wave.hasFluxFromConsCode = true
Wave.roeUseFluxFromCons = true

Wave.initStates = require 'init.euler'	 -- use rho as our initial condition


--[[
args:
	scalar = 'real' or 'cplx'

	usePressure = whether to use pressure or density from the initial conditions for the wave equation
					Sedov initial conditions uses pressure
					gaussian, etc use density
--]]
function Wave:init(args)
	if args and args.usePressure ~= nil then
		self.usePressure = not not args.usePressure
	else
		self.usePressure = false
	end
	
	self.scalar = (args and args.scalar) or 'real'


	self.vec3 = self.scalar..'3'
	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'
	
	self.numStates = 4 * self.numRealsInScalar

	if self.usePressure then
		self.consVars = table{
			{name='Pi', type=self.scalar, units='kg/(m*s^3)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^2*s^2)'},
		}
	else
		self.consVars = table{
			{name='Pi', type=self.scalar, units='kg/(m^3*s)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^4)'},
		}
	end

	local suffix = self.scalar == 'real' and '' or ' re'
	self.predefinedDisplayVars = {
		'U Pi'..suffix,
		'U Psi_l'..suffix,
		'U Psi_l x'..suffix,
		'U Psi_l y'..suffix,
		'U Psi_l z'..suffix,
		'U Psi_l mag',
	}

	self.init_alpha = args.alpha
	self.init_beta = args.beta
	self.init_K = args.K
	
	Wave.super.init(self, args)
end

function Wave:compile(expr)
	return self.solver.coord:compile(expr)
end

function Wave:createInitState()
	Wave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},
	}
end

Wave.initStateCode = [[
<?
local common = require 'common'
local xNames = common.xNames

local scalar = eqn.scalar
local vec3 = eqn.vec3
?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;
	
	<?=code?>

	UBuf[index] = (<?=eqn.cons_t?>){
<? if eqn.usePressure then
?>		.Pi = <?=scalar?>_from_real(P),
<? else		
?>		.Pi = <?=scalar?>_from_real(rho),
<? end		
?>		.Psi_l = <?=vec3?>_from_real3(cartesianToCoord(v, x)),
	};
}
]]

Wave.solverCodeFile = 'eqn/wave.cl'

function Wave:getCommonFuncCode()
	
	local symmath = require 'symmath'
	local Tensor = symmath.Tensor
	local Constant = symmath.Constant
	local var = symmath.var
	local vars = symmath.vars
	local fromlua = require 'ext.fromlua'
	local coords = Tensor.coords()[1].variables
	local t = var't'
	self.metric = {
		coords = coords,
		t = t,
		
		alpha = Constant(1),
		beta_u = {
			Constant(0),
			Constant(0),
			Constant(0),
		},
		K = Constant(0),
	}

	local x = self.solver.coord.vars.x
	local y = self.solver.coord.vars.y
	local z = self.solver.coord.vars.z
	local r = self.solver.coord.vars.r
	local function readarg(s)
		return symmath.clone(assert(load('local x,y,z,t,r = ... return '..s))(x,y,z,t,r))()
	end

	if self.init_alpha then
		self.metric.alpha = readarg(self.init_alpha)
	end
	if self.init_beta then
		for i=1,3 do
			if self.init_beta[i] then
				self.metric.beta_u[i] = readarg(self.init_beta[i])
			end
		end
	end
	if self.init_K then
		self.metric.K = readarg(self.init_K)
	end



	
	
	return template([[
<?
local scalar = eqn.scalar
local vec3 = eqn.vec3
local common = require 'common'
local xNames = common.xNames
?>
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

real metric_dt_alpha(real3 pt) { 
	return <?=eqn:compile(eqn.metric.alpha:diff(eqn.metric.t)())?>;
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

<?=scalar?> eqn_source(real3 pt) { return <?=scalar?>_zero; }

]], {
		eqn = self,
	})
end

Wave.eigenVars = {
	{name='unused', type='real'},
}

function Wave:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real wavespeed = solver->wavespeed / unit_m_per_s;
	real nLen = coord_sqrt_g_uu<?=side..side?>(<?=x?>);
	real alpha_nLen = metric_alpha(<?=x?>) * nLen;
	real3 beta_u = metric_beta_u(<?=x?>);
]], {
		side = side,
		x = x,
	})
end

function Wave:consWaveCodePrefix(side, U, x, W)
	return self:eigenWaveCodePrefix(side, nil, x)
end

function Wave:consWaveCode(side, U, x, waveIndex)
	waveIndex = math.floor(waveIndex / self.numRealsInScalar)
	local xside = xNames[side+1]
	if waveIndex == 0 then
		return 'wavespeed * (-beta_u.'..xside..' - alpha_nLen)' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return 'wavespeed * -beta_u.'..xside
	elseif waveIndex == 3 then
		return 'wavespeed * (-beta_u.'..xside..' + alpha_nLen)' 
	end
	error'got a bad waveIndex'
end

Wave.eigenWaveCode = Wave.consWaveCode

return Wave
