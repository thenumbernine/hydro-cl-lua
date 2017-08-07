--[[
based on Marti 1998, Marti & Muller 2008, and maybe some of Font 2008 (but that's grhd)

honestly I developed a Marti & Muller SRHD solver
then I bumped it up to GRHD by incorporating (fixed) alphas betas and gammas
and then I thought "why not keep the old SRHD solver around"
so viola, here it is, unnecessarily available.
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local SRHD = class(Equation)
SRHD.name = 'SRHD'
SRHD.numStates = 6
SRHD.numWaves = 5
SRHD.numIntStates = 5

SRHD.mirrorVars = {{'S.x'}, {'S.y'}, {'S.z'}}

SRHD.hasEigenCode = true 

-- SRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--SRHD.hasFluxFromCons = true

SRHD.hasCalcDT = true

SRHD.useConstrainU = true

SRHD.initStates = require 'init.euler'

local GuiFloat = require 'guivar.float'
local GuiInt = require 'guivar.int'
function SRHD:init(solver)

	-- hmm, setting the higher thresholds using double precision is less stable
	local double = false --solver.app.real == 'double'
	
	self.guiVars = table{
		GuiFloat{name='heatCapacityRatio', value=7/5},
		
		-- setting max iter to 100+ makes it freeze initially 
		-- but setting it to 100 after the first iteration is fine ...
		-- meaning the initial cons to prim is taking too long ...
		GuiInt{name='solvePrimMaxIter', value=10},	-- value=1000},

		GuiFloat{name='solvePrimStopEpsilon', value=1e-7},
	
		-- used by pressure solver
		-- velocity epsilon is how close we can get to the speed of light
		-- set ylabel "Lorentz factor"; set xlabel "velocity epsilon -log10"; set log xy; plot [1:10] 1/sqrt(1-(1-10**(-x))**2);
		--velEpsilon = 1e-5	-- <=> handles up to W = 500
		--velEpsilon = 1e-6	-- <=> handles up to W = 600
		--velEpsilon = 1e-7	-- <=> handles up to W = 2,000
		--velEpsilon = 1e-10	-- <=> handles up to W = 100,000
		-- <=> smaller than 1e-15 gnuplot x11 terminal breaks down past W = 1e+7 ...
		GuiFloat{name='solvePrimVelEpsilon', value=double and 1e-15 or 1e-7},
		
		GuiFloat{name='solvePrimPMinEpsilon', value=double and 1e-16 or 1e-7},
		
		GuiFloat{name='rhoMin', value=double and 1e-15 or 1e-7},
		GuiFloat{name='rhoMax', value=1e+20},
		GuiFloat{name='eIntMax', value=1e+20},
		GuiFloat{name='DMin', value=double and 1e-15 or 1e-7},
		GuiFloat{name='DMax', value=1e+20},
		GuiFloat{name='tauMin', value=double and 1e-15 or 1e-7},
		GuiFloat{name='tauMax', value=1e+20},
	}
	SRHD.super.init(self, solver)
end

function SRHD:getTypeCode()
	return template([[
typedef struct {
	real rho;
	real3 v;
	real eInt;
	real ePot;
} <?=eqn.prim_t?>;

typedef union {
	real ptr[6];
	struct {
		real D;
		real3 S;
		real tau;
		
		// TODO fix this.
		// it is here because prim_t is expected to be the same size as cons_t
		real unused;
	};
} <?=eqn.cons_t?>;
]], {
	eqn = self,
})
end

function SRHD:getCodePrefix()
	return table{
		SRHD.super.getCodePrefix(self),
		template([[

//pressure function for ideal gas
real calc_P(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * rho * eInt;
}	

//chi in most papers
real calc_dP_drho(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * eInt;
}

//kappa in most papers
real calc_dP_deInt(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * rho;
}

real calc_eInt_from_P(real rho, real P) {
	return P / ((heatCapacityRatio - 1.) * rho);
}

real calc_h(real rho, real P, real eInt) {
	return 1. + eInt + P / rho;
}

<?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> prim, real3 x) {
	real vSq = coordLenSq(prim.v, x);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 40-42:
	
	//rest-mass density = J^0 = rho u^0
	real D = prim.rho * W;	
	
	//momentum = T^0i = rho h u^0 u^i + P g^0i
	real3 S = real3_scale(prim.v, prim.rho * h * WSq);
	
	//energy = T^00 = rho h u^0 u^0 + P g^00
	real tau = prim.rho * h * WSq - D - P;
	
	return (<?=eqn.cons_t?>){.D=D, .S=S, .tau=tau};
}
]], {
	eqn = self,
}),
	}:concat'\n'
end

function SRHD:getInitStateCode()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)
	local code = initState.init(self.solver)
	return template([[

kernel void initState(
	global <?=eqn.cons_t?>* consBuf,
	global <?=eqn.prim_t?>* primBuf
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
	real P = 0;
	//ignored:
	real3 B = _real3(0,0,0);

]]..code..[[
	
	real eInt = calc_eInt_from_P(rho, P);
	real vSq = coordLenSq(v, x);
	real W = 1. / sqrt(1. - vSq);
	real h = calc_h(rho, P, eInt);

	<?=eqn.prim_t?> prim = {.rho=rho, .v=v, .eInt=eInt};
	primBuf[index] = prim;
	consBuf[index] = consFromPrim(prim, x);
}
]], {
	eqn = self,
})
end

function SRHD:getSolverCode()
	return template(file['eqn/srhd.cl'], {
		eqn = self,
		solver = self.solver,
	})
end

function SRHD:getDisplayVarCodePrefix()
	return template([[
	<?=eqn.cons_t?> U = buf[index];
	<?=eqn.prim_t?> prim = primBuf[index];
]], {
	eqn = self,
})
end

-- TODO put in common parent of Euler, SRHD, GRHD
-- k is 0,1,2
local function vorticity(eqn,k)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {['vorticity '..xs[k+1]] = template([[
	global const <?=eqn.prim_t?>* prim_im = primBuf + index - stepsize.s<?=i?>;
	global const <?=eqn.prim_t?>* prim_ip = primBuf + index + stepsize.s<?=i?>;
	global const <?=eqn.prim_t?>* prim_jm = primBuf + index - stepsize.s<?=j?>;
	global const <?=eqn.prim_t?>* prim_jp = primBuf + index + stepsize.s<?=j?>;

	//TODO incorporate metric
	//TODO 3-vorticity vs 4-vorticity?

	real vim_j = prim_im->v.s<?=j?>;
	real vip_j = prim_ip->v.s<?=j?>;
	
	real vjm_i = prim_jm->v.s<?=i?>;
	real vjp_i = prim_jp->v.s<?=i?>;
	
	*value = (vjp_i - vjm_i) / (2. * grid_dx<?=i?>)
			- (vip_j - vim_j) / (2. * grid_dx<?=j?>);
]], {
		i = i,
		j = j,
		eqn = eqn,
	})}
end

function SRHD:getDisplayVars()
	return table{
		{D = '*value = U.D;'},
		{Sx = '*value = U.S.x;'},
		{Sy = '*value = U.S.y;'},
		{Sz = '*value = U.S.z;'},
		{S = '*value = coordLen(U.S, x);'},
		{tau = '*value = U.tau;'},
		{['W based on D'] = '*value = U.D / prim.rho;'},
		{['W based on v'] = '*value = 1. / sqrt(1. - coordLenSq(prim.v, x));'},
		{['primitive reconstruction error'] = template([[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	{
		<?=eqn.cons_t?> U2 = consFromPrim(prim, x);
		*value = 0;
		for (int j = 0; j < numStates; ++j) {
			*value += fabs(U.ptr[j] - U2.ptr[j]);
		}
	}
	]], {eqn=self})},
		{['W error'] = [[
	real W1 = U.D / prim.rho;
	real W2 = 1. / sqrt(1. - coordLenSq(prim.v, x));
	*value = fabs(W1 - W2);
		]]},
	}:append( ({
	-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
			[1] = {},
			[2] = {vorticity(self,2)},
			[3] = range(0,2):map(function(i) return vorticity(self,i) end),
	})[self.solver.dim] )
end

function SRHD:getPrimDisplayVarCodePrefix()
	return template([[
	<?=eqn.prim_t?> prim = buf[index];
]], {
		eqn = self,
	})
end

SRHD.primDisplayVars = {
	{rho = '*value = prim.rho;'},
	{vx = '*value = prim.v.x;'},
	{vy = '*value = prim.v.y;'},
	{vz = '*value = prim.v.z;'},
	{v = '*value = coordLen(prim.v, x);'},
	{eInt = '*value = prim.eInt;'},
	{ePot = '*value = prim.ePot;'},
	{P = '*value = calc_P(prim.rho, prim.eInt);'},
	{h = '*value = calc_h(prim.rho, calc_P(prim.rho, prim.eInt), prim.eInt);'},
}

function SRHD:getVecDisplayVars()
	local vars = table{
		{v = 'valuevec = prim.v;'},
		{S = 'valuevec = U.S;'},
	}
	if self.solver.dim == 3 then
		local v = range(0,2):map(function(i) return vorticity(self,i) end)
		vars:insert{vorticityVec = template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
		++value;
	}<? end ?>
	value -= 3;
]], {v=v})}
	end
	return vars
end

SRHD.eigenStructFields = {
	{rho = 'real'},
	{v = 'real3'},
	{h = 'real'},
	{W = 'real'},
	{ATildeMinus = 'real'},
	{ATildePlus = 'real'},
	{VMinus = 'real'},
	{VPlus = 'real'},
	{CMinus = 'real'},
	{CPlus = 'real'},
	{Kappa = 'real'},
}

function SRHD:getEigenTypeCode()
	return 'typedef struct {\n'
		..table.map(self.eigenStructFields, function(field)
			local name, ctype = next(field)
			return '\t'..ctype..' '..name..';\n'
		end):concat'\n'
		..'} '..self.eigen_t..';\n'
end

function SRHD:getEigenDisplayVars()
	return {}
end

return SRHD
