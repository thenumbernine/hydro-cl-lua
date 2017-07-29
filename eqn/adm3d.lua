--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local Equation = require 'eqn.eqn'

local xNames = table{'x', 'y', 'z'}

-- symmetric indexes: xx xy xz yy yz zz
local symNames = table()
for i,xi in ipairs(xNames) do
	for j=i,3 do
		local xj = xNames[j]
		symNames:insert(xi..xj)
	end
end

local ADM_BonaMasso_3D = class(Equation)
ADM_BonaMasso_3D.name = 'ADM_BonaMasso_3D'

ADM_BonaMasso_3D.consVars = {
--[[ 0	]]	'alpha',
--[[ 1	]]	'gamma_xx', 'gamma_xy', 'gamma_xz', 'gamma_yy', 'gamma_yz', 'gamma_zz',
--[[ 7	]]	'a_x', 'a_y', 'a_z',
--[[ 10	]]	'd_xxx', 'd_xxy', 'd_xxz', 'd_xyy', 'd_xyz', 'd_xzz',
--[[ 16	]]	'd_yxx', 'd_yxy', 'd_yxz', 'd_yyy', 'd_yyz', 'd_yzz',
--[[ 22	]]	'd_zxx', 'd_zxy', 'd_zxz', 'd_zyy', 'd_zyz', 'd_zzz',
--[[ 28	]]	'K_xx', 'K_xy', 'K_xz', 'K_yy', 'K_yz', 'K_zz',
--[[ 34	]]	'V_x', 'V_y', 'V_z',
--[[ 37	]]
}
ADM_BonaMasso_3D.numStates = 37
ADM_BonaMasso_3D.numWaves = 30	-- skip alpha and gamma_ij

assert(ADM_BonaMasso_3D.numStates == #ADM_BonaMasso_3D.consVars)

ADM_BonaMasso_3D.hasCalcDT = true
ADM_BonaMasso_3D.hasEigenCode = true
ADM_BonaMasso_3D.useSourceTerm = true
ADM_BonaMasso_3D.useConstrainU = true

ADM_BonaMasso_3D.initStates = require 'init.adm'

function ADM_BonaMasso_3D:getTypeCode()
	return template([[
typedef union {
	real ptr[37];
	struct {
		real alpha;
		sym3 gamma;
		real3 a;
		sym3 d[3];
		sym3 K;
		real3 V;
	};
} <?=eqn.cons_t?>;
]], {
	eqn = self,
})
end

function ADM_BonaMasso_3D:init(...)
	self.guiVars = table{
		require 'guivar.combo'{
			name = 'f',
			options = {
				'1 + 1/alpha^2',
				'2/alpha', 
				'1', '.49', '.5', '1.5', '1.69', 
			},
		}
	}
	ADM_BonaMasso_3D.super.init(self, ...)
end

local symmath = require 'symmath'
function ADM_BonaMasso_3D:getCodePrefix()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)	
	
	return initState.init(self.solver, function(exprs, vars, args)
		print('building lapse partials...')
		exprs.a = table.map(vars, function(var)
			return (exprs.alpha:diff(var) / exprs.alpha)()
		end)
		print('...done building lapse partials')
		
		print('building metric partials...')
		exprs.d = table.map(vars, function(xk)
			return table.map(exprs.gamma, function(gamma_ij)
				print('differentiating '..gamma_ij)
				return (gamma_ij:diff(xk)/2)()
			end)
		end)
		print('...done building metric partials')

		--[[
		local gammaU = table{mat33.inv(gamma:unpack())} or calc.gammaU:map(function(gammaUij) return gammaUij(x,y,z) end) 
		exprs.V = table.map(vars, function(var)
			local s = 0
			for j=1,3 do
				for k=1,3 do
					local d_ijk = sym3x3(exprs.d[i],j,k)
					local d_kji = sym3x3(exprs.d[k],j,i)
					local gammaUjk = sym3x3(gammaU,j,k)
					local dg = (d_ijk - d_kji) * gammaUjk
					s = s + dg
				end
			end
			return s
		end)
		--]]

		return table(
			{alpha = exprs.alpha},
			symNames:map(function(xij,ij)
				return exprs.gamma[ij], 'gamma_'..xij
			end),
			xNames:map(function(xi,i)
				return exprs.a[i], 'a_'..xi
			end),
			table(xNames:map(function(xk,k,t)
				return symNames:map(function(xij,ij)
					return exprs.d[k][ij], 'd_'..xk..xij
				end), #t+1
			end):unpack()),
			symNames:map(function(xij,ij)
				return exprs.K[ij], 'K_'..xij
			end)
			--[[
			xNames:map(function(xi,i)
				return exprs.V[i], 'V_'..xi
			end),
			--]]
		)
	end)
end

local xNames = table{'x', 'y', 'z'}
local symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

function ADM_BonaMasso_3D:getInitStateCode()
	local lines = table{
		template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
]], {
	eqn = self,
}),
	}

	local function build(var)
		local prefix, suffix = var:match'(.*)_(.*)'
		local field = var
		if prefix then
			if #suffix == 3 then
				field = prefix..'['..(('xyz'):find(suffix:sub(1,1))-1)..'].'..suffix:sub(2)
			else
				field = prefix..'.'..suffix
			end
		end
		lines:insert('\tU->'..field..' = calc_'..var..'(x.x, x.y, x.z);')
	end
	
	build'alpha'
	symNames:map(function(xij) build('gamma_'..xij) end)
	xNames:map(function(xi) build('a_'..xi) end)	
	xNames:map(function(xk)
		symNames:map(function(xij) build('d_'..xk..xij) end)
	end)
	symNames:map(function(xij) build('K_'..xij) end)

	--[[ symbolic
	xNames:map(function(xi) build('V_'..xi) end)
	--]]
	-- [[
	lines:insert'	U->V = _real3(0,0,0);'
	--]]

	lines:insert'}'
	
	local code = lines:concat'\n'
	return code
end

function ADM_BonaMasso_3D:getSolverCode()
	return template(file['eqn/adm3d.cl'], {eqn=self, solver=self.solver})
end

function ADM_BonaMasso_3D:getDisplayVars()
	return table()
	:append(table.map(ADM_BonaMasso_3D.consVars, function(var)
		local code = var:gsub('_', '.')
		local dx, rest = code:match('^d%.([xyz])([xyz][xyz])$')
		if dx then
			local i = xNames:find(dx)
			code = 'd['..(i-1)..'].'..rest
		end
		return {[var] = '*value = U->'..code..';'}
	end))
	:append{
		{det_gamma = '*value = sym3_det(U->gamma);'},
		{volume = '*value = U->alpha * sqrt(sym3_det(U->gamma));'},
		{f = '*value = calc_f(U->alpha);'},
		{K = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = sym3_dot(gammaU, U->K);
]]},
		{expansion = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = -U->alpha * sym3_dot(gammaU, U->K);
]]},	
		-- TODO needs shift influence (which is lengthy)
		{gravityMagn = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = real3_len(sym3_real3_mul(gammaU, U->a));
]]},
		-- TODO V^i constraint codes

-- ... and here's the dt code ...
		{lambdaLight_x = template([[
	real det_gamma = sym3_det(U->gamma);
	real f = calc_f(U->alpha);
	real sqrt_f = sqrt(f);
	real gammaUxx = (U->gamma.yy * U->gamma.zz - U->gamma.yz * U->gamma.yz) / det_gamma;
	*value = U->alpha * sqrt(gammaUxx);
]], {solver=self.solver})},

		{lambdaLight_y = template([[
	real det_gamma = sym3_det(U->gamma);
	real f = calc_f(U->alpha);
	real sqrt_f = sqrt(f);
	real gammaUyy = (U->gamma.xx * U->gamma.zz - U->gamma.xz * U->gamma.xz) / det_gamma;
	*value = U->alpha * sqrt(gammaUyy);
]], {solver=self.solver})},
		
		{lambdaLight_z = template([[
	real det_gamma = sym3_det(U->gamma);
	real f = calc_f(U->alpha);
	real sqrt_f = sqrt(f);
	real gammaUzz = (U->gamma.xx * U->gamma.yy - U->gamma.xy * U->gamma.xy) / det_gamma;
	*value = U->alpha * sqrt(gammaUzz);
]], {solver=self.solver})},

		{dt = template([[
	real det_gamma = sym3_det(U->gamma);
	real f = calc_f(U->alpha);
	real sqrt_f = sqrt(f);

	*value = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side==0 then ?>
		real gammaUxx = (U->gamma.yy * U->gamma.zz - U->gamma.yz * U->gamma.yz) / det_gamma;
		real lambdaLight = U->alpha * sqrt(gammaUxx);
		<? elseif side==1 then ?>
		real gammaUyy = (U->gamma.xx * U->gamma.zz - U->gamma.xz * U->gamma.xz) / det_gamma;
		real lambdaLight = U->alpha * sqrt(gammaUyy);
		<? elseif side==2 then ?>
		real gammaUzz = (U->gamma.xx * U->gamma.yy - U->gamma.xy * U->gamma.xy) / det_gamma;
		real lambdaLight = U->alpha * sqrt(gammaUzz);
		<? end ?>	
		
		real lambdaGauge = lambdaLight * sqrt_f;
		real lambda = (real)max(lambdaGauge, lambdaLight);
		
		real lambdaMin = (real)min((real)0., -lambda);
		real lambdaMax = (real)max((real)0., lambda);
		*value = (real)min((real)*value, (real)(dx<?=side?>_at(i) / ((real)fabs(lambdaMax - lambdaMin) + (real)1e-9)));
	}<? end ?>
]], {solver=self.solver})},
	}
end

function ADM_BonaMasso_3D:getEigenTypeCode()
	return template([[
typedef struct {
	real alpha;	//used only by eigen_calcWaves ... makes me think eigen_forCell / eigen_forSide should both calculate waves and basis variables in the same go
	real sqrt_f;
	sym3 gammaU;
	
	//sqrt(gamma^jj) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
	real3 sqrt_gammaUjj;	
} <?=eqn.eigen_t?>;
]], {
	eqn = self,
})
end

function ADM_BonaMasso_3D:getVecDisplayVars()
	local vars = table()
	local function add(field)
		vars:insert{[field] = 'valuevec = U->'..field..';'}
	end
	local function addSym3(field)
		for i,xi in ipairs(xNames) do
			vars:insert{[field..'_'..xi] = 'valuevec = sym3_'..xi..'(U->'..field..');'}
		end
	end
	addSym3'gamma'
	add'a'
	addSym3'K'
	-- shift-less gravity only
	-- gravity with shift is much more complex
	vars:insert{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	valuevec = real3_scale(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]]}
	return vars
end

function ADM_BonaMasso_3D:getEigenDisplayVars()
	return {
		{sqrt_f = '*value = eigen->sqrt_f;'},
		{gammaUxx = '*value = eigen->gammaU.xx;'},
		{gammaUxy = '*value = eigen->gammaU.xy;'},
		{gammaUxz = '*value = eigen->gammaU.xz;'},
		{gammaUyy = '*value = eigen->gammaU.yy;'},
		{gammaUyz = '*value = eigen->gammaU.yz;'},
		{gammaUzz = '*value = eigen->gammaU.zz;'},
	}
end

return ADM_BonaMasso_3D
