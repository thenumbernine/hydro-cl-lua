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
	'alpha',
	'gamma_xx', 'gamma_xy', 'gamma_xz', 'gamma_yy', 'gamma_yz', 'gamma_zz',
	'a_x', 'a_y', 'a_z',
	'd_xxx', 'd_xxy', 'd_xxz', 'd_xyy', 'd_xyz', 'd_xzz',
	'd_yxx', 'd_yxy', 'd_yxz', 'd_yyy', 'd_yyz', 'd_yzz',
	'd_zxx', 'd_zxy', 'd_zxz', 'd_zyy', 'd_zyz', 'd_zzz',
	'K_xx', 'K_xy', 'K_xz', 'K_yy', 'K_yz', 'K_zz',
	'V_x', 'V_y', 'V_z',
}
ADM_BonaMasso_3D.numStates = 7 + 30	-- should equal # consVars
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

ADM_BonaMasso_3D.guiVars = table{
	require 'guivar.combo'{
		name = 'f',
		options = {
			'1 + 1/alpha^2',
			'2/alpha', 
			'1', '.49', '.5', '1.5', '1.69', 
		},
	}
}

function ADM_BonaMasso_3D:getCodePrefix()
	local initState = self.initStates[self.solver.initStatePtr[0]+1]
	assert(initState, "couldn't find initState "..self.solver.initStatePtr[0])	
	
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
	
	-- TODO V^i
	
	lines:insert'}'
	
	local code = lines:concat'\n'
	return code
end

function ADM_BonaMasso_3D:getSolverCode()
	return require 'template'(file['eqn/adm3d.cl'], {eqn=self, solver=self.solver})
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
		return {[var] = 'value = U->'..code..';'}
	end))
	:append{
		{det_gamma = 'value = sym3_det(U->gamma);'},
		{volume = 'value = U->alpha * sqrt(sym3_det(U->gamma));'},
		{f = 'value = calc_f(U->alpha);'},
		{K = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	value = sym3_dot(gammaU, U->K);
]]},
		{expansion = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	value = -U->alpha * sym3_dot(gammaU, U->K);
]]},	
		-- TODO needs shift influence (which is lengthy)
		{gravityMagn = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	value = real3_len(sym3_real3_mul(gammaU, U->a));
]]},
		-- TODO V^i constraint codes
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

function ADM_BonaMasso_3D:getEigenDisplayVars()
	return {
		{sqrt_f = 'value = eigen->sqrt_f;'},
		{gammaUxx = 'value = eigen->gammaU.xx;'},
		{gammaUxy = 'value = eigen->gammaU.xy;'},
		{gammaUxz = 'value = eigen->gammaU.xz;'},
		{gammaUyy = 'value = eigen->gammaU.yy;'},
		{gammaUyz = 'value = eigen->gammaU.yz;'},
		{gammaUzz = 'value = eigen->gammaU.zz;'},
	}
end

return ADM_BonaMasso_3D
