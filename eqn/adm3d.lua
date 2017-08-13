--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local NumRelEqn = require 'eqn.numrel'

local xNames = table{'x', 'y', 'z'}
local symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

local from3x3to6_table = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6},}
local function from3x3to6(i,j) return from3x3to6_table[i][j] end

local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
local function from6to3x3(i) return table.unpack(from6to3x3_table[i]) end

local function sym(a,b)
	assert(a >= 1 and a <= 3, "tried to index sym with "..tostring(a)..", "..tostring(b))
	assert(b >= 1 and b <= 3, "tried to index sym with "..tostring(a)..", "..tostring(b))
	if a > b then a,b = b,a end
	return xNames[a]..xNames[b]
end

local ADM_BonaMasso_3D = class(NumRelEqn)
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

local symmath = require 'symmath'
function ADM_BonaMasso_3D:getCodePrefix()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)	
	
	local lines = table()
		
	-- don't call super because it generates the guivar code
	-- which is already being generated in initState
	--lines:insert(ADM_BonaMasso_3D.super.getCodePrefix(self))
	
	lines:insert(template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1.;
	U->gamma = _sym3(1,0,0,1,0,1);
	U->a = _real3(0,0,0);
	U->d[0] = _sym3(0,0,0,0,0,0);
	U->d[1] = _sym3(0,0,0,0,0,0);
	U->d[2] = _sym3(0,0,0,0,0,0);
	U->K = _sym3(0,0,0,0,0,0);
	U->V = _real3(0,0,0);
}
]], {eqn=self}))
	
	lines:insert(initState.init(self.solver, function(exprs, vars, args)
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
	end))

	return lines:concat()
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
	return template(file['eqn/adm3d.cl'], {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
		sym = sym,
	})
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
]]		},
		{expansion = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = -U->alpha * sym3_dot(gammaU, U->K);
]]		},
		-- TODO needs shift influence (which is lengthy)
		{gravityMagn = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = real3_len(sym3_real3_mul(gammaU, U->a));
]]		},
	}:append( xNames:map(function(xi,i)
		-- V_i - (d_im^m - d^m_mi)
		return {['constraint_V_'..xi] = template([[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	real d1 = sym3_dot(U->d[<?=i-1?>], gammaU);
	real d2 = 0.<?
for j=1,3 do
	for k,xk in ipairs(xNames) do
?> + U->d[<?=j-1?>].<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
	end
end ?>;
	*value = U->V.<?=xi?> - d1 + d2;
]], {i=i, xi=xi, sym=sym, xNames=xNames})}
	
	end) ):append{
		{constraint_V_magn = template([[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = 0.;
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d[<?=i-1?>], gammaU);
		real d2 = 0.<?
	for j=1,3 do
		for k,xk in ipairs(xNames) do
?> + U->d[<?=j-1?>].<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
		end
	end ?>;
		real d3 = U->V.<?=xi?> - (d1 - d2);
		*value += d3 * d3;
	}<? end ?>
	*value = sqrt(*value);
]], {sym=sym, xNames=xNames})}

--[=[
	-- 1998 Bona et al
--[[
H = 1/2 ( R + K^2 - K_ij K^ij ) - alpha^2 8 pi rho
for 8 pi rho = G^00

momentum constraints
--]]
		{H = [[
	.5 * 
]]		},
--]=]
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
	
	-- shift-less gravity only
	-- gravity with shift is much more complex
	vars:insert{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	valuevec = real3_scale(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]]}
	
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

	vars:insert{constraint_V = template([[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d[<?=i-1?>], gammaU);
		real d2 = 0.<?
	for j=1,3 do
		for k,xk in ipairs(xNames) do
?> + U->d[<?=j-1?>].<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
		end
	end ?>;
		valuevec.<?=xi?> = U->V.<?=xi?> - (d1 - d2);
	}<? end ?>
]], {sym=sym, xNames=xNames})}

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
