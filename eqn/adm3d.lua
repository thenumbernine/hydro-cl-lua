--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local NumRelEqn = require 'eqn.numrel'
local symmath = require 'symmath'
local makeStruct = require 'eqn.makestruct'

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

local fluxVars = table{
	{a = 'real3'},
	{d = '_3sym3'},
	{K = 'sym3'},
	{V = 'real3'},
}

ADM_BonaMasso_3D.consVars = table{
	{alpha = 'real'},
	{gamma = 'sym3'},
}:append(fluxVars)

-- skip alpha and gamma
ADM_BonaMasso_3D.numWaves = makeStruct.countReals(fluxVars)
assert(ADM_BonaMasso_3D.numWaves == 30)

ADM_BonaMasso_3D.hasCalcDT = true
ADM_BonaMasso_3D.hasEigenCode = true
ADM_BonaMasso_3D.useSourceTerm = true
ADM_BonaMasso_3D.useConstrainU = true

function ADM_BonaMasso_3D:createInitState()
	ADM_BonaMasso_3D.super.createInitState(self)
	self:addGuiVar{
		type = 'combo',
		name = 'constrain V',
		options = {
			'none',
			'replace V',
			'average',	-- TODO add averaging weights, from 100% V (which works) to 100% d (which doesn't yet work)
		}
	}
end

function ADM_BonaMasso_3D:getCodePrefix()
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
	
	lines:insert(self.initState:getCodePrefix(self.solver, function(exprs, vars, args)
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
	local vars = table()
	for _,var in ipairs(self.consVars) do
		local varname, vartype = next(var)
		
		if vartype == 'real' then
			vars:insert{[varname] = '*value = U->'..varname..';'}
		elseif vartype == 'real3' then
			-- use vec display vars, it'll auto-create scalars
		elseif vartype == 'sym3' then
			for ij,xij in ipairs(symNames) do
				vars:insert{[varname..'_'..xij] = '*value = U->'..varname..'.'..xij..';'}
			end
		elseif vartype == '_3sym3' then
			for i,xi in ipairs(xNames) do
				for jk,xjk in ipairs(symNames) do
					vars:insert{[varname..'_'..xi..xjk] = '*value = U->'..varname..'['..(i-1)..'].'..xjk..';'}
				end
			end
		end
	end
	
	vars:append{
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

	return vars
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

	for _,var in ipairs(self.consVars) do
		local varname, vartype = next(var)
		
		if vartype == 'real3' then
			vars:insert{[varname] = 'valuevec = U->'..varname..';'}
		end
	end

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
	vars:insert{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	valuevec = real3_scale(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]]}
	
	local function addSym3(field)
		for i,xi in ipairs(xNames) do
			vars:insert{[field..'_'..xi] = 'valuevec = sym3_'..xi..'(U->'..field..');'}
		end
	end
	addSym3'gamma'
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

local sym3_delta = {1,0,0,1,0,1}

local ffi = require 'ffi'
local function crand() return 2 * math.random() - 1 end
function ADM_BonaMasso_3D:fillRandom(epsilon)
	local solver = self.solver
	local ptr = ffi.new(self.cons_t..'[?]', solver.volume)
	for i=0,solver.volume-1 do
		ptr[i].alpha = epsilon * crand()
		for jk=0,5 do
			ptr[i].gamma.s[jk] = sym3_delta[jk+1] + epsilon * crand()
		end
		for j=0,2 do
			ptr[i].a.s[j] = epsilon * crand()
		end
		for j=0,2 do
			for mn=0,5 do
				ptr[i].d[j].s[mn] = epsilon * crand()
			end
		end
		for jk=0,5 do
			ptr[i].K.s[jk] = epsilon * crand()
		end
		for j=0,2 do
			ptr[i].V.s[j] = epsilon * crand()
		end
	end
	solver.UBufObj:fromCPU(ptr)
end

return ADM_BonaMasso_3D
