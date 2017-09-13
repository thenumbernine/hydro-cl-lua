--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
With some hints from a few of the Z4 papers
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local NumRelEqn = require 'eqn.numrel'
local symmath = require 'symmath'
local makeStruct = require 'eqn.makestruct'
require 'common'(_G)

local Z4 = class(NumRelEqn)
Z4.name = 'Z4'

local fluxVars = table{
	{a = 'real3'},
	{d = '_3sym3'},
	{K = 'sym3'},
	{Theta = 'real'},
	{Z = 'real3'},
}

Z4.consVars = table{
	{alpha = 'real'},
	{gamma = 'sym3'},
}:append(fluxVars)

Z4.numWaves = makeStruct.countReals(fluxVars)

Z4.hasCalcDT = true
Z4.hasEigenCode = true
Z4.useSourceTerm = true

function Z4:getCodePrefix()
	local lines = table()
	lines:insert(template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1;
	U->gamma = _sym3(1,0,0,1,0,1);
	U->a = _real3(0,0,0);
	U->d[0] = _sym3(0,0,0,0,0,0);
	U->d[1] = _sym3(0,0,0,0,0,0);
	U->d[2] = _sym3(0,0,0,0,0,0);
	U->K = _sym3(0,0,0,0,0,0);
	U->Theta = 0;
	U->Z = _real3(0,0,0);
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
		
			-- TODO Theta and Z_i, or can I just leave them as zero?
		)
	end))

	return lines:concat()
end

function Z4:getInitStateCode()
	local lines = table{
		template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

	setFlatSpace(U);
	
]], 	{
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

	lines:insert'	U->Theta = 0;'
	lines:insert'	U->Z = _real3(0,0,0);'

	lines:insert'}'
	
	local code = lines:concat'\n'
	return code
end

function Z4:getSolverCode()
	return template(file['eqn/z4.cl'], {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
		sym = sym,
	})
end

function Z4:getDisplayVars()
	local vars = Z4.super.getDisplayVars(self)
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

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
		{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*valuevec = real3_scale(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]], type='real3'},
	}
	
	return vars
end

Z4.eigenVars = table{
	{alpha = 'real'},
	{sqrt_f = 'real'},
	{gammaU = 'sym3'},
	{sqrt_gammaUjj = 'real3'},
}

function Z4:fillRandom(epsilon)
	local ptr = ADM_BonaMasso_3D.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma.xx = ptr[i].gamma.xx + 1
		ptr[i].gamma.yy = ptr[i].gamma.yy + 1
		ptr[i].gamma.zz = ptr[i].gamma.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return Z4 
