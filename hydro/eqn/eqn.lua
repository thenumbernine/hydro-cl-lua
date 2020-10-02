local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local Equation = class()

-- this is passed on to hydro/solver/calcDerivFV.cl
-- it has the effect of adding the connection terms Conn^k_jk u^I_,j (for the I'th conserved quantity u^I)
-- TODO get rid of this altogether
Equation.weightFluxByGridVolume = true

--[[
This flag determines whether the evR . lambda . evL is factored outside the other flux limiter computations.

I found this was especially useful in providing with the ideal MHD and using in the Roe flux computation
which, when using evR . lambda . evL, developed numerical errors that the flux didn't.

I think other equations were better performing without this, like Euler.
--]]
Equation.roeUseFluxFromCons = nil

-- whether to use the 'addSource' kernel
Equation.useSourceTerm = nil


-- used by bssnok-fd-num and bssnok-fd-sym
-- useful with spherical grids
-- (which no other eqn has attempted to implement yet)
-- signs = array of 1 or -1 based on the index parity wrt reflection about the boundary condition
-- see table III in 2017 Ruchlin
function Equation:getParityVars(...)
	local sign = {...}
	local vars = table()
	for _,var in ipairs(self.consStruct.vars) do
		if var.type == 'real' then
		elseif var.type == 'cplx' then
		elseif var.type == 'real3' then
			for i,xi in ipairs(xNames) do
				if sign[i] == -1 then
					vars:insert(var.name..'.'..xi)
				end
			end
		elseif var.type == 'cplx3' then
			for i,xi in ipairs(xNames) do
				if sign[i] == -1 then
					-- should reals be inserted in and - used?
					-- or should scalars be inserted and <scalar>_neg be used?
					vars:insert(var.name..'.'..xi..'.re')
					vars:insert(var.name..'.'..xi..'.im')
				end
			end	
		elseif var.type == 'sym3' then
			for ij,xij in ipairs(symNames) do
				local i,j = from6to3x3(ij)
				if sign[i] * sign[j] == -1 then
					vars:insert(var.name..'.'..xij)
				end
			end
		elseif var.type == '_3sym3' then
			for k,xk in ipairs(xNames) do
				for ij,xij in ipairs(symNames) do
					local i,j = from6to3x3(ij)
					if sign[i] * sign[j] * sign[k] == -1 then
						vars:insert(var.name..'.'..xk..'.'..xij)
					end
				end		
			end
		else
			error("don't know how to handle type "..var.type)
		end
	end
	return vars
end


--[[
args:
	solver = required

	make sure self.consStruct and self.primStruct is defined beforehand
--]]
function Equation:init(args)
	local solver = assert(args.solver)
	self.solver = solver
	local app = solver.app

	
	-- build consStruct and primStruct from self.consVars and .primVars (then erase them -- don't use them anymore)
	-- TODO think of a better way for the eqn to pass this info along, maybe a callback instead of a field?
	if not self.consStruct then
		assert(self.consVars)
		self.consStruct = Struct{
			solver = solver,
			name = 'cons_t',
			vars = self.consVars,
		}
	end
	
	-- make sure all math types are cdef'd
	self:cdefAllVarTypes(solver, self.consStruct.vars)

	self.consStruct:makeType()	-- create consStruct.typename
	self.cons_t = self.consStruct.typename
	-- don't use consVars anymore ... use consStruct.vars instead
	self.consVars = nil
	-- if you have multiple eqns then the class needs to keep the field
	--getmetatable(self).consVars = nil

	if not self.primStruct and self.primVars then
		self.primStruct = Struct{
			solver = solver,
			name = 'prim_t',
			vars = self.primVars,
		}
	end
	if self.primStruct then
		self.primStruct:makeType()
		self.prim_t = self.primStruct.typename
		-- don't use primVars anymore ... use primStruct.vars instead
		self.primVars = nil
		-- if you have multiple eqns then the class needs to keep the field
		--getmetatable(self).primVars = nil
	else
		--self.prim_t = self.cons_t
		-- or you could typedef this ...
		self.prim_t = app:uniqueName'prim_t'
	end

	if not self.eigenVars then
		self.eigenStruct = Struct{
			solver = solver, 
			name = 'eigen_t', 
			dontUnion = true,
			vars = {
				{name='unused', type='char'},
			},
		}
	else
		self.eigenStruct = Struct{solver=solver, name='eigen_t', vars=self.eigenVars}
	end
	self.eigenStruct:makeType()
	self.eigen_t = assert(self.eigenStruct.typename)

	self.consLR_t = app:uniqueName'consLR_t'
	self.waves_t = app:uniqueName'waves_t'
	
	local numReals
	if self.consStruct.vars then
		numReals = self.consStruct:countScalars()
		if self.primStruct then
			local numPrimReals = self.primStruct:countScalars()
			assert(numPrimReals <= numReals, "hmm, this is awkward")
		end
	end
	
	if not self.numStates then 
		self.numStates = numReals
		if not self.numStates then
			error("you either need to define numStates or consVars or consStruct")
		end
	else
		if numReals then
			assert(self.numStates == numReals, "found numStates="..self.numStates.." but found numReals="..numReals)
		end
	end
	-- default # waves is the # of states
	if not self.numWaves then self.numWaves = self.numStates end 
	
	-- how many states are integratable
	-- (put static states at the end of your cons_t structures)
	if not self.numIntStates then self.numIntStates = self.numStates end
	
	self.initCondNames = table.map(self.initConds, function(info) return info.name end)

	
	self.reflectVars = self.reflectVars or {}

	-- r min, for spherical coordinates
	-- what variables to mirror at sphere center
	-- 2013 Baumgarte et al, "Numerical Relativity in Spherical Polar Coordinates...", IIIB
	-- 2017 Ruchlin et al, section E.1
	self.reflectVars.sphereRMin = self.reflectVars.sphereRMin or {
		self:getParityVars(-1, 1, -1),
		{},
		{},
	} 

	-- theta min/max, for spherical coordinates
	self.reflectVars.sphereTheta = self.reflectVars.sphereTheta or {
		{},
		self:getParityVars(1, -1, -1),
		{},
	}

	-- phi min/max, for cylindrical or for spherical coordinates
	self.reflectVars.cylinderRMin = self.reflectVars.cylinderRMin or {
		self:getParityVars(-1, -1, 1),
		{},
		{},
	}

	-- x,y,z min/max for cartesian coordinates
	self.reflectVars.mirror = self.reflectVars.mirror or {
		self:getParityVars(-1, 1, 1),
		self:getParityVars(1, -1, 1),
		self:getParityVars(1, 1, -1),
	}
end

function Equation:cdefAllVarTypes(solver, vars)
	require 'hydro.code.safecdef'(
		solver.app.modules:getTypeHeader(
			table.mapi(vars, function(var,i,t)
				return true, var.type
			end):keys():unpack()
		)
	)
end

function Equation:addGuiVars(args)
	for _,arg in ipairs(args) do
		self:addGuiVar(arg)
	end
end

function Equation:addGuiVar(args)
	local vartype = args.type
	if not vartype then
		vartype = type(args.value)
	end
	local cl = require('hydro.guivar.'..vartype)

	-- no non-strings allowed as names, especially not numbers,
	-- because I'm going to map names into guiVars' keys, and I don't want to overwrite any integer-indexed guiVars
	assert(args.name and type(args.name) == 'string')

	local var = cl(args)
	self.guiVars:insert(var)
	self.guiVars[var.name] = var
end

-- TODO this creates and finalizes initCond.guiVars within its lifetime (via member function calls)
-- which means Equation subclasses cannot modify initCond guiVars
-- maybe that's alright.  i don't think I use that functionality.
-- Equation subclasses are modifying solver_t guiVars anyways.
function Equation:createInitState()
	self.guiVars = table()

	-- start with units
	self:addGuiVars{
		{name='meter', value=1},
		{name='second', value=1},
		{name='kilogram', value=1},
		{name='coulomb', value=1},
		{name='kelvin', value=1},
	}

	local mt = getmetatable(self)
	if mt.guiVars then
		self:addGuiVars(mt.guiVars)
	end
	
	-- first create the init state
	assert(self.initConds, "expected Eqn.initConds")
	self.initCond = self.initConds[self.solver.initCondIndex](self.solver, self.solver.initCondArgs)
	assert(self.initCond, "couldn't find initCond "..self.solver.initCondIndex)
	self.initCond:createInitStruct(self.solver)
	self.initCond:finalizeInitStruct(self.solver)
	
	-- should ops add vars to initCond_t or solver_t?
	-- or should there be a new eqn_t?
	-- I would like init cond stuff in initCond_t so changing the init cond only recompiles initCond.cl
	-- so let's continue to put op guiVars in solver_t
	for _,op in ipairs(self.solver.ops) do
		if op.guiVars then
			self:addGuiVars(op.guiVars)
		end
	end
end

-- shorthand
function Equation:template(code, args)
	if args then
		args = table(self:getEnv(), args)
	else
		args = self:getEnv()
	end
	return template(code, args)
end

-- add to self.solver.modules, or add to self.modules and have solver add later?
function Equation:initCodeModules()
	local solver = self.solver

	assert(self.consStruct)
	solver.modules:add{
		name = 'eqn.cons_t',
		structs = {self.consStruct},
	}
	-- boundary, initCond, solver ... everyone needs this
	solver.sharedModulesEnabled['eqn.cons_t'] = true

	solver.modules:add{
		name = 'eqn.prim_t',
		structs = {self.primStruct},
		typecode = not self.primStruct and ('typedef '..self.cons_t..' '..self.prim_t..';') or nil,
	}
	
	solver.modules:add{
		name = 'eqn.waves_t',
		depends = {'real'},
		typecode = self:template[[
typedef union { 
	real ptr[<?=eqn.numWaves?>]; 
} <?=eqn.waves_t?>;
]],
	}
	
	assert(self.eigenStruct)
	solver.modules:add{
		name = 'eqn.eigen_t',
		structs = {self.eigenStruct},
	}

	-- only require this if we're a fvsolver

	-- parallel propagate autogen code 
	-- only used for finite-volume solvers
	-- also NOTICE there seems to be a bug where the CL compiler stalls when compiling the parallel-propagate code with bssnok-fd
	-- so hopefully the module system can help that out
	--
	-- TODO hmm, can't add the module unless it's being used
	--  because some that don't use it (bssnok-fd) use custom suffixes on var names
	if require 'hydro.solver.fvsolver'.is(solver) 
	-- TODO only if it's a mesh solver using a flux integrator ... which is currently all mesh solvers
	or require 'hydro.solver.meshsolver'.is(solver) 
	then
		local degreeForType = {
			real = 0,
			real3 = 1,
			sym3 = 2,
			real3x3 = 2,
			_3sym3 = 3,
			real3x3x3 = 3,
		}
		solver.modules:add{
			name = 'eqn.cons_parallelPropagate',
			depends = {
				'eqn.cons_t',
				'coord_parallelPropagate',
			},
			code = self:template([[<? 
for side=0,solver.dim-1 do
	if coord.vectorComponent == 'cartesian'
	or require 'hydro.coord.cartesian'.is(coord) 
	then
?>#define cons_parallelPropagate<?=side?>(U, x, dx) (U)
<?	else
?><?=eqn.cons_t?> cons_parallelPropagate<?=side?>(<?=eqn.cons_t?> U, real3 x, real dx) {
<?		for _,var in ipairs(eqn.consStruct.vars) do
			local variance = var.variance or var.name:match'_([^_]*)$' or ''
			local degree = degreeForType[var.type]
			if #variance ~= degree  then
				error("variable variance "..('%q'):format(variance).." does not match variable type "..var.type.." degree "..degree)
			end
--print(var.name, var.type, variance, degree)			
			if variance == '' then
			elseif variance == 'u' then
-- P_u^a T^u
?>	U.<?=var.name?> = coord_parallelPropagateU<?=side?>(U.<?=var.name?>, x, dx);
<?
			elseif variance == 'l' then
-- T_u (P^-T)_a^u
?>	U.<?=var.name?> = coord_parallelPropagateL<?=side?>(U.<?=var.name?>, x, dx);
<?
			elseif variance == 'll' then
-- (P^-T)_a^u (P^-T)_b^v T_uv
?>				real3x3 t = real3x3_from_<?=var.type?>(U.<?=var.name?>);
				t.x = coord_parallelPropagateL<?=side?>(t.x, x, dx);
				t.y = coord_parallelPropagateL<?=side?>(t.y, x, dx);
				t.z = coord_parallelPropagateL<?=side?>(t.z, x, dx);
				t = real3x3_transpose(t);
				t.x = coord_parallelPropagateL<?=side?>(t.x, x, dx);
				t.y = coord_parallelPropagateL<?=side?>(t.y, x, dx);
				t.z = coord_parallelPropagateL<?=side?>(t.z, x, dx);
				U.<?=var.name?> = <?=var.type?>_from_real3x3(t);
<?			elseif variance == 'lll' then
?>				real3x3x3 t = real3x3x3_from_<?=var.type?>(U.<?=var.name?>);
				real3 tmp;
<?				local is = table()
				for e=1,3 do
					for i,xi in ipairs(xNames) do
						is[e] = xi
						for j,xj in ipairs(xNames) do
							is[e%3+1] = xj
							for k,xk in ipairs(xNames) do
								is[(e+1)%3+1] = xk
?>					tmp.<?=xk?> = t.<?=is:concat'.'?>;
<?							end
?>					tmp = coord_parallelPropagateL<?=side?>(tmp, x, dx);
<?							for k,xk in ipairs(xNames) do
								is[(e+1)%3+1] = xk
?>					t.<?=is:concat'.'?> = tmp.<?=xk?>;
<?							end
						end
					end
				end
?>				U.<?=var.name?> = real3x3x3_from_<?=var.type?>(t);
<?			else
				error("don't know how to handle variance for "..('%q'):format(variance))
			end
		end
?>	return U;
}
<?	end
end
?>]], 		{
				coord = solver.coord,
				degreeForType = degreeForType,
			}),
		}
	end

	-- Put this here or in SolverBase?
	-- This calls initCond:getCodePrefix, which adds to eqn.guiVars.
	-- That means call this after eqn:createInitState (solverbase.lua:524 in :initObjs)
	-- eqn:initCodeModules is called after ... hmm
	self.initCond:initCodeModules(solver)

	self:initCodeModuleCommon()	-- eqn.common

	-- init eqn.prim-cons
	-- prim-cons should have access to all ... prefix stuff?
	-- but initstate has access to it
	self:initCodeModulePrimCons()

	-- these are only the compile-time gui vars
	-- the runtime ones are stored in solver_t
	solver.modules:add{
		name = 'eqn.guiVars.compileTime',
		headercode = table.mapi(self.guiVars or {}, function(var,i,t) 
			return (var.compileTime and var:getCode() or nil), #t+1
		end):concat'\n',
	}

	self:initCodeModuleSolver()
	self:initCodeModuleCalcDT()

	self:initCodeModule_fluxFromCons()
end

function Equation:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'eqn.solvercode',	-- eigen_fluxTransform, eigen_forCell
			'eqn.cons_t',
			'solver.solver_t',
			'normal_t',		-- normal_t
		},
		code = self:template[[
<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	return eigen_fluxTransform(solver, eigen_forCell(solver, U, x, n), U, x, n);
}
]],
	}
end

function Equation:getModuleDependsCommon() end	-- eqn.common, used by init and solver

function Equation:initCodeModuleCommon()
	-- TODO don't even use this, just subclass initModule
	self.solver.modules:add{
		name = 'eqn.common',
		depends = self:getModuleDependsCommon(),
		code = self.getCommonFuncCode and self:getCommonFuncCode() or nil,
	}
end

function Equation:initCodeModuleSolver()
	self.solver.modules:add{
		name = 'eqn.solvercode',
		depends = table{
			'eqn.cons_t',
			'eqn.prim_t',
			'eqn.waves_t',
			'eqn.eigen_t',
			'eqn.guiVars.compileTime',
			'coordLenSq',
			-- commonly used by most solvers ... and anything that uses kernels for that matter
			'INDEX', 'INDEXV', 'OOB', 'SETBOUNDS', 'SETBOUNDS_NOGHOST',
		}:append(self:getModuleDependsSolver()),
		code = self:template(file[self.solverCodeFile]),
	}
end

function Equation:getModuleDependsSolver() end	-- eqn.solver, used by solver

function Equation:getModuleDependsApplyInitCond() 
	return {
		'solver.solver_t',
		'eqn.cons_t',
	}
end	-- get'd in hydro/init/init.lua initCodeModules

-- Really really used by maxwell, glm-maxwell, and other things that vary their scalar type between real and cplx.  but it fits here just as well.
function Equation:getEnv()
	return {
		-- most have this
		eqn = self,
		solver = self.solver,
		coord = self.solver.coord,
	
		-- common 
		xNames = xNames,
		symNames = symNames,
		from3x3to6 = from3x3to6,
		from6to3x3 = from6to3x3,
		sym = sym,
	}
end

-- this only goes to hydro/init/init.lua
-- and this is influenced by the initCond object
-- changing initCond should only change this and not the solver program
function Equation:getInitCondCode()
	assert(self.initCondCode, "expected solver.eqn.initCondCode")
	return self:template(self.initCondCode, {
		code = self.initCond:getInitCondCode(self.solver),
	})
end

Equation.displayVarCodeUsesPrims = false
function Equation:getDisplayVarCodePrefix()
	return self:template[[
	const global <?=eqn.cons_t?>* U = buf + index;
<? if eqn.displayVarCodeUsesPrims then 
?>	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
<? end 
?>]]
end

-- accepts a list of struct var info {name=..., [type=..., units=...]}
-- returns a list of display var construction info
function Equation:getDisplayVarsForStructVars(structVarInfos, ptrName)
	ptrName = ptrName or 'U'
	ptrName = ptrName .. '->'
	local displayVarInfos = table()
	for _,structVarInfo in ipairs(structVarInfos) do
		
		local function addvar(name, varname, vartype)
			local assignvar = 'value.v'..vartype
			displayVarInfos:insert{
				name = name,
				code = assignvar..' = '..ptrName..varname..';', 
				type = vartype,
				units = structVarInfo.units,
				
				-- if a display var has 'field' set then use a predefined calcDisplayVar function to just read the field directly (without any computations required)
				-- ... unless it has units too ... in which case ... I'll be scaling the units
				-- ... of course I could do the scaling after reading the value ...
				field = varname,
			}
		end

		local varname = structVarInfo.name
		local vartype = structVarInfo.type
		if vartype == '_3sym3' then
			for i,xi in ipairs(xNames) do
				addvar(varname..' '..xi, varname..'.'..xi, 'sym3')
			end
		else
			addvar(varname, varname, vartype)
		end
	end
	return displayVarInfos	
end

function Equation:getDisplayVars()
	return self:getDisplayVarsForStructVars(self.consStruct.vars)
end

-- I would make this a component, but then the component code would need to access the entire previous buffer befure making its computation
-- that could be done in GLSL using the dx operators ... maybe I'll look into that later ...
-- TODO use the automatic arbitrary finite difference generator in bssnok
function Equation:createDivDisplayVar(args)
	if require 'hydro.solver.meshsolver'.is(self.solver) then return end
	
	local field = assert(args.field)
	local getField = args.getField
	local scalar = args.scalar or 'real'
	local units = args.units
	
	return {
		name = 'div '..field, 
		code = self:template([[
	if (OOB(1,1)) {
		value.v<?=scalar?> = 0./0.;
	} else {
		<?=scalar?> v = <?=scalar?>_zero;
		<? for j=0,solver.dim-1 do ?>{
			global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
			global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
			v = <?=scalar?>_add(v, <?=scalar?>_real_mul(
				<?=scalar?>_sub(
					<?=getField('Ujp', j)?>,
					<?=getField('Ujm', j)?>
				), .5 / solver->grid_dx.s<?=j?>));
		}<? end ?>
		value.v<?=scalar?> = v;
	}
]], 	{
			getField = getField or function(U, j)
				return U..'->'..field..'.s'..j	
			end,
			scalar = scalar,
		}),
		type = scalar, 
		units = units,
	}
end

-- curl = [,x ,y ,z] [v.x, v.y, v.z]
-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
-- TODO use the automatic arbitrary finite difference generator in bssnok
function Equation:createCurlDisplayVar(args)
	if require 'hydro.solver.meshsolver'.is(self.solver) then return end

	local field = assert(args.field)
	local getField = args.getField
	local units = args.units

	-- k is 0-based
	local function curl(k, result)
		local i = (k+1)%3	-- 0-based
		local j = (i+1)%3	-- 0-based
		return {
			name = 'curl '..field..' '..xNames[k+1],
			code = self:template([[
	if (OOB(1,1)) {
		<?=result?> = 0./0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = U - solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;

		//TODO incorporate metric

		real vim_j = <?=getField('Uim', j)?>;
		real vip_j = <?=getField('Uip', j)?>;
		
		real vjm_i = <?=getField('Ujm', i)?>;
		real vjp_i = <?=getField('Ujp', i)?>;
		
		<?=result?> = (vjp_i - vjm_i) / (2. * solver->grid_dx.s<?=i?>)
					- (vip_j - vim_j) / (2. * solver->grid_dx.s<?=j?>);
	}
]], 	{
			i = i,
			j = j,
			result = result,
			getField = getField or function(U, j)
				return U..'->'..field..'.s'..j
			end,
		})}
	end

	-- TODO just use LIC and the mag will be ... the abs of this
	if self.solver.dim == 2 then
		local var = curl(2,'value.vreal')
		var.name = 'curl '..field
		return var
	elseif self.solver.dim == 3 then
		local v = xNames:mapi(function(xi,i)
			return curl(i-1,'value.vreal3.'..xi)
		end)
		return {
			name = 'curl '..field,
			code = self:template([[
	<? for i,vi in ipairs(v) do ?>{
		<?=vi.code?>
	}<? end ?>
]], 			{
					v = v,
				}), 
			type = 'real3',
			units = units,
		}
	end
end

function Equation:getEigenDisplayVars()
	-- use the automatic codegen for display vars
	if not self.eigenVars then return end
	return self:getDisplayVarsForStructVars(self.eigenVars, '(&eig)')
end

function Equation:eigenWaveCodePrefix(n, eig, x)
	return ''
end
function Equation:eigenWaveCode(n, eig, x, waveIndex)
	return '\n#error :eigenWaveCode() not implemented'
end

function Equation:consWaveCodePrefix(n, U, x)
	return '\n#error :consWaveCodePrefix() not implemented'
end
function Equation:consWaveCode(n, U, x)
	return '\n#error :consWaveCode() not implemented'
end

-- default implementation -- the first is the min and the last is the max
-- however some don't do this, like GLM
function Equation:eigenMinWaveCode(n, eig, x)
	return self:eigenWaveCode(n, eig, x, 0)
end
function Equation:eigenMaxWaveCode(n, eig, x)
	return self:eigenWaveCode(n, eig, x, self.numWaves-1)
end
function Equation:consMinWaveCode(n, U, x)
	return self:consWaveCode(n, U, x, 0)
end
function Equation:consMaxWaveCode(n, U, x)
	return self:consWaveCode(n, U, x, self.numWaves-1)
end

-- Whether the eqn has its own calcDT.  Otherwise hydro/eqn/cl/calcDT.cl is used. 
Equation.hasCalcDTCode = nil

function Equation:initCodeModuleCalcDT()
	-- hmm can't do this anymore since even if calcDT is in cl code it won't be in the module system
	if self.hasCalcDTCode then return end
	
	self.solver.modules:add{
		name = 'eqn.calcDT',
		depends = {
			'eqn.cons_t',
			'eqn.prim_t',
			'eqn.waves_t',
			'eqn.eigen_t',
			'normal_t',
			-- so these dependencies are going to vary based on the eigen code of each eqn 
			'eqn.common',		-- used by eqn/wave
			'eqn.prim-cons',	-- used by eqn/shallow-water
		},
		code = self:template(file['hydro/eqn/cl/calcDT.cl']),
	}
end

--[[
Default code for the following:
	primFromCons
	consFromPrim
	apply_dU_dW : prim_t -> cons_t 
	apply_dW_dU : cons_t -> prim_t

The default assumes prim_t == cons_t and this transformation is identity
--]]
function Equation:initCodeModulePrimCons()
	assert(not self.primStruct, "if you're using the default prim<->cons code then you shouldn't have any primStruct")

	self.solver.modules:add{
		name = 'eqn.prim-cons',
		depends = {
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
		},
		code = self:template[[
#define primFromCons(solver, U, x)	U
/*
<?=eqn.prim_t?> primFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U, 
	real3 x
) { 
	return U; 
}
*/

#define consFromPrim(solver, W, x)	W
/*
<?=eqn.cons_t?> consFromPrim(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> W, 
	real3 x
) { 
	return W; 
}
*/
]],
	}

	-- only used by PLM
	self.solver.modules:add{
		name = 'eqn.dU-dW',
		depends = {
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
		},
		code = self:template[[
/*
WA = W components that make up the jacobian matrix
W = input vector
x = coordinate location
returns output vector
*/
#define apply_dU_dW(solver, WA, W, x)	W
/*
<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) { 
	return W; 
}
*/

/*
WA = W components that make up the jacobian matrix
U = input vector
x = coordinate location
returns output vector
*/
#define apply_dW_dU(solver, WA, U, x)	U
/*
<?=eqn.prim_t?> apply_dW_dU(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.cons_t?> U, 
	real3 x
) { 
	return U; 
}
*/
]],
	}
end

return Equation
