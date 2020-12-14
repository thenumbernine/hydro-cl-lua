local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local Struct = require 'hydro.code.struct'
local makePartials = require 'hydro.eqn.makepartial'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local Equation = class()

-- this is passed on to hydro/solver/fvsolver.cl
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

-- singleton
Equation.parityVarsGetters = table{
	real = function(sign, parityVars, field) end,
	cplx = function(sign, parityVars, field) end,
	real3 = function(sign, parityVars, field)
		for i,xi in ipairs(xNames) do
			if sign[i] == -1 then
				parityVars:insert(field..'.'..xi)
			end
		end
	end,
	cplx3 = function(sign, parityVars, field)
		for i,xi in ipairs(xNames) do
			if sign[i] == -1 then
				-- should reals be inserted in and - used?
				-- or should scalars be inserted and <scalar>_neg be used?
				parityVars:insert(field..'.'..xi..'.re')
				parityVars:insert(field..'.'..xi..'.im')
			end
		end	
	end,
	sym3 = function(sign, parityVars, field)
		for ij,xij in ipairs(symNames) do
			local i,j = from6to3x3(ij)
			if sign[i] * sign[j] == -1 then
				parityVars:insert(field..'.'..xij)
			end
		end
	end,
	real3x3 = function(sign, parityVars, field)
		for i,xi in ipairs(xNames) do
			for j,xj in ipairs(xNames) do
				if sign[i] * sign[j] == -1 then
					parityVars:insert(field..'.'..xi..'.'..xj)
				end
			end
		end	
	end,
	_3sym3 = function(sign, parityVars, field)
		for k,xk in ipairs(xNames) do
			for ij,xij in ipairs(symNames) do
				local i,j = from6to3x3(ij)
				if sign[i] * sign[j] * sign[k] == -1 then
					parityVars:insert(field..'.'..xk..'.'..xij)
				end
			end		
		end
	end,
}

local function getParityVars(structVars, sign, parityVars, field)
	field = field and (field..'.') or ''
	for _,var in ipairs(structVars) do
		local getter = Equation.parityVarsGetters[var.type]
		if not getter then
			error("don't know how to handle type "..var.type)
			-- TODO if it's a struct made up of these types, then handle it recursively
		end
		getter(sign, parityVars, field..var.name)
	end
end

-- used by bssnok-fd-num and bssnok-fd-sym
-- useful with spherical grids
-- (which no other eqn has attempted to implement yet)
-- signs = array of 1 or -1 based on the index parity wrt reflection about the boundary condition
-- see table III in 2017 Ruchlin
function Equation:getParityVars(...)
	local sign = {...}
	local parityVars = table()
	getParityVars(self.consStruct.vars, sign, parityVars)
	return parityVars
end


--[[
args:
	solver = required

	make sure self.consStruct and self.primStruct is defined beforehand
--]]
function Equation:init(args)
	local mt = getmetatable(self)
	setmetatable(self, nil)
	self.uid = assert(tostring(self):match'table: 0x(.*)')
	setmetatable(self, mt)
	self.symbolPrefix = self.name..'_'..self.uid..'_'

	-- [[ make C symbols unique to the eqn class.
	-- TODO later,unique to the eqn object, in case I make a composite of two like classes
	-- like euler + euler to do two-fluid simulations
	-- TODO this should be modular
	local symbolFields = self:getSymbolFields()

	-- TODO add the typedef names into this
	-- TODO don't use uniqueName() instead just use the object uid suffix (and use it for functions too)
	self.symbols = {}
	for _,field in ipairs(symbolFields) do
		self.symbols[field] = self.symbolPrefix..field
	end
--]]

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

-- make sure all math types are cdef'd,
-- for sizeof's for calculating the union ptr size 
-- when we makeType the consStruct
-- we can't do this after makeType unless we also put it after initCodeModule
-- (but eqn:initCodeModule is called after the consStruct type is defined)	
	self:cdefAllVarTypes(solver, self.consStruct.vars)

	self.consStruct:makeType()	-- create consStruct.typename
	-- TODO replace the cdef uniqueName with a unique eqn object name
	self.symbols.cons_t = self.consStruct.typename
	
	-- don't use consVars anymore ... use consStruct.vars instead
	self.consVars = nil
	-- if you have multiple eqns then the class needs to keep the field
	--getmetatable(self).consVars = nil
	
	self.consStruct.eqn = self	-- hack
	solver.structForType[self.consStruct.typename] = self.consStruct	-- hack


	if not self.primStruct and self.primVars then
		self.primStruct = Struct{
			solver = solver,
			name = 'prim_t',
			vars = self.primVars,
		}
	end
	if self.primStruct then

		-- cdef the used types before :makeType the Struct.  see consStruct for more notes.
		self:cdefAllVarTypes(solver, self.primStruct.vars)
		
		local res, err = xpcall(function()
			self.primStruct:makeType()
		end, function(err)
			return "eqn "..self.name.." primStruct:makeType() failed for type "..self.primStruct.name..'\n'
				..require 'ext.tolua'(table(
					self.primStruct,
					{app=false}))..'\n'
				..tostring(err)..'\n'
				..debug.traceback()
		end)
		if not res then error(err) end

		-- TODO replace the cdef uniqueName with a unique eqn object name
		self.symbols.prim_t = self.primStruct.typename
		
		-- don't use primVars anymore ... use primStruct.vars instead
		self.primVars = nil
		-- if you have multiple eqns then the class needs to keep the field
		--getmetatable(self).primVars = nil
		
		self.primStruct.eqn = self	-- hack
		solver.structForType[self.primStruct.typename] = self.primStruct
	else
		--self.symbols.prim_t = self.symbols.cons_t
		-- or you could typedef this ...
		-- TODO replace the cdef uniqueName with a unique eqn object name
		self.symbols.prim_t = app:uniqueName'prim_t'
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
	self.symbols.eigen_t = assert(self.eigenStruct.typename)
	
	self.eigenStruct.eqn = self	-- hack
	solver.structForType[self.eigenStruct.typename] = self.eigenStruct

	self.symbols.consLR_t = app:uniqueName'consLR_t'
	

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


	if not self.wavesVars then
		self.wavesVars = range(self.numWaves):mapi(function(i)
			return {name='wave'..(i-1), type='real'}
		end)
	end
	self.wavesStruct = Struct{
		solver = solver,
		name = 'waves_t',
		vars = self.wavesVars,
	}
	
	self.wavesStruct:makeType()
	self.symbols.waves_t = assert(self.wavesStruct.typename)

	self.wavesStruct.eqn = self	-- hack
	solver.structForType[self.wavesStruct.typename] = self.wavesStruct



	self.initCondNames = table.mapi(
		assert(self.initConds, "you forgot to specify your initConds in your hydro/eqn/* file"),
		function(info) return info.name end)

	
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

	-- now add our own prim_t, cons_t to the parityVarsGetters for recursive application
	self.parityVarsGetters[assert(self.symbols.cons_t)] = function(sign, parityVars, field)
		getParityVars(self.consStruct.vars, sign, parityVars, field)
	end
	self.parityVarsGetters[assert(self.symbols.prim_t)] = function(sign, vars, var)
		getParityVars(self.primStruct.vars, sign, parityVars, field)
	end
end

function Equation:getSymbolFields()
	return table{
		-- functions:
		'primFromCons',
		'consFromPrim',
		'apply_dU_dW',
		'apply_dW_dU',
		'fluxFromCons',
		'calcCellMinMaxEigenvalues',
		'eigen_forCell',
		'eigen_forInterface',
		'eigen_leftTransform',
		'eigen_rightTransform',
		'eigen_fluxTransform',
		'cons_parallelPropagate',
		'applyInitCondCell',
		'calcDTCell',
	
		-- kernels:
		'applyInitCond',	-- TODO don't uniquely identify this, only do the Cell version 
		'calcDT',			-- TODO don't uniquely identify this, only do the Cell version 
		'initDerivs',
		'addSource',
		'constrainU',
	
		-- placeholder modules for dependencies
		'eqn_guiVars_compileTime',	-- module of code for compile-time #defines of gui vars
		'eqn_waveCode_depends',		-- module of dependencies used by the eigen/cons wave code
		'eqn_common',				-- module of functions that are commonly used ... not required.
	}
end

function Equation:cdefAllVarTypes(solver, vars)
	-- TODO not just math, but also cons_t
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

	self:createInitState_createInitState()
	
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

function Equation:createInitState_createInitState()
	-- first create the init state
	assert(self.initConds, "expected Eqn.initConds")
	self.initCond = self.initConds[self.solver.initCondIndex](self.solver, self.solver.initCondArgs)
	assert(self.initCond, "couldn't find initCond "..self.solver.initCondIndex)
	self.initCond:createInitStruct(self.solver)
	self.initCond:finalizeInitStruct(self.solver)
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

function Equation:getEnv()
	local env = {
		-- most have this
		eqn = self,
		solver = self.solver,
		coord = self.solver.coord,
		initCond = self.initCond,
	
		-- type names
		solver_t = self.solver.solver_t,
		initCond_t = self.solver.initCond_t,
		cell_t = self.solver.coord.cell_t,
		face_t = self.solver.coord.face_t,
		
		-- macro numbers 
		numWaves = self.numWaves,

		-- common 
		xNames = xNames,
		symNames = symNames,
		from3x3to6 = from3x3to6,
		from6to3x3 = from6to3x3,
		sym = sym,
	
		-- really only used by applyInitCond
		initCode = function() 
			-- calls initCond:getInitCondCode
			return self.initCond:getInitCondCode(self.solver)
		end,
	}

	for k,v in pairs(self.symbols) do
		env[k] = v
	end

	return env
end

-- add to self.solver.modules, or add to self.modules and have solver add later?
function Equation:initCodeModules()
	local solver = self.solver

	self:initCodeModule_cons_prim_eigen_waves()

	self:initCodeModule_cons_parallelPropagate()

	-- Put this here or in SolverBase?
	-- This calls initCond:getCodePrefix, which adds to eqn.guiVars.
	-- That means call this after eqn:createInitState (solverbase.lua:524 in :initObjs)
	-- eqn:initCodeModules is called after ... hmm
	self.initCond:initCodeModules(solver)

	-- init primFromCons and consFromPrim
	-- prim-cons should have access to all ... prefix stuff?
	-- but initstate has access to it
	self:initCodeModule_consFromPrim_primFromCons()

	-- these are only the compile-time gui vars
	-- the runtime ones are stored in solver_t
	solver.modules:add{
		name = self.symbols.eqn_guiVars_compileTime,
		headercode = table.mapi(self.guiVars or {}, function(var,i,t) 
			return (var.compileTime and var:getCode() or nil), #t+1
		end):concat'\n',
	}

	-- this contains calcDTCell, which varies per-equation
	self:initCodeModule_calcDTCell()
	-- and here is calcDT, which is always the same
	solver.modules:add{
		name = self.symbols.calcDT,
		depends = {self.symbols.calcDTCell},
		code = self:template[[
kernel void <?=calcDT?>(
	constant <?=solver_t?> const * const solver,
	global real * const dtBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf<?
if require "hydro.solver.meshsolver".is(solver) then 
?>,
	global <?=face_t?> const * const faces,
	global int const * const cellFaceIndexes<?
end
?>
) {
	SETBOUNDS(0,0);
	
	//write inf to boundary cells
	//TODO why not write inf to boundary cells upon init,
	// and then give this the domain SETBOUNDS_NOGHOST?
	//that would work except that it is writing to reduceBuf, which is reused for any reduce operation
	// like display var min/max ranges
	global real * const dt = dtBuf + index;
	*dt = INFINITY;
	
	if (OOB(numGhost,numGhost)) return;
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	<?=calcDTCell?>(
		dt,
		solver,
		U,
		cell<?
if require "hydro.solver.meshsolver".is(solver) then 
?>,
		faces,
		cellFaceIndexes<?
end
?>
	);
}
]],
	}

	self:initCodeModule_fluxFromCons()
	self:initCodeModule_waveCode()

	solver.modules:addFromMarkup{
		code = self:template(file[self.solverCodeFile]),
		onAdd = function(args)
			-- special case for applyInitCond ...
			if args.name == self.symbols.applyInitCond then
				args.depends:append(self.initCond:getBaseDepends(solver))
				args.depends:append(self.initCond.depends)
				-- only used by hydro/eqn/bssnok-fd.lua:
				if self.getModuleDependsApplyInitCond then
					args.depends:append(self:getModuleDependsApplyInitCond())
				end
			end
		end,
	}

	solver.modules:add{
		name = self.symbols.applyInitCond,
		depends = {
			self.symbols.cell_t,
			self.symbols.initCond_t,
			self.symbols.applyInitCondCell,
			'SETBOUNDS',
		},
		code = self:template[[
kernel void <?=applyInitCond?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?>* UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	global <?=cons_t?> * const U = UBuf + index;
	global <?=cons_t?> const * const cell = cellBuf + index;
	<?=applyInitCondCell?>(solver, initCond, U, cell);
}
]],
	}
end

function Equation:initCodeModule_cons_parallelPropagate()
	local solver = self.solver

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
		
		for _,var in ipairs(self.consStruct.vars) do
			-- guess the variance if it isn't specified
			if not var.variance then
				var.variance = var.name:match'_([^_]*)$' or ''
			end
			-- also assert it is the right one for the struct
			local degree = degreeForType[var.type] or 0
			if #var.variance ~= degree  then
				local tolua = require 'ext.tolua'
				error("variable "..tolua(var.name).." variance "..tolua(var.variance).." does not match variable type "..tolua(var.type).." degree "..tolua(degree))
			end
		end

--[[			
cons_parallelPropagate is a macro
First argument is the name of the resulting local var
that will hold a ptr to the results.

In the event that propagation is identity, a pointer with name 'resultName' pointing to the original variable is created. 

In the event that a transformation is necessary, then a temp var is created, and 'resultName' points to it.
--]]
		solver.modules:add{
			name = self.symbols.cons_parallelPropagate,
			depends = table{
				self.symbols.cons_t,
				'coord_parallelPropagate',
			}:append(
				-- rank-2 always use real3x3 for transformation
				self.consStruct.vars:find(nil, function(var)
					return degreeForType[var.type] == 2
				end) and {'real3x3'} or nil
			):append(
				-- rank-3 always use real3x3x3 for transformation
				self.consStruct.vars:find(nil, function(var)
					return degreeForType[var.type] == 3
				end) and {'real3x3x3'} or nil
			),
			code = self:template([[<? 
for side=0,solver.dim-1 do
	if coord.vectorComponent == 'cartesian'
	or require 'hydro.coord.cartesian'.is(coord) 
	then
?>#define <?=cons_parallelPropagate?><?=side?>(resultName, U, pt, dx)\
	global <?=cons_t?> const * const resultName = U;
<?	else
?>#define <?=cons_parallelPropagate?><?=side?>(\
	resultName,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*real const */dx\
)\
/* TODO don't assign here, instead assign all fields and just don't propagate the scalars */\
	<?=cons_t?> resultName##base = *(U);\
<?		for _,var in ipairs(eqn.consStruct.vars) do
			local variance = assert(var.variance)
			local degree = degreeForType[var.type]
			if variance == '' then
			elseif variance == 'u' then
-- P_u^a T^u
?>	resultName##base.<?=var.name?> = coord_parallelPropagateU<?=side?>(resultName##base.<?=var.name?>, pt, dx);\
<?
			elseif variance == 'l' then
-- T_u (P^-T)_a^u
?>	resultName##base.<?=var.name?> = coord_parallelPropagateL<?=side?>(resultName##base.<?=var.name?>, pt, dx);\
<?
			elseif variance == 'll' then
-- (P^-T)_a^u (P^-T)_b^v T_uv
?>				{\
					real3x3 t = real3x3_from_<?=var.type?>(resultName##base.<?=var.name?>);\
					t.x = coord_parallelPropagateL<?=side?>(t.x, pt, dx);\
					t.y = coord_parallelPropagateL<?=side?>(t.y, pt, dx);\
					t.z = coord_parallelPropagateL<?=side?>(t.z, pt, dx);\
					t = real3x3_transpose(t);\
					t.x = coord_parallelPropagateL<?=side?>(t.x, pt, dx);\
					t.y = coord_parallelPropagateL<?=side?>(t.y, pt, dx);\
					t.z = coord_parallelPropagateL<?=side?>(t.z, pt, dx);\
					resultName##base.<?=var.name?> = <?=var.type?>_from_real3x3(t);\
				}\
<?			elseif variance == 'lll' then
?>				{\
					real3x3x3 t = real3x3x3_from_<?=var.type?>(resultName##base.<?=var.name?>);\
					real3 tmp;\
<?				local table = require 'ext.table'
				local is = table()
				for e=1,3 do
					for i,xi in ipairs(xNames) do
						is[e] = xi
						for j,xj in ipairs(xNames) do
							is[e%3+1] = xj
							for k,xk in ipairs(xNames) do
								is[(e+1)%3+1] = xk
?>						tmp.<?=xk?> = t.<?=is:concat'.'?>;\
<?							end
?>						tmp = coord_parallelPropagateL<?=side?>(tmp, pt, dx);\
<?							for k,xk in ipairs(xNames) do
								is[(e+1)%3+1] = xk
?>						t.<?=is:concat'.'?> = tmp.<?=xk?>;\
<?							end
						end
					end
				end
?>					resultName##base.<?=var.name?> = <?=var.type?>_from_real3x3x3(t);\
				}\
<?			else
				error("don't know how to handle variance for "..('%q'):format(variance))
			end
		end
?>	<?=cons_t?> const * const resultName = &resultName##base;
<?	end
end
?>]], 		{
				coord = solver.coord,
				degreeForType = degreeForType,
			}),
		}
	end
end

-- separate this out so composite solvers can init this and this alone -- without anything else (which would cause module name collisions)
function Equation:initCodeModule_cons_prim_eigen_waves()
	local solver = self.solver

	assert(self.consStruct)
	solver.modules:add{
		name = self.symbols.cons_t,
		structs = {self.consStruct},
		depends = self.getModuleDepends_cons_t and self:getModuleDepends_cons_t() or nil,
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.cons_t..' cons_t;',
	}

	if self.primStruct then
		solver.modules:add{
			name = self.symbols.prim_t,
			structs = {self.primStruct},
			depends = self.getModuleDepends_prim_t and self:getModuleDepends_prim_t() or nil,
			-- only generated for cl, not for ffi cdef
			headercode = 'typedef '..self.symbols.prim_t..' prim_t;',
		}
	else
		solver.modules:add{
			name = self.symbols.prim_t,
			depends = {self.symbols.cons_t},
			typecode = 'typedef '..self.symbols.cons_t..' '..self.symbols.prim_t..';',
			-- only generated for cl, not for ffi cdef
			headercode = 'typedef '..self.symbols.prim_t..' prim_t;',
		}
	end

	solver.modules:add{
		name = self.symbols.eigen_t,
		structs = {assert(self.eigenStruct)},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.eigen_t..' eigen_t;',
	}

	solver.modules:add{
		name = self.symbols.waves_t,
		structs = {assert(self.wavesStruct)},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.waves_t..' waves_t;',
	}
end

function Equation:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = self.symbols.fluxFromCons,
		depends = {
			self.solver.solver_t,
			self.symbols.cons_t,
			self.symbols.eigen_fluxTransform,
			self.symbols.eigen_forCell,
			'normal_t',
		},
		code = self:template[[
#define <?=fluxFromCons?>(\
	/*<?=cons_t?> const * const */flux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	<?=eigen_t?> eig;\
	<?=eigen_forCell?>(&eig, solver, U, x, n)\
	<?=eigen_fluxTransform?>(flux, solver, &eig, U, x, n);\
}
]],
	}
end

-- put your eigenWaveCode / consWaveCode dependencies here
function Equation:getModuleDepends_waveCode() end

function Equation:initCodeModule_waveCode()
	self.solver.modules:add{
		name = self.symbols.eqn_waveCode_depends,
		depends = self:getModuleDepends_waveCode(),
	}
end

-- put your getDisplayVars code dependencies here
-- called from solverbase
function Equation:getModuleDepends_displayCode() 
	return {
		'SETBOUNDS',
	}
end

Equation.displayVarCodeUsesPrims = false
function Equation:getDisplayVarCodePrefix()
	return self:template[[
	global <?=cons_t?> const * const U = buf + index;
<? if eqn.displayVarCodeUsesPrims then 
?>	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
<? end 
?>]]
end

function Equation:getDisplayVars()
	return self.solver:createDisplayVarArgsForStructVars(self.consStruct.vars)
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
			global <?=cons_t?> const * const Ujm = U - solver->stepsize.s<?=j?>;
			global <?=cons_t?> const * const Ujp = U + solver->stepsize.s<?=j?>;
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
		global <?=cons_t?> const * const Uim = U - solver->stepsize.s<?=i?>;
		global <?=cons_t?> const * const Uip = U + solver->stepsize.s<?=i?>;
		global <?=cons_t?> const * const Ujm = U - solver->stepsize.s<?=j?>;
		global <?=cons_t?> const * const Ujp = U + solver->stepsize.s<?=j?>;

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
	return self.solver:createDisplayVarArgsForStructVars(self.eigenVars, '(&eig)')
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

-- By default calcDT is taken from hydro/eqn/cl/calcDT.cl
-- Override to provide your own.
function Equation:initCodeModule_calcDTCell()
	self.solver.modules:addFromMarkup(self:template(file['hydro/eqn/cl/calcDT.cl']))
end

--[[
Default code for the following:
	primFromCons
	consFromPrim
	apply_dU_dW : prim_t -> cons_t 
	apply_dW_dU : cons_t -> prim_t

The default assumes prim_t == cons_t and this transformation is identity
--]]
function Equation:initCodeModule_consFromPrim_primFromCons()
	assert(not self.primStruct, "if you're using the default prim<->cons code then you shouldn't have any primStruct")

	self.solver.modules:add{
		name = self.symbols.primFromCons,
		depends = {self.solver.solver_t, self.symbols.prim_t, self.symbols.cons_t},
		code = self:template[[
#define <?=primFromCons?>(W, solver, U, x)	(*(W) = *(U))
/*
void <?=primFromCons?>(
	<?=prim_t?> * const W,
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const * const U, 
	real3 const x
) { 
	return U; 
}
*/
]],
	}
	
	self.solver.modules:add{
		name = self.symbols.consFromPrim,
		depends = {
			self.solver.solver_t,
			self.symbols.prim_t,
			self.symbols.cons_t,
		},
		code = self:template[[
#define <?=consFromPrim?>(U, solver, W, x)	(*(U) = *(W))
/*
void <?=consFromPrim?>(
	<?=cons_t?> * const U,
	constant solver_t const * const solver,
	<?=prim_t?> const * const W, 
	real3 const x
) { 
	return W; 
}
*/
]],
	}

	-- only used by PLM
	self.solver.modules:add{
		name = self.symbols.apply_dU_dW,
		depends = {self.solver.solver_t, self.symbols.prim_t, self.symbols.cons_t},
		code = self:template[[
/*
WA = W components that make up the jacobian matrix
W = input vector
x = coordinate location
returns output vector
*/
#define <?=apply_dU_dW?>(result, solver, WA, W, x)	(*(result) = *(W))
/*
void <?=apply_dU_dW?>(
	<?=cons_t?> * const result,
	constant solver_t* solver,
	<?=prim_t?> const * const WA, 
	<?=prim_t?> const * const W, 
	real3 const x
) { 
	return W; 
}
*/
]],
	}

	-- only used by PLM
	self.solver.modules:add{
		name = self.symbols.apply_dW_dU,
		depends = {self.solver.solver_t, self.symbols.prim_t, self.symbols.cons_t},
		code = self:template[[
/*
WA = W components that make up the jacobian matrix
U = input vector
x = coordinate location
returns output vector
*/
#define <?=apply_dW_dU?>(solver, WA, U, x)	(*(result) = (*U))
/*
void <?=apply_dW_dU?>(
	<?=prim_t?> const * W,
	constant solver_t const * const solver,
	<?=prim_t?> const * const WA, 
	<?=cons_t?> const * const U,
	real3 const x
) { 
	return U; 
}
*/
]],
	}
end

-- especially used by the num rel stuff
-- especially the finite-difference num rel (bsnsok-fd), but sometimes by the constrain equations of the finite-volume num rel (adm3d, z4, etc)
-- but anyone can use it.
function Equation:fieldTypeForVar(varname)
	local _, var = self.consStruct.vars:find(nil, function(v) return v.name == varname end)
	if not var then
		error("couldn't find var "..varname)
	end
	return var.type
end

function Equation:makePartial1(field, fieldType, nameOverride)
	-- order = 4 = 2 * 2 = 2 * (3 - 1), so numGhost == 3
	local derivOrder = 2 * (self.solver.numGhost - 1)
	fieldType = fieldType or self:fieldTypeForVar(field)
	return makePartials.makePartial1(derivOrder, self.solver, field, fieldType, nameOverride)
end

function Equation:makePartial2(field, fieldType, nameOverride)
	local derivOrder = 2 * (self.solver.numGhost - 1)
	fieldType = fieldType or self:fieldTypeForVar(field)
	return makePartials.makePartial2(derivOrder, self.solver, field, fieldType, nameOverride)
end

return Equation
