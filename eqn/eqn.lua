local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local makestruct = require 'eqn..makestruct'

local common = require 'common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local Equation = class()

-- this is passed on to solver/calcDerivFV.cl
-- it has the effect of adding the connection terms Conn^k_jk u^I_,j (for the I'th conserved quantity u^I)
Equation.weightFluxByGridVolume = true

-- Whether the eqn has its own eigen_*** code.
-- Otherwise eqn/cl/eigen.cl is used, which depends on the default eigen_t structures.
Equation.hasEigenCode = nil

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
	for _,var in ipairs(self.consVars) do
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
		else
			error"you are here"
		end
	end
	return vars
end

function Equation:init(args)
	self.solver = assert(args.solver)

	local app = args.solver.app
	self.prim_t = app:uniqueName'prim_t'
	self.cons_t = app:uniqueName'cons_t'
	self.consLR_t = app:uniqueName'consLR_t'
	self.eigen_t = app:uniqueName'eigen_t'
	self.waves_t = app:uniqueName'waves_t'
	
	local numReals
	if self.consVars then
		numReals = makestruct.countScalars(self.consVars)
		if self.primVars then
			local numPrimReals = makestruct.countScalars(self.primVars)
			assert(numPrimReals <= numReals, "hmm, this is awkward")
		end
	end
	
	if not self.numStates then 
		self.numStates = numReals
		if not self.numStates then
			error("you either need to define numStates or consVars")
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
	
	self.initStateNames = table.map(self.initStates, function(info) return info.name end)


	-- r min, for spherical coordinates
	-- what variables to mirror at sphere center
	-- 2013 Baumgarte et al, "Numerical Relativity in Spherical Polar Coordinates...", IIIB
	-- 2017 Ruchlin et al, section E.1
	self.boundarySphereRMinMirrorVars = {
		self:getParityVars(-1, 1, 1),
		{},
		{},
	} 

	-- theta min/max, for spherical coordinates
	self.boundarySphereThetaMirrorVars = {
		{},
		self:getParityVars(1, -1, 1),
		{},
	}

	-- phi min/max, for spherical coordinates
	self.boundarySpherePhiMirrorVars = {
		self:getParityVars(-1, -1, 1),
		{},
		{},
	}

	-- phi min/max, for cylindrical or for spherical coordinates
	self.boundaryCylinderCenterMirrorVars = {
		self:getParityVars(-1, -1, 1),
		{},
		{},
	}

	-- x,y,z min/max for cartesian coordinates
	self.boundaryCartesianMirrorVars = {
		self:getParityVars(-1, 1, 1),
		self:getParityVars(1, -1, 1),
		self:getParityVars(1, 1, -1),
	}
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
	local cl = require('guivar.'..vartype)

	-- no non-strings allowed as names, especially not numbers,
	-- because I'm going to map names into guiVars' keys, and I don't want to overwrite any integer-indexed guiVars
	assert(args.name and type(args.name) == 'string')

	local var = cl(args)
	self.guiVars:insert(var)
	self.guiVars[var.name] = var
end

-- always call super first
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
	assert(self.initStates, "expected Eqn.initStates")
	self.initState = self.initStates[self.solver.initStateIndex](self.solver, self.solver.initStateArgs)
	assert(self.initState, "couldn't find initState "..self.solver.initStateIndex)

	-- then setup the gui vars
	if self.initState.guiVars then
		self:addGuiVars(self.initState.guiVars)
	end
	for _,op in ipairs(self.solver.ops) do
		if op.guiVars then
			self:addGuiVars(op.guiVars)
		end
	end
end

function Equation:getCodePrefix()
	return (self.guiVars and table.mapi(self.guiVars, function(var,i,t) 
		return (var.compileTime and var:getCode() or nil), #t+1
	end) or table()):append{
		self.initState.getCodePrefix 
			and self.initState:getCodePrefix(self.solver)
			or '',
		
		-- functions that prim-cons code will use, but which use macros:
		self.getCommonFuncCode and self:getCommonFuncCode() or '',
		
		-- prim-cons goes here
		-- it goes last so it has access to everything above it
		-- but it must be in codeprefix so initstate has access to it
		self:getPrimConsCode() or '',
	}:concat'\n'
end

function Equation:getExtraTypeCode()
	return template([[
typedef union { 
	real ptr[<?=eqn.numWaves?>]; 
} <?=eqn.waves_t?>;
]],		{
			eqn = self,
		})
end

function Equation:getTypeCode()
	assert(self.consVars)
	local lines = table{
		makestruct.makeStruct(self.cons_t, self.consVars),
	}
	if self.primVars then
		lines:insert(makestruct.makeStruct(self.prim_t, self.primVars))
	else
		lines:insert('typedef '..self.cons_t..' '..self.prim_t..';')
	end
	return lines:concat'\n'
end

function Equation:getSolverCode()
	return template(file[self.solverCodeFile], {eqn=self, solver=self.solver})
end

function Equation:getInitStateCode()
	assert(self.initStateCode, "expected solver.eqn.initStateCode")
	return template(self.initStateCode, {
		eqn = self,
		code = self.initState:initState(self.solver),
		solver = self.solver,
	})
end

Equation.displayVarCodeUsesPrims = false
function Equation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
<? if eqn.displayVarCodeUsesPrims then 
?>	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
<? end 
?>]], {
		eqn = self,
	})
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
	return self:getDisplayVarsForStructVars(self.consVars)
end

-- does anyone even use this anymore?  nobody should...
function Equation:getEigenTypeCode()
	if self.eigenVars then
		return makestruct.makeStruct(self.eigen_t, self.eigenVars)
	
	-- use the default matrix structures	
	-- whose code is in eqn/cl/eigen.cl (included below)
	else
error("this is deprecated")
		return template([[
typedef struct {
	real evL[<?=numIntStates * numWaves?>];
	real evR[<?=numIntStates * numWaves?>];
	real A[<?=numIntStates * numIntStates?>];
} <?=eqn.eigen_t?>;
]], 	{
			numIntStates = self.numIntStates,
			numWaves = self.numWaves,
			solver = self.solver,
			eqn = self,
		})
	end
end

function Equation:getEigenCode()
	if self.hasEigenCode then return end
	return template(file['eqn/cl/eigen.cl'], {
		solver = self.solver,
		eqn = self,
	})
end

function Equation:getEigenDisplayVars()
	-- use the automatic codegen for display vars
	if self.eigenVars then
		return self:getDisplayVarsForStructVars(self.eigenVars, '(&eig)')

	-- use the autogen left & right eigenvector matrices
	else
		return range(self.numIntStates * self.numWaves):map(function(i)
			local row = (i-1)%self.numWaves
			local col = (i-1-row)/self.numWaves
			return {name='evL_'..row..'_'..col, code='value.vreal = eig.evL['..i..'];'}
		end):append(range(self.numIntStates * self.numWaves):map(function(i)
			local row = (i-1)%self.numIntStates
			local col = (i-1-row)/self.numIntStates
			return {name='evR_'..row..'_'..col, code='value.vreal = eig.evR['..i..'];'}
		end)):append(range(self.numIntStates * self.numIntStates):map(function(i)
			local row = (i-1)%self.numIntStates
			local col = (i-1-row)/self.numIntStates
			return {name='A_'..row..'_'..col, code='value.vreal = eig.A['..i..'];'}
		end))
	end
end

function Equation:eigenWaveCodePrefix(side, eig, x)
	return ''
end
function Equation:eigenWaveCode(side, eig, x, waveIndex)
	return '\n#error :eigenWaveCode() not implemented'
end

function Equation:consWaveCodePrefix(side, U, x)
	return '\n#error :consWaveCodePrefix() not implemented'
end
function Equation:consWaveCode(side, U, x)
	return '\n#error :consWaveCode() not implemented'
end

-- default implementation -- the first is the min and the last is the max
-- however some don't do this, like GLM
function Equation:eigenMinWaveCode(side, eig, x)
	return self:eigenWaveCode(side, eig, x, 0)
end
function Equation:eigenMaxWaveCode(side, eig, x)
	return self:eigenWaveCode(side, eig, x, self.numWaves-1)
end
function Equation:consMinWaveCode(side, U, x)
	return self:consWaveCode(side, U, x, 0)
end
function Equation:consMaxWaveCode(side, U, x)
	return self:consWaveCode(side, U, x, self.numWaves-1)
end

-- Whether the eqn has its own calcDT.  Otherwise eqn/cl/calcDT.cl is used. 
Equation.hasCalcDTCode = nil

function Equation:getCalcDTCode()
	if self.hasCalcDTCode then return end
	return template(file['eqn/cl/calcDT.cl'], {
		solver = self.solver, 
		eqn = self,
	})
end

Equation.hasFluxFromConsCode = nil

function Equation:getFluxFromConsCode()
	if self.hasFluxFromConsCode then return end
	return template([[
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	return eigen_fluxTransform_<?=side?>(solver, eigen_forCell_<?=side?>(solver, U, x), U, x);
}
<? end ?>
]], {
		solver = self.solver, 
		eqn = self,
	})
end

--[[
Default code for the following:
	primFromCons
	consFromPrim
	apply_dU_dW : prim_t -> cons_t 
	apply_dW_dU : cons_t -> prim_t

The default assumes prim_t == cons_t and this transformation is identity
--]]
function Equation:getPrimConsCode()
	assert(not self.primVars, "if you're using the default prim<->cons code then you shouldn't have any primVars")

	return template([[

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
]], {
		eqn = self,
		solver = self.solver,
	})
end

return Equation
