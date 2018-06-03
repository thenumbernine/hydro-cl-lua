-- TODO make this solver/solver.lua, and make the old solver.lua something like structured-grid-solver

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local vec3 = require 'vec.vec3'
local tooltip = require 'tooltip'
local time, getTime = table.unpack(require 'time')
require 'common'(_G)	-- xNames, symNames


local integrators = require 'int.all'
local integratorNames = integrators:map(function(integrator) return integrator.name end)


local SolverBase = class()

SolverBase.name = 'Solver'

-- override to specify which eqn/*.lua to use as the equation
SolverBase.eqnName = nil

--[[
args:
	app
	dim
	eqn = name of eqn/<name>.lua to use
	eqnArgs (optional) = args of eqn ctor
	initState = init state name
	initStateArgs (optional) = args of init state ctor
	integrator = name of integrator in int/all.lua to use
--]]
function SolverBase:init(args)
	assert(args)
	self.app = assert(args.app)
	self.dim = assert(args.dim)
	
	self.color = vec3(math.random(), math.random(), math.random()):normalize()


	-- operators for this solver
	self.ops = table()


	-- hmm, do some eqns create ops need to know the grid size? 
	self.eqnName = args.eqn
	self.eqnArgs = args.eqnArgs
	self:createEqn()


	self.name = self.eqn.name..' '..self.name


	self.initStateIndex = table.find(self.eqn.initStateNames, args.initState) or 1
	self.initStateArgs = args.initStateArgs


	self.integratorIndex = integratorNames:find(args.integrator) or 1
	
	self:createDisplayVars()	-- depends on eqn
end

function SolverBase:postInit()
	self:refreshEqnInitState()
	self:resetState()
end

-- call this when the solver initializes or changes the codePrefix (or changes initState)
-- it will build the code prefix and refresh everything related to it
-- TODO if you change cons_t then call resetState etc (below the refreshEqnInitState() call a few lines above) in addition to this -- or else your values will get messed up
function SolverBase:refreshEqnInitState()
	--[[
	circular dependency
	... I want to createInitState when the initState changes
	but if a gui var changes
	then I want to call everything other than initState
	--]]
	
	-- this influences createCodePrefix (via its call of eqn:getCodePrefix)
	--  and refreshInitStateProgram()
	self.eqn:createInitState()
	
	-- Right now within eqn:createInitState I'm adding any subclass-specific gui vars
	-- so only after it finishes and all gui vars are created, ask the eqn.initState object if it wants to modify anything.
	-- Don't do this during Solver:refreshInitStateProgram()->InitCond:initState() or the changes won't get into the header.
	-- Hmm... should the initState even have control over the eqn's vars?
	if self.eqn.initState.overrideGuiVars then
		for k,v in pairs(self.eqn.initState.overrideGuiVars) do
			if self.eqn.guiVars[k] then
				self.eqn.guiVars[k].value = v
			end
		end
	end

	self:refreshCodePrefix()
end

-- this is the general function - which just assigns the eqn provided by the arg
-- but it can be overridden for specific equations
function SolverBase:createEqn()
	self.eqn = require('eqn.'..assert(self.eqnName, "expected solver.eqnName"))(table(
		self.eqnArgs or {},
		{solver = self}
	))
end



function SolverBase:refreshCodePrefix()
	self:createCodePrefix()		-- depends on eqn, gridSize, displayVars
	self:refreshIntegrator()	-- depends on eqn & gridSize ... & ffi.cdef cons_t
end

function SolverBase:createCodePrefix()
	local lines = table()
	
	if self.dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end

	-- real3
	lines:insert(file['math.h'])
	lines:insert(template(file['math.cl'], {
		xNames = xNames,
		symNames = symNames,
	}))

	lines:append{
		'#define geometry_'..self.geometry.name..' 1',
		'#define dim '..self.dim,
		'#define numStates '..self.eqn.numStates,
		'#define numIntStates '..self.eqn.numIntStates,
		'#define numWaves '..self.eqn.numWaves,
	}
	
	lines:insert(self.geometry:getCode(self))

	-- this can use the coord_raise or coord_lower code
	-- which is associated with the coordinate system,
	-- TODO rename 'geom' to 'coord'
	lines:append{
		self.eqn:getTypeCode(),
		self.eqn:getExtraTypeCode(),
		self.eqn:getEigenTypeCode() or '',
		self.eqn:getCodePrefix() or '',
	}

	return lines:concat'\n'
end


function SolverBase:refreshIntegrator()
	self.integrator = integrators[self.integratorIndex](self)
end

function SolverBase:resetState()
	self.t = 0
	self.app.cmds:finish()

	self.eqn:resetState()

	for _,op in ipairs(self.ops) do
		if op.resetState then
			op:resetState()
		end
	end
end


local DisplayVarGroup = class()

function DisplayVarGroup:init(args)
	self.name = assert(args.name)
	self.vars = table(args.vars)
end


local DisplayVar = class()

-- this is the default DisplayVar
SolverBase.DisplayVar = DisplayVar



function SolverBase:newDisplayVarGroup(args)
	local displayVarGroup = DisplayVarGroup(args)
	self.displayVarGroups:insert(displayVarGroup)
	return displayVarGroup
end


-- still used by gr-hd-separate to add 'extraArgs'
function SolverBase:getUBufDisplayVarsArgs()
	return {
		type = self.eqn.cons_t,
		codePrefix = self.eqn:getDisplayVarCodePrefix(),
		bufferField = 'UBuf',
	}
end

function SolverBase:addUBufDisplayVars()
	local group = self:newDisplayVarGroup{name='U'}
	
	local args = self:getUBufDisplayVarsArgs()
	args.group = group
	
	-- TODO rename to 'getDisplayVarDescs()'
	-- gets var descriptions, which is {name=code, [type=type]}
	args.vars = self.eqn:getDisplayVars()
	
	self:addDisplayVarGroup(args, self.DisplayVar_U)
end

function SolverBase:getDisplayInfosForType()
	return {
		real3 = {
			{name = ' x', code = '	*valuevec = _real3(valuevec->x,0,0);'},
			{name = ' y', code = '	*valuevec = _real3(valuevec->y,0,0);'},
			{name = ' z', code = '	*valuevec = _real3(valuevec->z,0,0);'},
			{name = ' mag', code = '	*valuevec = _real3(real3_len(*valuevec),0,0);', magn=true},
		},
		
		-- hmm, valuevec has to be bigger for this to work
		-- but does that mean I have to store 6 components in valuevec?
		-- I suppose it does if I want a sym3-specific visualization
		sym3 = {
			{name = ' x', code = '	*valuevec = sym3_x(*valuesym3); *valuevec_hi = _real3(0,0,0);', vartype='real3'},
			{name = ' y', code = '	*valuevec = sym3_y(*valuesym3); *valuevec_hi = _real3(0,0,0);', vartype='real3'},
			{name = ' z', code = '	*valuevec = sym3_z(*valuesym3); *valuevec_hi = _real3(0,0,0);', vartype='real3'},
	
			--[[ these are already added through real3 x_i real x_j
			{name = ' xx', code = '	*valuesym3 = _sym3(valuesym3->xx, 0,0,0,0,0);'},
			{name = ' xy', code = '	*valuesym3 = _sym3(valuesym3->xy, 0,0,0,0,0);'},
			{name = ' xz', code = '	*valuesym3 = _sym3(valuesym3->xz, 0,0,0,0,0);'},
			{name = ' yy', code = '	*valuesym3 = _sym3(valuesym3->yy, 0,0,0,0,0);'},
			{name = ' yz', code = '	*valuesym3 = _sym3(valuesym3->yz, 0,0,0,0,0);'},
			{name = ' zz', code = '	*valuesym3 = _sym3(valuesym3->zz, 0,0,0,0,0);'},
			--]]

			{name = ' norm', code = '	*valuesym3 = _sym3( sqrt(sym3_dot(*valuesym3, *valuesym3)), 0,0,0,0,0);'},
			{name = ' tr', code = '	*valuesym3 = _sym3( sym3_trace(*valuesym3), 0,0,0,0,0);'},
		}
	}
end



function SolverBase:addDisplayVarGroup(args, cl)
	cl = cl or self.DisplayVar

	if not args.group then
		args.group = self:newDisplayVarGroup{name = args.name}
	end

	local group = args.group
	local varInfos = args.vars

	local enableScalar = true
	local enableVector = true

	for i,varInfo in ipairs(varInfos) do
	
		local name, code, vartype
		for k,v in pairs(varInfo) do
			if k == 'type' then
				vartype = v
			else
				assert(not name and not code)
				name = k
				code = v
			end
		end

		-- this is a quick fix for the magn vars associated with the axis real3 vars of the sym3s
		-- but it has lots of redundant symmetries for the sym3 real3 elements 
		local function addvar(args)
			-- enable the first scalar field
			-- also enable the first vector field on non-1D simulations
			local enabled
			if group.name == 'U' 
			or (group.name:sub(1,5) == 'error' and self.dim == 1) 
			then
				if args.vartype ~= 'real3' then
					enabled = enableScalar
					
					if self.dim ~= 1 then
						enableScalar = nil
					end
				else
					if self.dim ~= 1 then
						enabled = enableVector
						enableVector = nil
					end
				end
			end

			local var = cl(table(args, {
				vectorField = args.vartype == 'real3',
				enabled = enabled,
			}))
			group.vars:insert(var)

			local infosForType = self:getDisplayInfosForType()

			local infos = infosForType[args.vartype]
			if infos then
				for _,info in ipairs(infos) do
					local scalarVar = addvar(table(args, {
						name = args.name .. info.name,
						code = args.code .. info.code,
						vartype = info.vartype or 'real',
						magn = info.magn or false,
						vectorField = info.vartype == 'real3',
						enabled = group.name == 'U' and self.dim == 1 and info.vartype ~= 'real3',
					}))
				
					-- tie together vectors and magnitudes,
					-- since reduceMin and Max applied to vectors is gonna reduce their magnitude
					-- so I need to always compile the magnitude kernels, even if they are not enabled
					if info.magn then
						var.magVar = scalarVar
						scalarVar.vecVar = var
					end
				end
			end
		
			return var
		end
	
		addvar(table(args, {
			solver = self,
			name = group.name .. ' ' .. name,
			code = code,
			vartype = vartype or 'real',
		}))
	end

	return args.group
end



function SolverBase:addDisplayVars()
	self:addUBufDisplayVars()
	
	-- might contain nonsense :-p
	self:addDisplayVarGroup{
		name = 'reduce', 
		vars = {{['0'] = '*value = buf[index];'}},
	}
end

function SolverBase:createDisplayVars()
	self.displayVarGroups = table()
	self:addDisplayVars()
	self.displayVars = table()
	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		self.displayVars:append(displayVarGroup.vars)
	end

	-- make lookup by name
	self.displayVarForName = self.displayVars:map(function(var)
		return var, var.name
	end)
end

function SolverBase:initDraw()
end

SolverBase.fpsNumSamples = 30

function SolverBase:update()
	--[[
	Here's an update-based FPS counter.
	This isn't as precise as profiling events for all of my OpenCL calls
	but it will give a general idea while running simulations continuously.
	Because pauses will mess with the numbers, I'll only look at the last n many frames.  
	Maybe just the last 1.
	--]]
	local thisTime = getTime()
	if not self.fpsSamples then
		self.fpsIndex = 0
		self.fpsSamples = table()
	end
	if self.lastFrameTime then
		local deltaTime = thisTime - self.lastFrameTime
		local fps = 1 / deltaTime
		self.fpsIndex = (self.fpsIndex % self.fpsNumSamples) + 1
		self.fpsSamples[self.fpsIndex] = fps
		self.fps = self.fpsSamples:sum() / #self.fpsSamples
	end
	self.lastFrameTime = thisTime
end

function SolverBase:updateGUIParams()
	if tooltip.comboTable('integrator', self, 'integratorIndex', integratorNames) then
		self:refreshIntegrator()
	end

	-- I think I'll display my GMRES # steps to converge / epsilon error ... 
	if self.integrator.updateGUI then
		ig.igSameLine()
		ig.igPushIDStr'integrator'
		if ig.igCollapsingHeader':' then
			self.integrator:updateGUI()
		end
		ig.igPopID()
	end
	
	for i,op in ipairs(self.ops) do
		ig.igPushIDInt(i)
		op:updateGUI()
		ig.igPopID()
	end
end

return SolverBase
