local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local math = require 'ext.math'
local vec3 = require 'vec.vec3'
local CLImageGL = require 'cl.imagegl'
local CLBuffer = require 'cl.obj.buffer'
local GLTex2D = require 'gl.tex2d'
local GLTex3D = require 'gl.tex3d'
local glreport = require 'gl.report'
local clnumber = require 'cl.obj.number'
local template = require 'template'
local vec3sz = require 'ffi.vec.vec3sz'
local tooltip = require 'tooltip'
local roundup = require 'roundup'
local time, getTime = table.unpack(require 'time')

--local tryingAMR = 'dt vs 2dt'
--local tryingAMR = 'gradient'

local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}


local Solver = class()

Solver.name = 'Solver'
Solver.numGhost = 2

Solver.integrators = require 'int.all'
Solver.integratorNames = Solver.integrators:map(function(integrator) return integrator.name end)

--[[
args:
	app
	gridSize
	dim (optional, required if gridSize is a vec3sz)
--]]
function Solver:init(args)
	assert(args)
	self.app = assert(args.app)

	self.color = vec3(math.random(), math.random(), math.random()):normalize()

	self.mins = vec3(table.unpack(args.mins or {-1, -1, -1}))
	self.maxs = vec3(table.unpack(args.maxs or {1, 1, 1}))


	self.dim = assert(args.dim)
	
	-- TODO OK this is a little ambiguous ...
	-- gridSize is the desired grid size
	-- self.gridSize is gonna be that, plus numGhost on either side ...
	-- so maybe I should rename 'self.gridSize' to 'self.gridSizeWithBorder'
	-- and 'self.gridSizeWithoutBorder' to 'self.gridSize'

	local gridSize = assert(args.gridSize)
	if type(gridSize) == 'number' then 
		self.gridSize = vec3sz(gridSize,1,1)
	elseif type(gridSize) == 'cdata' 
	and ffi.typeof(gridSize) == vec3sz 
	then
		self.gridSize = vec3sz(gridSize)
	elseif type(gridSize) == 'table' then
		self.gridSize = vec3sz(table.unpack(gridSize))
	else
		error("can't understand args.gridSize type "..type(args.gridSize).." value "..tostring(args.gridSize))
	end

	for i=0,self.dim-1 do self.gridSize:ptr()[i] = self.gridSize:ptr()[i] + 2 * self.numGhost end
	for i=self.dim,2 do self.gridSize:ptr()[i] = 1 end

		-- operators for this solver
	self.ops = table()
	self:createEqn(args.eqn)

	self:createBoundaryOptions()
	self:finalizeBoundaryOptions()

	self.boundaryMethods = {}
	for i=1,3 do
		for _,minmax in ipairs(minmaxs) do
			local var = xs[i]..minmax
			self.boundaryMethods[var] = self.boundaryOptions:find(
				(args.boundary or {})[var] or 'freeflow',
				function(option, search)
					return search == next(option)
				end)
		end
	end

	self.name = self.eqn.name..' '..self.name

	self.checkNaNs = false
	self.useFixedDT = not not args.fixedDT
	self.fixedDT = args.fixedDT or self.fixedDT or .001
	self.cfl = args.cfl or .5	--/self.dim
	self.initStateIndex = table.find(self.eqn.initStateNames, args.initState) or 1
	self.integratorIndex = self.integratorNames:find(args.integrator) or 1
	self.fluxLimiter = ffi.new('int[1]', (self.app.limiterNames:find(args.fluxLimiter) or 1)-1)


	self.geometry = require('geom.'..args.geometry){solver=self}

	self.usePLM = args.usePLM
	assert(not self.usePLM or self.fluxLimiter[0] == 0, "are you sure you want to use flux and slope limiters at the same time?")
	self.slopeLimiter = ffi.new('int[1]', (self.app.limiterNames:find(args.slopeLimiter) or 1)-1)

	local solver = self

	-- my kernel objs are going to need workgroup info based on domain.size-2*noGhost as well as domain.size ... 
	-- ...and rather than require an extra argument, I think I'll just take advantage of a closure
	local Program = class(require 'cl.obj.program')
	function Program:init(args)
		args.env = solver.app.env
		args.domain = solver.domain
		Program.super.init(self, args)
	end
	self.Program = Program  

	self:refreshGridSize()
end

--[[
boundaryOptions is a table of {name = args => assign code}
	args of the boundary function are:
		index
		assign = function(dst,src) 
			assignment operator
			default is dst = src
			poisson uses dist.field = src.field
		rhs
		gridSizeSide = solver.gridSize[side]
		side = 1,2,3
		minmax = 'min' or 'max'
		mirrorVars = which vars to reflect.
			for solver.UBuf this is taken from eqn
			for poisson this is nothing
--]]
function Solver:createBoundaryOptions()
	self.boundaryOptions = table{
		{periodic = function(args)
			local gridSizeSide = 'gridSize_'..xs[args.side]
			if args.minmax == 'min' then
				return '\t\t\t'..args.assign(
					'buf['..args.index'j'..']', 
					'buf['..args.index(gridSizeSide..'-2*numGhost+j')..']'
				)..';'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return '\t\t\t'..args.assign(
					'buf['..args.index(rhs)..']', 
					'buf['..args.index'numGhost+j'..']'
				)..';'
			end
		end},
		{mirror = function(args)
			local gridSizeSide = 'gridSize_'..xs[args.side]
			if args.minmax == 'min' then
				return table{
					'\t\t\t'..args.assign(
						'buf['..args.index'j'..']',
						' buf['..args.index'2*numGhost-1-j'..']'
					)..';'
				}:append(table.map((args.mirrorVars or {})[args.side] or {}, function(var)
					return '\t\t\t'..'buf['..args.index'j'..'].'..var..' = -buf['..args.index'j'..'].'..var..';'
				end)):concat'\n'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return table{
					'\t\t\t'..args.assign(
						'buf['..args.index(rhs)..']',
						'buf['..args.index(gridSizeSide..'-numGhost-1-j')..']'
					)..';'
				}:append(table.map((args.mirrorVars or {})[args.side] or {}, function(var)
					return '\t\t\t'..'buf['..args.index(rhs)..'].'..var..' = -buf['..args.index(rhs)..'].'..var..';'
				end)):concat'\n'
			end
		end},
		{freeflow = function(args)
			local gridSizeSide = 'gridSize_'..xs[args.side]
			if args.minmax == 'min' then
				return '\t\t\t'..args.assign(
					'buf['..args.index'j'..']',
					'buf['..args.index'numGhost'..']'
				)..';'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return '\t\t\t'..args.assign(
					'buf['..args.index(rhs)..']',
					'buf['..args.index(gridSizeSide..'-numGhost-1'
				)..']')..';'
			end
		end},
	}

	if self.eqn.createBoundaryOptions then
		self.eqn:createBoundaryOptions()
	end
end
function Solver:finalizeBoundaryOptions()
	self.boundaryOptionNames = self.boundaryOptions:map(function(option) return (next(option)) end)
	-- hmm here's the one time that table.map using k,v comes in handy
	self.boundaryOptionForName = self.boundaryOptions:map(function(option) local k,v = next(option) return v,k end)
end

--[[
this is shorthand for assignment
you can also assign each solver.boundaryMethod.[x|y|z][min|max] to an index in solver.boundaryOptions
args:
	a string = sets all boundaryMethods to the boundaryOptions index with name of the string in args
	a table = sets each xmin...zmax boundaryMethod with the associated name
--]]
function Solver:setBoundaryMethods(args)
	for _,x in ipairs(xs) do
		for _,minmax in ipairs(minmaxs) do
			local k = x..minmax
			local name
			if type(args) == 'string' then
				name = args	
			elseif type(args) == 'table' then
				name = args[k]
			else
				error("don't know how to deal with these args") 
			end
			local i = self.boundaryOptions:find(nil, function(option) return next(option) == name end)
			if not i then
				io.stderr:write("unable to find boundary method "..name.." -- can't assign it to side "..k.."\n")
				io.stderr:flush()
			else
				self.boundaryMethods[k] = i
			end
		end
	end
end

-- this is the general function - which just assigns the eqn provided by the arg
-- but it can be overridden for specific equations
function Solver:createEqn(eqn)
	self.eqn = require('eqn.'..assert(eqn))(self)
end


-- this is only used by Solver, but each kernel gets its own,
-- so TODO get rid of this and just use kernel's sizes
function Solver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	local localSize1d = math.min(maxWorkGroupSize, tonumber(self.gridSize:volume()))

	local localSize2d
	if self.dim == 3 then
		local localSizeX = math.min(tonumber(self.gridSize.x), 2^math.ceil(math.log(maxWorkGroupSize,2)/2))
		local localSizeY = maxWorkGroupSize / localSizeX
		localSize2d = {localSizeX, localSizeY}
	end

--	self.localSize = self.dim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
	-- TODO better than constraining by math.min(gridSize),
	-- look at which gridSizes have the most room, and double them accordingly, until all of maxWorkGroupSize is taken up
	local localSize = vec3sz(1,1,1)
	local rest = maxWorkGroupSize
	local localSizeX = math.min(tonumber(self.gridSize.x), 2^math.ceil(math.log(rest,2)/self.dim))
	localSize.x = localSizeX
	if self.dim > 1 then
		rest = rest / localSizeX
		if self.dim == 2 then
			localSize.y = math.min(tonumber(self.gridSize.y), rest)
		elseif self.dim == 3 then
			local localSizeY = math.min(tonumber(self.gridSize.y), 2^math.ceil(math.log(math.sqrt(rest),2)))
			localSize.y = localSizeY
			localSize.z = math.min(tonumber(self.gridSize.z), rest / localSizeY)
		end
	end

	-- this is grid size, but rounded up to the next localSize
	local globalSize = vec3sz(
		roundup(self.gridSize.x, localSize.x),
		roundup(self.gridSize.y, localSize.y),
		roundup(self.gridSize.z, localSize.z))

	-- don't include the ghost cells as a part of the grid coordinate space
	local sizeWithoutBorder = vec3sz(self.gridSize:unpack())
	for i=0,self.dim-1 do
		sizeWithoutBorder:ptr()[i] = sizeWithoutBorder:ptr()[i] - 2 * self.numGhost
	end
	local volumeWithoutBorder = tonumber(sizeWithoutBorder:volume())
	
	local volume = tonumber(self.gridSize:volume())

	return {
		localSize1d = localSize1d,
		localSize2d = localSize2d,
		localSize = localSize,
		globalSize = globalSize,
		sizeWithoutBorder = sizeWithoutBorder,
		volumeWithoutBorder = volumeWithoutBorder,
		volume = volume,
	}
end

function Solver:refreshGridSize()
	self.stepSize = vec3sz()
	self.stepSize.x = 1
	for i=1,self.dim-1 do
		self.stepSize:ptr()[i] = self.stepSize:ptr()[i-1] * self.gridSize:ptr()[i-1]
	end
print('self.stepSize', self.stepSize)

	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	self.maxWorkGroupSize = tonumber(self.app.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')
print('maxWorkGroupSize', self.maxWorkGroupSize)

	self.offset = vec3sz(0,0,0)
	
	self.domain = self.app.env:domain{
		size = {self.gridSize:unpack()},
		dim = self.dim,
	}

	local sizeProps = self:getSizePropsForWorkGroupSize(self.maxWorkGroupSize)
	for k,v in pairs(sizeProps) do
print(k,v)	
		self[k] = v
	end

	self.domainWithoutBorder = self.app.env:domain{
		size = {self.sizeWithoutBorder:unpack()},
		dim = self.dim,
	}


	self:createDisplayVars()	-- depends on eqn
	
	-- depends on eqn & gridSize
	self.buffers = table()
	self:createBuffers()
	self:finalizeCLAllocs()

	-- create the code prefix, reflect changes
	self:refreshEqnInitState()
	-- initialize things dependent on cons_t alone
	self:refreshCommonProgram()
	self:resetState()
end
	
-- call this when the solver initializes or changes the codePrefix (or changes initState)
-- it will build the code prefix and refresh everything related to it
-- TODO if you change cons_t then call resetState etc (below the refreshEqnInitState() call a few lines above) in addition to this -- or else your values will get messed up
function Solver:refreshEqnInitState()	
	--[[
	circular dependency
	... I want to createInitState when the initState changes
	but if a gui var changes
	then I want to call everything other than initState
	--]]
	
	-- this influences createCodePrefix (via its call of eqn:getCodePrefix)
	--  and refreshInitStateProgram()
	self.eqn:createInitState()

	self:refreshCodePrefix()
end

-- call this when a gui var changes
-- it rebuilds the code prefix, but doesn't reset the initState
function Solver:refreshCodePrefix()
	self:createCodePrefix()		-- depends on eqn, gridSize, displayVars
	self:refreshIntegrator()	-- depends on eqn & gridSize ... & ffi.cdef cons_t
	self:refreshInitStateProgram()
	self:refreshSolverProgram()
	self:refreshDisplayProgram()
	-- changing initState calls this, and could change boundary programs, so I'm putting this here
	-- bad excuse, I know
	self:refreshBoundaryProgram()
end

function Solver:refreshIntegrator()
	self.integrator = self.integrators[self.integratorIndex](self)
end

--[[
it'd be convenient to do this for all properties of all types ...
	real -> scalar
	real3 -> x, y, z, magn scalar displays, real3 vector display
	sym3 -> x, y, z column and row vector displays, xx..zz, det, trace, |x|..|z| col & row scalar displays
--]]
local DisplayVar = class()

-- this is the default DisplayVar
Solver.DisplayVar = DisplayVar

-- and this is the DisplayVar used for UBuf
Solver.DisplayVar_U = DisplayVar

DisplayVar.type = 'real'	-- default

-- TODO buf (dest) shouldn't have ghost cells
-- and dstIndex should be based on the size without ghost cells
DisplayVar.displayCode = [[
kernel void <?=name?>(
	<?=input?>,
	const global <?= var.type ?>* buf<?= 
	var.extraArgs and #var.extraArgs > 0 
		and ',\n\t'..table.concat(var.extraArgs, ',\n\t')
		or '' 
?>
) {
	SETBOUNDS(0,0);

	int4 dsti = i;
	int dstindex = index;
	
	real3 x = cell_x(i);
	real3 xInt[<?=solver.dim?>];
<? for i=0,solver.dim-1 do
?>	xInt[<?=i?>] = x;
	xInt[<?=i?>].s<?=i?> -= .5 * grid_dx<?=i?>;
<? end
?>
	//now constrain
	if (i.x < numGhost) i.x = numGhost;
	if (i.x >= gridSize_x - numGhost) i.x = gridSize_x - numGhost-1;
<? 
if solver.dim >= 2 then
?>	if (i.y < numGhost) i.y = numGhost;
	if (i.y >= gridSize_y - numGhost) i.y = gridSize_y - numGhost-1;
<? 
end
if solver.dim >= 3 then
?>	if (i.z < numGhost) i.z = numGhost;
	if (i.z >= gridSize_z - numGhost) i.z = gridSize_z - numGhost-1;
<? end 
?>
	//and recalculate read index
	index = INDEXV(i);

	
	real value[6] = {0,0,0,0,0,0};	//size of largest struct
	sym3* valuesym3 = (sym3*)value;
	real3* valuevec = (real3*)value;
	real3* valuevec_hi = (real3*)(value+3);

<?= var.codePrefix or '' ?>
<?= var.code ?>

<?= output ?>
}
]]

function DisplayVar:init(args)
	self.code = assert(args.code)
	self.name = assert(args.name)
	self.solver = assert(args.solver)
	self.type = args.type	-- or self.type
	self.displayCode = args.displayCode 	-- or self.displayCode
	self.codePrefix = args.codePrefix
	
	-- display stuff
	self.enabled = not not args.enabled 
	self.useLog = args.useLog or false
	self.color = vec3(math.random(), math.random(), math.random()):normalize()
	self.heatMapFixedRange = false	-- args.name ~= 'error'
	self.heatMapValueMin = 0
	self.heatMapValueMax = 1

	-- is it a vector or a scalar?
	self.vectorField = args.vectorField

	-- maybe this should be in args too?
	-- or - instead of buffer - how about all the kernel's args?
	-- but the reason I have to store the field here is that the buffer isn't made yet 
	-- TODO? make display vars after buffers so I can store the buffer here?
	self.bufferField = args.bufferField
	self.extraArgs = args.extraArgs
end

function DisplayVar:setArgs(kernel)
	kernel:setArg(1, self.solver[self.bufferField])
end

function DisplayVar:setToTexArgs()
	self:setArgs(self.calcDisplayVarToTexKernelObj.obj)
end

function DisplayVar:setToBufferArgs(var)
	self:setArgs(self.calcDisplayVarToBufferKernelObj.obj)
end


local DisplayVarGroup = class()

function DisplayVarGroup:init(args)
	self.name = assert(args.name)
	self.vars = table(args.vars)
end

function Solver:newDisplayVarGroup(args)
	local displayVarGroup = DisplayVarGroup(args)
	self.displayVarGroups:insert(displayVarGroup)
	return displayVarGroup
end


function Solver:addDisplayVarGroup(args, cl)
	cl = cl or self.DisplayVar

	local displayVarGroup = args.group or self:newDisplayVarGroup{name=args.name}
	local groupName = displayVarGroup.name

	for i,var in ipairs(args.vars) do
		assert(type(var) == 'table', "failed on var "..self.name)
		local name, code = next(var)
		displayVarGroup.vars:insert(cl{
			solver = self,
			
			-- right now I use var.name to determine unique variables when overlapping the different graphs in the 1D simulations
			-- so if I don't prefix it correctly then eig_0 and flux_0 will draw overlapped 
			-- TODO have that check for var name and group name, and have this name=name only,
			name = groupName..'_'..name,
			code = code,
					
			type = args.type, -- or self.type
			displayCode = args.displayCode,	-- or self.displayCode
			codePrefix = args.codePrefix,
			
			enabled = groupName == 'U' and (
					(self.dim ~= 1 and i == 1 
						--and not args.vectorField
						and not args.vectorFieldMagn
					)
					or (self.dim == 1 and not args.vectorField)
				)
				or (groupName == 'error' and self.dim==1),
			useLog = args.useLog,
			vectorField = args.vectorField,
			bufferField = groupName..'Buf',
		})
	end

	return displayVarGroup
end

function Solver:createDisplayVars()
	self.displayVarGroups = table()
	self:addDisplayVars()
	self.displayVars = table()
	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		self.displayVars:append(displayVarGroup.vars)
	end
end

-- still used by gr-hd-separate to add 'extraArgs'
function Solver:getUBufDisplayVarsArgs()
	return {
		type = self.eqn.cons_t,
		codePrefix = self.eqn:getDisplayVarCodePrefix(),
		bufferField = 'UBuf',
	}
end

function Solver:getDisplayInfosForType()
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

function Solver:addUBufDisplayVars()
	local group = self:newDisplayVarGroup{name='U'}
	
	local args = self:getUBufDisplayVarsArgs()
	args.group = group

	local varInfos = self.eqn:getDisplayVars()

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

			local var = self.DisplayVar(table(args, {
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
						enabled = self.dim == 1 and info.vartype ~= 'real3',
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
			name = group.name .. '_' .. name,
			code = code,
			vartype = vartype or 'real',
		}))
	end
end

function Solver:addDisplayVars()
	self:addUBufDisplayVars()
	
	-- might contain nonsense :-p
	self:addDisplayVarGroup{
		name = 'reduce', 
		vars = {{['0'] = '*value = buf[index];'}},
	}

if tryingAMR == 'dt vs 2dt' then
	self:addDisplayVarGroup{
		name = 'U2',
		type = self.eqn.cons_t,
		vars = {
			{[0] = '*value = buf[index].ptr[0];'},
		}
	}
elseif tryingAMR == 'gradient' then
	self:addDisplayVarGroup{
		name = 'amrError',
		type = 'real',
		vars = {
			{[0] = '*value = buf[index];'},
		}
	}
end
end

-- my best idea to work around the stupid 8-arg max kernel restriction
-- this is almost as bad of an idea as using OpenCL was to begin with
Solver.allocateOneBigStructure = false

function Solver:clalloc(name, size)
	self.buffers:insert{name=name, size=size}
end

function Solver:finalizeCLAllocs()
	local total = 0
	for _,buffer in ipairs(self.buffers) do
		buffer.offset = total
		local name = buffer.name
		local size = buffer.size
		if not self.allocateOneBigStructure then
			local mod = size % ffi.sizeof(self.app.env.real)
			if mod ~= 0 then
				-- WARNING?
				size = size - mod + ffi.sizeof(self.app.env.real)
				buffer.size = size
			end
		end
		total = total + size
		if not self.allocateOneBigStructure then
			local bufObj = CLBuffer{
				env = self.app.env,
				name = name,
				type = 'real',
				size = size / ffi.sizeof(self.app.env.real),
			}
			self[name..'Obj'] = bufObj
			self[name] = bufObj.obj
		end
	end
	if self.allocateOneBigStructure then
		self.oneBigBuf = self.app.ctx:buffer{rw=true, size=total}
	end
end

function Solver:getConsLRTypeCode()
	return template([[
typedef union {
	<?=eqn.cons_t?> LR[2];
	struct {
		<?=eqn.cons_t?> L, R;
	};
} <?=eqn.consLR_t?>;
]], {
	eqn = self.eqn,
})
end

function Solver:createBuffers()
	local realSize = ffi.sizeof(self.app.real)

	-- to get sizeof
	ffi.cdef(self.eqn:getTypeCode())
	ffi.cdef(self:getConsLRTypeCode())

	-- for twofluid, cons_t has been renamed to euler_maxwell_t and maxwell_cons_t
	if ffi.sizeof(self.eqn.cons_t) ~= self.eqn.numStates * realSize then
	   error('Expected sizeof('..self.eqn.cons_t..') to be '
		   ..self.eqn.numStates..' * sizeof(real) = '..(self.eqn.numStates * realSize)
		   ..' but found '..ffi.sizeof(self.eqn.cons_t)..' = '..(ffi.sizeof(self.eqn.cons_t) / realSize)..' * sizeof(real). '
		   ..'Maybe you need to update Eqn.numStates?')
	end

	-- should I put these all in one AoS?

if tryingAMR == 'gradient' then
	--[[
	ok here's my thoughts on the size ...
	I'm gonna try to do AMR
	that means storing the nodes in the same buffer
	I think I'll pad the end of the UBuf with the leaves
	and then store a tree of index information somewhere else that says what leaf goes where
	(maybe that will go at the end)
	
	how should memory breakdown look?
	how big should the leaves be?

	how about do like reduce ...
	leaves can be 16x16 blocks
	a kernel can cycle through the root leaves and sum te amrError 
	-- in the same kernel as it is calculated
	-- then I don't need to store so much memory, only one value per leaf, not per grid ...
	-- then ... split or merge leaves, based on their error ...

	how to modify that?  in anothe kernel of its own ...
	bit array of whether each cell is divided ...
	-- then in one kernel we update all leaves
	-- in another kernel, populate leaf ghost cells
	-- in another kernel, 
	-- 		decide what should be merged and, based on that, copy from parent into this
	--		if something should be split then look for unflagged children and copy into it
	
	how many bits will we need?
	volume / leafSize bits for the root level
	
	then we have parameters of 
	- how big each node is
	- what level of refinement it is
	
	ex:
	nodes in the root level are 2^4 x 2^4 = 2^8 = 256 cells = 256 bits = 2^5 = 32 bytes
	

	leafs are 2^4 x 2^4 = 2^8 = 256 cells
	... but ghost cells are 2 border, so we need to allocate (2^n+2*2)^2 cells ... for n=4 this is 400 cells ... 
		so we lose (2^n+2*2)^2 - 2^(2n) = 2^(2n) + 2^(n+3) + 2^4 - 2^(2n) = 2^(n+3) + 2^4) cells are lost
	
	so leafs multipy at a factor of 2^2 x 2^2 = 2^4 = 16
	so the next level has 2^(8+4) = 2^12 = 4096 bits = 2^9 = 512 bytes

	--]]

	-- the size, in cells, which a node replaces
	self.amrNodeFromSize = ({
		vec3sz(8, 1, 1),
		vec3sz(8, 8, 1),
		vec3sz(8, 8, 8),
	})[self.dim]
print('self.amrNodeFromSize', self.amrNodeFromSize)

	-- the size, in cells, of each node, excluding border, for each dimension
	self.amrNodeSizeWithoutBorder = ({
		vec3sz(16, 1, 1),
		vec3sz(16, 16, 1),
		vec3sz(16, 16, 16),
	})[self.dim]
print('self.amrNodeSizeWithoutBorder', self.amrNodeSizeWithoutBorder)

	-- size of the root level, in terms of nodes ('from' size)
	self.amrRootSizeInFromSize = vec3sz(1,1,1)
	for i=0,self.dim-1 do
		self.amrRootSizeInFromSize:ptr()[i] = 
			roundup(self.sizeWithoutBorder:ptr()[i], self.amrNodeFromSize:ptr()[i]) 
				/ self.amrNodeFromSize:ptr()[i]
	end
print('self.amrRootSizeInFromSize', self.amrRootSizeInFromSize)

	-- how big each node is
	self.amrNodeSize = self.amrNodeSizeWithoutBorder + 2 * self.numGhost

	-- how many nodes to allocate and use
	-- here's the next dilemma in terms of memory layout
	-- specifically in terms of rendering
	-- if I want to easily copy and render the texture information then I will need to package the leafs into the same texture as the root
	-- which means extending the texture buffer in some particular direction.
	-- since i'm already adding the leaf information to the end of the buffer
	-- and since appending to the end of a texture buffer coincides with adding extra rows to the texture
	-- why not just put our leafs in extra rows of -- both our display texture and of  
	self.amrMaxNodes = 1
	
	-- this will hold info on what leafs have yet been used
	self.amrLeafs = table()

	-- hmm, this is the info for the root node ...
	-- do I want to keep the root level data separate?
	-- or do I just want to represent everything as a collection of leaf nodes?
	-- I'll keep the root structure separate for now
	-- so I can keep the original non-amr solver untouched
	self.amrLayers = table()
	self.amrLayers[1] = table()	-- here's the root
end

	-- this much for the first level U data
	local UBufSize = self.volume * ffi.sizeof(self.eqn.cons_t)
if tryingARM == 'gradient' then	
	UBufSize = UBufSize + self.amrMaxNodes * self.amrNodeSize:volume()
end	
	self:clalloc('UBuf', UBufSize)


if tryingAMR == 'dt vs 2dt' then
	-- here's my start at AMR, using the 1989 Berger, Collela two-small-steps vs one-big-step method
	self:clalloc('lastUBuf', self.volume * ffi.sizeof(self.eqn.cons_t))
	self:clalloc('U2Buf', self.volume * ffi.sizeof(self.eqn.cons_t))
elseif tryingAMR == 'gradient' then
	
	-- this is going to be a single value for each leaf
	-- that means the destination will be the number of nodes it takes to cover the grid (excluding the border)
	-- however, do I want this to be a larger buffer, and then perform reduce on it?
	self:clalloc('amrErrorBuf', 
		-- self.volume 
		tonumber(self.amrRootSizeInFromSize:volume())
		* ffi.sizeof(self.app.real))
end

	if self.usePLM then
		self:clalloc('ULRBuf', self.volume * self.dim * ffi.sizeof(self.eqn.consLR_t))
	end

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar 
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', self.volume * realSize * 3)
	local reduceSwapBufSize = roundup(self.volume * realSize / self.localSize1d, realSize)
	self:clalloc('reduceSwapBuf', reduceSwapBufSize)
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- CL/GL interop

	-- hmm, notice I'm still keeping the numGhost border on my texture 
	-- if I remove the border altogether then I get wrap-around
	-- maybe I should just keep a border of 1?
	-- for now i'll leave it as it is
	local cl = self.dim < 3 and GLTex2D or GLTex3D
	self.tex = cl{
		width = tonumber(self.gridSize.x),
		height = tonumber(self.gridSize.y),
		depth = tonumber(self.gridSize.z),
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		--magFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
	}

	if self.app.useGLSharing then
		self.texCLMem = CLImageGL{context=self.app.ctx, tex=self.tex, write=true}
	else
		self.calcDisplayVarToTexPtr = ffi.new(self.app.real..'[?]', self.volume * 3)
		
		--[[ PBOs?
		self.calcDisplayVarToTexPBO = ffi.new('gl_int[1]', 0)
		gl.glGenBuffers(1, self.calcDisplayVarToTexPBO)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, self.calcDisplayVarToTexPBO[0])
		gl.glBufferData(gl.GL_PIXEL_UNPACK_BUFFER, self.tex.width * self.tex.height * ffi.sizeof(self.app.real) * 4, nil, gl.GL_STREAM_READ)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
		--]]
	end
end

function Solver:createCodePrefix()
	local lines = table()

	if self.dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end

	-- real3
	lines:insert(file['math.h'])
	lines:insert(template(file['math.cl']))

	lines:append{
		'#define geometry_'..self.geometry.name..' 1',
		'#define dim '..self.dim,
		'#define numGhost '..self.numGhost,
		'#define numStates '..self.eqn.numStates,
		'#define numIntStates '..self.eqn.numIntStates,
		'#define numWaves '..self.eqn.numWaves,
	}:append(xs:map(function(x,i)
	-- coordinate space = u,v,w
	-- cartesian space = x,y,z
	-- min and max in coordinate space
		return '#define mins_'..x..' '..clnumber(self.mins[i])..'\n'
			.. '#define maxs_'..x..' '..clnumber(self.maxs[i])..'\n'
	end)):append{
		'constant real3 mins = _real3(mins_x, '..(self.dim<2 and '0' or 'mins_y')..', '..(self.dim<3 and '0' or 'mins_z')..');', 
		'constant real3 maxs = _real3(maxs_x, '..(self.dim<2 and '0' or 'maxs_y')..', '..(self.dim<3 and '0' or 'maxs_z')..');', 
	}:append(xs:map(function(name,i)
	-- grid size
		return '#define gridSize_'..name..' '..tonumber(self.gridSize[name])
	end)):append{
		'constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);',
		'constant int4 stepsize = (int4)(1, gridSize_x, gridSize_x * gridSize_y, gridSize_x * gridSize_y * gridSize_z);',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		indexForInt4ForSize(i, gridSize_x, gridSize_y, gridSize_z)',
	
	--[[
	naming conventions ...
	* the grid indexes i_1..i_n that span 1 through gridSize_1..gridSize_n
	(between the index and the coordinate space:)
		- grid_dx? is the change in coordinate space wrt the change in grid space
		- cell_x(i) calculates the coordinates at index i
	* the coordinates space x_1..x_n that spans mins.s1..maxs.sn
	(between the coordinate and the embedded space:)
		- vectors can be calculated from Cartesian by cartesianToCoord
		- the length of the basis vectors wrt the change in indexes is given by dx?_at(x)
		- the Cartesian length of the holonomic basis vectors is given by coordHolBasisLen.  
			This is like dx?_at except not scaled by grid_dx?
			This is just the change in embedded wrt the change in coordinate, not wrt the change in grid
		- volume_at(x) gives the volume between indexes at the coordinate x
		- the Cartesian length of a vector in coordinate space is given by coordLen and coordLenSq
	* the embedded Cartesian space ... idk what letters I should use for this.  
		Some literature uses x^I vs coordinate space x^a, but that's still an 'x', no good for programming.
		Maybe I'll use 'xc' vs 'x', like I've already started to do in the initial conditions.
	--]]
	
	}:append(range(3):map(function(i)
	-- this is the change in coordinate wrt the change in code
	-- delta in coordinate space along one grid cell
		return (('#define grid_dx{i} ((maxs_{x} - mins_{x}) / (real)(gridSize_{x} - '..(2*self.numGhost)..'))')
			:gsub('{i}', i-1)
			:gsub('{x}', xs[i]))
	end)):append(range(3):map(function(i)
	-- mapping from index to coordinate 
		return (('#define cell_x{i}(i) ((real)(i + '..clnumber(.5-self.numGhost)..') * grid_dx{i} + mins_'..xs[i]..')')
			:gsub('{i}', i-1))
	end)):append{
		'#define cell_x(i) _real3(cell_x0(i.x), cell_x1(i.y), cell_x2(i.z))',

		self.eqn.getTypeCode and self.eqn:getTypeCode() or '',

		-- run here for the code, and in buffer for the sizeof()
		self.eqn:getEigenTypeCode() or '',
		
		-- bounds-check macro
		'#define OOB(lhs,rhs) (i.x < (lhs) || i.x >= gridSize_x - (rhs)'
			.. (self.dim < 2 and '' or ' || i.y < (lhs) || i.y >= gridSize_y - (rhs)')
			.. (self.dim < 3 and '' or ' || i.z < (lhs) || i.z >= gridSize_z - (rhs)')
			.. ')',
		
		template([[

// define i, index, and bounds-check
#define SETBOUNDS(lhs,rhs)	\
	int4 i = globalInt4(); \
	if (OOB(lhs,rhs)) return; \
	int index = INDEXV(i);
		
// same as above, except for kernels that don't use the boundary
// index operates on buffers of 'gridSize' (with border)
// but the kernel must be invoked across sizeWithoutBorder
<? local range = require'ext.range'
?>#define SETBOUNDS_NOGHOST() \
	int4 i = globalInt4(); \
	if (OOB(0,2*numGhost)) return; \
	i += (int4)(<?=range(4):map(function(i) 
		return i <= solver.dim and '2' or '0' 
	end):concat','?>); \
	int index = INDEXV(i);
]], {solver = self}),
	}

	lines:insert(self.geometry:getCode(self))

	lines:append{
		-- not messing with this one yet
		self.allocateOneBigStructure and '#define allocateOneBigStructure' or '',
		
		-- this is dependent on coord map / length code
		self.eqn:getCodePrefix() or '',
		self:getConsLRTypeCode(),
	}

	self.codePrefix = lines:concat'\n'

print'done building solver.codePrefix'
--print(self.codePrefix)
end

function Solver:refreshInitStateProgram()
	self.eqn.initState:refreshInitStateProgram(self)
end

function Solver:resetState()
	self.t = 0
	self.app.cmds:finish()
		
	-- start off by filling all buffers with zero, just in case initState forgets something ...
	for _,bufferInfo in ipairs(self.buffers) do
		self.app.cmds:enqueueFillBuffer{buffer=self[bufferInfo.name], size=bufferInfo.size}
	end

	self.eqn:resetState()
	
	self:boundary()
	if self.eqn.useConstrainU then
		self.constrainUKernelObj()
	end
	self.app.cmds:finish()

	for _,op in ipairs(self.ops) do
		op:resetState()
	end
end

function Solver:getCalcDTCode()
	if self.eqn.hasCalcDT then return end
	return template(file['solver/calcDT.cl'], {solver=self, eqn=self.eqn})
end

function Solver:refreshCommonProgram()
	-- code that depend on real and nothing else
	-- TODO move to app, along with reduceBuf

	local commonCode = table():append{
		self.codePrefix,
	}:append{
		template([[
kernel void multAdd(
	global <?=eqn.cons_t?>* a,
	const global <?=eqn.cons_t?>* b,
	const global <?=eqn.cons_t?>* c,
	real d
) {
	SETBOUNDS_NOGHOST();	
<? for i=0,eqn.numIntStates-1 do
?>	a[index].ptr[<?=i?>] = b[index].ptr[<?=i?>] + c[index].ptr[<?=i?>] * d;
<? end
?>}
]], 	{
			solver = self,
			eqn = self.eqn,
		})
	}:concat'\n'

	time('compiling common program', function()
		-- TODO rename :compile() to :build() to be like cl.program?
		self.commonProgramObj = self.Program{code=commonCode}
		self.commonProgramObj:compile()
	end)

	-- used by the integrators
	-- needs the same globalSize and localSize as the typical simulation kernels
	-- TODO exclude states which are not supposed to be integrated
	self.multAddKernelObj = self.commonProgramObj:kernel{name='multAdd', domain=self.domainWithoutBorder}

	self.reduceMin = self.app.env:reduce{
		size = self.volume,
		op = function(x,y) return 'min('..x..', '..y..')' end,
		initValue = 'INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceMax = self.app.env:reduce{
		size = self.volume,
		op = function(x,y) return 'max('..x..', '..y..')' end,
		initValue = '-INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceSum = self.app.env:reduce{
		size = self.volume,
		op = function(x,y) return x..' + '..y end,
		initValue = '0.',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
end

function Solver:getSolverCode()
	local fluxLimiterCode = 'real fluxLimiter(real r) {'
		.. self.app.limiters[1+self.fluxLimiter[0]].code 
		.. '}'

	local slopeLimiterCode = 'real slopeLimiter(real r) {'
		.. self.app.limiters[1+self.slopeLimiter[0]].code 
		.. '}'
	
	return table{
		self.codePrefix,
		
		-- TODO move to Roe, or FiniteVolumeSolver as a parent of Roe and HLL?
		self.eqn:getEigenCode() or '',
		
		fluxLimiterCode,
		slopeLimiterCode,
		
		'typedef struct { real min, max; } range_t;',
		self.eqn:getSolverCode() or '',

		self:getCalcDTCode() or '',
	
		-- messing with this ...
		self.usePLM and template(file['solver/plm.cl'], {solver=self, eqn=self.eqn}) or '',

		template(
			({
				['dt vs 2dt'] = [[
kernel void compareUvsU2(
	global <?=eqn.cons_t?>* U2Buf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?> *U2 = U2Buf + index;
	const global <?=eqn.cons_t?> *U = UBuf + index;
	
	//what to use to compare values ...
	//if we combine all primitives, they'll have to be appropriately weighted ...
	real sum = 0.;
	real tmp;
<? for i=0,eqn.numStates-1 do
?>	tmp = U2->ptr[<?=i?>] - U->ptr[<?=i?>]; sum += tmp * tmp;
<? end
?>	U2->ptr[0] = sum * 1e+5;
}
]],
				gradient = [==[
kernel void calcAMRError(
	global real* amrErrorBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	int4 nodei = globalInt4();
	if (nodei.x >= <?=solver.amrRootSizeInFromSize.x?> || 
		nodei.y >= <?=solver.amrRootSizeInFromSize.y?>) 
	{
		return;
	}

	int nodeIndex = nodei.x + <?=solver.amrRootSizeInFromSize.x?> * nodei.y;

	real dV_dx;	
	real sum = 0.;
	
	//hmm, it's less memory, but it's probably slower to iterate across all values as I build them here
	for (int nx = 0; nx < <?=solver.amrNodeFromSize.x?>; ++nx) {
		for (int ny = 0; ny < <?=solver.amrNodeFromSize.y?>; ++ny) {
			int4 Ui = (int4)(0,0,0,0);
			
			Ui.x = nodei.x * <?=solver.amrNodeFromSize.x?> + nx + numGhost;
			Ui.y = nodei.y * <?=solver.amrNodeFromSize.y?> + ny + numGhost;
			
			int Uindex = INDEXV(Ui);
			const global <?=eqn.cons_t?>* U = UBuf + Uindex;
				
	//TODO this wasn't the exact formula ...
	// and TODO make this modular.  some papers use velocity vector instead of density.  
	// why not total energy -- that incorporates everything?
<? for i=0,solver.dim-1 do
?>			dV_dx = (U[stepsize.s<?=i?>].rho - U[-stepsize.s<?=i?>].rho) / (2. * grid_dx<?=i?>);
			sum += dV_dx * dV_dx;
<? end
?>	
		}
	}
	amrErrorBuf[nodeIndex] = sum * 1e-2 * <?=clnumber(1/tonumber( solver.amrNodeFromSize:volume() ))?>;
}

//from is the position on the root level to read from
//to is which node to copy into
kernel void initNodeFromRoot(
	global <?=eqn.cons_t?>* UBuf,
	int4 from,
	int toNodeIndex
) {
	int4 i = (int4)(0,0,0,0);
	i.x = get_global_id(0);
	i.y = get_global_id(1);
	int dstIndex = i.x + numGhost + <?=solver.amrNodeSize.x?> * (i.y + numGhost);
	int srcIndex = from.x + (i.x>>1) + numGhost + gridSize_x * (from.y + (i.y>>1) + numGhost);

	global <?=eqn.cons_t?>* dstU = UBuf + <?=solver.volume?> + toNodeIndex * <?=solver.amrNodeSize:volume()?>;
	
	//blitter srcU sized solver.amrNodeFromSize (in a patch of size solver.gridSize)
	// to dstU sized solver.amrNodeSize (in a patch of solver.amrNodeSize)
	
	dstU[dstIndex] = UBuf[srcIndex];
}
]==],
			})[tryingAMR] or '', {
				solver = self,
				eqn = self.eqn,
				clnumber = clnumber,
			}),
	}:append(self.ops:map(function(op)
		return op:getSolverCode()
	end)):concat'\n'
end

-- depends on buffers
function Solver:refreshSolverProgram()
	-- set pointer to the buffer holding the LR state information
	-- for piecewise-constant that is the original UBuf
	-- for piecewise-linear that is the ULRBuf
	self.getULRBuf = self.usePLM and self.ULRBuf or self.UBuf

	self.getULRArg = self.usePLM 
		and ('const global '..self.eqn.consLR_t..'* ULRBuf')
		or ('const global '..self.eqn.cons_t..'* UBuf')

	-- this code creates the const global cons_t* UL, UR variables
	-- it assumes that indexL, indexR, and side are already defined
	self.getULRCode = self.usePLM and ([[
	const global ]]..self.eqn.cons_t..[[* UL = &ULRBuf[side + dim * indexL].R;
	const global ]]..self.eqn.cons_t..[[* UR = &ULRBuf[side + dim * indexR].L;
]]) or ([[
	const global ]]..self.eqn.cons_t..[[* UL = UBuf + indexL;
	const global ]]..self.eqn.cons_t..[[* UR = UBuf + indexR;
]])

	local code = self:getSolverCode()

	time('compiling solver program', function()
		self.solverProgramObj = self.Program{code=code}
		self.solverProgramObj:compile()
	end)

	self:refreshCalcDTKernel()

	if self.eqn.useConstrainU then
		self.constrainUKernelObj = self.solverProgramObj:kernel('constrainU', self.UBuf)
	end

	if self.usePLM then
		self.calcLRKernelObj = self.solverProgramObj:kernel('calcLR', self.ULRBuf, self.UBuf)
	end

if tryingAMR == 'dt vs 2dt' then
	self.compareUvsU2KernelObj = self.solverProgramObj:kernel('compareUvsU2', self.U2Buf, self.UBuf)
elseif tryingAMR == 'gradient' then
	self.calcAMRErrorKernelObj = self.solverProgramObj:kernel('calcAMRError', self.amrErrorBuf, self.UBuf)
	self.initNodeFromRootKernelObj = self.solverProgramObj:kernel('initNodeFromRoot', self.UBuf)
end

	for _,op in ipairs(self.ops) do
		op:refreshSolverProgram()
	end
end

-- for solvers who don't rely on calcDT
function Solver:refreshCalcDTKernel()
	self.calcDTKernelObj = self.solverProgramObj:kernel('calcDT', self.reduceBuf, self.UBuf)
end

function Solver:refreshDisplayProgram()

	local lines = table{
		self.codePrefix,
	}
	
	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			var.id = tostring(var):sub(10)
		end
	end

	if self.app.useGLSharing then
		for _,displayVarGroup in ipairs(self.displayVarGroups) do
			for _,var in ipairs(displayVarGroup.vars) do
				if var.enabled 
				or (var.vecVar and var.vecVar.enabled)
				then
					lines:append{
						template(var.displayCode, {
							solver = self,
							var = var,
							name = 'calcDisplayVarToTex_'..var.id,
							input = 'write_only '
								..(self.dim == 3 
									and 'image3d_t' 
									or 'image2d_t'
								)..' tex',
							output = '	write_imagef(tex, '
								..(self.dim == 3 
									and '(int4)(dsti.x, dsti.y, dsti.z, 0)' 
									or '(int2)(dsti.x, dsti.y)'
								)..', (float4)(value[0], value[1], value[2], 0.));',
						})
					}
				end
			end
		end
	end

	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			if var.enabled 
			or (var.vecVar and var.vecVar.enabled)
			then
				lines:append{
					template(var.displayCode, {
						solver = self,
						var = var,
						name = 'calcDisplayVarToBuffer_'..var.id,
						input = 'global real* dest',
						output = var.vectorField and [[
	dest[0+3*dstindex] = valuevec->x;
	dest[1+3*dstindex] = valuevec->y;
	dest[2+3*dstindex] = valuevec->z;
]] or [[
	dest[dstindex] = value[0];
]],
					})
				}
			end
		end
	end
	
	local code = lines:concat'\n'
	time('compiling display program', function()
		self.displayProgramObj = self.Program{code=code}
		self.displayProgramObj:compile()
	end)

	if self.app.useGLSharing then
		for _,displayVarGroup in ipairs(self.displayVarGroups) do
			for _,var in ipairs(displayVarGroup.vars) do
				if var.enabled 
				or (var.vecVar and var.vecVar.enabled)
				then
					var.calcDisplayVarToTexKernelObj = self.displayProgramObj:kernel('calcDisplayVarToTex_'..var.id, self.texCLMem)
				end
			end
		end
	end

	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			if var.enabled 
			or (var.vecVar and var.vecVar.enabled)
			then
				var.calcDisplayVarToBufferKernelObj = self.displayProgramObj:kernel('calcDisplayVarToBuffer_'..var.id, self.reduceBuf)
			end
		end
	end
end

--[[
args:
	type = type of buffer to create boundary information
	extraArgs = any extra args to add to the kernel (other than a global of the buffer type)
	methods = {
		[x|y|z..min|max] = 
			either 'peroidic', 'mirror', 'freeflow',
			or a function of the name of variable pointer 
			that returns the code to apply to that particular pointer 
	}
	mirrorVars = {x vars, y vars, z vars} = table of tables of strings of what fields should be negative'd on mirror condition
		this tends to be vectors' components that point into the boundary.
--]]
function Solver:createBoundaryProgramAndKernel(args)
	local assign = args.assign or function(a, b) return a .. ' = ' .. b end
	
	local lines = table()
	lines:insert(self.codePrefix)
	lines:insert(template([[
kernel void boundary(
	global <?=args.type?>* buf
<?= args.extraArgs and #args.extraArgs > 0 
	and ','..table.concat(args.extraArgs, ',\n\t')
	or '' ?>
) {
]], {args = args}))
	if self.dim == 2 then
		lines:insert[[
	int i = get_global_id(0);
]]
	elseif self.dim == 3 then
		lines:insert[[
	int2 i = (int2)(get_global_id(0), get_global_id(1));
]]
	end
	
	-- 1D: use a small 1D kernel and just run once 
	-- 2D: use a 1D kernel the size of the max dim 
		
	lines:insert[[
	for (int j = 0; j < numGhost; ++j) {
]]

	for side=1,self.dim do

		if self.dim == 2 then
			if side == 1 then
				lines:insert[[
		if (i < gridSize_y) {
]]
			elseif side == 2 then
				lines:insert[[
		if (i < gridSize_x) {
]]
			end
		elseif self.dim == 3 then
			if side == 1 then
				lines:insert[[
		if (i.x < gridSize_y && i.y < gridSize_z) {
]]
			elseif side == 2 then
				lines:insert[[
		if (i.x < gridSize_x && i.y < gridSize_z) {
]]
			elseif side == 3 then
				lines:insert[[
		if (i.x < gridSize_x && i.y < gridSize_y) {
]]
			end
		end

		local function index(j)
			if self.dim == 1 then
				return j
			elseif self.dim == 2 then
				if side == 1 then
					return 'INDEX('..j..',i,0)'
				elseif side == 2 then
					return 'INDEX(i,'..j..',0)'
				end
			elseif self.dim == 3 then
				if side == 1 then
					return 'INDEX('..j..',i.x,i.y)'
				elseif side == 2 then
					return 'INDEX(i.x,'..j..',i.y)'
				elseif side == 3 then
					return 'INDEX(i.x,i.y,'..j..')'
				end
			end
		end
	
		for _,minmax in ipairs{'min', 'max'} do
			local method = args.methods[xs[side]..minmax]
			lines:insert(method{
				index = index,
				assign = assign,
				side = side,
				mirrorVars = args.mirrorVars,
				minmax = minmax,
			})
		end
	
		if self.dim > 1 then
			lines:insert'		}'
		end
	end
	
	lines:insert'	}'
	lines:insert'}'

	local code = lines:concat'\n'

	local boundaryProgramObj
	time('compiling boundary program', function()
		boundaryProgramObj = self.Program{code=code}
		boundaryProgramObj:compile()
	end)
	local boundaryKernelObj = boundaryProgramObj:kernel'boundary'
	
	-- TODO switch these over to obj
	return boundaryProgramObj, boundaryKernelObj
end

function Solver:getBoundaryProgramArgs()
	return {
		type = self.eqn.cons_t,
		-- remap from enum/combobox int values to functions from the solver.boundaryOptions table
		methods = table.map(self.boundaryMethods, function(v)
			if type(v) == 'function' then return v end
			return (select(2, next(self.boundaryOptions[v])))
		end),
		mirrorVars = self.eqn.mirrorVars,
	}
end

function Solver:refreshBoundaryProgram()
	self.boundaryProgramObj, self.boundaryKernelObj = 
		self:createBoundaryProgramAndKernel(self:getBoundaryProgramArgs())
	self.boundaryKernelObj.obj:setArg(0, self.UBuf)
	for _,op in ipairs(self.ops) do
		op:refreshBoundaryProgram()
	end
end

-- assumes the buffer is already in the kernel's arg
function Solver:applyBoundaryToBuffer(kernelObj)
	-- 1D:
	if self.dim == 1 then
		-- if you do change this size from anything but 1, make sure to add conditions to the boundary kernel code
		self.app.cmds:enqueueNDRangeKernel{
			kernel = kernelObj.obj,
			globalSize = 1,
			localSize = 1,
		}
	elseif self.dim == 2 then
		-- TODO what happens when kernelObj.maxWorkGroupSize 
		-- is lower than the grid size?
		local maxSize = math.min(
				roundup(
					math.max(
						tonumber(self.gridSize.x),
						tonumber(self.gridSize.y)),
				self.localSize1d),
			kernelObj.maxWorkGroupSize)
		local localSize = math.min(self.localSize1d, maxSize)
		self.app.cmds:enqueueNDRangeKernel{
			kernel = kernelObj.obj,
			globalSize = maxSize, 
			localSize = localSize,
		}
	elseif self.dim == 3 then
		-- xy xz yz
		local maxSizeX = roundup(
			math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y)),
			self.localSize2d[1])
		local maxSizeY = roundup(
			math.max(tonumber(self.gridSize.y), tonumber(self.gridSize.z)),
			self.localSize2d[2])
		self.app.cmds:enqueueNDRangeKernel{
			kernel = kernelObj.obj,
			globalSize = {maxSizeX, maxSizeY},
			localSize = {
				math.min(kernelObj.localSize2d[1], maxSizeX),
				math.min(kernelObj.localSize2d[2], maxSizeY),
			},
		}
	else
		error("can't run boundary for dim "..tonumber(self.dim))
	end
end

function Solver:boundary()
	self:applyBoundaryToBuffer(self.boundaryKernelObj)
end

function Solver:calcDT()
	local dt
	-- calc cell wavespeeds -> dts
	if self.useFixedDT then
		dt = self.fixedDT
	else
		-- TODO this without the border, but that means changing reduce *and display*
		self.calcDTKernelObj()
		dt = self.cfl * self.reduceMin()
		if not math.isfinite(dt) then
			print("got a bad dt!") -- TODO dump all buffers
		end
		self.fixedDT = dt
	end
	return dt
end

Solver.fpsNumSamples = 30

function Solver:update()
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


	if self.checkNaNs then
		if self:checkFinite(self.UBufObj, self.volume) then return end
	end

--local before = self.UBufObj:toCPU()
--print'\nself.UBufObj before boundary:' self:printBuf(self.UBufObj) print(debug.traceback(),'\n\n')	
	self:boundary()
--[[
local after = self.UBufObj:toCPU()
print'\nself.UBufObj after boundary:' self:printBuf(self.UBufObj) print(debug.traceback(),'\n\n')	
if self.app.dim == 2 then
	for i=0,tonumber(self.gridSize.x)-1 do
		for j=0,self.numGhost-1 do
			for _,s in ipairs{0,4} do
				for _,offset in ipairs{
					s + self.eqn.numStates * (i + self.gridSize.x * j),
					s + self.eqn.numStates * (j + self.gridSize.x * i),
					s + self.eqn.numStates * (i + self.gridSize.x * (self.gridSize.y - 1 - j)),
					s + self.eqn.numStates * ((self.gridSize.x - 1 - j) + self.gridSize.x * i),
				} do
					if after[offset] < 0 then
						local msg = "negative found in boundary after :boundary() was called"
						msg = msg .. "\nat position " .. tostring(offset / self.eqn.numStates)
						msg = msg .. '\nit was '..(before[offset] < 0 and 'negative' or 'positive')..' before the boundary update'
						error(msg)
					end
				end
			end
		end
	end
end
--]]	
	local dt = self:calcDT()

if tryingAMR == 'dt vs 2dt' then
	-- back up the last buffer
	self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.lastUBuf, size=self.volume * self.eqn.numStates * ffi.sizeof(self.app.real)}
end
	
	-- first do a big step
	self:step(dt)

	-- now copy it to the backup buffer
if tryingAMR == 'dt vs 2dt' then
	-- TODO have step() provide a target, and just update directly into U2Buf?
	self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.U2Buf, size=self.volume * self.eqn.numStates * ffi.sizeof(self.app.real)}
	self.app.cmds:enqueueCopyBuffer{src=self.lastUBuf, dst=self.UBuf, size=self.volume * self.eqn.numStates * ffi.sizeof(self.app.real)}

	local t = self.t
	self:step(.5 * dt)
	self.t = t + .5 * dt
	self:step(.5 * dt)
	self.t = t + dt

	-- now compare UBuf and U2Buf, store in U2Buf in the first real of cons_t
	self.compareUvsU2KernelObj()
elseif tryingAMR == 'gradient' then
	
	-- 1) compute errors from gradient, sum up errors in each root node, and output on a per-node basis
	local amrRootSizeInFromGlobalSize = vec3sz(
		roundup(self.amrRootSizeInFromSize.x, self.localSize.x),
		roundup(self.amrRootSizeInFromSize.y, self.localSize.y),
		roundup(self.amrRootSizeInFromSize.z, self.localSize.z))
	
	self.app.cmds:enqueueNDRangeKernel{
		kernel = self.calcAMRErrorKernelObj.obj, 
		dim = self.dim, 
		globalSize = amrRootSizeInFromGlobalSize:ptr(), 
		localSize = self.localSize:ptr(),
	}

	-- 2) based on what nodes' errors are past some value, split or merge...
	--[[
	1) initial tree will have nothing flagged as split
	2) then we get some split data - gradients -> errors -> thresholds -> flags 
		... which are lined up with the layout of the patches ...
		... which doesn't necessarily match the tree structure ...
	3) look through all used patches' error thresholds, min and max
		if it says to split ... 
			then look and see if we have room for any more free leafs in our state buffer
		
			the first iteration will request to split on some cells
			so go through the error buffer for each (root?) node,
			see if the error is bigger than some threshold then this node needs to be split
				then we have to add a new leaf node
			
			so i have to hold a table of what in the U extra leaf buffer is used
			which means looking
		
		if it says to merge ...
			clear the 'used' flag in the overall tree / in the layout of leafs in our state buffer
	--]]
	local vol = tonumber(self.amrRootSizeInFromSize:volume())
	local ptr = ffi.new('real[?]', vol)
	self.app.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * vol, ptr=ptr}
-- [[
print'armErrors:'
for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
	for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
		local i = nx + self.amrRootSizeInFromSize.x * ny
		io.write('\t', ('%.5f'):format(ptr[i]))
	end
	print()
end
--]]

for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
	for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
		local i = nx + self.amrRootSizeInFromSize.x * ny
		local nodeErr = ptr[i]
		if nodeErr > .2 then
			print('root node '..tostring(i)..' needs to be split')
			
			-- flag for a split
			-- look for a free node to allocate in the buffer
			-- if there's one available then ...
			-- store it in a map

			local amrLayer = self.amrLayers[1]
			
			-- see if there's an entry in this layer
			-- if there's not then ...
			-- allocate a new patch and make an entry
			if not amrLayer[i+1] then
			
				-- next: find a new unused leaf
				-- for now, just this one node
				local leafIndex = 0
				
				if not self.amrLeafs[leafIndex+1] then
			
					print('splitting root node '..tostring(i)..' and putting in leaf node '..leafIndex)
		
					-- create info about the leaf
					self.amrLeafs[leafIndex+1] = {
						level = 0,	-- root
						layer = amrLayer,
						layerX = nx,		-- node x and y in the root
						layerY = ny,
						leafIndex = leafIndex,	-- which leaf we are using
					}
			
					-- tell the root layer table which node is used
					--  by pointing it back to the table of the leaf nodes 
					amrLayer[i+1] = self.amrLeafs[1]
				
					-- copy data from the root node location into the new node
					-- upsample as we go ... by nearest?
				
					-- TODO setup kernel args
					self.initNodeFromRootKernelObj.obj:setArg(1, ffi.new('int[4]', {nx, ny, 0, 0}))
					self.initNodeFromRootKernelObj.obj:setArg(2, ffi.new('int[1]', 0))
					self.app.cmds:enqueueNDRangeKernel{
						kernel = self.initNodeFromRootKernelObj.obj,
						dim = self.dim, 
						globalSize = self.amrNodeSizeWithoutBorder:ptr(),
						localSize = self.amrNodeSizeWithoutBorder:ptr(),
					}
				end
			end
		end
	end
end

os.exit()
	
	self.t = self.t + dt
	self.dt = dt
else
	self.t = self.t + dt
	self.dt = dt
end

end

function Solver:step(dt)
	self.integrator:integrate(dt, function(derivBuf)
		self:calcDeriv(derivBuf, dt)
	end)
	
	if self.eqn.useConstrainU then
		self:boundary()
		self.constrainUKernelObj()
	end

	for _,op in ipairs(self.ops) do
		if op.step then
			op:step(dt)
		end
	end
end

function Solver:printBuf(buf, ptr)
	ptr = ptr or buf:toCPU()
	local max = #tostring(self.volume-1)
	for i=0,self.volume-1 do
		io.write((' '):rep(max-#tostring(i)), i,':')
		for j=0,self.eqn.numStates-1 do
			io.write(' ', ptr[j + self.eqn.numStates * i])
		end 
		print()
	end
end

-- check for nans
-- expects buf to be of type cons_t, made up of numStates real variables
function Solver:checkFinite(buf, volume)
	local ptr = buf:toCPU()
	local found
	for i=0,buf.size-1 do
		if not math.isfinite(ptr[i]) then
			found = found or table()
			found:insert(i)
		end
	end
	if not found then return end
--	self:printBuf(nil, ptr)
--	print(found:map(tostring):concat', ')
--	error'found non-finite numbers'
	return true
end

function Solver:getTex(var) 
	return self.tex
end

function Solver:calcDisplayVarToTex(var)
	local app = self.app
	local displayVarGroup = var.displayVarGroup
	if app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		app.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	
		var:setToTexArgs()
		var.calcDisplayVarToTexKernelObj()
		
		app.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		app.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		-- calcDisplayVarToTexPtr is sized without ghost cells 
		-- so is the GL texture
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex

		local channels = var.vectorField and 3 or 1
		local format = var.vectorField and gl.GL_RGB or gl.GL_RED

		var:setToBufferArgs()
		var.calcDisplayVarToBufferKernelObj()
		app.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * self.volume * channels, ptr=ptr}
		local destPtr = ptr
		if app.is64bit then
			-- can this run in place?
			destPtr = ffi.cast('float*', ptr)
			for i=0,self.volume*channels-1 do
				destPtr[i] = ptr[i]
			end
		end
		tex:bind()
		if self.dim < 3 then
			gl.glTexSubImage2D(gl.GL_TEXTURE_2D, 0, 0, 0, tex.width, tex.height, format, gl.GL_FLOAT, destPtr)
		else
			for z=0,tex.depth-1 do
				gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, tex.width, tex.height, 1, format, gl.GL_FLOAT, destPtr + channels * tex.width * tex.height * z)
			end
		end
		tex:unbind()
		glreport'here'
	end
end

-- used by the display code to dynamically adjust ranges
function Solver:calcDisplayVarRange(var)
	if var.lastTime == self.t then
		return var.lastMin, var.lastMax
	end
	var.lastTime = self.t

	var:setToBufferArgs()

	var.calcDisplayVarToBufferKernelObj()
	local min = self.reduceMin()
	var.calcDisplayVarToBufferKernelObj()
	local max = self.reduceMax()
	
	var.lastMin = min
	var.lastMax = max
	
	return min, max
end

-- used by the output to print out avg, min, max
function Solver:calcDisplayVarRangeAndAvg(var)
	local needsUpdate = var.lastTime ~= self.t

	-- this will update lastTime if necessary
	local min, max = self:calcDisplayVarRange(var)
	-- displayVarGroup has already set up the appropriate args

	local avg
	if needsUpdate then
		var.calcDisplayVarToBufferKernelObj()
		avg = self.reduceSum(nil, self.volume) / tonumber(self.volume)
	else
		avg = var.lastAvg
	end

	return min, max, avg
end

function Solver:updateGUIParams()
	ig.igText('t: '..self.t)

	-- hmm put fps somewhere else, or put ms update here
	ig.igText('fps: '..(self.fps and tostring(self.fps) or ''))
	
	if ig.igCollapsingHeader'parameters:' then

		tooltip.checkboxTable('check NaNs', self, 'checkNaNs')

		tooltip.checkboxTable('use fixed dt', self, 'useFixedDT')
		ig.igSameLine()
		
		tooltip.numberTable('fixed dt', self, 'fixedDT')
		tooltip.numberTable('CFL', self, 'cfl')

		for i=1,self.dim do
			for _,minmax in ipairs(minmaxs) do
				local k = xs[i]..minmax
				if tooltip.numberTable(k, self[minmax..'s'], i, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
					-- if the domain changes
					-- then the dx has to change
					-- and all the stuff based on codePrefix has to change
					self:refreshCodePrefix()
				end
			end
		end

		if tooltip.comboTable('integrator', self, 'integratorIndex', self.integratorNames) then
			self:refreshIntegrator()
		end

		-- I think I'll display my GMRES # steps to converge / epsilon error ... 
		if self.integrator.updateGUI then
			ig.igSameLine()
			ig.igPushIdStr'integrator'
			if ig.igCollapsingHeader'' then
				self.integrator:updateGUI()
			end
			ig.igPopId()
		end

		for _,op in ipairs(self.ops) do
			op:updateGUI()
		end

		if tooltip.combo('slope limiter', self.fluxLimiter, self.app.limiterNames) then
			self:refreshSolverProgram()
		end

		for i=1,self.dim do
			for _,minmax in ipairs(minmaxs) do
				local var = xs[i]..minmax
				if tooltip.comboTable(var, self.boundaryMethods, var, self.boundaryOptionNames) then
					self:refreshBoundaryProgram()
				end
			end
		end
	end
end

function Solver:updateGUIEqnSpecific()
	if ig.igCollapsingHeader'equation:' then
		if tooltip.comboTable('init state', self, 'initStateIndex', self.eqn.initStateNames) then
			-- TODO hmm ... the whole point of making a separate initStateProgram was to be able to refresh it without rebuilding all of the solver ...
			self:refreshEqnInitState()
		end	
		
		for _,var in ipairs(self.eqn.guiVars) do
			var:updateGUI(self)
		end
	end
end

do
	-- display vars: TODO graph vars
	local function handle(var, title)
		ig.igPushIdStr(title)
		
		local enableChanged = tooltip.checkboxTable('enabled', var, 'enabled') 
		ig.igSameLine()
		
		tooltip.checkboxTable('log', var, 'useLog')
		ig.igSameLine()

		tooltip.checkboxTable('fixed range', var, 'heatMapFixedRange')
		ig.igSameLine()
		
		if ig.igCollapsingHeader(var.name) then
			tooltip.numberTable('value min', var, 'heatMapValueMin')
			tooltip.numberTable('value max', var, 'heatMapValueMax')
		end
		
		ig.igPopId()
	
		return enableChanged
	end

	-- do one for 'all'
	local function _and(a,b) return a and b end
	local fields = {'enabled', 'useLog', 'heatMapFixedRange', 'heatMapValueMin', 'heatMapValueMax'}
	local defaults = {true, true, true, math.huge, -math.huge}
	local combines = {_and, _and, _and, math.min, math.max}
	local all = {name='all'}
	for i=1,#fields do
		all[fields[i]] = defaults[i]
	end
	local original = {}
	for i,field in ipairs(fields) do
		original[field] = defaults[i]
	end

	function Solver:updateGUIDisplay()
		local refresh 
		if ig.igCollapsingHeader'display:' then
			for i,displayVarGroup in ipairs(self.displayVarGroups) do
				ig.igPushIdStr('display '..i)
				if ig.igCollapsingHeader(displayVarGroup.name) then				
					for i=1,#fields do
						all[fields[i]] = defaults[i]
					end
					for _,var in ipairs(displayVarGroup.vars) do
						for i,field in ipairs(fields) do
							all[field] = combines[i](all[field], var[field])
						end
					end
					for _,field in ipairs(fields) do
						original[field] = all[field]
					end
					local enableChanged = handle(all, 'all')
					refresh = refresh or enableChanged
					for _,field in ipairs(fields) do
						if all[field] ~= original[field] then
							for _,var in ipairs(displayVarGroup.vars) do
								var[field] = all[field]
							end
						end
					end

					for _,var in ipairs(displayVarGroup.vars) do
						local enableChanged = handle(var, displayVarGroup.name..' '..var.name)
						refresh = refresh or enableChanged
					end
				end
				ig.igPopId()
			end
		end
		if refresh then
			self:refreshDisplayProgram()
		end
	end
end

function Solver:updateGUI()
	self:updateGUIParams()
	self:updateGUIEqnSpecific()
	self:updateGUIDisplay()

	-- heat map var

	-- TODO volumetric var
end

local Image = require 'image'
function Solver:save(prefix)
	-- TODO add planes to image, then have the FITS module use planes and not channels
	-- so the dimension layout of the buffer is [channels][width][height][planes]
	local width = tonumber(self.gridSize.x)
	local height = tonumber(self.gridSize.y)
	local depth = tonumber(self.gridSize.z)

	for _,bufferInfo in ipairs(self.buffers) do
		local name = bufferInfo.name
		local channels = bufferInfo.size / self.volume / ffi.sizeof(self.app.real)
		if channels ~= math.floor(channels) then
			print("can't save buffer "..name.." due to its size not being divisible by the solver volume")
		else
			local buffer = self[name]

			local numReals = self.volume * channels
--[[ 3D interleave the depth and the channels ... causes FV to break 
			local image = Image(width, height, depth * channels, assert(self.app.real))
			local function getIndex(ch,i,j,k) return i + width * (j + height * (ch + channels * k)) end
--]]
-- [[ 3D interleave the depth and the width
			local image = Image(width * depth, height, channels, assert(self.app.real))
			local function getIndex(ch,i,j,k) return i + width * (k + depth * (j + height * ch)) end
--]]
			self.app.cmds:enqueueReadBuffer{buffer=buffer, block=true, size=ffi.sizeof(self.app.real) * numReals, ptr=image.buffer}
			local src = image.buffer
			
			-- now convert from interleaved to planar
			-- *OR* add planes to the FITS output
			local tmp = ffi.new(self.app.real..'[?]', numReals)
			for ch=0,channels-1 do
				for k=0,depth-1 do
					for j=0,height-1 do
						for i=0,width-1 do
							tmp[getIndex(ch,i,j,k)] = src[ch + channels * (i + width * (j + height * k))]
						end
					end
				end
			end
			image.buffer = tmp
		
			local filename = prefix..'_'..name..'.fits'
			print('saving '..filename)
			image:save(filename)
		end
	end
end

return Solver
