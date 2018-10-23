-- TODO make this solver/solver.lua, and make the old solver.lua something like structured-grid-solver

local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local math = require 'ext.math'
local CLBuffer = require 'cl.obj.buffer'
local template = require 'template'
local vec3 = require 'vec.vec3'
local roundup = require 'roundup'
local tooltip = require 'tooltip'
local makestruct = require 'eqn.makestruct'
local time, getTime = table.unpack(require 'time')


local common = require 'common'()	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6
local from6to3x3 = common.from6to3x3
local sym = common.sym


local integrators = require 'int.all'
local integratorNames = integrators:map(function(integrator) return integrator.name end)


local SolverBase = class()

SolverBase.name = 'Solver'

-- override to specify which eqn/*.lua to use as the equation
SolverBase.eqnName = nil

-- enable for us to let the user accum/not.  requires an extra buffer allocation
SolverBase.allowAccum = true
SolverBase.displayVarAccumFunc = false

--[[
args:
	app
	dim
	eqn = name of eqn/<name>.lua to use
	eqnArgs (optional) = args of eqn ctor
	coord = coordinate system name to use, associated with coord/<name>.lua
	initState = init state name
	initStateArgs (optional) = args of init state ctor
	integrator = name of integrator in int/all.lua to use
--]]
function SolverBase:init(args)
	self:initL1(args)
	self:preInit(args)
	self:postInit()
end

function SolverBase:initL1(args)
	assert(args)
	self.app = assert(args.app)
	self.dim = assert(args.dim)
	
	self.color = vec3(math.random(), math.random(), math.random()):normalize()

	-- operators for this solver
	self.ops = table()
end

function SolverBase:getSolverTypeCode()
	error("not implemented in base class")
end


function SolverBase:preInit(args)

-- despite the name, this doesn't have anything to do with the grid size ...
-- ... except in the GridSolver class
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	self.maxWorkGroupSize = tonumber(self.app.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')
print('maxWorkGroupSize', self.maxWorkGroupSize)
	
	local sizeProps = self:getSizePropsForWorkGroupSize(self.maxWorkGroupSize)
	for k,v in pairs(sizeProps) do
print(k,v)
		self[k] = v
	end
	
	-- hmm, do some eqns create ops need to know the grid size?
	self.eqnName = args.eqn
	self.eqnArgs = args.eqnArgs
	self:createEqn()

	self.name = self.eqn.name..' '..self.name

	self.initStateIndex = table.find(self.eqn.initStateNames, args.initState) or 1
	self.initStateArgs = args.initStateArgs

	self.integratorIndex = integratorNames:find(args.integrator) or 1
	
	if require 'coord.coord'.is(args.coord) then
		self.coord = args.coord	-- ptr copy expected by AMR
	else
		self.coord = require('coord.'..args.coord)(table({solver=self}, args.coordArgs))
	end

	self.checkNaNs = false
	self.useFixedDT = not not args.fixedDT
	self.fixedDT = args.fixedDT or self.fixedDT or .001
	self.cfl = args.cfl or .5	--/self.dim
	self.fluxLimiter = self.app.limiterNames:find(args.fluxLimiter) or 1


	self:createDisplayVars()	-- depends on eqn


	-- do this before any call to createBuffers or createCodePrefix
	self.solver_t = require 'eqn.makestruct'.uniqueName'solver_t'
	require'eqn.makestruct'.safeFFICDef(self:getSolverTypeCode())
	self.solverPtr = ffi.new(self.solver_t)
end

function SolverBase:postInit()
	self:refreshGridSize()
end

function SolverBase:createSolverBuf()
	-- doesn't use clalloc ...
	-- should happen before any other buffer allocs
	self.solverBuf = CLBuffer{
		env = self.app.env,
		name = 'solver',
		type = self.solver_t,
		size = 1,	-- should be 'count'
		readwrite = 'read',
	}
end

function SolverBase:refreshGridSize()

	self:createSolverBuf()
	
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

function SolverBase:refreshCommonProgram()
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
	realparam d
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
		self.commonProgramObj = self.Program{name='common', code=commonCode}
		self.commonProgramObj:compile()
	end)

	-- used by the integrators
	-- needs the same globalSize and localSize as the typical simulation kernels
	-- TODO exclude states which are not supposed to be integrated
	self.multAddKernelObj = self.commonProgramObj:kernel{name='multAdd', domain=self.domainWithoutBorder}

	self.reduceMin = self.app.env:reduce{
		size = self.numCells,
		op = function(x,y) return 'min('..x..', '..y..')' end,
		initValue = 'INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceMax = self.app.env:reduce{
		size = self.numCells,
		op = function(x,y) return 'max('..x..', '..y..')' end,
		initValue = '-INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceSum = self.app.env:reduce{
		size = self.numCells,
		op = function(x,y) return x..' + '..y end,
		initValue = '0.',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
end

-- here's another function that needs to be renamed ...
function SolverBase:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	error'abstract'
end

function SolverBase:clalloc(name, size, sizevec)
	size = math.ceil(size / 4) * 4
	self.buffers:insert{name=name, size=size, sizevec=sizevec}
end

function SolverBase:finalizeCLAllocs()
	local total = 0
	for _,info in ipairs(self.buffers) do
		info.offset = total
		local name = info.name
		local size = info.size
		if not self.allocateOneBigStructure then
			local mod = size % ffi.sizeof(self.app.env.real)
			if mod ~= 0 then
				-- WARNING?
				size = size - mod + ffi.sizeof(self.app.env.real)
				info.size = size
			end
		end
		total = total + size
		if not self.allocateOneBigStructure then
			local bufObj = CLBuffer{
				env = self.app.env,
				name = name,
				type = 'real',
				size = size / ffi.sizeof(self.app.real),
			}
			self[name..'Obj'] = bufObj
			self[name] = bufObj.obj
			
			-- I don't know where else to put this ...
			self[name].sizevec = info.sizevec
		end
	end
	if self.allocateOneBigStructure then
		self.oneBigBuf = self.app.ctx:buffer{rw=true, size=total}
	end
end


-- TODO this matches GridSolver very closely.  merge somehow?
function SolverBase:createBuffers()
	local realSize = ffi.sizeof(self.app.real)
	
	-- for twofluid, cons_t has been renamed to euler_maxwell_t and maxwell_cons_t
	if ffi.sizeof(self.eqn.cons_t) ~= self.eqn.numStates * realSize then
	   error('Expected sizeof('..self.eqn.cons_t..') to be '
		   ..self.eqn.numStates..' * sizeof(real) = '..(self.eqn.numStates * realSize)
		   ..' but found '..ffi.sizeof(self.eqn.cons_t)..' = '..(ffi.sizeof(self.eqn.cons_t) / realSize)..' * sizeof(real). '
		   ..'Maybe you need to update Eqn.numStates?')
	end
	if ffi.sizeof(self.eqn.cons_t) < ffi.sizeof(self.eqn.prim_t) then
		error("for PLM's sake I might need sizeof(prim_t) <= sizeof(cons_t)")
	end
	
	local UBufSize = self.numCells * ffi.sizeof(self.eqn.cons_t)
	self:clalloc('UBuf', UBufSize)

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', self.numCells * realSize * 3)
	local reduceSwapBufSize = roundup(self.numCells * realSize / self.localSize1d, realSize)
	self:clalloc('reduceSwapBuf', reduceSwapBufSize)
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- TODO CLImageGL ?
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
	
	ffi.cdef(self.eqn:getTypeCode())
	ffi.cdef(self.eqn:getExtraTypeCode())
end


function SolverBase:refreshCodePrefix()
	self:createCodePrefix()		-- depends on eqn, gridSize, displayVars
	self:refreshIntegrator()	-- depends on eqn & gridSize ... & ffi.cdef cons_t
	self:refreshInitStateProgram()
	self:refreshSolverProgram()
--[[ 
TODO here -- refresh init state
but what does that look like for mesh solvers?
where should the initial state be stored?
in an external file?
or should it be calculated?
or should overriding the state be allowed?
--]]
end

-- TODO this has a lot in common with GridSolver
function SolverBase:refreshSolverProgram()
	-- [[ from GritSolver:refreshSolverProgram
	
	self.getULRArg = 'const global '..self.eqn.cons_t..'* UBuf'
	
	self.getULRCode = function(self, args)
		args = args or {}
		local suffix = args.suffix or ''
		return template([[
	const global <?=eqn.cons_t?>* UL<?=suffix?> = UBuf + <?=indexL?>;
	const global <?=eqn.cons_t?>* UR<?=suffix?> = UBuf + <?=indexR?>;
]],		{
			eqn = self.eqn,
			suffix = suffix,
			indexL = args.indexL or 'indexL'..suffix,
			indexR = args.indexR or 'indexR'..suffix,
		})
	end
	--]]

	local code = self:getSolverCode()

	time('compiling solver program', function()
		self.solverProgramObj = self.Program{name='solver', code=code}
		self.solverProgramObj:compile()
	end)

	self:refreshCalcDTKernel()

	-- this is created in the parent class, however it isn't called by the parent class.
	--  instead it has to be called by the individual implementation classes
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name='addSource', domain=self.domainWithoutBorder}
	end

	if self.eqn.useConstrainU then
		self.constrainUKernelObj = self.solverProgramObj:kernel'constrainU'
	end

	for _,op in ipairs(self.ops) do
		op:refreshSolverProgram()
	end

	-- display stuff:

	if self.app.useGLSharing then
		for _,displayVarGroup in ipairs(self.displayVarGroups) do
			for _,var in ipairs(displayVarGroup.vars) do
				--[[
				if var.enabled
				or (var.vecVar and var.vecVar.enabled)
				then
				--]]do
					var.calcDisplayVarToTexKernelObj = self.solverProgramObj:kernel(var.toTexKernelName)
					var.calcDisplayVarToTexKernelObj.obj:setArg(1, self.texCLMem)
				end
			end
		end
	end

	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			--[[
			if var.enabled
			or (var.vecVar and var.vecVar.enabled)
			then
			--]]do
				var.calcDisplayVarToBufferKernelObj = self.solverProgramObj:kernel(var.toBufferKernelName)
				var.calcDisplayVarToBufferKernelObj.obj:setArg(1, self.reduceBuf)
			end
		end
	end
end

-- for solvers who don't rely on calcDT
function SolverBase:refreshCalcDTKernel()
	self.calcDTKernelObj = self.solverProgramObj:kernel'calcDT'
	self.calcDTKernelObj.obj:setArg(1, self.reduceBuf)
end


function SolverBase:getSolverCode()
	local fluxLimiterCode = 'real fluxLimiter(real r) {'
		.. self.app.limiters[self.fluxLimiter].code
		.. '}'

	return table{
		self.codePrefix,
		
		fluxLimiterCode,
		
		'typedef struct { real min, max; } range_t;',
		
		-- TODO move to Roe, or FiniteVolumeSolver as a parent of Roe and HLL?
		self.eqn:getEigenCode() or '',
		self.eqn:getSolverCode() or '',
		self.eqn:getCalcDTCode() or '',
		self.eqn:getFluxFromConsCode() or '',
		self.eqn:getCalcEigenBasisCode() or '',
	
	}:append(self.ops:map(function(op)
		return op:getSolverCode()
	end)):append{
		self:getDisplayCode(),
	}:concat'\n'
end

function SolverBase:getDisplayCode()
	local lines = table()
	
	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			var.id = tostring(var):sub(10)
		end
	end

	if self.app.useGLSharing then
		for _,displayVarGroup in ipairs(self.displayVarGroups) do
			for _,var in ipairs(displayVarGroup.vars) do
				var.toTexKernelName = makestruct.uniqueName'calcDisplayVarToTex'
				--[[
				if var.enabled
				or (var.vecVar and var.vecVar.enabled)
				then
				--]]do
					lines:append{
						template(var.displayCode, {
							solver = self,
							var = var,
							name = var.toTexKernelName,
							
							input = 
							
							-- nvidia needed this, but I don't want to write only -- I want to accumulate and do other operations
							'write_only '	..
							-- if I do accumulate, then I will need to ensure the buffer is initialized to zero ...
								
								(self.dim == 3
									and 'image3d_t'
									or 'image2d_t'
								)..' tex',
							output = template([[
#warning no accum function support yet
write_imagef(
	tex,
	<? if solver.dim == 3 then ?> i <? else ?> i.xy <? end ?>,
	(float4)(value[0], value[1], value[2], 0.));
]], {
		solver = self,
		accumFunc = self.displayVarAccumFunc and 'max' or nil,
	}),
						})
					}
				end
			end
		end
	end

	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			var.toBufferKernelName = makestruct.uniqueName'calcDisplayVarToBuffer'
			--[[
			if var.enabled
			or (var.vecVar and var.vecVar.enabled)
			then
			--]]do
				lines:append{
					template(var.displayCode, {
						solver = self,
						var = var,
						name = var.toBufferKernelName,
						input = 'global real* dest',
						output = var.vectorField and [[
	dest[0+3*dstindex] = value[0];
	dest[1+3*dstindex] = value[1];
	dest[2+3*dstindex] = value[2];
]] or template([[
	dest[dstindex] = <?
if accumFunc then
	?>
	<?=accumFunc?>(value[0], dest[dstindex])
	<?
else
	?> value[0] <?
end
?>
;
]], {
		accumFunc = self.displayVarAccumFunc and 'max' or nil,
	}),
					})
				}
			end
		end
	end
	
	local code = lines:concat'\n'
	
	return code
end


function SolverBase:refreshInitStateProgram()
	self.eqn.initState:refreshInitStateProgram(self)
end

function SolverBase:createCodePrefix()
	local lines = table()
	
	if self.dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end

	-- real3
	lines:insert(template(file['math.h'], {app=self.app}))
	
	lines:append{
		'#ifndef M_PI',
		'#define M_PI '..('%.50f'):format(math.pi),
		'#endif',
	}

	lines:append{
		'#define dim '..self.dim,
		'#define numStates '..self.eqn.numStates,
		'#define numIntStates '..self.eqn.numIntStates,
		'#define numWaves '..self.eqn.numWaves,
	}
	
	-- real3_rotateFrom, real3_rotateTo depend on 'dim'
	lines:insert(template(file['math.cl']))
	
	lines:insert(self.coord:getCode(self))

	-- this can use the coord_raise or coord_lower code
	-- which is associated with the coordinate system,
	lines:append{
		self.eqn:getTypeCode(),
		self.eqn:getExtraTypeCode(),
		self.eqn:getEigenTypeCode() or '',
		self.eqn:getCodePrefix() or '',
	}

	self.codePrefix = lines:concat'\n'
end


function SolverBase:refreshIntegrator()
	self.integrator = integrators[self.integratorIndex](self)
end

SolverBase.t = 0
function SolverBase:resetState()
	self.t = 0
	self.app.cmds:finish()

	self:applyInitCond()
	self:resetOps()
end

-- override this by the mesh solver ... since I don't know what it will be doing
function SolverBase:applyInitCond()
	self.eqn.initState:resetState(self)
	if self.allowAccum then
		self.app.cmds:enqueueFillBuffer{buffer=self.accumBuf, size=ffi.sizeof(self.app.real) * self.numCells * 3}
	end
end

function SolverBase:resetOps()
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

DisplayVar.type = 'real'	-- default

--[[
TODO buf (dest) shouldn't have ghost cells
and dstIndex should be based on the size without ghost cells

why would I bother write to the ghost cells?
the only reason I can think of is for good subtexel lookup when rendering
--]]
DisplayVar.displayCode = [[
kernel void <?=name?>(
	constant <?=solver.solver_t?>* solver,
	<?=input?>,
	const global <?= var.type ?>* buf
<? if require 'solver.meshsolver'.is(solver) then ?>
	,const global cell_t* cells			//[numCells]
	,const global iface_t* ifaces		//[numInterfaces]
<? end ?>
	<?=
	var.extraArgs and #var.extraArgs > 0
		and ',\n\t'..table.concat(var.extraArgs, ',\n\t')
		or ''
?>
) {
	SETBOUNDS(0,0);
<? if not require 'solver.meshsolver'.is(solver) then ?>
	int4 dsti = i;
	int dstindex = index;
	
	real3 x = cell_x(i);

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
<? else	-- mesh ?>
	int dstindex = index;
	real3 x = cells[index].x;
<? end 		-- mesh vs grid ?>

	//TODO rename to value_real
	real value[6] = {0,0,0,0,0,0};	//size of largest struct
	sym3* valuesym3 = (sym3*)value;

	real* value_real = value;
	cplx* value_cplx = (cplx*)value;
	real3* value_real3 = (real3*)value;
	real3* value_real3_hi = (real3*)(value+3);
	cplx3* value_cplx3 = (cplx3*)value;

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
	local buffer = assert(self.solver[self.bufferField], "failed to find buffer "..tostring(self.bufferField))
	kernel:setArg(0, self.solver.solverBuf)
	kernel:setArg(2, buffer)
end

function DisplayVar:setToTexArgs()
	self:setArgs(self.calcDisplayVarToTexKernelObj.obj)
end

function DisplayVar:setToBufferArgs(var)
	self:setArgs(self.calcDisplayVarToBufferKernelObj.obj)
end


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
		cplx = {
			{name = ' re', code = '	*value_cplx = cplx_from_real(value_cplx->re);'},
			{name = ' im', code = '	*value_cplx = cplx_from_real(value_cplx->im);'},
			{name = ' abs', code = '	*value_cplx = cplx_from_real(cplx_abs(*value_cplx));'},
			{name = ' arg', code = '	*value_cplx = cplx_from_real(cplx_arg(*value_cplx));'},
		},
		
		real3 = {
			{name = ' x', code = '	*value_real3 = _real3(value_real3->x,0,0);'},
			{name = ' y', code = '	*value_real3 = _real3(value_real3->y,0,0);'},
			{name = ' z', code = '	*value_real3 = _real3(value_real3->z,0,0);'},
			{name = ' mag', code = '	*value_real3 = _real3(real3_len(*value_real3),0,0);', magn=true},
		},


		cplx3 = {
			{name = ' mag', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx3_len(*value_cplx3)), cplx_zero, cplx_zero);', magn=true},
			
			-- TODO these two are crashing
			{name = ' re', code = '	*value_real3 = cplx3_re(*value_cplx3); *value_real3_hi = real3_zero;', vartype='real3'},
			{name = ' im', code = '	*value_real3 = cplx3_im(*value_cplx3); *value_real3_hi = real3_zero;', vartype='real3'},
			
			-- re and im will include re len, im len, re xyz, im xyz
			-- but will skip the x,y,z cplx abs and arg:
			{name = ' x abs', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx_abs(value_cplx3->x)), cplx_zero, cplx_zero);'},
			{name = ' y abs', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx_abs(value_cplx3->y)), cplx_zero, cplx_zero);'},
			{name = ' z abs', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx_abs(value_cplx3->z)), cplx_zero, cplx_zero);'},
			{name = ' x arg', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx_arg(value_cplx3->x)), cplx_zero, cplx_zero);'},
			{name = ' y arg', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx_arg(value_cplx3->y)), cplx_zero, cplx_zero);'},
			{name = ' z arg', code = '	*value_cplx3 = _cplx3(cplx_from_real(cplx_arg(value_cplx3->z)), cplx_zero, cplx_zero);'},
		},
		
		-- hmm, value_real3 has to be bigger for this to work
		-- but does that mean I have to store 6 components in value_real3?
		-- I suppose it does if I want a sym3-specific visualization
		sym3 = {
			{name = ' x', code = '	*value_real3 = sym3_x(*valuesym3); *value_real3_hi = real3_zero;', vartype='real3'},
			{name = ' y', code = '	*value_real3 = sym3_y(*valuesym3); *value_real3_hi = real3_zero;', vartype='real3'},
			{name = ' z', code = '	*value_real3 = sym3_z(*valuesym3); *value_real3_hi = real3_zero;', vartype='real3'},
	
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
enableVector = false

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
			
			-- TODO how about saving somewhere what should be enabled by default?
			-- TODO pick predefined somewhere?
			if self.predefinedDisplayVars then
				enabled = not not table.find(self.predefinedDisplayVars, args.name)
			else
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

-- used by the display code to dynamically adjust ranges
function SolverBase:calcDisplayVarRange(var)
	if var.lastTime == self.t then
		return var.lastMin, var.lastMax
	end
	var.lastTime = self.t
	
	var:setToBufferArgs()

	local channels = var.vectorField and 3 or 1
	-- this size stuff is very GridSolver-based
	local volume = self.numCells
	local sizevec = self[var.bufferField].sizevec
	if sizevec then
		volume = tonumber(sizevec:volume())
	end
	
	self:calcDisplayVarToBuffer(var)

-- why (when sizevec is explicitly set) does cpu reduce work, but gpu reduce not work?
--[[
if var.name == 'amrError 0' then
local volume = tonumber(self.amrRootSizeInFromSize:volume())
print('self.amrRootSizeInFromSize',self.amrRootSizeInFromSize)
local ptr = ffi.new('real[?]', volume*channels)
self.app.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * volume * channels, ptr=ptr}
print'buffer:'
local min = math.huge
local max = -math.huge
for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
local i = nx + self.amrRootSizeInFromSize.x * ny
io.write('\t', ('%.5f'):format(ptr[i]))
min = math.min(min, ptr[i])
max = math.max(max, ptr[i])
end
print()
end
print('reduce min',min,'max',max,'volume',volume,'name',var.name,'channels',channels)
var.lastMin = min
var.lastMax = max
return min, max
else
--]] do

	-- this is failing with 1D for channels == 3, for vectors (but works for 2D etc)
	local min = self.reduceMin(nil, volume*channels)
	self:calcDisplayVarToBuffer(var)
	local max = self.reduceMax(nil, volume*channels)

--print('reduce min',min,'max',max,'volume',volume,'name',var.name,'channels',channels)
	var.lastMin = min
	var.lastMax = max

	return min, max
end

end

-- used by the output to print out avg, min, max
function SolverBase:calcDisplayVarRangeAndAvg(var)
	local needsUpdate = var.lastTime ~= self.t

	-- this will update lastTime if necessary
	local min, max = self:calcDisplayVarRange(var)
	-- displayVarGroup has already set up the appropriate args
	
	-- duplicated in calcDisplayVarRange
	local size = self.numCells
	local sizevec = self[var.bufferField].sizevec
	if sizevec then
		size = tonumber(sizevec:volume())
	end
	
	local avg
	if needsUpdate then
		self:calcDisplayVarToBuffer(var)
		avg = self.reduceSum(nil, size) / tonumber(size)
	else
		avg = var.lastAvg
	end

	return min, max, avg
end





function SolverBase:initDraw()
end

SolverBase.fpsNumSamples = 30


function SolverBase:calcDT()
	local dt
	-- calc cell wavespeeds -> dts
	if self.useFixedDT then
		dt = self.fixedDT
	else
		-- TODO this without the border, but that means changing reduce *and display*
		self.calcDTKernelObj.obj:setArg(0, self.solverBuf)
		self.calcDTKernelObj.obj:setArg(2, self.UBuf)
		self.calcDTKernelObj()
		dt = self.cfl * self.reduceMin()
		if not math.isfinite(dt) then
			print("got a bad dt at time "..self.t) -- TODO dump all buffers
		end
		self.fixedDT = dt
	end
	return dt
end



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

	if self.checkNaNs then
		if self:checkFinite(self.UBufObj, self.numCells) then return end
	end
end


-- this is abstracted because accumBuf might want to be used ...
function SolverBase:calcDisplayVarToBuffer(var)
	local channels = var.vectorField and 3 or 1
	local app = self.app

	-- duplicated in calcDisplayVarRange
	local volume = self.numCells
	local sizevec = self[var.bufferField].sizevec
	if sizevec then
		volume = tonumber(sizevec:volume())
	end
	
	if self.displayVarAccumFunc	then
		app.cmds:enqueueCopyBuffer{src=self.accumBuf, dst=self.reduceBuf, size=ffi.sizeof(app.real) * volume * channels}
	end
	var.calcDisplayVarToBufferKernelObj()
	if self.displayVarAccumFunc then
		app.cmds:enqueueCopyBuffer{src=self.reduceBuf, dst=self.accumBuf, size=ffi.sizeof(app.real) * volume * channels}
	end
end


function SolverBase:updateGUIParams()
	ig.igText('t: '..self.t)
	
	-- hmm put fps somewhere else, or put ms update here
	ig.igText('fps: '..(self.fps and tostring(self.fps) or ''))
		
	tooltip.checkboxTable('check NaNs', self, 'checkNaNs')

	tooltip.checkboxTable('use fixed dt', self, 'useFixedDT')
	ig.igSameLine()
	
	tooltip.numberTable('fixed dt', self, 'fixedDT')
	tooltip.numberTable('CFL', self, 'cfl')


	if self.allowAccum then
		if tooltip.checkboxTable('accum', self, 'displayVarAccumFunc') then
			self:refreshSolverProgram()	-- I guess getDisplayCode is now in getSolverCode
			self.app.cmds:enqueueFillBuffer{buffer=self.accumBuf, size=ffi.sizeof(self.app.real) * self.numCells * 3}
		end
	end

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

	if tooltip.comboTable('flux limiter', self, 'fluxLimiter', self.app.limiterNames) then
		self:refreshSolverProgram()
	end
end

require 'draw.vectorfield'.applyToSolver(SolverBase)

return SolverBase
