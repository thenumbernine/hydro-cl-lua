-- TODO make this solver/solver.lua, and make the old solver.lua something like structured-grid-solver

local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local file = require 'ext.file'
local math = require 'ext.math'
local CLBuffer = require 'cl.obj.buffer'
local template = require 'template'
local vec3 = require 'vec.vec3'
local tooltip = require 'tooltip'
local makestruct = require'eqn.makestruct'
local roundup = require 'util.roundup'
local time, getTime = table.unpack(require 'util.time')


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

-- whether to use separate linked binaries.  would it save on compile time?
-- this does work, however I don't have caching for linked libraries set up yet
SolverBase.useCLLinkLibraries = false

-- whether to check for NaNs
SolverBase.checkNaNs = cmdline.checknans or false

SolverBase.showFPS = cmdline.showfps or false

-- enable for us to let the user accum/not.  requires an extra buffer allocation
SolverBase.allowAccum = true
SolverBase.displayVarAccumFunc = false

SolverBase.solverVars = table{
	-- [[ right now the mesh initial conditions use these, but otherwise they can be GridSolver-specific
	{name='mins', type='real3'},
	{name='maxs', type='real3'},
	--]]
}

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
	integratorArgs = integrator args
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
	-- this is moved from #define to constant, so AMR leaf nodes can change it.
	-- coordinate space = u,v,w
	-- cartesian space = x,y,z
	-- min and max in coordinate space
	local code = makestruct.makeStruct(self.solver_t, self.solverVars, nil, true)
	if self.lastSolverTypeCode then
		assert(code == self.lastSolverTypeCode)
	end
	self.lastSolverTypeCode = code
	return code
end

function SolverBase:preInit(args)

-- despite the name, this doesn't have anything to do with the grid size ...
-- ... except in the GridSolver class
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	self.maxWorkGroupSize = tonumber(self.app.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')
--print('maxWorkGroupSize', self.maxWorkGroupSize)
	
	local sizeProps = self:getSizePropsForWorkGroupSize(self.maxWorkGroupSize)
	for k,v in pairs(sizeProps) do
--print(k,v)
		self[k] = v
	end
	
	-- hmm, do some eqns create ops need to know the grid size?
	self.eqnName = args.eqn
	self.eqnArgs = args.eqnArgs
	self:createEqn()

	self.name = self.eqn.name..' '..self.name

	self.initStateIndex = table.find(self.eqn.initStateNames, args.initState) or 1
	self.initStateArgs = args.initStateArgs

	self.integratorArgs = args.integratorArgs
	self.integratorIndex = integratorNames:find(args.integrator) or 1
	
	if require 'coord.coord'.is(args.coord) then
		self.coord = args.coord	-- ptr copy expected by AMR
	else
		self.coord = require('coord.'..args.coord)(table({solver=self}, args.coordArgs))
	end

	self.checkNaNs = self.checkNaNs
	self.useFixedDT = not not args.fixedDT
	self.fixedDT = args.fixedDT or self.fixedDT or .001
	self.cfl = args.cfl or .5	--/self.dim
	self.fluxLimiter = self.app.limiterNames:find(args.fluxLimiter) or 1



	-- this influences createCodePrefix (via its call of eqn:getCodePrefix)
	--  and refreshInitStateProgram()
	self.eqn:createInitState()

	-- add eqn vars to solver_t
	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			var:addToSolver(self)
		end
	end

	self:refreshGetULR()

	-- do this before any call to createBuffers or createCodePrefix
	-- make sure it's done after createEqn (for the solver_t struct to be filled out by the eqn)
	-- actually this has to go after self.eqn:createInitState, 
	--  which is called in refreshEqnInitState
	--  which is called in refreshGridSize
	--  which is called in postInit
	-- the createInitState also creates kernels and runs them on the solver buffers
	-- so I probably need a separate call to eqn.initState, earlier, which constructs the object and the guiVars, but runs no kernels
	self.solver_t = self.app:uniqueName'solver_t'
	makestruct.safeFFICDef(self:getSolverTypeCode())
	self.solverPtr = ffi.new(self.solver_t)
end

function SolverBase:refreshGetULR()
	self.getULRBufType = self.eqn.cons_t
	self.getULRBufName = 'UBuf'
	self.getULRArg = self.getULRBufType..'* '..self.getULRBufName
	self.getULRCode = function(self, args)
		args = args or {}
		local suffix = args.suffix or ''
		return template([[
	const global <?=eqn.cons_t?>* UL<?=suffix?> = <?=bufName?> + <?=indexL?>;
	const global <?=eqn.cons_t?>* UR<?=suffix?> = <?=bufName?> + <?=indexR?>;
]],		{
			solver = self,
			eqn = self.eqn,
			suffix = suffix,
			indexL = args.indexL or 'indexL'..suffix,
			indexR = args.indexR or 'indexR'..suffix,
			bufName = args.bufName or self.getULRBufName,	-- for displayVars the variable name is 'buf', so I need to override it either in displayCode or here
		})
	end
end

function SolverBase:postInit()
	self:createDisplayVars()	-- depends on self.eqn
	
	self:refreshGridSize()		-- depends on createDisplayVars
	-- refreshGridSize calls refreshCodePrefix 
	-- ... calls refreshEqnInitState
	-- ... calls refreshSolverProgram 
	-- ... which needs display vars
	-- TODO get rid of this 'refresh' stuff.  
	-- it was designed for callbacks when changing grid resolution, integrator, etc while the simulation was live.
	-- doing this now is unreasonable, considering how solver_t is tightly wound with initState, eqn, and the scheme

	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			var:setToSolver(self)
		end
	end

	self:refreshSolverBuf()
end

function SolverBase:createSolverBuf()
	-- doesn't use clalloc ...
	-- should happen before any other buffer allocs
	self.solverBuf = CLBuffer{
		env = self.app.env,
		name = 'solver',
		type = self.solver_t,
		count = 1,
		readwrite = 'read',
		constant = true,
	}
end

function SolverBase:refreshSolverBuf()
	self.solverBuf:fromCPU(self.solverPtr)
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
	constant <?=solver.solver_t?>* solver,	// TODO just 'n'?
	global <?=eqn.cons_t?>* a,
	const global <?=eqn.cons_t?>* b,
	const global <?=eqn.cons_t?>* c,
	realparam d
) {
	SETBOUNDS_NOGHOST();
<? 
-- hmm, I only need numIntStates integrated
-- but the remaining variables I need initialized at least
-- and in RK, the U's are initialized to zero ...
-- how to get around this?
for i=0,eqn.numStates-1 do
?>	a[index].ptr[<?=i?>] = b[index].ptr[<?=i?>] + c[index].ptr[<?=i?>] * d;
<? 
end
?>}
]], 	{
			solver = self,
			eqn = self.eqn,
		})
	}:concat'\n'

if self.useCLLinkLibraries then 
	time('compiling common program', function()
		self.commonUnlinkedObj = self.Program{name='common', code=commonCode}
		self.commonUnlinkedObj:compile{dontLink=true}
	end)
	time('linking common program', function()
		self.commonProgramObj = self.Program{
			programs = {
				self.mathUnlinkedObj, 
				self.commonUnlinkedObj,
			},
		}
	end)
else
	time('building common program', function()
		self.commonProgramObj = self.Program{name='common', code=commonCode}
		self.commonProgramObj:compile()
	end)
end

	-- used by the integrators
	-- needs the same globalSize and localSize as the typical simulation kernels
	-- TODO exclude states which are not supposed to be integrated
	self.multAddKernelObj = self.commonProgramObj:kernel{name='multAdd', domain=self.domainWithoutBorder}
	self.multAddKernelObj.obj:setArg(0, self.solverBuf)

	-- TODO vectors won't reduce anymore unless the reduceMin is constructed with 3* the # of count
	-- but this breaks reduce for scalars (it includes those extra zeros)
	-- which makes calcDT fail
	self.reduceMin = self.app.env:reduce{
		count = self.numCells,
		op = function(x,y) return 'min('..x..', '..y..')' end,
		initValue = 'INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceMax = self.app.env:reduce{
		count = self.numCells,
		op = function(x,y) return 'max('..x..', '..y..')' end,
		initValue = '-INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceSum = self.app.env:reduce{
		count = self.numCells,
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

--[[
name = name of the buffer.  The buffer is assigned to self[name].
ctype = buffer type
count = the number of the elements in the buffer
sizevec = used to specify the grid size of the buffer, used with overriding the default solver gridsize, was going to be used with AMR
--]]
function SolverBase:clalloc(name, ctype, count, sizevec)
	self.buffers:insert{
		name = name,
		count = count,
		type = ctype,
		sizevec = sizevec,
	}
end

function SolverBase:finalizeCLAllocs()
	for _,info in ipairs(self.buffers) do
		local name = info.name
		local ctype = info.type
		local count = info.count
		local size = count * ffi.sizeof(ctype)
		local mod = size % ffi.sizeof(self.app.real)
		if mod ~= 0 then
			-- WARNING?
			print('resizing buffer to be aligned with sizeof('..self.app.real..')')
			size = size - mod + ffi.sizeof(self.app.real)
			info.size = size
			count = count + math.ceil(ffi.sizeof(self.app.real) / ffi.sizeof(ctype))
			info.count = count
		end
		local bufObj = CLBuffer{
			env = self.app.env,
			name = name,
			type = ctype,
			count = count,
		}
		self[name..'Obj'] = bufObj
		self[name] = bufObj.obj
		
		-- I don't know where else to put this ...
		self[name].sizevec = info.sizevec
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
	
	self:clalloc('UBuf', self.eqn.cons_t, self.numCells)

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', self.app.real, self.numCells * 3)
	self:clalloc('reduceSwapBuf', self.app.real, math.ceil(self.numCells / self.localSize1d))
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- TODO CLImageGL ?
end

-- call this when the solver initializes or changes the codePrefix (or changes initState)
-- it will build the code prefix and refresh everything related to it
-- TODO if you change cons_t then call resetState etc (below the refreshEqnInitState() call a few lines above) in addition to this -- or else your values will get messed up
function SolverBase:refreshEqnInitState()	
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
	self.eqn = require('eqn.'..assert(self.eqnName, "expected solver.eqnName or solver args.eqn"))(table(
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
	self:refreshSolverProgram()	-- depends on createDisplayVars
--[[ 
TODO here -- refresh init state
but what does that look like for mesh solvers?
where should the initial state be stored?
in an external file?
or should it be calculated?
or should overriding the state be allowed?
--]]
end

function SolverBase:buildMathCLUnlinked()
if not SolverBase.useCLLinkLibraries then return end
	if self.mathUnlinkedObj then return end
	-- build math cl binary obj
	time('compiling math program', function()
		self.mathUnlinkedObj = self.Program{
			name = 'math',
			code = template(table{
				file['math.types.h'],
				file['math.h'],
				file['math.cl'],
			}:concat'\n', {app=self.app}),
		}
		self.mathUnlinkedObj:compile{
			dontLink = true,
			buildOptions = '-create-library',
		}
	end)
end

-- depends on buffers
function SolverBase:refreshSolverProgram()
	self:refreshGetULR()	

	self:buildMathCLUnlinked()

	local code
	time('generating solver code', function()
		code = self:getSolverCode()	-- depends on createDisplayVars
	end)

if SolverBase.useCLLinkLibraries then 
	time('compiling solver program', function()
		self.solverUnlinkedObj = self.Program{
			name = 'solver',
			code = code,
		}
		self.solverUnlinkedObj:compile{dontLink=true}
	end)

	time('linking solver program', function()
		self.solverProgramObj = self.Program{
			programs = {self.mathUnlinkedObj, self.solverUnlinkedObj},
		}
	end)
else
	time('building solver program', function()
		self.solverProgramObj = self.Program{
			name = 'solver',
			code = code,
		}
		self.solverProgramObj:compile()
	end)
end

	self:refreshCalcDTKernel()

	-- this is created in the parent class, however it isn't called by the parent class.
	--  instead it has to be called by the individual implementation classes
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name='addSource', domain=self.domainWithoutBorder}
	end

	if self.eqn.useConstrainU then
		self.constrainUKernelObj = self.solverProgramObj:kernel'constrainU'
	end

	if self.usePLM then
		self.calcLRKernelObj = self.solverProgramObj:kernel'calcLR'
	end
	if self.useCTU then
		-- currently implemented in solver/roe.cl
		-- not available for any other flux method
		assert(self.fluxBuf)
		self.updateCTUKernelObj = self.solverProgramObj:kernel'updateCTU'
	end


	for _,op in ipairs(self.ops) do
		op:refreshSolverProgram()
	end

	-- display stuff:
	if self.app.targetSystem ~= 'console' then
	
		if self.app.useGLSharing then
			for _,group in ipairs(self.displayVarGroups) do
				for _,var in ipairs(group.vars) do
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

		for _,group in ipairs(self.displayVarGroups) do
			for _,var in ipairs(group.vars) do
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
end


-- for solvers who don't rely on calcDT
function SolverBase:refreshCalcDTKernel()
	self.calcDTKernelObj = self.solverProgramObj:kernel'calcDT'
	self.calcDTKernelObj.obj:setArg(1, self.reduceBuf)
end


function SolverBase:getSolverCode()
	local fluxLimiterCode = 'real fluxLimiter(real r) {\n'
		.. '\t' ..self.app.limiters[self.fluxLimiter].code..'\n'
		.. '}\n'

	return table{
		self.codePrefix,
		
		fluxLimiterCode,
		
		'typedef struct { real min, max; } range_t;',
		
		-- TODO move to Roe, or FiniteVolumeSolver as a parent of Roe and HLL?
		self.eqn:getEigenCode() or '',
		self.eqn:getSolverCode() or '',
		self.eqn:getCalcDTCode() or '',
		self.eqn:getFluxFromConsCode() or '',
	
	}:append(self.ops:mapi(function(op)
		return op:getCode() or ''
	end)):append{
		self:getDisplayCode(),
	}:concat'\n'
end

local function shouldDeferCode(code)
	do return false end
	return #string.split(string.trim(code), '\n') > 3
end

function SolverBase:getDisplayCode()
	if self.app.targetSystem == 'console' then return '' end
	local texVsBufInfo = {
		Tex = {
			outputArg = function(var)
				-- nvidia needed this, but I don't want to write only -- I want to accumulate and do other operations
				return 'write_only '..
				-- if I do accumulate, then I will need to ensure the buffer is initialized to zero ...
				(self.dim == 3
					and 'image3d_t'
					or 'image2d_t'
				)..' tex'
			end,
			output = function(var)
				return template([[
<? if accumFunc then ?>
	float4 texel = read_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>);
	texel.x = <?=accumFunc?>(texel.x, value[0]);
	if (vectorField) {
		texel.y = <?=accumFunc?>(texel.y, value[1]);
		texel.z = <?=accumFunc?>(texel.z, value[2]);
	}
	write_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>, texel);
<? else ?>
	write_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>, (float4)(value[0], value[1], value[2], 0.));
<? end ?>
]], 			{
					solver = self,
					accumFunc = self.displayVarAccumFunc and 'max' or nil,
				})
			end,
		},
		Buffer = {
			outputArg = function(var) 
				return 'global real* dest' 
			end,
			output = function(var)
				return template([[
if (vectorField) {	
	dest[0+3*dstindex] = <? if accumFunc then ?><?=accumFunc?>(value[0], dest[0+3*dstindex])?><? else ?>value[0]<? end ?>;
	dest[1+3*dstindex] = <? if accumFunc then ?><?=accumFunc?>(value[1], dest[1+3*dstindex])?><? else ?>value[1]<? end ?>;
	dest[2+3*dstindex] = <? if accumFunc then ?><?=accumFunc?>(value[2], dest[2+3*dstindex])?><? else ?>value[2]<? end ?>;
} else {
	dest[dstindex] = <? if accumFunc then ?><?=accumFunc?>(value[0], dest[dstindex])<? else ?>value[0]<? end ?>;
}
]], 			{
					accumFunc = self.displayVarAccumFunc and 'max' or nil,
				})
			end,
		},
	}
	
	local lines = table()

	local alreadyAddedComponentForGroup = {}
	local function addPickComponetForGroup(var)
		local name = self:getPickComponentNameForGroup(var)
		if alreadyAddedComponentForGroup[name] then return end
		alreadyAddedComponentForGroup[name] = true
		lines:insert((template([[
void <?=name?>(
	constant <?=solver.solver_t?>* solver,
	global const <?=var and var.bufferType or 'real'?>* buf,
	int component,
	int* vectorField,
	real* value,
	int4 i
) {
	real* value_real = value;
	sym3* value_sym3 = (sym3*)value;
	cplx* value_cplx = (cplx*)value;
	real3* value_real3 = (real3*)value;
	cplx3* value_cplx3 = (cplx3*)value;
	real3x3* value_real3x3 = (real3x3*)value;

	switch (component) {
<? 
for i,component in ipairs(solver.displayComponentFlatList) do
	if not component.onlyFor 
	or (var and var.group.name == component.onlyFor)
	then
?>	case <?=i?>:	//<?=component.type or 'real'?> <?=component.name?>
		{
<?=component.code?>
			*vectorField = <?= solver:isVarTypeAVectorField(component.type) and '1' or '0' ?>;
			break;
		}
<? 
	end
end
?>	}
}
]], 	{
			name = name,
			solver = self,
			var = var,
		})))
	end
	addPickComponetForGroup()

--[=[
-- ok here's another idea for saving lines of code
-- add some predefined functions for getters for real, real3, sym3, cplx, cplx3, real3x3
-- and, if your variable is just a member of a struct of one of these types, use that.
-- (this means adding in extra params for display: the offset and the struct size)
	
	for _,ctype in ipairs{'real', 'real3', 'sym3', 'cplx', 'cplx3', 'real3x3'} do
		for _, texVsBuf in ipairs{'Tex', 'Buffer'} do	
			local tempvar = {
				code = template([[
	*(<?=ctype?>*)value = *(global const <?=ctype?> *)((global const byte*)buf + index * structSize + structOffset);
]], 			{
					ctype = ctype,
				}),
			}
			lines:insert(template(self.displayCode, {
				solver = self,
				var = tempvar,
				name = 'calcDisplayVarTo'..texVsBuf..'_Simple_'..ctype,
				outputArg = texVsBufInfo[texVsBuf].outputArg(var),
				output = texVsBufInfo[texVsBuf].output(var),
				extraArgs = {
					'int structSize',
					'int structOffset',
				},
			}))
		end
	end
--]=]

	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			if var.originalVar then
				if self.app.useGLSharing then
					var.toTexKernelName = assert(var.originalVar.toTexKernelName)
				end
				var.toBufferKernelName = assert(var.originalVar.toBufferKernelName)
			else
				lines:insert('//'..var.name..'\n')

				addPickComponetForGroup(var)

				if shouldDeferCode(var.code) then
					local clFuncName = self.app:uniqueName'calcDisplayVar'
					lines:insert(template([[
void <?=clFuncName?>(
	constant <?=solver.solver_t?>* solver,
	const global <?= var.bufferType ?>* buf,
	int4 i,
	int index, 
	int4 dsti,
	int dstindex,
	real3 x,
	real* value<?=
	var.extraArgs and #var.extraArgs > 0
		and ',\n\t'..table.concat(var.extraArgs, ',\n\t')
		or ''
?>
) {
	real* value_real = value;
	sym3* value_sym3 = (sym3*)value;
	cplx* value_cplx = (cplx*)value;
	real3* value_real3 = (real3*)value;
	cplx3* value_cplx3 = (cplx3*)value;
	real3x3* value_real3x3 = (real3x3*)value;

<?=var.codePrefix or ''?>
<?=var.code?>
}
]], 				{
						solver = self,
						var = var,
						clFuncName = clFuncName,
					}))
				
					var.code = template([[
	<?=clFuncName?>(solver, buf, i, index, dsti, dstindex, x, value<?=
		var.extraArgNames 
		and #var.extraArgNames > 0 
		and ', '..var.extraArgNames:concat', '
		or ''
	?>);
]], 				{
						var = var,
						clFuncName = clFuncName,
					})
				else
					-- hmm, this will change its success the next time through this test
					-- so this is destructive, only run it once per display var?
					-- or better yet TODO somewhere earlier, maybe before the 'prefix func' stuff,
					-- prepend 'var.codePrefix' onto 'var.code'
					var.code = template([[
	real* value_real = value;
	sym3* value_sym3 = (sym3*)value;
	cplx* value_cplx = (cplx*)value;
	real3* value_real3 = (real3*)value;
	cplx3* value_cplx3 = (cplx3*)value;
	real3x3* value_real3x3 = (real3x3*)value;

<?=var.codePrefix or ''?>
<?=var.code?>
]],					{
						var = var,
					})
				end

				if self.app.useGLSharing then
					var.toTexKernelName = self.app:uniqueName'calcDisplayVarToTex'
					--[[
					if var.enabled
					or (var.vecVar and var.vecVar.enabled)
					then
					--]]do
						lines:insert(
							template(var.displayCode, {
								solver = self,
								var = var,
								name = var.toTexKernelName,
								outputArg = texVsBufInfo.Tex.outputArg(var),
								output = texVsBufInfo.Tex.output(var),
							})
						)
					end
				end

				var.toBufferKernelName = self.app:uniqueName'calcDisplayVarToBuffer'
				--[[
				if var.enabled
				or (var.vecVar and var.vecVar.enabled)
				then
				--]]do
					lines:insert(
						template(var.displayCode, {
							solver = self,
							var = var,
							name = var.toBufferKernelName,
							outputArg = texVsBufInfo.Buffer.outputArg(var),
							output = texVsBufInfo.Buffer.output(var),
						})
					)
				end
			end
		end
	end
	
	local code = lines:concat'\n'
	
	return code
end


function SolverBase:refreshInitStateProgram()
	self.eqn.initState:refreshInitStateProgram(self)
end

-- TODO compile the code of CommonCode into a bin of its own
--  and link against it instead of recopying and recompiling
function SolverBase:createCodePrefixHeader()
	-- header
	
	local lines = table()
	
	-- real3
	lines:insert(template(file['math.types.h'], {app=self.app}))
	lines:insert(template(file['math.h']))

	if self.dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end

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

	-- this can use the coord_raise or coord_lower code
	-- which is associated with the coordinate system,
	lines:append{
		'//SolverBase:createCodePrefix() begin',
		
		'//SolverBase.eqn:getTypeCode()',
		self.eqn:getTypeCode(),
		
		'//SolverBase:getSolverTypeCode()',
		self:getSolverTypeCode(),
		
		'//SolverBase.eqn:getExtraTypeCode()',
		self.eqn:getExtraTypeCode(),
		
		'//SolverBase.eqn:getEigenTypeCode()',
		self.eqn:getEigenTypeCode() or '',
	}
	
	return lines:concat'\n'
end

-- TODO split 'codePrefix' into the common header and common source

function SolverBase:createCodePrefixSource()
	lines = table()
if not SolverBase.useCLLinkLibraries then 
	lines:append{
		'//math.ch',
		template(file['math.cl']),
	}
end	
	lines:append{
		'//solver.coord:getCode()',
		self.coord:getCode(self) or '',
		
		'//solver.eqn:getCodePrefix()',
		self.eqn:getCodePrefix() or '',
	}
	return lines:concat'\n'
end

function SolverBase:createCodePrefix()
	self.codePrefix = table{
		self:createCodePrefixHeader(),
		self:createCodePrefixSource(),
	}:concat'\n'
end


function SolverBase:refreshIntegrator()
	self.integrator = integrators[self.integratorIndex](self, self.integratorArgs)
end

SolverBase.t = 0
function SolverBase:resetState()
	self.t = 0
	self.app.cmds:finish()

	self:applyInitCond()
	self:boundary()
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

function SolverBase:convertToSIUnitsCode(units)
	self.unitCodeCache = self.unitCodeCache or {}
	if self.unitCodeCache[units] then 
		return self.unitCodeCache[units]
	end
	local symmath = require 'symmath'
	local vars = symmath.vars
	local m, s, kg, C, K = vars('m', 's', 'kg', 'C', 'K')
	local expr = assert(load([[
local m, s, kg, C, K = ...
return ]]..units), "failed to compile unit expression "..units)(m, s, kg, C, K)
	expr = expr()
	expr = expr:map(function(ex)
		if symmath.op.pow.is(ex) then
			local power = ex[2].value
			assert(type(power) == 'number')
			if power == math.floor(power) and power > 0 then
				if power == 1 then return ex[1] end
				return symmath.op.mul(table{ex[1]}:rep(power):unpack())
			end
		end
	end)
	
	local Ccode = expr:compile({
		{__solver_meter = m}, 
		{__solver_second = s},
		{__solver_kilogram = kg},
		{__solver_coulomb = C},
		{__solver_kelvin = K},
	}, 'C')
	Ccode = Ccode:gsub('__solver_', 'solver->')
	Ccode = assert((Ccode:match'{ return (.*); }'), "failed to find code block")

	local luaFunc, luaCode = expr:compile{m, s, kg, C, K}

	local code = {
		C = Ccode,
		lua = luaCode,
		func = function()
			return luaFunc(
				self.solverPtr.meter,
				self.solverPtr.second,
				self.solverPtr.kilogram,
				self.solverPtr.coulomb,
				self.solverPtr.kelvin)
		end,
	}
	
	self.unitCodeCache[units] = code
	return code
end


-------------------------------------------------------------------------------
--                                 display vars                              --
-------------------------------------------------------------------------------


local DisplayVarGroup = class()

function DisplayVarGroup:init(args)
	self.name = assert(args.name)
	self.vars = table(args.vars)
end


local DisplayVar = class()

-- this is the default DisplayVar
SolverBase.DisplayVar = DisplayVar

DisplayVar.bufferType = 'real'

--[[
component is going to be one of several pathways to modify the data
I'm doing this in hopes to reduce the number of display kernels
The downside is that I now need to consider all permutations up front, rather than recursively create display kernels

TODO buf (dest) shouldn't have ghost cells
and dstIndex should be based on the size without ghost cells

why would I bother write to the ghost cells?
the only reason I can think of is for good subtexel lookup when rendering
--]]
DisplayVar.displayCode = [[
kernel void <?=name?>(
	constant <?=solver.solver_t?>* solver,
	<?=outputArg?>,
	global const <?=var.bufferType?>* buf,
	int component<?
if require 'solver.meshsolver'.is(solver) then ?>
	,const global cell_t* cells			//[numCells]
	,const global iface_t* ifaces		//[numInterfaces]<?
end ?><?= var.extraArgs and #var.extraArgs > 0
		and ',\n\t'..table.concat(var.extraArgs, ',\n\t')
		or '' ?>
) {
	SETBOUNDS(0,0);
<? if not require 'solver.meshsolver'.is(solver) then 
?>	int4 dsti = i;
	int dstindex = index;
	real3 x = cell_x(i);
<? for j=0,solver.dim-1 do 
?>	i.s<?=j?> = clamp(i.s<?=j?>, numGhost, solver->gridSize.s<?=j?> - numGhost - 1);
<? end
?>	index = INDEXV(i);
<? else	-- mesh 
?>	int dstindex = index;
	real3 x = cells[index].x;
<? end 		-- mesh vs grid 
?>	
	real value[9] = {0,0,0,0,0,0,0,0,0};	<? -- size of largest struct.  TODO how about a union? ?>
	int vectorField = <?=solver:isVarTypeAVectorField(var.type) and '1' or '0'?>;
<?=var.code?>
	<?=solver:getPickComponentNameForGroup(var)?>(solver, buf, component, &vectorField, value, i);
<?=output
?>}
]]

function DisplayVar:init(args)
	self.code = assert(args.code)
	self.name = assert(args.name)
	self.solver = assert(args.solver)
	self.bufferType = args.bufferType	-- or self.bufferType
	self.displayCode = args.displayCode 	-- or self.displayCode
	self.codePrefix = args.codePrefix
	
	-- used in rescaling after calcDisplayVar, and in displaying.  toggled by a flag
	self.units = args.units
	self.showInUnits = true	-- set to false to display raw values, true to display in units.

	self.group = args.group

	self.type = assert(args.type)

	self.component = assert((self.solver.displayComponentFlatList:find(nil, function(component)
		return component.base == self.type
	end)), "failed to find component for base type "..self.type)

	-- display stuff
	self.enabled = not not args.enabled
	self.useLog = args.useLog or false
	self.color = vec3(math.random(), math.random(), math.random()):normalize()
	self.heatMapFixedRange = false	-- args.name ~= 'error'
	self.heatMapValueMin = 0
	self.heatMapValueMax = 1

	-- maybe this should be in args too?
	-- or - instead of buffer - how about all the kernel's args?
	-- but the reason I have to store the field here is that the buffer isn't made yet
	-- TODO? make display vars after buffers so I can store the buffer here?
	self.getBuffer = args.getBuffer
	if not self.getBuffer then
	local bufferField = args.bufferField
		self.getBuffer = function()
			return self.solver[bufferField]
		end
	end
	self.extraArgs = args.extraArgs
end

local intptr = ffi.new('int[1]', 0)
local function int(x)
	intptr[0] = x
	return intptr
end
function DisplayVar:setArgs(kernel)
	local buffer = assert(self.getBuffer(), "failed to find buffer for var "..tostring(self.name))
	kernel:setArg(0, self.solver.solverBuf)
	kernel:setArg(2, buffer)
	kernel:setArg(3, int(self.component))
end

function DisplayVar:setToTexArgs()
	self:setArgs(self.calcDisplayVarToTexKernelObj.obj)
end

function DisplayVar:setToBufferArgs()
	self:setArgs(self.calcDisplayVarToBufferKernelObj.obj)
end



--[[
properties:
	base = base type of the variable.
		if var.type == component.base then this component applies to them
		(the old table was nested: info[base] = {...}, maybe I should still do that?
	name = componet name
	code = component code.  nil means empty string.
	type = result type of the component.  nil means real.  only other option is 'real3' or 'cplx' for vector-field display.
	magn = optional, for type=='real3' or type=='cplx', this is the name of the associated variable that gives the vector magnitude.

	onlyFor = the component is only applied to the specified display group.  hack fo bssnok solver.
--]]
function SolverBase:createDisplayComponents()
	self.displayComponents = table()
	self:addDisplayComponents('real', {
		{name = 'default'},
	})
	self:addDisplayComponents('real3', {
		{name = 'default', type = 'real3', magn='mag'},
		{name = 'mag', code = '*value_real3 = _real3(real3_len(*value_real3),0,0);'},
		{name = 'x', code = '*value_real3 = _real3(value_real3->x,0,0);'},
		{name = 'y', code = '*value_real3 = _real3(value_real3->y,0,0);'},
		{name = 'z', code = '*value_real3 = _real3(value_real3->z,0,0);'},
	})
	self:addDisplayComponents('sym3', {
		{name = 'xx', code = '*value_sym3 = _sym3(value_sym3->xx,0,0,0,0,0);'},
		{name = 'xy', code = '*value_sym3 = _sym3(value_sym3->xy,0,0,0,0,0);'},
		{name = 'xz', code = '*value_sym3 = _sym3(value_sym3->xz,0,0,0,0,0);'},
		{name = 'yy', code = '*value_sym3 = _sym3(value_sym3->yy,0,0,0,0,0);'},
		{name = 'yz', code = '*value_sym3 = _sym3(value_sym3->yz,0,0,0,0,0);'},
		{name = 'zz', code = '*value_sym3 = _sym3(value_sym3->zz,0,0,0,0,0);'},
		{name = 'norm', code = '*value_sym3 = _sym3(sqrt(sym3_dot(*value_sym3, *value_sym3)), 0,0,0,0,0);'},
		{name = 'tr', code = '*value_sym3 = _sym3(sym3_trace(*value_sym3), 0,0,0,0,0);'},
		
		{name = 'x', code = '*value_real3 = _real3(value_sym3->xx, value_sym3->xy, value_sym3->xz); value_real3[1] = real3_zero;', type = 'real3', magn='x mag'},
		{name = 'y', code = '*value_real3 = _real3(value_sym3->xy, value_sym3->yy, value_sym3->yz); value_real3[1] = real3_zero;', type = 'real3', magn='y mag'},
		{name = 'z', code = '*value_real3 = _real3(value_sym3->xz, value_sym3->yz, value_sym3->zz); value_real3[1] = real3_zero;', type = 'real3', magn='z mag'},
		{name = 'x mag', code = '*value_sym3 = _sym3(real3_len(sym3_x(*value_sym3)), 0,0,0,0,0);'},
		{name = 'y mag', code = '*value_sym3 = _sym3(real3_len(sym3_y(*value_sym3)), 0,0,0,0,0);'},
		{name = 'z mag', code = '*value_sym3 = _sym3(real3_len(sym3_z(*value_sym3)), 0,0,0,0,0);'},
	})
	self:addDisplayComponents('cplx', {
		{name = 'default', type='cplx', magn='abs'},
		{name = 're', code = '*value_cplx = cplx_from_real(value_cplx->re);'},
		{name = 'im', code = '*value_cplx = cplx_from_real(value_cplx->im);'},
		{name = 'abs', code = '*value_cplx = cplx_from_real(cplx_abs(*value_cplx));'},
		{name = 'arg', code = '*value_cplx = cplx_from_real(cplx_arg(*value_cplx));'},
	})
	self:addDisplayComponents('cplx3', {
		{name = 'mag', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx3_len(*value_cplx3)), cplx_zero, cplx_zero);'},
		
		{name = 'x', code='*value_cplx3 = _cplx3(value_cplx3->x, cplx_zero, cplx_zero);', type='cplx', magn='x abs'},
		{name = 'x abs', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx_abs(value_cplx3->x)), cplx_zero, cplx_zero);'},
	
		{name = 'y', code='*value_cplx3 = _cplx3(value_cplx3->y, cplx_zero, cplx_zero);', type='cplx', magn='y abs'},
		{name = 'y abs', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx_abs(value_cplx3->y)), cplx_zero, cplx_zero);'},
		
		{name = 'z', code='*value_cplx3 = _cplx3(value_cplx3->z, cplx_zero, cplx_zero);', type='cplx', magn='z abs'},
		{name = 'z abs', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx_abs(value_cplx3->z)), cplx_zero, cplx_zero);'},
		
		{name = 'x arg', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx_arg(value_cplx3->x)), cplx_zero, cplx_zero);'},
		{name = 'y arg', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx_arg(value_cplx3->y)), cplx_zero, cplx_zero);'},
		{name = 'z arg', code = '*value_cplx3 = _cplx3(cplx_from_real(cplx_arg(value_cplx3->z)), cplx_zero, cplx_zero);'},
		{name = 'x re', code = '*value_cplx3 = _cplx3(cplx_from_real(value_cplx3->x.re), cplx_zero, cplx_zero);'},
		{name = 'x im', code = '*value_cplx3 = _cplx3(cplx_from_real(value_cplx3->x.im), cplx_zero, cplx_zero);'},
		{name = 'y re', code = '*value_cplx3 = _cplx3(cplx_from_real(value_cplx3->y.re), cplx_zero, cplx_zero);'},
		{name = 'y im', code = '*value_cplx3 = _cplx3(cplx_from_real(value_cplx3->y.im), cplx_zero, cplx_zero);'},
		{name = 'z re', code = '*value_cplx3 = _cplx3(cplx_from_real(value_cplx3->z.re), cplx_zero, cplx_zero);'},
		{name = 'z im', code = '*value_cplx3 = _cplx3(cplx_from_real(value_cplx3->z.im), cplx_zero, cplx_zero);'},
		
		{name = 're', code = '*value_real3 = cplx3_re(*value_cplx3); *(real3*)(value+3) = real3_zero;', type = 'real3', magn='re mag'},
		{name = 're mag', code = '*value_cplx3 = _cplx3(cplx_from_real(real3_len(cplx3_re(*value_cplx3))), cplx_zero, cplx_zero);'},
		{name = 'im', code = '*value_real3 = cplx3_im(*value_cplx3); *(real3*)(value+3) = real3_zero;', type = 'real3', magn='im mag'},
		{name = 'im mag', code = '*value_cplx3 = _cplx3(cplx_from_real(real3_len(cplx3_im(*value_cplx3))), cplx_zero, cplx_zero);'},
	})
	self:addDisplayComponents('real3x3', {
		{name = 'xx', code = '*value_real3x3 = _real3x3(value_real3x3->x.x, 0,0,0,0,0,0,0,0);'},
		{name = 'xy', code = '*value_real3x3 = _real3x3(value_real3x3->x.y, 0,0,0,0,0,0,0,0);'},
		{name = 'xz', code = '*value_real3x3 = _real3x3(value_real3x3->x.z, 0,0,0,0,0,0,0,0);'},
		
		{name = 'yx', code = '*value_real3x3 = _real3x3(value_real3x3->y.x, 0,0,0,0,0,0,0,0);'},
		{name = 'yy', code = '*value_real3x3 = _real3x3(value_real3x3->y.y, 0,0,0,0,0,0,0,0);'},
		{name = 'yz', code = '*value_real3x3 = _real3x3(value_real3x3->y.z, 0,0,0,0,0,0,0,0);'},
		
		{name = 'zx', code = '*value_real3x3 = _real3x3(value_real3x3->z.x, 0,0,0,0,0,0,0,0);'},
		{name = 'zy', code = '*value_real3x3 = _real3x3(value_real3x3->z.y, 0,0,0,0,0,0,0,0);'},
		{name = 'zz', code = '*value_real3x3 = _real3x3(value_real3x3->z.z, 0,0,0,0,0,0,0,0);'},
		
		{name = 'norm', code = '*value_real3x3 = _real3x3(sqrt(real3x3_dot(*value_real3x3, *value_real3x3)), 0,0,0,0,0,0,0,0);'},
		{name = 'tr', code = '*value_real3x3 = _real3x3(real3x3_trace(*value_real3x3), 0,0,0,0,0,0,0,0);'},
		
		{name = 'x', code = '*value_real3 = value_real3x3->x; value_real3[1] = real3_zero; value_real3[2] = real3_zero;', type = 'real3', magn='x mag'},
		{name = 'y', code = '*value_real3 = value_real3x3->y; value_real3[1] = real3_zero; value_real3[2] = real3_zero;', type = 'real3', magn='y mag'},
		{name = 'z', code = '*value_real3 = value_real3x3->z; value_real3[1] = real3_zero; value_real3[2] = real3_zero;', type = 'real3', magn='z mag'},
		{name = 'x mag', code = '*value_real3 = _real3(real3_len(value_real3x3->x), 0,0); value_real3[1] = real3_zero; value_real3[2] = real3_zero;'},
		{name = 'y mag', code = '*value_real3 = _real3(real3_len(value_real3x3->y), 0,0); value_real3[1] = real3_zero; value_real3[2] = real3_zero;'},
		{name = 'z mag', code = '*value_real3 = _real3(real3_len(value_real3x3->z), 0,0); value_real3[1] = real3_zero; value_real3[2] = real3_zero;'},
		
		{name = 'T x', code = '*value_real3 = _real3(value_real3x3->x.x, value_real3x3->y.x, value_real3x3->z.x); value_real3[1] = real3_zero; value_real3[2] = real3_zero;', type = 'real3', magn='T x mag'},
		{name = 'T y', code = '*value_real3 = _real3(value_real3x3->x.y, value_real3x3->y.y, value_real3x3->z.y); value_real3[1] = real3_zero; value_real3[2] = real3_zero;', type = 'real3', magn='T y mag'},
		{name = 'T z', code = '*value_real3 = _real3(value_real3x3->x.z, value_real3x3->y.z, value_real3x3->z.z); value_real3[1] = real3_zero; value_real3[2] = real3_zero;', type = 'real3', magn='T z mag'},
		{name = 'T x mag', code = '*value_real3 = _real3(real3_len(_real3(value_real3x3->x.x, value_real3x3->y.x, value_real3x3->z.x)), 0,0); value_real3[1] = real3_zero; value_real3[2] = real3_zero;'},
		{name = 'T y mag', code = '*value_real3 = _real3(real3_len(_real3(value_real3x3->x.y, value_real3x3->y.y, value_real3x3->z.y)), 0,0); value_real3[1] = real3_zero; value_real3[2] = real3_zero;'},
		{name = 'T z mag', code = '*value_real3 = _real3(real3_len(_real3(value_real3x3->x.z, value_real3x3->y.z, value_real3x3->z.z)), 0,0); value_real3[1] = real3_zero; value_real3[2] = real3_zero;'},
	})

end
function SolverBase:addDisplayComponents(basetype, components)
	for _,component in ipairs(components) do
		self:addDisplayComponent(basetype, component)
	end
end
function SolverBase:addDisplayComponent(basetype, component)
	if not self.displayComponents[basetype] then self.displayComponents[basetype] = table() end
	-- add defaults
	component.base = basetype
	-- TODO evaluate here or later?
	component.code = template(component.code or '', {
		solver = self,
		eqn = self.eqn,
	})
	component.type = component.type or 'real'
	self.displayComponents[basetype]:insert(component)
end
function SolverBase:finalizeDisplayComponents()
	-- build a 1-based enum of all components
	self.displayComponentFlatList = table()
	self.displayComponentNames = table()
	for basetype,components in pairs(self.displayComponents) do
		for _,component in ipairs(components) do
			self.displayComponentFlatList:insert(component)
			-- TODO one name per gruop and only show that list to the dropdown
			self.displayComponentNames:insert(basetype..' '..component.name)
			component.globalIndex = #self.displayComponentFlatList
		end
	end
	-- now find the magnitude components
	for basetype,components in pairs(self.displayComponents) do
		for _,component in ipairs(components) do
			assert(self:isVarTypeAVectorField(component.type) == (not not component.magn), "expected all vector field variables to have an associated magnitude variable.  missed "..basetype..' '..component.name)
			if component.magn then
				local _, magn = components:find(nil, function(component2) return component2.name == component.magn end)
				assert(magn, "couldn't find magn component "..component.magn)
				component.magn = magn.globalIndex
			end	
		end
	end
end
function SolverBase:getPickComponentNameForGroup(var)
	local name = 'pickComponent'
	if var 
	and var.group 
	and self.displayComponentFlatList:find(nil, function(component)
		return component.onlyFor
	end)
	then 
		name = name..'_'
			-- TODO further sanitization?
			..var.group.name:gsub(' ', '_')
	end
	return name
end


function SolverBase:isVarTypeAVectorField(vartype)
	return vartype == 'real3' or vartype == 'cplx'
	-- sym3 and cplx3 are too complex to merely be vector fields.
	--  maybe I'll add another display for them later.
end





function SolverBase:newDisplayVarGroup(args)
	local displayVarGroup = DisplayVarGroup(args)
	self.displayVarGroups:insert(displayVarGroup)
	return displayVarGroup
end


-- still used by gr-hd-separate to add 'extraArgs'
function SolverBase:getUBufDisplayVarsArgs()
	return {
		codePrefix = self.eqn:getDisplayVarCodePrefix(),
		bufferType = self.eqn.cons_t,
		bufferField = 'UBuf',
	}
end

function SolverBase:addUBufDisplayVars()
	local group = self:newDisplayVarGroup{name='U'}
	
	local args = self:getUBufDisplayVarsArgs()
	args.group = group
	args.vars = self.eqn:getDisplayVars()
	
	self:addDisplayVarGroup(args, self.DisplayVar_U)
end

function SolverBase:addDisplayVarGroup(args, cl)
	cl = cl or self.DisplayVar

	if not args.group then
		args.group = self:newDisplayVarGroup{name = args.name}
	end

	local group = args.group

	local enableScalar = true
	local enableVector = true
enableVector = false

	for i,var in ipairs(args.vars) do

		local units = var.units
		var = table(var)
		var.units = nil
	
		local name = assert(var.name, "expected to find name in "..require'ext.tolua'(var))
		local code = assert(var.code, "expected to find code")

		args = table(args, {
			solver = self,
			name = group.name .. ' ' .. name,
			code = code,
			units = units,
			group = group,
			type = var.type or 'real',
		})

		-- add units suffix
		-- don't add it to args.name in case subvars are created
		if args.units then
			args.name = args.name .. ' ('..args.units:gsub('%*', ' ')..')'
		end

		-- enable the first scalar field
		-- also enable the first vector field on non-1D simulations
		local enabled
		
		-- TODO how about saving somewhere what should be enabled by default?
		-- TODO pick predefined somewhere?
		if self.eqn.predefinedDisplayVars then
			-- moved this to after the duplicate var creation
			--enabled = not not table.find(self.eqn.predefinedDisplayVars, args.name)
		else
			if group.name == 'U'
			or (group.name:sub(1,5) == 'error' and self.dim == 1)
			then
				if not self:isVarTypeAVectorField(args.type) then
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
		
		args.enabled = enabled
		local varForGroup = cl(args)
		group.vars:insert(varForGroup)
	end

	return args.group
end

function SolverBase:addDisplayVars()
	self:addUBufDisplayVars()
	
	-- might contain nonsense :-p
	-- TODO it also might contain vector components
	self:addDisplayVarGroup{
		name = 'reduce', 
		getBuffer = function() return self.reduceBuf end,
		vars = {{name='0', code='*value = buf[index];'}},
	}

-- [[ use for debugging only for the time being
	-- TODO make this flexible for our integrator
	-- if I put this here then integrator isn't created yet
	-- but if I put this after integrator then display variable init has already happened
	-- also TODO - either only use the UBuf's state variables, 
	-- or don't regen all the display var code somehow and just bind derivbuf to the ubuf functions
	do	--if self.integrator.derivBufObj then
		local group = self:newDisplayVarGroup{name='deriv'}
		local args = self:getUBufDisplayVarsArgs()
		args.getBuffer = function() 
			local int = self.integrator
			-- Euler's deriv buffer
			if int.derivBufObj then return int.derivBufObj.obj end
			-- RK4's first deriv buffer
			if int.derivBufObjs and int.derivBufObjs[1] then return int.derivBufObjs[1].obj end
			-- BE's deriv buffer
			if int.krylov_dUdtObj then return int.krylov_dUdtObj.obj end
			print"HERE"
		end
		args.group = group
		args.vars = self.eqn:getDisplayVarsForStructVars(self.eqn.consVars)
		-- why in addUBufDisplayVars() do I make a new group and assign args.group to it?
		self:addDisplayVarGroup(args, self.DisplayVar_U)
	end
--]]
end

function SolverBase:finalizeDisplayVars()
	self.displayVars = table()
	
	for _,group in ipairs(self.displayVarGroups) do
		local newGroupVars = table()
		for _,var in ipairs(group.vars) do
			--[[ don't insert the original vars
			-- instead use the first duplicated var in the group as the original
			-- this way they can use the default component as their component
			-- (and I don't have to build a global default component and then make exceptions in all the component lookup code in the draw/*)
			newGroupVars:insert(var)
			self.displayVars:insert(var)
			--]]
			-- [[ otherwise...
			local originalVarForGroup
			--]]
			-- while we're here, add duplicates for all components
			for _,component in ipairs(self.displayComponents[var.type]) do
				-- shallow copy the var, so it contains matching kernel references
				-- and then change only the component index
				-- TODO ... but the kernel names haven't been assigned yet ...
				-- so I have to do this after they are assigned, which is after refreshSolverProgram
				-- or I can flag the copied vars and have refreshSolverProgram not make duplicate names for them
				local dupvar = setmetatable(table(var), getmetatable(var))
				dupvar.component = component.globalIndex
				if not originalVarForGroup then
					originalVarForGroup = dupvar
				else
					dupvar.originalVar = originalVarForGroup
				end
				if component.name ~= 'default' then
					dupvar.name = dupvar.name..' '..component.name
				end
				self.displayVars:insert(dupvar)
				newGroupVars:insert(dupvar)
			end
		end
		group.vars = newGroupVars
	end

	-- make lookup by name
	self.displayVarForName = self.displayVars:map(function(var)
		return var, var.name
	end)

	-- only now that we're here, enable predefined vars
	if self.eqn.predefinedDisplayVars then
		for _,name in ipairs(self.eqn.predefinedDisplayVars) do
			local var = self.displayVarForName[name]
			if not var then
				print('predefined var '..name..' not found')
			else
				var.enabled = true
			end
		end
	end
end

-- depends on self.eqn
function SolverBase:createDisplayVars()
	-- create component map for each display var
	-- 'component' should be 'varproperty' or something else
	self:createDisplayComponents()
	self:finalizeDisplayComponents()

	self.displayVarGroups = table()
	self:addDisplayVars()		
	self:finalizeDisplayVars()
end

-- used by the display code to dynamically adjust ranges
-- this returns raw values, not scaled by units
function SolverBase:calcDisplayVarRange(var, componentIndex)
	componentIndex = componentIndex or var.component
	if var.lastTime == self.t then
		return var.lastMin, var.lastMax
	end
	var.lastTime = self.t
	
	var:setToBufferArgs(componentIndex)

	local channels = 1
	-- this size stuff is very GridSolver-based
	local volume = self.numCells
	local sizevec = var.getBuffer().sizevec
	if sizevec then
		volume = tonumber(sizevec:volume())
	end
	
	self:calcDisplayVarToBuffer(var, componentIndex)

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

	-- this fails with 1D for channels == 3, for vectors (but works for 2D etc)
	local min = self.reduceMin(nil, volume*channels)
	self:calcDisplayVarToBuffer(var, componentIndex)
	local max = self.reduceMax(nil, volume*channels)

--print('reduce min',min,'max',max,'volume',volume,'name',var.name,'channels',channels)
	var.lastMin = min
	var.lastMax = max

	return min, max
end

end

-- used by the output to print out avg, min, max
function SolverBase:calcDisplayVarRangeAndAvg(var, componentIndex)
	componentIndex = componentIndex or var.component
	if var.lastTime == self.t then
		return var.lastMin, var.lastMax, var.lastAvg
	end
	--don't do this -- instead let calcDisplayVarRange update and do this:
	-- var.lastTime = self.t

	-- this will update lastTime if necessary
	local min, max = self:calcDisplayVarRange(var)
	-- displayVarGroup has already set up the appropriate args
	
	-- duplicated in calcDisplayVarRange
	local size = self.numCells
	local sizevec = var.getBuffer().sizevec
	if sizevec then
		size = tonumber(sizevec:volume())
	end
	
	self:calcDisplayVarToBuffer(var, componentIndex)
	local avg = self.reduceSum(nil, size) / tonumber(size)
	var.lastAvg = avg
	
	return min, max, avg
end


-------------------------------------------------------------------------------
--                              gui                                          --
-------------------------------------------------------------------------------



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
		if self.showFPS and math.floor(thisTime) ~= math.floor(self.lastFrameTime) then 
			print('fps='..self.fps) 
		end
	end
	self.lastFrameTime = thisTime

if self.checkNaNs then assert(self:checkFinite(self.UBufObj, self.numCells)) end
end

-- check for nans
-- expects buf to be of type cons_t, made up of numStates real variables
function SolverBase:checkFinite(buf)
	local ptr = buf:toCPU()
	local ptr0size = tonumber(ffi.sizeof(buf.type))
	local realSize = tonumber(ffi.sizeof'real')
	local ptrsPerReal = ptr0size / realSize
	assert(ptrsPerReal == math.floor(ptrsPerReal))
	ptr = ffi.cast('real*', ptr)
	local size = buf.count * ptrsPerReal
	local found
	for i=0,size-1 do
		local x = tonumber(ptr[i])
		if not math.isfinite(x) then
			found = found or table()
			found:insert{i, x}
		end
	end
	if not found then return true end
--	self:printBuf(nil, ptr)
	return false, 'found non-finite offsets and numbers: '..require'ext.tolua'(found)..' at t='..self.t
end

function SolverBase:printBuf(buf, ptr)
	ptr = ptr or buf:toCPU()
	local ptr0size = tonumber(ffi.sizeof(buf.type))
	local realSize = tonumber(ffi.sizeof'real')
	local ptrsPerReal = ptr0size / realSize
	assert(ptrsPerReal == math.floor(ptrsPerReal))
	ptr = ffi.cast('real*', ptr)
	local size = buf.size * ptrsPerReal
	local realsPerCell = math.floor(size / self.numCells)
	local max = #tostring(size-1)
	for i=0,self.numCells-1 do
		io.write((' '):rep(max-#tostring(i)), i,':')
		for j=0,realsPerCell-1 do
			print('\t'
				..(j==0 and '[' or '')
				..('%.50f'):format(ptr[j + realsPerCell * i])
				..(j==self.eqn.numStates-1 and ']' or ',')
			)
		end 
	end
end



-- this is abstracted because accumBuf might want to be used ...
function SolverBase:calcDisplayVarToBuffer(var, componentIndex)
	componentIndex = componentIndex or var.component
	local component = self.displayComponentFlatList[componentIndex]
	local vectorField = self:isVarTypeAVectorField(component.type)
	local channels = vectorField and 3 or 1
	local app = self.app

	-- duplicated in calcDisplayVarRange
	local volume = self.numCells
	local sizevec = var.getBuffer().sizevec
	if sizevec then
		volume = tonumber(sizevec:volume())
	end
	
	if self.displayVarAccumFunc	then
		app.cmds:enqueueCopyBuffer{src=self.accumBuf, dst=self.reduceBuf, size=ffi.sizeof(app.real) * volume * channels}
	end
	var:setToBufferArgs()
	var.calcDisplayVarToBufferKernelObj.obj:setArg(1, self.reduceBuf)
	var.calcDisplayVarToBufferKernelObj.obj:setArg(3, int(componentIndex))
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
