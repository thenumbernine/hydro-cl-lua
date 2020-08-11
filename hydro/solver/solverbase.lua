--[[
TODO make this hydro/solver/solver.lua, and make the old solver.lua something like structured-grid-solver

TODO break up the construction of solver and all its members into the following:

0) init mesh size information / vars (numCells, maxWorkGroupSize, etc)
1) init objects
	- create flux
	- create integrator
	- create eqn
		- create initCond
2) init ffi.cdef's
3) init solver code
4) build programs
5) init buffers
6) init kernels


here's the current init structure for fvsolver:
SolverBase:init
	FVSolver:initMeshVars
		GridSolver:initMeshVars
			SolverBase:initMeshVars
				self.app = ...
				self.dim = ...
			
				-- hmm, meshsolver needs this before the numCells etc vars
				-- but coord.sphere-log-radial needs this after .maxs is finalized, after getInitCondCode() is called
				self.coord = ...		<- this object creation is wedged between the other mesh vars because meshsolver needs it early
				
				self.device = ...
				self.cmds = ...
				self.color = ...
				self.ops = ...
				self.solverStruct = ...
				self.solverStruct.vars:append(...)
			self.solverStruct.vars:append(...)	
			self.mins = ...
			self.maxs = ...
			self.initCondMins = ...
			self.initCondMaxs = ...
			self.gridSize = ...
	
	SolverBase:initCLDomainVars
		self.maxWorkGroupSize = ...
		SolverBase:getSizePropsForWorkGroupSize
		self.localSize1d = ...
		self.localSize2d = ...
		self.localSize = ...
		self.globalSize = ...
		self.sizeWithoutBorder = ...
		self.volumeWithoutBorder = ...
		self.numCells = ...
		self.stepSize = ...
		self.offset = ...
	
	FVSolver:initObjs
		GridSolver:initObjs
			SolverBase:initObjs
				SolverBase:createEqn
					self.eqn = ...
				self.fluxLimiter = ...
				self.eqn:createInitState
					self.eqn.initCond:createInitStruct
						self.eqn.initCond.initStruct = ...
					self.eqn.initCond:finalizeInitStruct
						self.eqn.initCond.initCond_t = ...
				
				self.solverStruct:makeType				-- now that initCond does not modify solver_t, this can be moved back -- but eqn still modifies it, so not back by far
				self.solver_t = ...
			
-- Here's where the solverPtr is created.
-- It's based on solver_t's creation so maybe it should go after the initCDef?
-- Creating the initCond_t / initCondPtr object is very parallel to this -- and both use guiVars -- so that should go around here somewhere.

				self.solverPtr = ...
				self.initCond_t = ...
				self.initCondPtr = ...
			GridSolver:createBoundaryOptions
				self.boundaryOptions = ...
				self.eqn:createBoundaryOptions
			self.boundaryMethods = ...
			self.usePLM = ...
			self.slopeLimiter = ...
			self.useCTU = ...
		FVSolver:createFlux
			self.flux = ...

	SolverBase:initCodeModules
		self.modules = ...
		self.coord:initCodeModules
		self.reqmdules = ...

--------- here is where the ffi.cdef is called --------- 

	SolverBase:initCDefs
		SolverBase:getTypeCode()
			self.eqn:getTypeCode
			self.eqn:getEigenTypeCode
			self.eqn:getExtraTypeCode
			self.coord:getCellTypeCode
			self.eqn.initCond.typecode
	
	SolverBase:postInit
		SolverBase:createDisplayVars
			SolverBase:createDisplayComponents
			SolverBase:finalizeDisplayComponents
			FiniteVolume:addDisplayVars
				SolverBase:addDisplayVars
			SolverBase:finalizeDisplayVars
		SolverBase:refreshGridSize
			GridSolver:createSolverBuf
				SolverBase:createSolverBuf
					self.solverBuf = ...

--------- here's where buffers are created ---------
			
			FiniteVolume:createBuffers
				GridSolver:createBuffers
					SolverBase:createBuffers
			
			SolverBase:finalizeCLAllocs
			SolverBase:refreshEqnInitState
				self.eqn.guiVars[k] = ...

--------- here is where the code header is created --------- 
				
				GridSolver:refreshCodePrefix
					SolverBase:refreshCodePrefix
						GridSolver:createCodePrefix
							SolverBase:createCodePrefix
								SolverBase:createCodePrefixHeader
									self.eqn:getTypeCode
									self.initCond.typecode
									self.coord:getTypeCode
									self.solverStruct.typecode
									self.eqn:getTypeCode
									self.eqn:getEigenTypeCode
								SolverBase:createCodePrefixSource
									self.eqn:getCodePrefix
						SolverBase:refreshIntegrator
							self.integrator = ...
						SolverBase:refreshInitStateProgram
							self.eqn.initCond:refreshInitStateProgram
								
--------- here's the last place where the grid mins/maxs can be changed --------- 

								self.eqn:getInitCondCode
						

						FiniteVolumeSolver:refreshSolverProgram
							SolverBase:refreshSolverProgram
								GridSolver:refreshGetULR
									self.getURLBufType = ...
									self.getULRBufName = ...
									self.getURLArg = ...
									self.getURLCode = ...
								SolverBase:buildMathCLUnlinked
									self.mathUnlinkedObj = ...
								
--------- here's where the program code is collected ---------

								FiniteVolumeSolver:getSolverCode
									GridSolver:getSolverCode
										SolverBase:getSolverCode
											self.eqn:getSolverCode
											self.eqn:getCalcDTCode
											self.eqn:getFluxFromConsCode
											self.ops[i]:getCode
											SolverBase:getDisplayCode
												self.displayVarGroups[i].vars[j].toTexKernelName = ...
									self.flux:getSolverCode

--------- here's where programs are created ---------

								self.solverProgramObj = ...
								
--------- here's where kernels are created ---------

								SolverBase:refreshCalcDTKernel
									self.calcDTKernelObj = ...
								self.addSourceKernelObj = ...
								self.constrainUKernelObj = ...
								self.calcLRKernelObj = ...
								self.updateCTUKernelObj = ...
								self.ops[i]:refreshSolverProgram
								self.displayVarGroups[i].vars[j].calcDisplayVarToTexKernelObj = ...
								self.displayVarGroups[i].vars[j].calcDisplayVarToBufferKernelObj = ...
							self.calcFluxKernelObj = ...
							self.calcDerivFromFluxKernelObj = ...
					
					GridSolver:refreshBoundaryProgram
						self.boundaryProgramObj = ...
						self.boundaryKernelObj = ...
				
				self:checkStructSizes
				self.solverPtr.mins|maxs[k] = ...
				self:refreshSolverBufMinsMaxs
				SolverBase:copyGuiVarsToBufs
					self.solverPtr[var.name] = ...
					self:refreshSolverBuf()
					self.initCondPtr[var.name] = ...
					self:refreshInitCondBuf()
			SolverBase:refreshCommonProgram
			GridSolver:resetState
				SolverBase:resetState
					GridSolver:applyInitCond
						self.coord:fillGridCellBuf	<- this depends on solver.gridSize, solver.dim, solver.numGhost, solver.mins, solver.maxs
						SolverBase:applyInitCond
							self.eqn.initCond:resetState
								self.eqn.initCond.applyInitCondKernelObj()
					self:boundary
					self:resetOps
					self:constrainU
		SolverBase:copyGuiVarsToBufs		<- is this needed?  does anything touch the guivars between here and the copyGuiVarsToBufs() call inside refreshEqnInitState() ?
			... same calls as above

so if I wanted to cache all code against its permutation - and not rebuild anything unnecessarily ...

cache/{app.real}/{solver.name and other options}/{eqn.name and other options}/

--]]

local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local io = require 'ext.io'
local file = require 'ext.file'
local math = require 'ext.math'
local gl = require 'gl'
local glreport = require 'gl.report'
local CLBuffer = require 'cl.obj.buffer'
local template = require 'template'
local vec3d = require 'vec-ffi.vec3d'
local tooltip = require 'hydro.tooltip'
local roundup = require 'hydro.util.roundup'
local time, getTime = table.unpack(require 'hydro.util.time')
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6
local from6to3x3 = common.from6to3x3
local sym = common.sym

local half = require 'hydro.half'
local toreal, fromreal = half.toreal, half.fromreal


-- whether to cache the opencl binaries
local useCache = cmdline.useCache
if useCache == nil then useCache = true end


local function addTab(s)
	s = tostring(s)
	if s:sub(1,1) ~= '\t' then s = '\t' .. s end
	if s:sub(-1) ~= '\n' then s = s .. '\n' end
	return s
end


local integrators = require 'hydro.int.all'
local integratorNames = integrators:map(function(integrator) return integrator.name end)


local SolverBase = class()

SolverBase.name = 'Solver'

-- override to specify which hydro/eqn/*.lua to use as the equation
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

--[[
args:
	app
	dim
	eqn = name of hydro/eqn/<name>.lua to use
	eqnArgs (optional) = args of eqn ctor
	coord = coordinate system name to use, associated with coord/<name>.lua
	initCond = init state name
	initCondArgs (optional) = args of init state ctor
	integrator = name of integrator in int/all.lua to use
	integratorArgs = integrator args
--]]
function SolverBase:init(args)
	time('SolverBase:init()', function()
		self:initMeshVars(args)
		self:initCLDomainVars(args)
		self:initObjs(args)
		self:initCodeModules()
		self:initCDefs()
		self:postInit()
	end)
end

function SolverBase:initMeshVars(args)
	assert(args, "expected named parameter table")
	self.app = assert(args.app, "expected app")
	self.dim = assert(args.dim, "expected dim")

	-- MeshSolver needs coord created early
	--  so that it can allocate cell_t's which are defined in coord
	--  (come to think of it, GridSolver allocates cellBuf also)
	-- so cell_t has to be finished before createBuffers.
	self:createCoord(args)


	self.device = args.device or self.app.env.devices[1]
	self.cmds = args.cmds or self.app.env.cmds[1]
	
	self.color = vec3d(math.random(), math.random(), math.random()):normalize()

	-- operators for this solver
	self.ops = table()

	self.solverStruct = Struct{
		solver = self,
		name = 'solver_t',
		dontUnion = true,
	}
	self.solverStruct.vars:append{
	-- [[ right now the mesh initial conditions use these, but otherwise they can be GridSolver-specific
		{name='mins', type='real3'},
		{name='maxs', type='real3'},
	--]]
	-- [[ the mins/maxs, or the super-solver's mins/maxs.  only needed because of the composite solvers. 
		{name='initCondMins', type='real3'},
		{name='initCondMaxs', type='real3'},
	--]]
	}


	-- in GridSolver this was 'initMeshVars' which comes first
	-- in MeshSolver this was 'preInit' which comes later
	
	-- my kernel objs are going to need workgroup info based on domain.size-2*noGhost as well as domain.size ... 
	-- ...and rather than require an extra argument, I think I'll just take advantage of a closure
	local solver = self
	local Program = class(require 'cl.obj.program')
	if ffi.os == 'Windows' then
		os.execute'mkdir cache-cl 2> nul'
	else
		os.execute'mkdir cache-cl 2> /dev/null'
	end
	function Program:init(args)
		self.name = args.name
		args.env = solver.app.env
		args.domain = solver.domain
		
		-- caching binaries, which doesn't write unless the program successfully compiles
		if not cmdline.usecachedcode then 
			if args.name then
				local path = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
				if useCache then args.cacheFile = path end
				file[path..'-attempt'] = args.code	-- write before attempt to compile.  I don't want to write .cl before so that .cl and .bin are always related.
			end
			Program.super.init(self, args)
		
		-- Write generated code the first time.  Subsequent times use the pre-existing code.  Useful for debugging things in the generated OpenCL.
		else
			local path = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
			local clfn = path..'.cl'
			if io.fileexists(clfn) then
				local cachedCode = file[clfn]
				assert(cachedCode:sub(1,#args.env.code) == args.env.code, "seems you have changed the cl env code")
				args.code = cachedCode:sub(#args.env.code+1)
				Program.super.init(self, args)
			else
				Program.super.init(self, args)	-- do this so getCode() works
				file[clfn] = self:getCode()
			end
		end
	end
	function Program:compile(args)
		args = args or {}
		args.buildOptions = '-w'	-- show warnings
		local results = Program.super.compile(self, args)
		if self.obj then	-- did compile
			print((self.name and self.name..' ' or '')..'log:')
			-- TODO log per device ...			
			print(string.trim(self.obj:getLog(solver.device)))
		end
		return results
	end
	self.Program = Program

	if self.app.targetSystem ~= 'console' then
		local GLProgram = class(require 'gl.program')
		function GLProgram:init(...)
			local args = ...
			
			-- Write generated code
			local fn = 'cache-cl/'..solver.app:uniqueName(args.name or 'shader')
			file[fn..'.vert'] = args.vertexCode
			file[fn..'.frag'] = args.fragmentCode

			GLProgram.super.init(self, ...)
		end
		self.GLProgram = GLProgram 
	end
end

function SolverBase:createCoord(args)
	-- not sure where to put this, but it must precede MeshSolver:initMeshVars
	-- also I'm pretty sure CoordinateSystem:init() doesn't require anything from SolverBase, only the later member functions
	if require 'hydro.coord.coord'.is(args.coord) then
		self.coord = args.coord	-- ptr copy expected by AMR
		self.coord.solver = self
	else
		self.coord = require('hydro.coord.'..(args.coord or 'cartesian'))(table({solver=self}, args.coordArgs))
	end
end

function SolverBase:initCLDomainVars(args)

-- despite the name, this doesn't have anything to do with the grid size ...
-- ... except in the GridSolver class
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	self.maxWorkGroupSize = tonumber(self.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')
	if self.app.verbose then print('maxWorkGroupSize', self.maxWorkGroupSize) end
	
	local sizeProps = self:getSizePropsForWorkGroupSize(self.maxWorkGroupSize)
	for k,v in pairs(sizeProps) do
		if self.app.verbose then print(k,v) end
		self[k] = v
	end
end

function SolverBase:initObjs(args)

	-- hmm, do some eqns create ops need to know the grid size?
	self.eqnName = args.eqn
	self.eqnArgs = args.eqnArgs
	self:createEqn()

	self.name = self.eqn.name..' '..self.name

	self.initCondIndex = table.find(self.eqn.initCondNames, args.initCond) or 1
	self.initCondArgs = args.initCondArgs

	self.integratorArgs = args.integratorArgs
	self.integratorIndex = integratorNames:find(args.integrator) or 1

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
			self.solverStruct.vars:insert{name=var.name, type=var.ctype}
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
	-- so I probably need a separate call to eqn.initCond, earlier, which constructs the object and the guiVars, but runs no kernels
	self.solverStruct:makeType()
	self.solver_t = self.solverStruct.typename
	self.solverPtr = ffi.new(self.solver_t)
	-- TODO I see room for a separation between solver_t and eqn_t
	-- but right now both structs are tied directly to solver.cl, so there's no need to separate them at the moment.

	-- not sure where this should go, but probably somewhere parallel to solverPtr
	-- initStruct:makeType() is already called in eqn:createInitState
	-- TODO if initCond is supposed to be modular then this would have to be created after initCond is changed
	self.initCond_t = self.eqn.initCond.initStruct.typename
	self.initCondPtr = ffi.new(self.initCond_t)
end

-- collect *all* modules from all sub-objects
-- do this after objs are created and before codegen
-- first the global ones from math and app, then coord, initCond, boundary, eqn, solver ...
-- then with all of them, specify which ones to target (for .h and .cl code) and they will trim the others
function SolverBase:initCodeModules()
	self.modules = require 'hydro.code.moduleset'(self.app.modules)
	self.coord:initCodeModules()
	-- what to compile?
	self.reqmodules = table{
		'math',
		'coord',
--		'conn',
		'metric',
	}
	self.eqn:initCodeModules()
end

-- TODO if you want to define ffi ctype metatable
-- then put them all in one spot here
-- if they are using the Struct:makeType function then don't bother use SolverBase:getTypeCode()
function SolverBase:initCDefs()
	require 'hydro.code.safecdef'(self:getTypeCode())
end

-- this is code that goes in the codePrefix header (for CL use), as well as ffi.cdef (for ffi C use)
function SolverBase:getTypeCode()
	local lines = table()

	lines:insert''
	lines:insert'//self.eqn:getTypeCode()'
	lines:insert(self.eqn:getTypeCode() or nil)

	lines:insert''
	lines:insert'//self.coord:getCellTypeCode'
	lines:insert(self.coord:getCellTypeCode() or nil)
		
	lines:insert''
	lines:insert'//self.solverStruct.typecode'
	lines:insert(assert(self.solverStruct.typecode))

	lines:insert''
	lines:insert'//self.eqn:getExtraTypeCode'
	lines:insert(self.eqn:getExtraTypeCode() or nil)

	if self.eqn.getEigenTypeCode then
		lines:insert''
		lines:insert'//self.eqn:getEigenTypeCode'
		lines:insert(self.eqn:getEigenTypeCode() or nil)
	end

	return lines:concat'\n'
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
	time('SolverBase:postInit()', function()
		self:createDisplayVars()	-- depends on self.eqn
		
		self:refreshGridSize()		-- depends on createDisplayVars
		-- refreshGridSize calls refreshCodePrefix 
		-- ... calls refreshEqnInitState
		-- ... calls refreshSolverProgram 
		-- ... which needs display vars
		-- TODO get rid of this 'refresh' stuff.  
		-- it was designed for callbacks when changing grid resolution, integrator, etc while the simulation was live.
		-- doing this now is unreasonable, considering how solver_t is tightly wound with initCond, eqn, and the scheme

		self:copyGuiVarsToBufs()
	end)
end

-- buffers that don't use clalloc ...
-- should happen before any other buffer allocs
function SolverBase:createSolverBuf()
	self.solverBuf = CLBuffer{
		env = self.app.env,
		name = 'solver',
		type = self.solver_t,
		count = 1,
		readwrite = 'read',
		constant = true,
	}
end
function SolverBase:createInitCondBuf()
	self.initCondBuf = CLBuffer{
		env = self.app.env,
		name = 'initCond',
		type = self.initCond_t,
		count = 1,
		readwrite = 'read',
		constant = true,
	}
end

function SolverBase:refreshSolverBuf()
	self.solverBuf:fromCPU(self.solverPtr)
end

function SolverBase:refreshInitCondBuf()
	self.initCondBuf:fromCPU(self.initCondPtr)
end

function SolverBase:copyGuiVarsToBufs()
	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			if var.ctype == 'real' then
				self.solverPtr[var.name] = toreal(var.value)
			else
				self.solverPtr[var.name] = var.value
			end
		end
	end
	self:refreshSolverBuf()

	for _,var in ipairs(self.eqn.initCond.guiVars) do
		if not var.compileTime then
			if var.ctype == 'real' then
				self.initCondPtr[var.name] = toreal(var.value)
			else
				self.initCondPtr[var.name] = var.value
			end
		end
	end
	self:refreshInitCondBuf()
end

function SolverBase:refreshGridSize()
	time('SolverBase:refreshGridSize()', function()
		self:createSolverBuf()
		self:createInitCondBuf()
		
		-- depends on eqn & gridSize
		self.buffers = table()
		self:createBuffers()
		self:finalizeCLAllocs()

		-- create the code prefix, reflect changes
		self:refreshEqnInitState()
		
		-- initialize things dependent on cons_t alone
		self:refreshCommonProgram()
		self:resetState()
	end)
end

function SolverBase:refreshCommonProgram()
	-- code that depend on real and nothing else
	-- TODO move to app, along with reduceBuf

	local commonCode = table():append{
		self.codePrefix,
	}:append{
		template([[
kernel void multAddInto(
	constant <?=solver.solver_t?>* solver,	// TODO just 'n'?
	global <?=eqn.cons_t?>* a,
	const global <?=eqn.cons_t?>* b,
	realparam c
) {
	SETBOUNDS_NOGHOST();
<? 
for i=0,eqn.numIntStates-1 do
?>	a[index].ptr[<?=i?>] += b[index].ptr[<?=i?>] * c;
<? 
end
?>}

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
-- Another thought, if I'm scaling *everything* in the struct, then just use reals and scale up the kernel size by numReals
for i=0,eqn.numIntStates-1 do
?>	a[index].ptr[<?=i?>] = b[index].ptr[<?=i?>] + c[index].ptr[<?=i?>] * d;
<? 
end
?>}

kernel void square(
	constant <?=solver.solver_t?>* solver,	// TODO just 'n'?
	global <?=eqn.cons_t?>* a
) {
	SETBOUNDS_NOGHOST();
<?	-- numStates or numIntStates?
for i=0,eqn.numIntStates-1 do
?>	a[index].ptr[<?=i?>] *= a[index].ptr[<?=i?>];
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
	
	self.multAddIntoKernelObj = self.commonProgramObj:kernel{name='multAddInto', domain=self.domainWithoutBorder}
	self.multAddIntoKernelObj.obj:setArg(0, self.solverBuf)
	
	self.squareKernelObj = self.commonProgramObj:kernel{name='square', domain=self.domainWithoutBorder}
	self.squareKernelObj.obj:setArg(0, self.solverBuf)

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
print(require 'ext.tolua'(info))
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
	local app = self.app
	local realSize = ffi.sizeof(app.real)
	
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
	
	-- should I put these all in one AoS?
	-- or finally make use of constant args ...
	
	-- this much for the first level U data
	self:clalloc('UBuf', self.eqn.cons_t, self.numCells)

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', app.real, self.numCells * 3)
	self:clalloc('reduceSwapBuf', app.real, math.ceil(self.numCells / self.localSize1d))
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- as big as reduceBuf, because it is a replacement for reduceBuf
	-- ... though I don't accum on vector fields yet, so it doesn't need the x3 really
	if self.allowAccum then
		self:clalloc('accumBuf', app.real, self.numCells * 3)
	end

	-- CL/GL interop buffers

	if app.targetSystem ~= 'console' then
		assert(self.texSize, "you forgot to define solver.texSize for the CL/GL interop buffers")

		-- hmm, notice I'm still keeping the numGhost border on my texture 
		-- if I remove the border altogether then I get wrap-around
		-- maybe I should just keep a border of 1?
		-- for now i'll leave it as it is
		local GLTex2D = require 'gl.tex2d'
		local GLTex3D = require 'gl.tex3d'
		-- use texSize.z instead of dim because meshsolver might not have dim defined correctly
		local cl = self.texSize.z == 1 and GLTex2D or GLTex3D
		-- TODO check for extension GL_ARB_half_float_pixel
		local gltype = app.real == 'half' and gl.GL_HALF_FLOAT_ARB or gl.GL_FLOAT
		self.tex = cl{
			width = tonumber(self.texSize.x),
			height = tonumber(self.texSize.y),
			depth = tonumber(self.texSize.z),
			internalFormat = gl.GL_RGBA32F,
			format = gl.GL_RGBA,
			type = gltype,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_LINEAR,
			wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
		}
		
		local CLImageGL = require 'cl.imagegl'
		if app.useGLSharing then
			self.texCLMem = CLImageGL{context=app.ctx, tex=self.tex, write=true}
		else
			-- use texSize:volume() so the glTexSubImage can use the whole buffer, in the event of meshsolver where texSize:volume can be > numCells
			self.calcDisplayVarToTexPtr = ffi.new(app.real..'[?]', self.texSize:volume() * 3)
		end
	end
end

function SolverBase:getTex(var) 
	return self.tex
end

-- call this when the solver initializes or changes the codePrefix (or changes initCond)
-- it will build the code prefix and refresh everything related to it
-- TODO if you change cons_t then call resetState etc (below the refreshEqnInitState() call a few lines above) in addition to this -- or else your values will get messed up
function SolverBase:refreshEqnInitState()
	-- Right now within eqn:createInitState I'm adding any subclass-specific gui vars
	-- so only after it finishes and all gui vars are created, ask the eqn.initCond object if it wants to modify anything.
	-- Don't do this during Solver:refreshInitStateProgram()->InitCond:getInitCondCode() or the changes won't get into the header.
	-- Hmm... should the initCond even have control over the eqn's vars?
	if self.eqn.initCond.solverVars then
		for k,v in pairs(self.eqn.initCond.solverVars) do
			if self.eqn.guiVars[k] then
				self.eqn.guiVars[k].value = v
			end
		end
	end
	
	self:refreshCodePrefix()

	-- do this after the codePrefix has been created
	if cmdline.checkStructSizes then
		self:checkStructSizes()
	end

	-- bounds don't get set until getInitCondCode() is called, but code prefix needs them ...
	-- TODO do a proper refresh so mins/maxs can be properly refreshed
	local initCond = self.eqn.initCond
	if initCond.mins then 
		self.mins = vec3d(unpack(initCond.mins)) 
		for j=1,3 do
			self.solverPtr.mins.s[j-1] = toreal(self.mins.s[j-1])
		end
	end
	if initCond.maxs then 
		self.maxs = vec3d(unpack(initCond.maxs)) 
		for j=1,3 do
			self.solverPtr.maxs.s[j-1] = toreal(self.maxs.s[j-1])
		end
	end

	-- there's a lot of overlap between this and the solverBuf creation... 
	self:refreshSolverBufMinsMaxs()

	-- while we're here, write all gui vars to the solver_t
	self:copyGuiVarsToBufs()
end

function SolverBase:refreshSolverBufMinsMaxs()
	for j=1,3 do
		self.solverPtr.mins.s[j-1] = toreal(self.mins.s[j-1])
		self.solverPtr.maxs.s[j-1] = toreal(self.maxs.s[j-1])
	end
	if self.app.verbose then
		print('mins = '..fromreal(self.solverPtr.mins.x)..', '..fromreal(self.solverPtr.mins.y)..', '..fromreal(self.solverPtr.mins.z))
		print('maxs = '..fromreal(self.solverPtr.maxs.x)..', '..fromreal(self.solverPtr.maxs.y)..', '..fromreal(self.solverPtr.maxs.z))
	end
end

-- this is the general function - which just assigns the eqn provided by the arg
-- but it can be overridden for specific equations
function SolverBase:createEqn()
	self.eqn = require('hydro.eqn.'..assert(self.eqnName, "expected solver.eqnName or solver args.eqn"))(table(
		self.eqnArgs or {},
		{solver = self}
	))
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
			code = table{
				self.modules:getHeader(self.reqmodules:unpack()),
				self.modules:getCode(self.reqmodules:unpack()),
			}:concat'\n',
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
		-- currently implemented in hydro/solver/roe.cl
		-- not available for any other flux method
		assert(self.fluxBuf)
		self.updateCTUKernelObj = self.solverProgramObj:kernel'updateCTU'
	end


	for _,op in ipairs(self.ops) do
		op:refreshSolverProgram()
	end

	-- display stuff.  build these just in case trackvars is used.
	do	--if self.app.targetSystem ~= 'console' then
	
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
	--if self.app.targetSystem == 'console' then return '' end

--[[
TODO better hasmodule, let eqn and solver etc request modules
and have this search through their requests
that might mean holding a state while adding modules
and that would allow for a templated function to require modules from within module code,
 although it is nice keeping the type, header and cl code separate, so why not keep the dependencies separated as well?
--]]
	local function hasmodule(name)
		return table.find(self.modules.set.math.depends, name)
	end

	local lines = table()

	lines:insert(template([[

typedef union {
	real	ptr[9];
	real	vreal;
<? if hasmodule'sym3' then ?>
	sym3	vsym3;
<? end ?>
<? if hasmodule'cplx' then ?>
	cplx	vcplx;
<? end ?>
<? if hasmodule'real3' then ?>
	real3	vreal3;
<? end ?>
<? if hasmodule'cplx3' then ?>
	cplx3	vcplx3;
<? end ?>
<? if hasmodule'real3x3' then ?>
	real3x3	vreal3x3;
<? end ?>
} displayValue_t;

#define INIT_DISPLAYFUNC()\
	SETBOUNDS(0,0);\
<? if not require 'hydro.solver.meshsolver'.is(solver) then 
?>	int4 dsti = i;\
	int dstindex = index;\
	real3 x = cellBuf[index].pos;\
<? for j=0,solver.dim-1 do 
?>	i.s<?=j?> = clamp(i.s<?=j?>, numGhost, solver->gridSize.s<?=j?> - numGhost - 1);\
<? end
?>	index = INDEXV(i);\
<? else	-- mesh 
?>	int dstindex = index;\
	real3 x = cellBuf[index].pos;\
<? end 		-- mesh vs grid 
?>	displayValue_t value = {.ptr={0}};

<?-- nvidia needed 'write_only', but I don't want to write only -- I want to accumulate and do other operations 
-- TODO if I do accumulate, then I will need to ensure the buffer is initialized to zero ...
?>
#define DISPLAYFUNC_OUTPUTARGS_TEX() write_only <?=
	solver.dim == 3 and 'image3d_t' or 'image2d_t'?> tex

#define DISPLAYFUNC_OUTPUTARGS_BUFFER() global real* dest
]], {
		solver = self,
		hasmodule = hasmodule,
	}))

	local accumFunc = self.displayVarAccumFunc and 'max' or nil
	if accumFunc then
		lines:insert(template([[
#define END_DISPLAYFUNC_TEX()\
	float4 texel = read_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>);\
	texel.x = <?=accumFunc?>(texel.x, value.ptr[0]);\
	if (vectorField) {\
		texel.y = <?=accumFunc?>(texel.y, value.ptr[1]);\
		texel.z = <?=accumFunc?>(texel.z, value.ptr[2]);\
	}\
	write_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>, texel);

#define END_DISPLAYFUNC_BUFFER()\
	if (vectorField) {\
		dest[0+3*dstindex] = <?=accumFunc?>(value.ptr[0], dest[0+3*dstindex]);\
		dest[1+3*dstindex] = <?=accumFunc?>(value.ptr[1], dest[1+3*dstindex]);\
		dest[2+3*dstindex] = <?=accumFunc?>(value.ptr[2], dest[2+3*dstindex]);\
	} else {\
		dest[dstindex] = <?=accumFunc?>(value.ptr[0], dest[dstindex]);\
	}

]],		{
			solver = self,
		}))
	else
		lines:insert(template([[
#define END_DISPLAYFUNC_TEX()\
	write_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>, (float4)(value.ptr[0], value.ptr[1], value.ptr[2], 0.));

#define END_DISPLAYFUNC_BUFFER()\
	if (vectorField) {\
		((global real3*)dest)[dstindex] = value.vreal3;\
	} else {\
		dest[dstindex] = value.vreal;\
	}

]], 	{
			solver = self,
		}))
	end


	local alreadyAddedComponentForGroup = {}
	local function addPickComponetForGroup(var)
		local name = self:getPickComponentNameForGroup(var)
		if alreadyAddedComponentForGroup[name] then return end
		alreadyAddedComponentForGroup[name] = true
		lines:insert((template([[
void <?=name?>(
	constant <?=solver.solver_t?>* solver,
	global const <?=var.bufferType?>* buf,
	int component,
	int* vectorField,
	displayValue_t* value,
	int4 i,
	int index,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	real3 x = cellBuf[index].pos;
	switch (component) {
<? 
for i,component in ipairs(solver.displayComponentFlatList) do
	if not component.onlyFor 
	or (var and var.group and var.group.name == component.onlyFor)
	then
		if hasmodule(component.base) then
?>	case <?=i?>:	//<?=component.base or 'real'?> <?=component.name?>
		{
			<?=component.code?>
			*vectorField = <?= solver:isVarTypeAVectorField(component.type) and '1' or '0' ?>;
			break;
		}
<? 
		end
	end
end
?>	}
}
]], 	{
			name = name,
			solver = self,
			var = var,
			hasmodule = hasmodule,
		})))
	end
	addPickComponetForGroup{
		bufferType = self.eqn.cons_t,
	}

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
				texVsBuf = texVsBuf,
				extraArgs = {
					'int structSize',
					'int structOffset',
				},
				addTab = addTab,
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
				lines:insert('//'..var.name)

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
	displayValue_t* value<?=
	var.extraArgs and #var.extraArgs > 0
		and ',\n\t'..table.concat(var.extraArgs, ',\n\t')
		or ''
?>
) {
<?=addTab(var.codePrefix or '')
?><?=addTab(var.code)
?>}
]], 				{
						solver = self,
						var = var,
						clFuncName = clFuncName,
						addTab = addTab,
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
<?=addTab(var.codePrefix or '')
?><?=addTab(var.code)
?>]],					{
						var = var,
						addTab = addTab,
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
								texVsBuf = 'Tex',
								addTab = addTab,
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
							texVsBuf = 'Buffer',
							addTab = addTab,
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
	self.eqn.initCond:refreshInitStateProgram(self)
end

-- TODO compile the code of CommonCode into a bin of its own
--  and link against it instead of recopying and recompiling
function SolverBase:createCodePrefixHeader()
	-- header
	
	local lines = table()
	
	-- real3
	lines:insert(self.modules:getHeader(self.reqmodules:unpack()))

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
		self:getTypeCode(),
	}
	
	return lines:concat'\n'
end

function SolverBase:createCodePrefixSource()
	lines = table()
if not SolverBase.useCLLinkLibraries then 
	lines:append{
		self.modules:getCode(self.reqmodules:unpack()),
	}
end	
	lines:append{
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
	self.cmds:finish()

	self:applyInitCond()
	self:boundary()
	
	self:resetOps()

	self:constrainU()
end

function SolverBase:constrainU()
	if self.eqn.useConstrainU then
		self.constrainUKernelObj(self.solverBuf, self.UBuf, self.cellBuf)
		if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
		self:boundary()
	end
end

-- override this by the mesh solver ... since I don't know what it will be doing
function SolverBase:applyInitCond()
	self.eqn.initCond:resetState(self)
	if self.allowAccum then
		self.cmds:enqueueFillBuffer{buffer=self.accumBuf, size=ffi.sizeof(self.app.real) * self.numCells * 3}
	end
end

function SolverBase:resetOps()
	for _,op in ipairs(self.ops) do
		if op.resetState then
			op:resetState()
			self:boundary()
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
				fromreal(self.solverPtr.meter),
				fromreal(self.solverPtr.second),
				fromreal(self.solverPtr.kilogram),
				fromreal(self.solverPtr.coulomb),
				fromreal(self.solverPtr.kelvin))
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
	DISPLAYFUNC_OUTPUTARGS_<?=texVsBuf:upper()?>(),
	global const <?=var.bufferType?>* buf,
	int component,
	const global <?=solver.coord.cell_t?>* cellBuf<? 
if require 'hydro.solver.meshsolver'.is(solver) then ?>
	,const global face_t* faces							//[numFaces]<?
end ?><?= var.extraArgs and #var.extraArgs > 0
		and ',\n\t'..table.concat(var.extraArgs, ',\n\t')
		or '' ?>
) {
	INIT_DISPLAYFUNC()
<?=addTab(var.code)
?>	int vectorField = <?=solver:isVarTypeAVectorField(var.type) and '1' or '0'?>;
	<?=solver:getPickComponentNameForGroup(var)?>(solver, buf, component, &vectorField, &value, i, index, cellBuf);
	END_DISPLAYFUNC_<?=texVsBuf:upper()?>()
}
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
	self.color = vec3d(math.random(), math.random(), math.random()):normalize()
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
	kernel:setArg(4, self.solver.cellBuf)
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
		{name = 'mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3),0,0);'},
		{name = 'mag metric', code = 'value->vreal3 = _real3(coordLen(value->vreal3, x),0,0);'},
		{name = 'x', code = 'value->vreal3 = _real3(value->vreal3.x,0,0);'},
		{name = 'y', code = 'value->vreal3 = _real3(value->vreal3.y,0,0);'},
		{name = 'z', code = 'value->vreal3 = _real3(value->vreal3.z,0,0);'},
	})
	self:addDisplayComponents('sym3', {
		{name = 'xx', code = 'value->vsym3 = _sym3(value->vsym3.xx,0,0,0,0,0);'},
		{name = 'xy', code = 'value->vsym3 = _sym3(value->vsym3.xy,0,0,0,0,0);'},
		{name = 'xz', code = 'value->vsym3 = _sym3(value->vsym3.xz,0,0,0,0,0);'},
		{name = 'yy', code = 'value->vsym3 = _sym3(value->vsym3.yy,0,0,0,0,0);'},
		{name = 'yz', code = 'value->vsym3 = _sym3(value->vsym3.yz,0,0,0,0,0);'},
		{name = 'zz', code = 'value->vsym3 = _sym3(value->vsym3.zz,0,0,0,0,0);'},
		{name = 'norm', code = 'value->vsym3 = _sym3(sqrt(sym3_dot(value->vsym3, value->vsym3)), 0,0,0,0,0);'},
		{name = 'tr', code = 'value->vsym3 = _sym3(sym3_trace(value->vsym3), 0,0,0,0,0);'},
		{name = 'det', code = 'value->vsym3 = _sym3(sym3_det(value->vsym3), 0,0,0,0,0);'},
		
		{name = 'x', code = 'value->vsym3 = _sym3(value->vsym3.xx, value->vsym3.xy, value->vsym3.xz, 0,0,0);', type = 'real3', magn='x mag'},
		{name = 'y', code = 'value->vsym3 = _sym3(value->vsym3.xy, value->vsym3.yy, value->vsym3.yz, 0,0,0);', type = 'real3', magn='y mag'},
		{name = 'z', code = 'value->vsym3 = _sym3(value->vsym3.xz, value->vsym3.yz, value->vsym3.zz, 0,0,0);', type = 'real3', magn='z mag'},
		{name = 'x mag', code = 'value->vsym3 = _sym3(real3_len(sym3_x(value->vsym3)), 0,0,0,0,0);'},
		{name = 'y mag', code = 'value->vsym3 = _sym3(real3_len(sym3_y(value->vsym3)), 0,0,0,0,0);'},
		{name = 'z mag', code = 'value->vsym3 = _sym3(real3_len(sym3_z(value->vsym3)), 0,0,0,0,0);'},
		{name = 'x mag metric', code = 'value->vsym3 = _sym3(coordLen(sym3_x(value->vsym3), x), 0,0,0,0,0);'},
		{name = 'y mag metric', code = 'value->vsym3 = _sym3(coordLen(sym3_y(value->vsym3), x), 0,0,0,0,0);'},
		{name = 'z mag metric', code = 'value->vsym3 = _sym3(coordLen(sym3_z(value->vsym3), x), 0,0,0,0,0);'},
	})
	self:addDisplayComponents('cplx', {
		{name = 'default', type='cplx', magn='abs'},
		{name = 're', code = 'value->vcplx = cplx_from_real(value->vcplx.re);'},
		{name = 'im', code = 'value->vcplx = cplx_from_real(value->vcplx.im);'},
		{name = 'abs', code = 'value->vcplx = cplx_from_real(cplx_abs(value->vcplx));'},
		{name = 'arg', code = 'value->vcplx = cplx_from_real(cplx_arg(value->vcplx));'},
	})
	self:addDisplayComponents('cplx3', {
		{name = 'mag', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx3_len(value->vcplx3)), cplx_zero, cplx_zero);'},
		{name = 'mag metric', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx3_weightedLenSq(value->vcplx3, coord_g_ll(x))), cplx_zero, cplx_zero);'},

		{name = 'x', code='value->vcplx3 = _cplx3(value->vcplx3.x, cplx_zero, cplx_zero);', type='cplx', magn='x abs'},
		{name = 'x abs', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx_abs(value->vcplx3.x)), cplx_zero, cplx_zero);'},
	
		{name = 'y', code='value->vcplx3 = _cplx3(value->vcplx3.y, cplx_zero, cplx_zero);', type='cplx', magn='y abs'},
		{name = 'y abs', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx_abs(value->vcplx3.y)), cplx_zero, cplx_zero);'},
		
		{name = 'z', code='value->vcplx3 = _cplx3(value->vcplx3.z, cplx_zero, cplx_zero);', type='cplx', magn='z abs'},
		{name = 'z abs', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx_abs(value->vcplx3.z)), cplx_zero, cplx_zero);'},
		
		{name = 'x arg', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx_arg(value->vcplx3.x)), cplx_zero, cplx_zero);'},
		{name = 'y arg', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx_arg(value->vcplx3.y)), cplx_zero, cplx_zero);'},
		{name = 'z arg', code = 'value->vcplx3 = _cplx3(cplx_from_real(cplx_arg(value->vcplx3.z)), cplx_zero, cplx_zero);'},
		{name = 'x re', code = 'value->vcplx3 = _cplx3(cplx_from_real(value->vcplx3.x.re), cplx_zero, cplx_zero);'},
		{name = 'x im', code = 'value->vcplx3 = _cplx3(cplx_from_real(value->vcplx3.x.im), cplx_zero, cplx_zero);'},
		{name = 'y re', code = 'value->vcplx3 = _cplx3(cplx_from_real(value->vcplx3.y.re), cplx_zero, cplx_zero);'},
		{name = 'y im', code = 'value->vcplx3 = _cplx3(cplx_from_real(value->vcplx3.y.im), cplx_zero, cplx_zero);'},
		{name = 'z re', code = 'value->vcplx3 = _cplx3(cplx_from_real(value->vcplx3.z.re), cplx_zero, cplx_zero);'},
		{name = 'z im', code = 'value->vcplx3 = _cplx3(cplx_from_real(value->vcplx3.z.im), cplx_zero, cplx_zero);'},
		
		{name = 're', code = 'value->vreal3 = cplx3_re(value->vcplx3); *(real3*)(value+3) = real3_zero;', type = 'real3', magn='re mag'},
		{name = 're mag', code = 'value->vcplx3 = _cplx3(cplx_from_real(real3_len(cplx3_re(value->vcplx3))), cplx_zero, cplx_zero);'},
		{name = 're mag metric', code = 'value->vcplx3 = _cplx3(cplx_from_real(coordLen(cplx3_re(value->vcplx3), x)), cplx_zero, cplx_zero);'},
		{name = 'im', code = 'value->vreal3 = cplx3_im(value->vcplx3); *(real3*)(value+3) = real3_zero;', type = 'real3', magn='im mag'},
		{name = 'im mag', code = 'value->vcplx3 = _cplx3(cplx_from_real(real3_len(cplx3_im(value->vcplx3))), cplx_zero, cplx_zero);'},
		{name = 'im mag metric', code = 'value->vcplx3 = _cplx3(cplx_from_real(coordLen(cplx3_im(value->vcplx3), x)), cplx_zero, cplx_zero);'},
	})
	self:addDisplayComponents('real3x3', {
		{name = 'xx', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.x.x, 0,0,0,0,0,0,0,0);'},
		{name = 'xy', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.x.y, 0,0,0,0,0,0,0,0);'},
		{name = 'xz', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.x.z, 0,0,0,0,0,0,0,0);'},
		
		{name = 'yx', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.y.x, 0,0,0,0,0,0,0,0);'},
		{name = 'yy', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.y.y, 0,0,0,0,0,0,0,0);'},
		{name = 'yz', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.y.z, 0,0,0,0,0,0,0,0);'},
		
		{name = 'zx', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.z.x, 0,0,0,0,0,0,0,0);'},
		{name = 'zy', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.z.y, 0,0,0,0,0,0,0,0);'},
		{name = 'zz', code = 'value->vreal3x3 = _real3x3(value->vreal3x3.z.z, 0,0,0,0,0,0,0,0);'},
		
		{name = 'norm', code = 'value->vreal3x3 = _real3x3(sqrt(real3x3_dot(value->vreal3x3, value->vreal3x3)), 0,0,0,0,0,0,0,0);'},
		{name = 'tr', code = 'value->vreal3x3 = _real3x3(real3x3_trace(value->vreal3x3), 0,0,0,0,0,0,0,0);'},
		{name = 'tr metric', code = 'value->vreal3x3 = _real3x3(real3x3_sym3_dot(value->vreal3x3, coord_g_ll(x)), 0,0,0,0,0,0,0,0);'},
		
		{name = 'x', code = 'value->vreal3 = value->vreal3x3.x; value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='x mag'},
		{name = 'y', code = 'value->vreal3 = value->vreal3x3.y; value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='y mag'},
		{name = 'z', code = 'value->vreal3 = value->vreal3x3.z; value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='z mag'},
		{name = 'x mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3x3.x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'y mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3x3.y), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'z mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3x3.z), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'x mag metrc', code = 'value->vreal3 = _real3(coordLen(value->vreal3x3.x, x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'y mag metrc', code = 'value->vreal3 = _real3(coordLen(value->vreal3x3.y, x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'z mag metrc', code = 'value->vreal3 = _real3(coordLen(value->vreal3x3.z, x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		
		{name = 'T x', code = 'value->vreal3 = _real3(value->vreal3x3.x.x, value->vreal3x3.y.x, value->vreal3x3.z.x); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='T x mag'},
		{name = 'T y', code = 'value->vreal3 = _real3(value->vreal3x3.x.y, value->vreal3x3.y.y, value->vreal3x3.z.y); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='T y mag'},
		{name = 'T z', code = 'value->vreal3 = _real3(value->vreal3x3.x.z, value->vreal3x3.y.z, value->vreal3x3.z.z); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='T z mag'},
		{name = 'T x mag', code = 'value->vreal3 = _real3(real3_len(_real3(value->vreal3x3.x.x, value->vreal3x3.y.x, value->vreal3x3.z.x)), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T y mag', code = 'value->vreal3 = _real3(real3_len(_real3(value->vreal3x3.x.y, value->vreal3x3.y.y, value->vreal3x3.z.y)), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T z mag', code = 'value->vreal3 = _real3(real3_len(_real3(value->vreal3x3.x.z, value->vreal3x3.y.z, value->vreal3x3.z.z)), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T x mag metric', code = 'value->vreal3 = _real3(coordLen(_real3(value->vreal3x3.x.x, value->vreal3x3.y.x, value->vreal3x3.z.x), x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T y mag metric', code = 'value->vreal3 = _real3(coordLen(_real3(value->vreal3x3.x.y, value->vreal3x3.y.y, value->vreal3x3.z.y), x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T z mag metric', code = 'value->vreal3 = _real3(coordLen(_real3(value->vreal3x3.x.z, value->vreal3x3.y.z, value->vreal3x3.z.z), x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
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
	
		local name = assert(var.name, "expected to find name in "..require 'ext.tolua'(var))
		local code = assert(var.code, "expected to find code")

		args = table(args, {
			solver = self,
			name = group.name .. ' ' .. name,
			code = code,
			units = units,
			group = group,
			type = var.type or 'real',
		})

		-- enable the first scalar field
		-- also enable the first vector field on non-1D simulations
		local enabled
		
		-- TODO how about saving somewhere what should be enabled by default?
		-- TODO pick predefined somewhere?
		if cmdline.displayvars or self.eqn.predefinedDisplayVars then
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
		vars = {{name='0', code='value.vreal = buf[index];'}},
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
		args.vars = self.eqn:getDisplayVarsForStructVars(self.eqn.consStruct.vars)
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
				dupvar.varIndex = #self.displayVars
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
	if cmdline.displayvars or self.eqn.predefinedDisplayVars then
		for _,name in ipairs(cmdline.displayvars and string.split(cmdline.displayvars, ',') or self.eqn.predefinedDisplayVars) do
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
self.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * volume * channels, ptr=ptr}
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
	local min = fromreal(self.reduceMin(nil, volume*channels))
	
	self:calcDisplayVarToBuffer(var, componentIndex)
	local max = fromreal(self.reduceMax(nil, volume*channels))
--print('reduce min',min,'max',max,'volume',volume,'name',var.name,'channels',channels)
	var.lastMin = min
	var.lastMax = max
	var.lastAvg = nil	-- invalidate

	return min, max
end

end

-- used by the output to print out avg, min, max
function SolverBase:calcDisplayVarRangeAndAvg(var, componentIndex)
	componentIndex = componentIndex or var.component

	if var.lastTime ~= self.t then
		--don't set lastTime yet -- instead let calcDisplayVarRange update and do this:
		-- var.lastTime = self.t

		-- this will update lastTime if necessary
		self:calcDisplayVarRange(var, componentIndex)
		-- displayVarGroup has already set up the appropriate args
	end

	if not var.lastAvg then
		-- duplicated in calcDisplayVarRange
		local size = self.numCells
		local sizevec = var.getBuffer().sizevec
		if sizevec then
			size = tonumber(sizevec:volume())
		end

		self:calcDisplayVarToBuffer(var, componentIndex)
		
		-- [[ avg
		local lastAvg = self.reduceSum(nil, size) / tonumber(size)
		var.lastAvg = fromreal(lastAvg)
		--]]
		--[[ rms
		self.squareKernelObj(self.solverBuf, fromreal(self.reduceBuf))
		var.lastAvg = math.sqrt(self.reduceSum(nil, size) / tonumber(size))
		--]]
	end

	return var.lastMin, var.lastMax, var.lastAvg
end


-------------------------------------------------------------------------------
--                              gui                                          --
-------------------------------------------------------------------------------



function SolverBase:initDraw()
end

-- some common draw code that everyone seems to use
-- can GLSL use #include ?
function SolverBase:getGradientGLSLCode()
	return template([[
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>

uniform bool useLog;
uniform float valueMin;
uniform float valueMax;
uniform sampler1D gradientTex;

float logmap(float x) {
	return log(1. + abs(x)) * _1_LN_10;
}

float getGradientFrac(float value) {
	if (useLog) {
		
		// TODO all of this on CPU before setting the uniforms
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		if ((valueMin < 0.) != (valueMax < 0.)) {
			logValueMin = 0.;	//logmap(0)
		} else {
			if (valueMin < 0.) {
				float tmp = logValueMin;
				logValueMin = logValueMax;
				logValueMax = tmp;
			}
		}
		
		return (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		return (value - valueMin) / (valueMax - valueMin);
	}
}

//TODO just change texel lookup in gradTex?
float getGradientTexCoord(float frac) {
	return (frac * <?=clnumber(app.gradientTex.width-1)?> + .5) * <?=clnumber(1 / app.gradientTex.width)?>;
}

vec4 getGradientColor(float value) {
	float frac = getGradientFrac(value);
	float tc = getGradientTexCoord(frac);
	return texture1D(gradientTex, tc);
}

]], {
		clnumber = require 'cl.obj.number',
		app = self.app,
	})
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
		self.calcDTKernelObj.obj:setArg(3, self.cellBuf)
		self.calcDTKernelObj()
		dt = self.cfl * fromreal(self.reduceMin())
		if not math.isfinite(dt) then
			print("got a bad dt at time "..self.t) -- TODO dump all buffers
		end
		self.fixedDT = dt
	end
	return dt
end

-- how often to update the console
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
		local tick = cmdline.tick or 1
		local deltaTime = thisTime - self.lastFrameTime
		local fps = 1 / deltaTime
		self.fpsIndex = (self.fpsIndex % self.fpsNumSamples) + 1
		self.fpsSamples[self.fpsIndex] = fps
		self.fps = self.fpsSamples:sum() / #self.fpsSamples
		if (self.showFPS or cmdline.trackvars)
		and math.floor(thisTime / tick) ~= math.floor(self.lastFrameTime / tick) then 
			--io.write(tostring(self), ' ')
			local sep = ''
			if self.showFPS then
				io.write(sep, 'fps=', self.fps)
				sep = '\t'
			end
			io.write(sep, 't=', self.t)
			sep = '\t'
			if cmdline.trackvars then
				local varnames = string.split(cmdline.trackvars, ','):map(string.trim)
				if varnames:find'dt' then
					io.write(sep, 'dt=', self.dt)
				end
				for _,varname in ipairs(varnames) do
					if varname ~= 'dt' then
						local var = assert(self.displayVarForName[varname], "couldn't find "..varname)
						local ymin, ymax, yavg = self:calcDisplayVarRangeAndAvg(var)
						io.write(sep, varname, '=[', ymin, '..', yavg, '..', ymax, ']')
					end
				end
			end
			print()
		end
	end
	self.lastFrameTime = thisTime

	if self.checkNaNs then 
		assert(self:checkFinite(self.UBufObj)) 
	end

	-- [[ stop condition: '|U H| > 1'
	-- math should be in namespace ... but setfenv seems to run really really slow
	-- what about min vs max vs avg?  dereference them as table fields
	-- use || as shorthand for max(abs(${x}.min), abs(${x}.max))
	if cmdline.stopcond then
		-- cache in global namespace
		if not stopcondfunc then
			local code = cmdline.stopcond
			local prefix = table()
			prefix:insert[[
local self = ...
]]
			local varindex = 1
			code = code:gsub('|([^|]*)|', function(name)
				local varmin = 'varmin'..varindex
				local varmax = 'varmax'..varindex
				local varabsmax = 'varabsmax'..varindex
				prefix:insert(template([[
local <?=varmin?>, <?=varmax?> = self:calcDisplayVarRange(self.displayVarForName['<?=name?>'])
local <?=varabsmax?> = math.max(math.abs(<?=varmin?>, <?=varmax?>))
]],				{
					varmin = varmin,
					varmax = varmax,
					varabsmax = varabsmax,
					name = name,
				}))
				varindex = varindex + 1
				return varabsmax
			end)
			code = prefix:concat'\n'..'\nreturn '..code
			stopcondfunc = assert(load(code))
		end
		if stopcondfunc(self) then
			error(cmdline.stopcond..' met at t='..self.t)
		end
	end
	--]]



	-- the rest of this used to be in gridsolver
	--  so if anything gets errors from this (like composite solvers, etc)
	-- then ... mve the stuff above this line into its own function, and only call that, and don't call SolverBase.update



	local dt = self:calcDT()

	-- first do a step
	self:step(dt)

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	-- why was this moved out of :step() ?
	self.t = self.t + dt
	self.dt = dt

	if cmdline.testAccuracy then
		local err = self:calcExactError()
		if #self.app.solvers > 1 then
			io.write(self.name,'\t')
		end
		print('t='..self.t..' L1-error='..err)
	end

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	self:boundary()

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
end


function SolverBase:step(dt)

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	self.integrator:integrate(dt, function(derivBufObj)

if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
	
		if self.calcDeriv then
			self:calcDeriv(derivBufObj, dt)
		end

if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

		if self.eqn.useSourceTerm then
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
			self.addSourceKernelObj(self.solverBuf, derivBufObj.obj, self.UBuf, self.cellBuf)
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

		for _,op in ipairs(self.ops) do
			if op.addSource then
				op:addSource(derivBufObj)
			end
		end
	end)

	-- I moved the constrainU() and boundary() functions from here to be within the integrator, 
	-- so that each sub-step that the RK integrator performs can apply these and be physically correct states.

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	for _,op in ipairs(self.ops) do
		if op.step then
			self:boundary()
			op:step(dt)
		end
	end

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
end


-- check for nans
-- expects buf to be of type cons_t, made up of numStates real variables
function SolverBase:checkFinite(buf)
	local ptrorig = buf:toCPU()
	local ptr0size = tonumber(ffi.sizeof(buf.type))
	local realSize = tonumber(ffi.sizeof'real')
	local ptrsPerReal = ptr0size / realSize
	assert(ptrsPerReal == math.floor(ptrsPerReal))
	-- don't free the original ptr too soon
	local ptr = ffi.cast('real*', ptrorig)
	local size = buf.count * ptrsPerReal
	local found
	for i=0,size-1 do
		local x = tonumber(ptr[i])
		if not math.isfinite(x) then
			found = found or table()
			local ins = {i, x}
			-- for certain bufs show the field
			-- TODO associate each type with the array of fields creating the struct, then reverse lookup on arbitrary types to find the field
			if buf == self.UBufObj then
				local vars = self.eqn.consStruct.vars
				local numScalars = self.eqn.consStruct:countScalars()
				local offset = (i % numScalars)
				local cellIndex = (i - offset) / numScalars
				local cellpos = vec3sz()
				cellpos.x = tonumber(cellIndex % self.gridSize.x)
				cellIndex = (cellIndex - cellpos.x) / self.gridSize.x
				cellpos.y = tonumber(cellIndex % self.gridSize.y)
				cellIndex = (cellIndex - cellpos.y) / self.gridSize.y
				cellpos.z = tonumber(cellIndex)
				assert(cellpos.z < self.gridSize.z)
				table.insert(ins, tostring(cellpos))
				
				offset = offset * ffi.sizeof'real'
				local field
				for _,var in ipairs(vars) do
					offset = offset - ffi.sizeof(var.type)
					if offset < 0 then
						field = var.name
						break
					end
				end
				-- assert(field, "shouldn't have got to this point")
				table.insert(ins, field)
			end
			found:insert(ins)
		end
	end
	if not found then return true end
--	self:printBuf(nil, ptr)
	return false, 'found non-finite offsets and numbers: '..require 'ext.tolua'(found)..' at t='..self.t
end

function SolverBase:printBuf(buf, ptrorig, colsize, colmax)
	ptrorig = ptrorig or buf:toCPU()
	local ptr0size = tonumber(ffi.sizeof(buf.type))
--print('ptr0size', ptr0size)
	local realSize = tonumber(ffi.sizeof'real')
--print('realSize', realSize)
	local ptrsPerReal = ptr0size / realSize
--print('ptrsPerReal', ptrsPerReal)
	assert(ptrsPerReal == math.floor(ptrsPerReal))
	-- I've seen this problem before, in gcmem ... if you assign a ptr to a cast of itself, luajit can segfault
	-- fix?  save the old pointer, and luajit doesn't try to free it too early
	local ptr = ffi.cast('real*', ptrorig)
	local size = buf.count * ptrsPerReal
--print('size', size)
	if buf.type == self.eqn.cons_t then
		local maxdigitlen = #tostring(self.numCells-1)
--print('maxdigitlen', maxdigitlen)
		local realsPerCell = math.floor(size / self.numCells)
--print('realsPerCell', realsPerCell)
		colmax = colmax or realsPerCell
--print('colmax', colmax)
		if colmax > realsPerCell then
			error("got too many realsPerCell\n"..require 'ext.tolua'{
				colmax = colmax,
				realsPerCell = realsPerCell,
				ptr0size = ptr0size,
				realSize = realSize,
				ptrsPerReal = ptrsPerReal,
				size = size,
			})
		end
		for i=0,self.numCells-1 do
			io.write((' '):rep(maxdigitlen-#tostring(i)), i,':')
			for j=0,colmax-1 do
				io.write('\t',
					--..(j==0 and '[' or '')..
					('%f'):format(fromreal(ptr[j + realsPerCell * i]))
					--..(j==self.eqn.numStates-1 and ']' or ',')
				)
				if self.app.real == 'half' then
					io.write(('/0x%x'):format(ptr[j + realsPerCell * i].i))
				end
			end 
			print()
		end
	else
		local maxdigitlen = #tostring(size-1)
		colsize = colsize or 1
		for i=0,size-1 do
			if i % colsize == 0 then
				io.write((' '):rep(maxdigitlen-#tostring(i)), i,':')
			end
			io.write(' ', ('%f'):format(ptr[i]))
			if i % colsize == colsize-1 then 
				print() 
			end
		end
		if size % colsize ~= 0 then print() end
	end
	return ptrorig 	-- in case someone wants it
end


function SolverBase:calcDisplayVarToTex(var, componentIndex)
	componentIndex = componentIndex or var.component
	local component = self.displayComponentFlatList[componentIndex]
	local vectorField = self:isVarTypeAVectorField(component.type)
	
	local app = self.app
	local displayVarGroup = var.displayVarGroup
	if app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		self.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	
		var:setToTexArgs()
		var.calcDisplayVarToTexKernelObj()
		
		self.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		self.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		-- calcDisplayVarToTexPtr is sized without ghost cells 
		-- so is the GL texture
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex
		
		local channels = vectorField and 3 or 1
		local format = vectorField and gl.GL_RGB or gl.GL_RED
		
		var:setToBufferArgs()
		
		self:calcDisplayVarToBuffer(var)

		local sizevec = var.getBuffer().sizevec or self.texSize
		assert(sizevec.x <= tex.width and sizevec.y <= tex.height)
		local volume = tonumber(sizevec:volume())
		
		-- use 'volume' here reads the whole texture, even if it exceeds the number of cells.  calcDisplayVarToTexPtr should be big enough.
		-- use 'solver.numCells' to just read what is needed
		self.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * self.numCells * channels, ptr=ptr}
		local destPtr = ptr
		-- TODO check for extension GL_ARB_half_float_pixel
		local gltype = app.real == 'half' and gl.GL_HALF_FLOAT_ARB or gl.GL_FLOAT
		
		if app.real == 'double' then
			-- can this run in place?  seems like it
			destPtr = ffi.cast('float*', ptr)
			for i=0,volume*channels-1 do
				destPtr[i] = ptr[i]
			end
		end
		tex:bind()
		
		-- use texSize because meshsolver dim isn't reliable
		if self.texSize.z == 1 then
			gl.glTexSubImage2D(gl.GL_TEXTURE_2D, 0, 0, 0, sizevec.x, sizevec.y, format, gltype, destPtr)
		else
			for z=0,tex.depth-1 do
				gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, sizevec.x, sizevec.y, 1, format, gltype, destPtr + channels * sizevec.x * sizevec.y * z)
			end
		end
		tex:unbind()
		glreport'here'
	end
end

-- this is abstracted because accumBuf might want to be used ...
function SolverBase:calcDisplayVarToBuffer(var, componentIndex)
	componentIndex = componentIndex or var.component
	local component = self.displayComponentFlatList[componentIndex]
	local vectorField = self:isVarTypeAVectorField(component.type)
	local channels = vectorField and 3 or 1

	-- duplicated in calcDisplayVarRange
	local volume = self.numCells
	local sizevec = var.getBuffer().sizevec
	if sizevec then
		volume = tonumber(sizevec:volume())
	end
	
	if self.displayVarAccumFunc	then
		self.cmds:enqueueCopyBuffer{src=self.accumBuf, dst=self.reduceBuf, size=ffi.sizeof(app.real) * volume * channels}
	end
	var:setToBufferArgs()
	var.calcDisplayVarToBufferKernelObj.obj:setArg(1, self.reduceBuf)
	var.calcDisplayVarToBufferKernelObj.obj:setArg(3, int(componentIndex))
	var.calcDisplayVarToBufferKernelObj()
	if self.displayVarAccumFunc then
		self.cmds:enqueueCopyBuffer{src=self.reduceBuf, dst=self.accumBuf, size=ffi.sizeof(app.real) * volume * channels}
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
			self.cmds:enqueueFillBuffer{buffer=self.accumBuf, size=ffi.sizeof(self.app.real) * self.numCells * 3}
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
		if op.updateGUI then
			ig.igPushIDInt(i)
			op:updateGUI()
			ig.igPopID()
		end
	end

	if tooltip.comboTable('flux limiter', self, 'fluxLimiter', self.app.limiterNames) then
		self:refreshSolverProgram()
	end
end

function SolverBase:updateGUIEqnSpecific()
--[[ TODO why is this crashing
	if tooltip.comboTable('init state', self, 'initCondIndex', self.eqn.initCondNames) then
		-- TODO hmm ... the whole point of making a separate initCondProgram was to be able to refresh it without rebuilding all of the solver ...
		-- TODO try again once initCond_t is separated from solver_t
		self:refreshEqnInitState()
	end	
--]]		
	for _,var in ipairs(self.eqn.guiVars) do
		var:updateGUI(self)
	end
end

do
	local function handle(self, var, title)
		local anyChanged = false
		ig.igPushIDStr(title)

		var.enabled = not not var.enabled
		local enableChanged = tooltip.checkboxTable('enabled', var, 'enabled') 
		anyChanged = anyChanged or enableChanged
		ig.igSameLine()
		
		anyChanged = anyChanged or tooltip.checkboxTable('log', var, 'useLog')
		ig.igSameLine()

		anyChanged = anyChanged or tooltip.checkboxTable('units', var, 'showInUnits') 
		ig.igSameLine()

		anyChanged = anyChanged or tooltip.checkboxTable('fixed range', var, 'heatMapFixedRange')
		ig.igSameLine()

		--tooltip.comboTable('component', var, 'component', self.displayComponentNames)
	
		local name = var.name
		if var.units then name = name..' '..var.units end
		if ig.igCollapsingHeader(name) then
			local unitScale = 1
			if var.units and var.showInUnits then -- convert our ranges from raw to units
				unitScale = self:convertToSIUnitsCode(var.units).func()
				var.heatMapValueMin = var.heatMapValueMin * unitScale
				var.heatMapValueMax = var.heatMapValueMax * unitScale
			end
			anyChanged = anyChanged or tooltip.numberTable('value min', var, 'heatMapValueMin')
			anyChanged = anyChanged or tooltip.numberTable('value max', var, 'heatMapValueMax')
			if var.units and var.showInUnits then -- convert our ranges from raw to units
				var.heatMapValueMin = var.heatMapValueMin / unitScale
				var.heatMapValueMax = var.heatMapValueMax / unitScale
			end
		end
		
		ig.igPopID()
	
		return enableChanged, anyChanged
	end

	-- do one for 'all'
	local op = require 'ext.op'
	local fields = {
		'enabled',
		'useLog',
		'showInUnits',
		'heatMapFixedRange',
		'heatMapValueMin',
		'heatMapValueMax',
	}
	local defaults = {
		true,
		true,
		true,
		true,
		math.huge,
		-math.huge,
	}
	local combines = {
		op.land,
		op.land,
		op.land,
		op.land,
		math.min,
		math.max,
	}
	local all = {name='all'}
	for i=1,#fields do
		all[fields[i]] = defaults[i]
	end
	local original = {}
	for i,field in ipairs(fields) do
		original[field] = defaults[i]
	end

	local function search(str, what)
		return str:match(what:gsub('%s+', '.*'))
	end

	function SolverBase:updateGUIDisplay()
		if self.guiDisplayFilterStr == nil then self.guiDisplayFilterStr = '' end
		if self.guiDisplayFilterEnabledVars == nil then self.guiDisplayFilterEnabledVars = false end
		local refresh 
			
		tooltip.textTable('filter', self, 'guiDisplayFilterStr')
		ig.igSameLine()
		tooltip.checkboxTable('filter enabled', self, 'guiDisplayFilterEnabledVars')
		
		for i,displayVarGroup in ipairs(self.displayVarGroups) do
			ig.igPushIDStr('display '..i)
			if ig.igCollapsingHeader(displayVarGroup.name) then
				for i=1,#fields do
					all[fields[i]] = defaults[i]
				end
				for _,var in ipairs(displayVarGroup.vars) do
					local title = displayVarGroup.name..' '..var.name
					if (#self.guiDisplayFilterStr == 0 or search(title, self.guiDisplayFilterStr))
					and (not self.guiDisplayFilterEnabledVars or var.enabled)
					then
						for i,field in ipairs(fields) do
							all[field] = combines[i](all[field], var[field])
						end
					end
				end
				for _,field in ipairs(fields) do
					original[field] = all[field]
				end
				local enableChanged, anyChanged = handle(self, all, 'all')
				--refresh = refresh or enableChanged
				if anyChanged then
					for _,field in ipairs(fields) do
						if all[field] ~= original[field] then
							for _,var in ipairs(displayVarGroup.vars) do
								local title = displayVarGroup.name..' '..var.name
								if (#self.guiDisplayFilterStr == 0 or search(title, self.guiDisplayFilterStr))
								and (not self.guiDisplayFilterEnabledVars or var.enabled)
								then
									var[field] = all[field]
								end
							end
						end
					end
				end

				for _,var in ipairs(displayVarGroup.vars) do
					local title = displayVarGroup.name..' '..var.name
					if (#self.guiDisplayFilterStr == 0 or search(title, self.guiDisplayFilterStr))
					and (not self.guiDisplayFilterEnabledVars or var.enabled)
					then
						local enableChanged = handle(self, var, title)
						--refresh = refresh or enableChanged
					end
				end
			end
			ig.igPopID()
		end
	end
	if refresh then
		-- solver and display are now tied together
		self:refreshSolverProgram()
	end
end

function SolverBase:updateGUI()
	if ig.igCollapsingHeader'parameters:' then
		self:updateGUIParams()
	end
	if ig.igCollapsingHeader'equation:' then
		self:updateGUIEqnSpecific()
	end
	if ig.igCollapsingHeader'display:' then
		self:updateGUIDisplay()
	end

	-- heat map var

	-- TODO volumetric var
end

-- [[ debugging -- determine sizeof
function SolverBase:checkStructSizes_getTypes()
	return 	table{
		'real',
		'real2',
		'real3',
		'real4',
		self.eqn.consStruct or self.eqn.cons_t,
		self.eqn.primStruct or self.eqn.prim_t,
		self.eqn.eigen_t,
		self.eqn.waves_t,
		self.solverStruct,
	}
end
function SolverBase:checkStructSizes()
	local typeinfos = self:checkStructSizes_getTypes()

	local varcount = 0
	for _,typeinfo in ipairs(typeinfos) do
		varcount = varcount + 1
		if Struct.is(typeinfo) then
			varcount = varcount + #typeinfo.vars
		end
	end

	local cmd = self.cmds
	local _1x1_domain = self.app.env:domain{size={1}, dim=1}
	local resultPtr = ffi.new('size_t[?]', varcount)
	local resultBuf = self.app.env:buffer{name='result', type='size_t', count=varcount, data=resultPtr}

	--print(self.codePrefix)
	require 'cl.obj.kernel'{
		env = self.app.env,
		domain = _1x1_domain,
		argsOut = {resultBuf},
		header = self.codePrefix,
		body = template([[
#define offsetof __builtin_offsetof

<? 
local index = 0
for i,typeinfo in ipairs(typeinfos) do 
	local typename
	if type(typeinfo) == 'string' then
?>	result[<?=index?>] = sizeof(<?=typeinfo?>);
<? 
		index = index + 1
	else
?>	result[<?=index?>] = sizeof(<?=typeinfo.typename?>);
<? 
		index = index + 1
		for _,var in ipairs(typeinfo.vars) do
?>	result[<?=index?>] = offsetof(<?=typeinfo.typename?>, <?=var.name?>);
<? 		
			index = index + 1
		end
	end
end 
?>
]], 	{
			typeinfos = typeinfos,
		}),
	}()
	resultBuf:toCPU(resultPtr)
	local index = 0
	for i,typeinfo in ipairs(typeinfos) do
		if type(typeinfo) == 'string' then
			local clsize = tostring(resultPtr[index]):match'%d+'
			index = index + 1
			local ffisize = tostring(ffi.sizeof(typeinfo))
			print('sizeof('..typeinfo..'): OpenCL='..clsize..', ffi='..ffisize..(clsize == ffisize and '' or ' -- !!!DANGER!!!'))
		else
			local clsize = tostring(resultPtr[index]):match'%d+'
			index = index + 1
			local ffisize = tostring(ffi.sizeof(typeinfo.typename))
			print('sizeof('..typeinfo.typename..'): OpenCL='..clsize..', ffi='..ffisize..(clsize == ffisize and '' or ' -- !!!DANGER!!!'))
			
			for _,var in ipairs(typeinfo.vars) do
				local cloffset = tostring(resultPtr[index]):match'%d+'
				index = index + 1
				local ffioffset = tostring(ffi.offsetof(typeinfo.typename, var.name))
				print('offsetof('..typeinfo.typename..', '..var.name..'): OpenCL='..cloffset..', ffi='..ffioffset..(cloffset == ffioffset and '' or ' -- !!!DANGER!!!'))
			end
		end
	end
	os.exit()
end
--]]

return SolverBase
