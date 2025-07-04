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
				-- but coord.sphere_sinh_radial needs this after .maxs is finalized, after getInitCondCode() is called
				self.coord = ...		<- this object creation is wedged between the other mesh vars because meshsolver needs it early

				self.device = ...
				self.cmds = ...
				self.color = ...
				self.ops = ...
				self.solverStructFields = ...
				self.solverStructFields:append(...)
			self.solverStructFields:append(...)
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

				self.solverStruct = ...				-- now that initCond does not modify solver_t, this can be moved back -- but eqn still modifies it, so not back by far
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

	-- code modules will collect the typedef info, the header info, and the code info
	-- subsequent cl programs will collect the required ones for building
	SolverBase:initCodeModules
		self.modules = ...
		self.sharedModulesEnabled = ...

		self.coord:initCodeModules
		self.eqn:initCodeModules
		self.solverStruct module
		self.ops[i]:initCodeModules

	SolverBase:initCodeModuleDisplay
		SolverBase:createDisplayVars
			SolverBase:createDisplayComponents
			SolverBase:finalizeDisplayComponents
			FiniteVolume:addDisplayVars
				SolverBase:addDisplayVars
			SolverBase:finalizeDisplayVars
		SolverBase:getDisplayCode
			self.displayVarGroups[i].vars[j].toTexKernelName = ...

--------- here is where the ffi.cdef is called ---------

	SolverBase:initCDefs

	SolverBase:postInit
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

--------- here's where programs are created ---------

								self.solverProgramObj = ...

--------- here's where kernels are created ---------

								SolverBase:refreshCalcDTKernel
									self.calcDTKernelObj = ...
								self.addSourceKernelObj = ...
								self.constrainUKernelObj = ...
								self.updateCTUKernelObj = ...
								self.ops[i]:refreshSolverProgram
								self.displayVarGroups[i].calcDisplayVarToTexKernelObj = ...
								self.displayVarGroups[i].calcDisplayVarToBufferKernelObj = ...
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
local ig = require 'imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local path = require 'ext.path'
local math = require 'ext.math'
local tolua = require 'ext.tolua'
local range = require 'ext.range'
local os = require 'ext.os'
local gl = require 'gl'
local glreport = require 'gl.report'
local CLBuffer = require 'cl.obj.buffer'
local template = require 'template'
local vec3d = require 'vec-ffi.vec3d'
local vec3sz = require 'vec-ffi.vec3sz'
local roundup = require 'hydro.util.roundup'
local timer = require 'ext.timer'.timer
local getTime = require 'ext.timer'.getTime
local Struct = require 'struct'
local HydroStruct = require 'hydro.code.struct'

local common = require 'hydro.common'	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6
local from6to3x3 = common.from6to3x3
local sym = common.sym

-- TODO call this 'torealparam' instead?
local real = require 'hydro.real'

local half = require 'cl.obj.half'
local toreal, fromreal = half.toreal, half.fromreal


-- whether to cache the opencl binaries
local useCache = cmdline.useCache
if useCache == nil then useCache = true end


local function addTab(s)
	s = tostring(s)
	if s:sub(-1) ~= '\n' then s = s .. '\n' end
	s = '\t'..s:gsub('\n([^\n])', '\n\t%1')
	return s
end


local integrators = require 'hydro.int.all'
local integratorNames = integrators:mapi(function(integrator) return integrator.name end)

local useClang = cmdline.useClang

local SolverBase = class()

SolverBase.name = 'solverbase'

-- override to specify which hydro/eqn/*.lua to use as the equation
SolverBase.eqnName = nil


-- singleton mapping from struct C name to Struct object
-- TODO move this to Struct?
-- but maybe don't do this until you fix / remove the uniqueName stuff, and replace it with forcing the names to be unique, and generate the names based on the eqn/solver objects themselves
SolverBase.structForType = {}


-- whether to use separate linked binaries.  would it save on compile time?
-- this does work, however I don't have caching for linked libraries set up yet
-- TODO with the new code module system, I'm removing this for now
--  not sure how I will design cl->obj->bin, whether it will be 1-1 with the modules (which are numerous) or 1-1 with the bins (which are few)
--[[
SolverBase.useCLLinkLibraries = false
-- for reference on how it was being used:

	if self.useCLLinkLibraries then
		timer('compiling common program', function()
			self.commonUnlinkedObj = self.Program{name='common', code=commonCode}
			self.commonUnlinkedObj:compile{dontLink=true}
		end)
		timer('linking common program', function()
			self.commonProgramObj = self.Program{
				programs = {
					self.mathUnlinkedObj,
					self.commonUnlinkedObj,
				},
			}
		end)
	else
		timer('building common program', function()
			self.commonProgramObj = self.Program{name='common', code=commonCode}
			self.commonProgramObj:compile()
		end)
	end

	...

	if not SolverBase.useCLLinkLibraries then return end
	if self.mathUnlinkedObj then return end
	-- build math cl binary obj
	timer('compiling math program', function()
		self.mathUnlinkedObj = self.Program{
			name = 'math',
			code = self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():unpack()),
		}
		self.mathUnlinkedObj:compile{
			dontLink = true,
			buildOptions = '-create-library',
		}
	end)

--]]

-- whether to check for NaNs
SolverBase.checkNaNs = cmdline.checknans or false

SolverBase.showFPS = cmdline.showfps or false

-- TODO this is made to be compile-time, but i don't support that right now (i think?)
-- enable for us to let the user accum/not.  requires an extra buffer allocation
SolverBase.allowAccum = true
SolverBase.displayVarAccumFunc = cmdline.accum or false

local solverUniqueIndex = 1

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

TODO this should be 'final' i.e. no child inherits it
so that SolverBase:init can run stuff after all child classes have initialized

TODO TODO put this in its own function: ":createIdentSerializationTable" or something and call it from :init
* have each class build it *only* from the same values that the class accepts as ctor arguments
* make sure the objects it gets are all either tables or strings
* no table circular references (tolua can enforce this)
* if you do come across an object, just call its own :createIdentSerializationTable() or whatever it's called, and add that into the str ser table
* if you add an extra arg in then you'll get two cache entries for the same solver.
* if you leave an arg out then you'll get recompilation for the same solver
* otehrwise that should restore cl bin cache functionality
--]]
function SolverBase:init(args)


-- [[ save this for later
	-- right now this is only used for serialization of the config of the solver,
	-- which used to be used for caching binaries for fast compiling
	-- but now that the code/obj names have their mem locs in them, the serialize str is dif every time
	-- so i need to fix that
	-- adn i need to fix this
	self.initArgsForSerialization = table(args)
	-- remove/replace object references
	self.initArgsForSerialization.app = nil
	self.initArgsForSerialization.solver = getmetatable(self).name	-- not in initArgsForSerialization but is unique
	if self.initArgsForSerialization.subsolverClass then
		self.initArgsForSerialization.subsolverClass = self.initArgsForSerialization.subsolverClass.name
	end
	-- remove runtime variables
	self.initArgsForSerialization.fixedDT = nil
	self.initArgsForSerialization.cfl = nil
	self.initArgsForSerialization.mins = nil
	self.initArgsForSerialization.maxs = nil
	self.initArgsForSerialization.gridSize = nil
	self.initArgsForSerialization.cmds = nil		-- in choppedup
	self.initArgsForSerialization.device = nil		-- in choppedup
	self.initArgsForSerialization.id = nil			-- in choppedup
-- also include # devices.  since, on the nvidia cluster, the binaries compiled for 1 device will segfault if they are loaded when using >1 device
	self.initArgsForSerialization.numDevices = #args.app.env.devices
	self.initArgsForSerialization.app_real = args.app.real
	-- hmm, this is defeating the whole purpose of this, ...
	-- but thanks to some function serialization in BoundaryFixed ... i'm getting rid of boundary as well
	self.initArgsForSerialization.mesh = nil
--]]

	timer('SolverBase:init()', function()
		require 'hydro.code.symbols'(self, self:getSymbolFields())	-- make unique symbols
		self:initMeshVars(args)
		self:initCLDomainVars(args)
		self:initObjs(args)

		self:initCodeModules()
		self:initCodeModuleDisplay()
		self:initCDefs()
		self:postInit()	--> refreshGridSize -> createBuffers
	end)
end

-- any module is going to need a symbol
-- another TODO could be to keep track of the modules added by an object, and then just change the names after-the-fact
-- then these symbols wouldn't need to be listed separately here ...
function SolverBase:getSymbolFields()
	return table{
		-- placeholder, used by solver
		-- it turns out all module names need to be unique in order to run more than one solver at a time.
		-- makes me think we should move all these symbols/generation of them into solver
		-- TODO put these' symbol generation in solver?
		'solver_macros',
		'solver_displayCode',
		'fluxLimiter',
		'range_t',
	}
end

-- an identifying string
-- used for filesystem
-- I could hash it and use that, for shorter names, but I'd rather things be identifiable
-- TODO: make use of this in tests/test-order" etc when saving files according to their parameters
function SolverBase:getIdent()
	if self.ident then return self.ident end

--[[ the original way
	if not self.uniqueIndex then
		self.uniqueIndex = solverUniqueIndex
		solverUniqueIndex = solverUniqueIndex + 1
	end
	self.ident = tostring(self.uniqueIndex)
--]]
--[[ TODO derive this from the solver's state
	self.ident = tolua{
		solver = getmetatable(self).name,	-- TODO ensure this matches the require('hydro/solver/$name')
		eqn = self.eqn.name,				-- TODO ensure this matches the require('hydro/eqn/$name')
		--eqnArgs = {}, -- TODO params here. anything beyond the default.
		coord = self.coord.name,			-- this already matches.
	}
--]]
-- [[ until then, trust the init args to be serialized already (and not objects ... which can be accepted as init args in some cases)
	--return tolua(self.initArgsForSerialization, {indent=false})
	-- here's how tests/test-order does it:
	-- (this is causing errors, because i think these filenames are coming out too long, because they include stuff that test-order didn't
	-- 	such as boundary condition info, eqn args, etc ... these need to be shortened)
	assert(self.initArgsForSerialization)
-- [=[ make sure there's no objects that are cdata, make sure there's no loops
	do
		local allTables = {}
		local function checkForLoops(o)
			for k,v in pairs(o) do
				if type(v) == 'cdata' then
					error("initArgsForSerialization field "..k.." is cdata")
				end
				if type(k) == 'cdata' then
					error("initArgsForSerialization key "..k.." is cdata ... you are really in trouble")
				end
				if allTables[v] then
					error("found a loop in the initArgsForSerialization ... field "..k.." ... try (a) replacing the field with an equally-unique string/table, or (b) removing loops")
				end
				if type(v) == 'table' then
					checkForLoops(v)
				end
			end
		end
		assert(type(self.initArgsForSerialization) == 'table')
		allTables[self.initArgsForSerialization] = true
		checkForLoops(self.initArgsForSerialization)
	end
--]=]
	local serStr = tolua(self.initArgsForSerialization)
	local destName = serStr:match('^{(.*)}$')
	assert(destName, "we must have got a circular reference in:\n"..serStr)
	destName = destName
		:gsub('%s+', ' ')
		:gsub('"', '')
	destName = string.trim(destName)
	--
	local configStr = destName
		:gsub('/', '')
		:gsub('{', '(')
		:gsub('}', ')')

--	self.ident = configStr
-- [=[
-- hmm, seems the filenames are too long
-- how about hashing them?
	local bit = require 'bit'
	local len = 8
	local chs = range(len):mapi(function() return 0 end)
	for i=1,#configStr do
		local j = (i-1)%len+1
		chs[j] = bit.bxor(chs[j], configStr:sub(i,i):byte())
	end
	self.ident = chs
		:sub(1, math.min(len, #configStr))
		:mapi(function(ch)
			--return string.char(ch)
			return ('%02x'):format(ch)
		end):concat()
	if self.app.verbose then
		print('ident config str: '..configStr)
	end
	print("writing files to cache/"..self.ident)
--]=]
--]]


	-- debugging?
	local dir = 'cache/'..self.ident
	path(dir):mkdir(true)
	path(dir..'/config'):write(configStr)

	return self.ident
end

function SolverBase:initMeshVars(args)
	assert(args, "expected named parameter table")
	self.app = assert(args.app, "expected app")
	self.dim = assert(args.dim, "expected dim")

-- [[ do this before creating coord so it can populate them optionally
-- (though in theory / philosophy of design, why would coord populate ops?)
	-- operators for this solver
	self.ops = table()

	-- struct for the solver
	self.solverStructFields = table()
--]]

	-- MeshSolver needs coord created early
	--  so that it can allocate cell_t's which are defined in coord
	--  (come to think of it, GridSolver allocates cellBuf also)
	-- so cell_t has to be finished before createBuffers.
	self:createCoord(args)


	self.device = args.device or self.app.env.devices[1]
	self.cmds = args.cmds or self.app.env.cmds[1]

	self.color = vec3d(math.random(), math.random(), math.random()):normalize()


	self.solverStructFields:append{
	-- [[ right now the mesh initial conditions use these, but otherwise they can be GridSolver-specific
		{name='mins', type='real3'},
		{name='maxs', type='real3'},
	--]]
	-- [[ the mins/maxs, or the super-solver's mins/maxs.  only needed because of the composite solvers.
		{name='initCondMins', type='real3'},
		{name='initCondMaxs', type='real3'},
	--]]
		{name='t', type='real'},
		{name='dt', type='real'},
	}

	-- in GridSolver this was 'initMeshVars' which comes first
	-- in MeshSolver this was 'preInit' which comes later

	local solver = self

	-- my kernel objs are going to need workgroup info based on domain.size-2*noGhost as well as domain.size ...
	-- ...and rather than require an extra argument, I think I'll just take advantage of a closure
	local Program = require 'cl.obj.program':subclass()

	-- I'm tempted to merge the two cl oop libs...
	require 'cl.program'.showCodeOnError = false

	function Program:init(args)
		self.name = args.name
		args.env = solver.app.env
		args.domain = solver.domain

		local cldir = 'cache/'..solver:getIdent()..'/src'
		local bindir = 'cache/'..solver:getIdent()..'/bin'
		path(cldir):mkdir(true)
		path(bindir):mkdir(true)

		--[[
		https://github.com/KhronosGroup/SPIR/tree/spirv-1.1
		how to compile 32-bit SPIR-V:
			clang -cc1 -emit-spirv -triple <triple> -cl-std=c++ -I <libclcxx dir> -x cl -o <output> <input> 				#For OpenCL C++
			clang -cc1 -emit-spirv -triple <triple> -cl-std=<CLversion> -include opencl.h -x cl -o <output> <input> 		#For OpenCL C
		how to compile 64-bit SPIR-V:
			clang -cc1 -emit-spirv -triple=spir-unknown-unknown -cl-std=c++ -I include kernel.cl -o kernel.spv 				#For OpenCL C++
			clang -cc1 -emit-spirv -triple=spir-unknown-unknown -cl-std=CL2.0 -include opencl.h kernel.cl -o kernel.spv 	#For OpenCL C
		--]]
		if useClang then
			assert(args.name, "clang needs a program to have a name")
			local cldir = 'cache/'..solver:getIdent()..'/src'
			local bcdir = 'cache/'..solver:getIdent()..'/bc'
			local spvdir = 'cache/'..solver:getIdent()..'/spv'
			path(cldir):mkdir(true)
			path(bcdir):mkdir(true)
			path(spvdir):mkdir(true)
			-- ok I don't want super to hit the args.code condition
			-- but I want to store the code for later, so
			self.code = args.code
			args.code = nil
			-- TODO hmm this dependency graph for code->cl->bin is hardcoded into cl/obj/program.lua
			--  but could I generalize this somehow, and integrate it into lua-make ...
			-- and then make a separate build-target chain for .clcpp -> .bc -> .spv ...
			self.cacheFileCL = cldir..'/'..args.name..'.cl'
			self.cacheFileBC = bcdir..'/'..args.name..'.bc'
			self.cacheFileSPV = spvdir..'/'..args.name..'.spv'
			Program.super.init(self, args)
			return
		end

		local clfn = cldir..'/'..args.name..'.cl'

		-- caching binaries, which doesn't write unless the program successfully compiles
		if not cmdline.usecachedcode
		and not useClang
		then
			if args.name then
				if useCache then
					args.cacheFileCL = clfn
					args.cacheFileBin = bindir..'/'..args.name..'.bin'
				end
			end
			Program.super.init(self, args)
			return
		end

		-- Write generated code the first time.  Subsequent times use the pre-existing code.  Useful for debugging things in the generated OpenCL.
		if path(clfn):exists() then
			local cachedCode = path(clfn):read()
			assert(cachedCode:sub(1,#args.env.code) == args.env.code, "seems you have changed the cl env code")
			args.code = cachedCode:sub(#args.env.code+1)	-- because the program will prepend env.code ... hmm, this could be done a better way.
			Program.super.init(self, args)
			return
		end

		Program.super.init(self, args)	-- do this so getCode() works
		path(clfn):write(self:getCode())
	end

	function Program:compile(args)
		if useClang then
			-- if cl file is out of date then regen bytecode

			local oldCode
			if path(self.cacheFileCL):exists() then
				oldCode = path(self.cacheFileCL):read()
			end
			local newCode = self:getCode()
			if oldCode ~= newCode then
				-- only write the new code if it is outdated -- to preserve file timestamps -- so the build system doesn't rebuild needlessly
				if oldCode then path(self.cacheFileCL..'.old'):write(oldCode) end
				path(self.cacheFileCL):write(newCode)
			end
			-- so cl.obj.program :compile using .code and .cacheFile basically does the same thing, but less flexible
			local exec = require 'make.exec'
			require 'make.targets'{
				verbose = true,
				{
					srcs = {self.cacheFileCL},
					dsts = {self.cacheFileBC},
					rule = function()
						exec(table{
							'clang',
							'-v',
							'-Xclang -finclude-default-header',
							'--target=spir64-unknown-unknown',
							'-emit-llvm',
							'-c',
							'-o', ('%q'):format(self.cacheFileBC),
							('%q'):format(self.cacheFileCL)
						}:concat' ')
					end,
				}, {
					srcs = {self.cacheFileBC},
					dsts = {self.cacheFileSPV},
					rule = function()
						exec(table{
							'llvm-spirv',
							('%q'):format(self.cacheFileBC),
							'-o', ('%q'):format(self.cacheFileSPV),
						}:concat' ')
					end,
				},
			}:run(self.cacheFileSPV)

			self.IL = path(self.cacheFileSPV):read()
			assert(self.IL, "failed to read the IL code")

			args = table(args):setmetatable(nil)
			args.verbose = solver.app.verbose
			local results = Program.super.compile(self, args)
			assert(self.obj, "there must have been an error in your error handler")	-- otherwise it would have thrown an error
			do--if self.obj then	-- did compile
				print((self.name and self.name..' ' or '')..'log:')
				-- TODO log per device ...
				print(string.trim(self.obj:getLog(solver.device)))
			end
			return results
		else
			args = args or {}
			args.verbose = solver.app.verbose
			local opts = table()
			opts:insert'-w'	-- show warnings
			--opts:insert'-cl-std=CLC++'	-- I don't have  cl_ext_cxx_for_opencl  so ... I can only do this with IL
			args.buildOptions = opts:concat' '
			local results = Program.super.compile(self, args)
			assert(self.obj, "there must have been an error in your error handler")	-- otherwise it would have thrown an error
			do--if self.obj then	-- did compile
				print((self.name and self.name..' ' or '')..'log:')
				-- TODO log per device ...
				print(string.trim(self.obj:getLog(solver.device)))
			end
			-- if we are using cached code then manually write binaries
			if cmdline.usecachedcode and useCache then
				local binfn = 'cache/'..solver:getIdent()..'/bin/'..self.name..'.bin'
				path(binfn):write(tolua(self.obj:getBinaries()))
			end
			return results
		end
	end

	self.Program = Program


	if self.app.targetSystem ~= 'console' then
		local GLProgram = require 'gl.program':subclass()
		function GLProgram:init(...)
			local args = ...
			args.version = require 'hydro.draw.draw'.glslVersion:match'^#version (.*)$'
			args.precision = 'best'

			local dir = 'cache/'..solver:getIdent()..'/shader'
			print('building '..dir..'/'..args.name..'.vert.glsl & .frag.glsl')

			path(dir):mkdir(true)
			local pathstr = dir..'/'..args.name
			-- Write generated code
			path(pathstr..'.vert.glsl'):write(args.vertexCode)
			path(pathstr..'.frag.glsl'):write(args.fragmentCode)

			GLProgram.super.init(self, ...)
		end
		self.GLProgram = GLProgram
	end
end

function SolverBase:createCoord(args)
	-- not sure where to put this, but it must precede MeshSolver:initMeshVars
	-- also I'm pretty sure CoordinateSystem:init() doesn't require anything from SolverBase, only the later member functions
	if require 'hydro.coord.coord':isa(args.coord) then
		self.coord = args.coord	-- ptr copy expected by AMR
		self.coord.solver = self
	else
		self.coord = require('hydro.coord.'..(args.coord or 'cartesian'))(table({solver=self}, args.coordArgs))
	end
end

-- this is called back from within coord creation *BEFORE* it is finished (and assigned)
-- hmm, here is a dilemma, eqn has not yet been created, but I would like eqn to be able to populate cells in this.
-- so that means ... delay this until after eqn has been created (which is also after solver.coord has been assigned)
-- and also delay use of any cellBuf's until after that point.
function SolverBase:createCellStruct(coord)

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

	self.initCondIndex = table.find(self.eqn.initCondNames, args.initCond)
	if not self.initCondIndex then
		assert(#self.eqn.initCondNames > 0, "couldn't find initCond "..('%q'):format(tostring(args.initCond)).." ... and there are no initConds to search.")
		self.initCondIndex = 1
		-- TODO do this fallback before storing the cache string
		-- this way the cached code / binaries don't get messed up if we are loading extra init conds and they are not present
		print("!!!!! couldn't find initCond "..('%q'):format(tostring(args.initCond)).." so falling back on initCond "..('%q'):format(tostring(self.eqn.initCondNames[self.initCondIndex])).." !!!!!")
	end
	self.initCondArgs = args.initCondArgs

	self.integratorArgs = args.integratorArgs
	self.integratorIndex = integratorNames:find(args.integrator)
	if not self.integratorIndex then
		assert(#integratorNames > 0, "couldn't find integrator "..('%q'):format(tostring(args.integrator)).." ... and there are no integrators to search.")
		self.integratorIndex = 1
		print("!!!!! couldn't find integrator "..('%q'):format(tostring(args.integrator)).." so falling back on integrator "..('%q'):format(tostring(integratorNames[self.integratorIndex])).." !!!!!")
		-- TODO shouldn't this be an error? if the explicitly-asked-for integrator isn't there...
	end

	self.checkNaNs = self.checkNaNs or args.checkNaNs
	self.useFixedDT = not not args.fixedDT
	self.fixedDT = args.fixedDT or self.fixedDT or .001
	self.cfl = args.cfl or .5	--/self.dim
	self.fluxLimiter = self.app.limiterNames:find(args.fluxLimiter)
	if not self.fluxLimiter then
		assert(#self.app.limiterNames > 0, "couldn't find fluxLimiter "..('%q'):format(tostring(args.fluxLimiter)).." ... and there are no limiters to search.")
		self.fluxLimiter = 1
		print("!!!! couldn't find fluxLimiter "..('%q'):format(tostring(args.fluxLimiter)).." so falling back on limiter "..('%q'):format(tostring(self.app.limiterNames[self.fluxLimiter])).." !!!!!")
	end

	-- this influences refreshInitStateProgram()
	self.eqn:createInitState()

	-- add eqn vars to solver_t
	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			self.solverStructFields:insert{name=var.name, type=var.ctype}
		end
	end

	self:refreshGetULR()

	-- do this before any call to createBuffers
	-- make sure it's done after createEqn (for the solver_t struct to be filled out by the eqn)
	-- actually this has to go after self.eqn:createInitState,
	--  which is called in refreshEqnInitState
	--  which is called in refreshGridSize
	--  which is called in postInit
	-- the createInitState also creates kernels and runs them on the solver buffers
	-- so I probably need a separate call to eqn.initCond, earlier, which constructs the object and the guiVars, but runs no kernels
	self.solver_t = self.app:uniqueName'solver_t'
	self.solverStruct = Struct{
		name = self.solver_t,
		fields = self.solverStructFields,
		-- without 'packed' I get misalignments between ffi and opencl
		packed = true,
	}.class
	self.solverStructFields = nil
	self.solverPtr = ffi.new(self.solver_t)
	-- TODO I see room for a separation between solver_t and eqn_t
	-- but right now both structs are tied directly to solver.cl, so there's no need to separate them at the moment.

	-- not sure where this should go, but probably somewhere parallel to solverPtr
	-- initStruct is already created in eqn:createInitState
	-- TODO if initCond is supposed to be modular then this would have to be created after initCond is changed
	self.initCond_t = self.eqn.initCond.initStruct.name
	self.initCondPtr = ffi.new(self.initCond_t)



	-- [[
	-- This code should go after all objs are created -- so they can each call their own callbacks to modify this
	-- and it should go before the next call, which is initCodeModules, which requires the finalizeCellStruct to be called

	-- this creates coord.cellStruct and coord.faceStruct, but doesn't yet create the types until coord:finalizeCellStruct()
	self.coord:createCellStruct()

	-- moving the solver-specific changes to cell_t into this function:
	-- especially MeshSolver uses this to add mesh fields to cell_t
	if self.createCellStruct then
		self:createCellStruct()
	end

	if self.eqn.createCellStruct then
		self.eqn:createCellStruct()
	end

	for _,op in ipairs(self.ops) do
		if op.createCellStruct then
			op:createCellStruct()
		end
	end

	self.coord:finalizeCellStruct()
	--]]

	-- [[ init boundary stuff
	self.boundaryOptions = table()
	self.boundaryOptionNames = table()
	self.boundaryOptionForName = {}

	-- This is either in GridSolver of MeshSolver
	-- each has dif Boundary objs, one for structured grids, another for the mesh/face/cell structs
	self:createBoundaryOptions()

	if self.eqn.createBoundaryOptions then
		self.eqn:createBoundaryOptions()
	end
	--]]
end

function SolverBase:addBoundaryOptions(args)
	for _,arg in ipairs(args) do
		self:addBoundaryOption(arg)
	end
end

function SolverBase:addBoundaryOption(boundaryClass)
	self.boundaryOptions:insert(assert(boundaryClass))
	self.boundaryOptionNames:insert(assert(boundaryClass.name))
	self.boundaryOptionForName[boundaryClass.name] = boundaryClass
end


-- collect *all* modules from all sub-objects
-- do this after objs are created and before codegen
-- first the global ones from math and app, then coord, initCond, boundary, eqn, solver ...
-- then with all of them, specify which ones to target (for .h and .cl code) and they will trim the others
function SolverBase:initCodeModules()
	--self.modules = require 'modules'(self.app.modules)
-- TODO just change the references
self.modules = self.app.modules

	-- what to compile?
	-- use keys of this
	-- solverModulesEnabled = modules for solvers
	-- initModulesEnabled = modules for init
	-- sharedModulesEnabled = modules for both.  previously 'codeprefix'.
	self.initModulesEnabled = table()
	self.solverModulesEnabled = table()
	self.sharedModulesEnabled = table()

	-- hmm, if solverStruct needs to be cdef'd immediately ...
	-- I get the impression that, previously, ffi.cdef was pulling nothing from solverStruct
	-- while opencl was getting the proper struct
	self.modules:add{
		name = assert(self.solver_t),
		-- I could do like in hydro/app, replace the .structs = ... with .typecode = function() ... end
		-- and then make that function conditional on the buildingOpenCL flag
		-- but then I'd be omitting the structs from the checkStructSizes() alignment-test
		-- and solver_t here is one of the first to fall out of alignment
		-- so it is probably better to just add/remove the struct object from the module's .structs
		-- TODO consider doing the same as this with the other .buildingOpenCL modules, and maybe make a function for switching the contents back and forth ...
		-- ... so I can include those in checkStructSizes as well ...
		structs = {self.solverStruct},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.solver_t..' solver_t;',
	}

	-- header
	self.modules:add{
		name = self.symbols.solver_macros,
		headercode = table{
			self.dim == 3 and '#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable' or '',
			'#ifndef M_PI',
			'#define M_PI '..('%.50f'):format(math.pi),
			'#endif',
			'#define dim '..self.dim,
			'#define numStates '..self.eqn.numStates,
			'#define numIntStates '..self.eqn.numIntStates,
			'#define numWaves '..self.eqn.numWaves,
		}:concat'\n',
	}

	self.coord:initCodeModules()
	self.eqn:initCodeModules()	-- calls eqn.initCond:initCodeModules()

	-------- solver modules --------

	self.modules:add{
		name = self.symbols.fluxLimiter,
		code = self.eqn:template[[
real <?=fluxLimiter?>(real r) {
	<?=solver.app.limiters[solver.fluxLimiter].code?>
}
]],
	}

	-- this could go into app's modules ...
	-- I'll keep it here, because as I make solvers more flexible, modules will be pushed solver-ward
	-- and in the case that i want to run mixed float/double solvers, this will be solver-specific
	self.modules:add{
		name = self.symbols.range_t,
		-- TODO use struct so I can verify cl/ffi struct alignment
		typecode = 'typedef struct range_t { real min, max; } range_t;',
	}

	-- when building modules for ops, only add them to solverModulesEnabled, not sharedModulesEnabled, since init doesn't need them
	for _,op in ipairs(self.ops) do
		if op.initCodeModules then
			op:initCodeModules()
		end
	end
end

-- Calling :getDisplayCode() will query other modules (esp type info) for what to produce
--  so add this last.
-- I can get around this if I go back to the previous functionality of only testing what dependencies are listed in the math module
--  but this won't work well with moving any math dependencies closer to the kernel functions
function SolverBase:initCodeModuleDisplay()
	-- depends on self.eqn
	self:createDisplayVars()

	-- this depends on :createDisplayVars()
	self.modules:addFromMarkup(
		'//// MODULE_NAME: '..self.symbols.solver_displayCode..'\n'
		..self:getDisplayCode()			-- this call constructs displayValue_t
	)
	self.solverModulesEnabled[self.symbols.solver_displayCode] = true
end

-- TODO if you want to define ffi ctype metatable then put them all in one spot here
-- TODO TODO since i'm switching to the modules,
-- and since eqn's init's ffi calls need the ctypes of the fields defined up-front in order to use them (like counting scalars in cons_t),
-- how about instead I call cdef() immediately upon request?
-- and then just get of this function completely.
function SolverBase:initCDefs()
	local moduleNames = table(
		self.initModulesEnabled,
		self.solverModulesEnabled,
		self.sharedModulesEnabled,
		{
			[self.eqn.symbols.cons_t] = true,
			[self.eqn.symbols.prim_t] = true,
		}
	):keys()
	if self.app.verbose then
		print("ffi.cdef'ing: "..moduleNames:concat', ')
	end

	-- opencl needs solver_t
	-- but ffi.cdef doesn't (cuz it was already cdef'd)
	-- (and I can't defer cdef, cuz I already need it asap for solverPtr)
	self.modules.set[self.solver_t].structs:remove()
	-- ... same with initCond_t?
	-- how come euler doesn't mind initCond_t being removed, but einstein-fd requires it to be removed?
	self.modules.set[self.initCond_t].structs:remove()

	require 'hydro.code.safecdef'(
		self.modules:getTypeHeader(moduleNames:unpack())
	)
	-- ... and re-add
	self.modules.set[self.solver_t].structs:insert(self.solverStruct)
	self.modules.set[self.initCond_t].structs:insert(self.eqn.initCond.initStruct)
end

function SolverBase:refreshGetULR()
	self.getULRBufType = self.eqn.symbols.cons_t
	self.getULRBufName = 'UBuf'
	self.getULRArg = self.getULRBufType..'* '..self.getULRBufName
	self.getULRCode = function(self, args)
		args = args or {}
		local suffix = args.suffix or ''
		return self.eqn:template([[
global <?=cons_t?> const * const UL<?=suffix?> = <?=bufName?> + <?=indexL?>;
global <?=cons_t?> const * const UR<?=suffix?> = <?=bufName?> + <?=indexR?>;
]],		{
			suffix = suffix,
			indexL = args.indexL or 'indexL'..suffix,
			indexR = args.indexR or 'indexR'..suffix,
			bufName = args.bufName or self.getULRBufName,	-- for displayVars the variable name is 'buf', so I need to override it either in displayCode or here
		})
	end
end

function SolverBase:postInit()
	timer('SolverBase:postInit()', function()
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

--[[
TODO NVIDIA driver, upon multiple writes, is - upon clEnqueueWriteBuffer, giving me CL_OUT_OF_RESOURCES.
This is documented nowhere (not in the official opencl docs),
however in some forums people complain of a similar thing with nvidia under clEnqueueReadBuffer.
	https://www.nvidia.com/en-us/geforce/forums/geforce-graphics-cards/5/229335/opencl-out-of-resources-error/
	https://forums.developer.nvidia.com/t/cl-out-of-resources-error-on-clenqueuereadbuffer-driver-crashes-and-i-get-an-cl-out-of-resources-err/16426
	https://stackoverflow.com/questions/34943647/c-opencl-return-cl-out-of-resources
To supplicant a similar clEnqueueReadBuffer problem, the solution is "just wait a few seconds for NVIDIA to settle its jimmies.  I tried waiting 30 seconds, no difference, so it is a state problem.
Another cause of this error on NVIDIA seems to be writing OOB ... which shouldn't be the case here, solverPtr, solverBuf, and solver_t are all the same size.
It seems to always happen in the same place, and only after several writes to the same buffer.
So my attempted fix: to minimize all writes to the solverBuf until just before the simulation starts.

Turns out that's not the problem.
The problem happens upon multiple writes ... but even only if ever writing to the same device.
If slurm allocates 2 GPUs, then OpenCL creates a context with 2 devices, and we only use cmd-queue #1 associated with device #1 ... we will get this mystery error.
But if slurm allocates 2 GPUs, then OpenCL creates a context with just 1 device, and we only use cmd-queue #1 associated with device #1 ... then things work fine.

Also don't forget: only test multi-gpu with the choppedup solver, which creates sub-solvers, each with a unique device for each.
But does choppedup itself ever allocate any memory?  I'm almost thinking no it doesn't.
In that case ... each sub-solver has its own cmd-queue, associated with a unique (and ideally separate) device ...
Next question: Should we create a unique context for each sub-solver?
	Will that solve our problem? (it should, from what I've seen so far)
	Will that just introduce more complexities with having to copy GPU data between different contexts? (in the GL world, copying between contexts is a headache)
--]]
function SolverBase:refreshSolverBuf()
--[[ bssnok-fd-senr is crashing here even though all the sizes seem to match up
print('self.solverPtr', self.solverPtr)
	assert(not rawequal(self.solverPtr, nil))
print('ffi.sizeof(self.solverPtr)', ffi.sizeof(self.solverPtr))
print('self.solverBuf.type', self.solverBuf.type)
print('ffi.sizeof(self.solverBuf.type)', ffi.sizeof(self.solverBuf.type))
print('self.solverBuf.count', self.solverBuf.count)
print('ffi.sizeof(self.solverBuf.type) * self.solverBuf.count', ffi.sizeof(self.solverBuf.type) * self.solverBuf.count)
	assert(ffi.sizeof(self.solverBuf.type) * self.solverBuf.count == ffi.sizeof(self.solverPtr))
--]]
	self.solverBuf:fromCPU(self.solverPtr)
end

function SolverBase:refreshInitCondBuf()
	self.initCondBuf:fromCPU(self.initCondPtr)
end

function SolverBase:copyGuiVarsToBufs()
	-- copy solverBuf vars
	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			if var.ctype == 'real' then
				self.solverPtr[var.name] = toreal(var.value)
			else
				self.solverPtr[var.name] = var.value
			end
		end
	end
	self:refreshSolverBuf()		-- solver=bssnok-fd, eqn=bssnok-fd-senr ... crashes here - can't access GPU memory

	-- copy initCondBuf vars
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
	timer('SolverBase:refreshGridSize()', function()
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

	--local moduleNames = self.sharedModulesEnabled:keys()
	-- what code does common use?
	-- in fact, what types does it use?
	-- it only seems to use solver_t and cons_t
	local moduleNames = table{
		'realparam',
		self.solver_t,
		self.eqn.symbols.cons_t,

		-- This is in GridSolver, a subclass.
		-- In fact, all the display stuff is pretty specific to cartesian grids.
		-- Not 100% though, since the MeshSolver stuff was working with it before I introduced the code module stuff.
		self.symbols.SETBOUNDS_NOGHOST,

		self.symbols.SETBOUNDS,
	}
	if self.app.verbose then
		print('common modules: '..moduleNames:concat', ')
	end

	self.app.buildingOpenCL = true
	local commonCode = table{
		-- just header, no function calls needed
		self.modules:getHeader(moduleNames:unpack()),
		-- and i need this. .. but I need its code ... but if it requests any of 'moduleNames' headers , then i'll have duplicates ... so ... ?
		self.modules:getCodeAndHeader'sqr',
		-- and here's the code
		self.eqn:template[[
kernel void multAddInto(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const a,
	global <?=cons_t?> const * const b,
	realparam const c
) {
	<?=SETBOUNDS_NOGHOST?>();
<?
for i=0,eqn.numIntStates-1 do
?>	a[index].ptr[<?=i?>] += b[index].ptr[<?=i?>] * c;
<?
end
?>}

kernel void multAdd(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const a,
	global <?=cons_t?> const * const b,
	global <?=cons_t?> const * const c,
	realparam const d
) {
	<?=SETBOUNDS_NOGHOST?>();
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

// cons x cons -> cons
// a -= b
// TODO cut out ghost cells?
kernel void subtractIntoCons(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const a,
	global <?=cons_t?> const * const b
) {
	<?=SETBOUNDS_NOGHOST?>();
<? for i=0,eqn.numIntStates-1 do
?>	a[index].ptr[<?=i?>] -= b[index].ptr[<?=i?>];
<? end
?>
}

// whereas the multAddInto and multAdd operate on cons_t,
// this operates on a buffer of reals.
// and this should only operate on non-ghost cells.
// This is used on reduceBuf / displayVars to calc their stddev.
kernel void subtractAndSquare(
	constant <?=solver_t?> const * const solver,
	global real * const a,
	realparam const mu
) {
	<?=SETBOUNDS_NOGHOST?>();
	global real * const ai = a + index;
	*ai -= mu;
	*ai *= *ai;
}

// cons -> real
// a = b^2
// TODO cut out ghost cells
kernel void squareCons(
	constant <?=solver_t?> const * const solver,
	global real * const a,
	global <?=cons_t?> const * const b
) {
	<?=SETBOUNDS_NOGHOST?>();
	real sum = 0.;
<? for i=0,eqn.numIntStates-1 do
?>	sum += sqr(b[index].ptr[<?=i?>]);
<? end
?>
	a[index] = sum;
}

<? if solver.checkNaNs then ?>
kernel void findNaNs(
	constant <?=solver_t?> const * const solver,
	global real * const dst,
	global <?=cons_t?> const * const src
) {
<? if solver.checkNaNs == 'gpu-noghost' then ?>
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
<? else ?>
	<?=SETBOUNDS?>(0, 0);
<? end ?>

	dst[index] = 0;
<? for i=0,eqn.numStates-1 do
?>	dst[index] += (real)(!isfinite(src[index].ptr[<?=i?>]));
<? end ?>
}
<? end ?>

]]
	}:concat'\n'
	self.app.buildingOpenCL = false

	timer('building program cache/'..self:getIdent()..'/src/common.cl ', function()
		self.commonProgramObj = self.Program{name='common', code=commonCode}
		self.commonProgramObj:compile()
	end)

	-- used by the integrators
	-- needs the same globalSize and localSize as the typical simulation kernels
	-- TODO exclude states which are not supposed to be integrated
	self.multAddKernelObj = self.commonProgramObj:kernel{name='multAdd', domain=self.domainWithoutBorder}
	self.multAddKernelObj.obj:setArg(0, self.solverBuf)

	self.multAddIntoKernelObj = self.commonProgramObj:kernel{name='multAddInto', domain=self.domainWithoutBorder}
	self.multAddIntoKernelObj.obj:setArg(0, self.solverBuf)

	self.subtractAndSquareKernelObj = self.commonProgramObj:kernel{name='subtractAndSquare', domain=self.domainWithoutBorder}
	self.subtractAndSquareKernelObj.obj:setArg(0, self.solverBuf)

	self.subtractIntoConsKernelObj = self.commonProgramObj:kernel{name='subtractIntoCons', domain=self.domainWithoutBorder}
	self.subtractIntoConsKernelObj.obj:setArgs(self.solverBuf)

	self.squareConsKernelObj = self.commonProgramObj:kernel{name='squareCons', domain=self.domainWithoutBorder}
	self.squareConsKernelObj.obj:setArgs(self.solverBuf, self.reduceBuf)

	if self.checkNaNs then
		self.findNaNsKernelObj = self.commonProgramObj:kernel('findNaNs', self.solverBuf, self.reduceBuf)
	end

	-- TODO vectors won't reduce anymore unless the reduceMin is constructed with 3* the # of count
	-- but this breaks reduce for scalars (it includes those extra zeros)
	-- which makes calcDT fail
	self.reduceMin = self.app.env:reduce{
		secondPassInCPU = cmdline.secondPassInCPU,
		count = self.numCells,
		op = function(x,y) return 'min('..x..', '..y..')' end,
		initValue = 'INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceMax = self.app.env:reduce{
		secondPassInCPU = cmdline.secondPassInCPU,
		count = self.numCells,
		op = function(x,y) return 'max('..x..', '..y..')' end,
		initValue = '-INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceSum = self.app.env:reduce{
		secondPassInCPU = cmdline.secondPassInCPU,
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
		if self.app.verbose then
			print(tolua(info))
		end
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
	if ffi.sizeof(self.eqn.symbols.cons_t) ~= self.eqn.numStates * realSize then
	   error('Expected sizeof('..self.eqn.symbols.cons_t..') to be '
		   ..self.eqn.numStates..' * sizeof(real) = '..(self.eqn.numStates * realSize)
		   ..' but found '..ffi.sizeof(self.eqn.symbols.cons_t)..' = '..(ffi.sizeof(self.eqn.symbols.cons_t) / realSize)..' * sizeof(real). '
		   ..'Maybe you need to update Eqn.numStates?')
	end
	if ffi.sizeof(self.eqn.symbols.cons_t) < ffi.sizeof(self.eqn.symbols.prim_t) then
		print("sizeof(cons_t) =", ffi.sizeof(self.eqn.symbols.cons_t))
		print("sizeof(prim_t) =", ffi.sizeof(self.eqn.symbols.prim_t))
		error("for PLM's sake I might need sizeof(prim_t) <= sizeof(cons_t)")
	end

	-- should I put these all in one AoS?
	-- or finally make use of constant args ...

	-- this much for the first level U data
	self:clalloc('UBuf', self.eqn.symbols.cons_t, self.numCells)

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', app.real, self.numCells * 3)
	self:clalloc('reduceSwapBuf', app.real, math.ceil(self.numCells * 3 / self.localSize1d))
	self.reduceResultPtr = ffi.new'real[1]'
	self.reduceResultPtr[0] = toreal(0)

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
			--[[ on AMD:
			0 gives CL_INVALID_VALUE
			CL_MEM_READ_ONLY gives CL_MEM_OBJECT_ALLOCATION_FAILURE
			CL_MEM_WRITE_ONLY gives CL_MEM_OBJECT_ALLOCATION_FAILURE
			CL_MEM_READ_WRITE gives CL_MEM_OBJECT_ALLOCATION_FAILURE
			so ... ?
			--]]
			self.texCLMem = CLImageGL{context=app.ctx, tex=self.tex, write=true}
		else
			-- use texSize:volume() so the glTexSubImage can use the whole buffer, in the event of meshsolver where texSize:volume can be > numCells
			if self.app.verbose then
				print('allocating gpu-cpu display copy buffer of '..app.real..'['..tonumber(self.texSize:volume() * 3)..']')
			end
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
		local mins = initCond.mins
		if type(mins) == 'function' then mins = assert(mins(initCond)) end
		self.mins = vec3d(table.unpack(mins))
		for j=1,3 do
			self.solverPtr.mins.s[j-1] = toreal(self.mins.s[j-1])
		end
		self.initCondMins = vec3d(self.mins:unpack())
	end
	if initCond.maxs then
		local maxs = initCond.maxs
		if type(maxs) == 'function' then maxs = assert(maxs(initCond)) end
		self.maxs = vec3d(table.unpack(maxs))
		for j=1,3 do
			self.solverPtr.maxs.s[j-1] = toreal(self.maxs.s[j-1])
		end
		self.initCondMaxs = vec3d(self.maxs:unpack())
	end

	-- there's a lot of overlap between this and the solverBuf creation...
	self:refreshSolverBufMinsMaxs()

	-- while we're here, write all gui vars to the solver_t
	self:copyGuiVarsToBufs()


-- [[
	-- get the symbolic function for the x, y, z
	-- and find its range for a domain of the solverMins/solverMaxs
	local coord = self.coord
	local symmath = require 'symmath'
	local var = symmath.var
	local u,v,w = table.unpack(coord.baseCoords)
	local chart = coord.chart()
	for i=1,#chart do
--print'before replace'
--print(chart[i])
		chart[i] = coord:applyReplVars(chart[i])
-- This is what applyReplVars is for, right?
-- but I wanted to keep these as code, and use CL #define's to replace their values
-- so they could be changed at runtime without regenerating the code
-- just with a kernel recompile.
-- so TODO if I change them to be based on applyReplVars instead then I don't need this extra code
-- or TODO otherwise I can just put these repls somewhere else ... ? idk
		chart[i] = coord:applyReplDefines(chart[i])
--print'after replcae'
--print(chart[i])
	end

-- ok for sphere_sinh_radial, at this point after replacement, we still have vars for AMPL and SINHW
-- and I wrote those to be macros

	local domains = range(3):mapi(function(i)
		return symmath.set.RealSubset(
			tonumber(self.mins.s[i-1]),
			tonumber(self.maxs.s[i-1]),
			true, true)
	end)

	-- find chart[1]'s range for domain
	chart = chart
		:replace(u, var('tmpU', nil, nil, domains[1]))
		:replace(v, var('tmpV', nil, nil, domains[2]))
		:replace(w, var('tmpW', nil, nil, domains[3]))

	local ranges = range(3):mapi(function(i)
		return chart[i]:getRealRange()
	end)

	self.cartesianMin = vec3d(
		ranges[1][1].start,
		ranges[2][1].start,
		ranges[3][1].start
	)
	self.cartesianMax = vec3d(
		table.last(ranges[1]).finish,
		table.last(ranges[2]).finish,
		table.last(ranges[3]).finish
	)
	if self.app.verbose then
		print('cartesianMin = '..self.cartesianMin)
		print('cartesianMin = '..self.cartesianMax)
	end

--]]
end

function SolverBase:refreshSolverBufMinsMaxs()
	for j=1,3 do
		self.solverPtr.mins.s[j-1] = toreal(self.mins.s[j-1])
		self.solverPtr.maxs.s[j-1] = toreal(self.maxs.s[j-1])
	end
	if self.app.verbose then
		print('coord min = '..fromreal(self.solverPtr.mins.x)..', '..fromreal(self.solverPtr.mins.y)..', '..fromreal(self.solverPtr.mins.z))
		print('coord max = '..fromreal(self.solverPtr.maxs.x)..', '..fromreal(self.solverPtr.maxs.y)..', '..fromreal(self.solverPtr.maxs.z))
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

-- depends on buffers
function SolverBase:refreshSolverProgram()
	local eqn = self.eqn

	-- enable here after all modules are provided
	if self:hasModule(eqn.symbols.calcDT) then
		self.solverModulesEnabled[eqn.symbols.calcDT] = true
	end
	if self:hasModule(eqn.symbols.addSource) then
		self.solverModulesEnabled[eqn.symbols.addSource] = true
	end
	if self:hasModule(eqn.symbols.constrainU) then
		self.solverModulesEnabled[eqn.symbols.constrainU] = true
	end

	local code
	timer('generating solver code', function()
		local moduleNames = table(self.sharedModulesEnabled, self.solverModulesEnabled):keys()
		if self.app.verbose then
			print('solver modules: '..moduleNames:concat', ')
		end

		self.app.buildingOpenCL = true
		code = self.modules:getCodeAndHeader(moduleNames:unpack())
		self.app.buildingOpenCL = false
	end)

	timer('building program cache/'..self:getIdent()..'/src/solver.cl ', function()
		self.solverProgramObj = self.Program{name='solver', code=code}
		self.solverProgramObj:compile()
	end)

	self:refreshCalcDTKernel()

	-- this is created in the parent class, however it isn't called by the parent class.
	--  instead it has to be called by the individual implementation classes
	if self:hasModule(eqn.symbols.addSource) then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name=eqn.symbols.addSource, domain=self.domainWithoutBorder}
	end

	if self:hasModule(eqn.symbols.constrainU) then
		self.constrainUKernelObj = self.solverProgramObj:kernel(eqn.symbols.constrainU)
	end

	for _,op in ipairs(self.ops) do
		if op.refreshSolverProgram then
			op:refreshSolverProgram()
		end
	end

	-- display stuff.  build these just in case trackvars is used.
	do	--if self.app.targetSystem ~= 'console' then

		if self.app.useGLSharing then
			for _,group in ipairs(self.displayVarGroups) do
				group.calcDisplayVarToTexKernelObj = self.solverProgramObj:kernel(group.toTexKernelName)
				group.calcDisplayVarToTexKernelObj.obj:setArg(1, self.texCLMem)
			end
		end

		for _,group in ipairs(self.displayVarGroups) do
			group.calcDisplayVarToBufferKernelObj = self.solverProgramObj:kernel(group.toBufferKernelName)
			group.calcDisplayVarToBufferKernelObj.obj:setArg(1, self.reduceBuf)
		end
	end
end


-- for solvers who don't rely on calcDT
function SolverBase:refreshCalcDTKernel()
	-- if the eqn doesn't have calcDT, let it still run if it provides a fixedDT
	if self:hasModule(self.eqn.symbols.calcDT) then
		self.calcDTKernelObj = self.solverProgramObj:kernel(self.eqn.symbols.calcDT)
		self.calcDTKernelObj.obj:setArg(1, self.reduceBuf)
	end
end

function SolverBase:hasModule(name)
	return self.modules.set[name]
end

function SolverBase:isModuleUsed(name)
--print('self.sharedModulesEnabled', table.keys(self.sharedModulesEnabled):unpack())
--print('self.solverModulesEnabled', table.keys(self.solverModulesEnabled):unpack())
	local moduleNames = table(self.sharedModulesEnabled, self.solverModulesEnabled):keys()
--print('moduleNames', moduleNames:unpack())
	local modulesEnabled = self.modules:getDependentModules(moduleNames:unpack())
		:mapi(function(module) return true, module.name end)
	return modulesEnabled[name]
end

function SolverBase:getDisplayCode()
	local lines = table()
	lines:insert(template([[

typedef union displayValue_t {
	real	ptr[9];
	real	vreal;
<? if solver:isModuleUsed'real3s3' then ?>
	real3s3	vreal3s3;
<? end ?>
<? if solver:isModuleUsed'cplx' then ?>
	cplx	vcplx;
<? end ?>
<? if solver:isModuleUsed'real3' then ?>
	real3	vreal3;
<? end ?>
<? if solver:isModuleUsed'cplx3' then ?>
	cplx3	vcplx3;
<? end ?>
<? if solver:isModuleUsed'real3x3' then ?>
	real3x3	vreal3x3;
<? end ?>
} displayValue_t;
]], {
		solver = self,
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
<? if solver:isModuleUsed'real3' then --\
?>	if (vectorField) {\
		dest[0+3*dstindex] = <?=accumFunc?>(value.ptr[0], dest[0+3*dstindex]);\
		dest[1+3*dstindex] = <?=accumFunc?>(value.ptr[1], dest[1+3*dstindex]);\
		dest[2+3*dstindex] = <?=accumFunc?>(value.ptr[2], dest[2+3*dstindex]);\
	} else\
<? end --\
?>	{\
		dest[dstindex] = <?=accumFunc?>(value.ptr[0], dest[dstindex]);\
	}

]],		{
			solver = self,
			accumFunc = accumFunc,
		}))
	else
		lines:insert(template([[
#define END_DISPLAYFUNC_TEX()\
	write_imagef(tex, <?= solver.dim == 3 and 'i' or 'i.xy'?>, (float4)(value.ptr[0], value.ptr[1], value.ptr[2], 0.));

#define END_DISPLAYFUNC_BUFFER()\
<? if solver:isModuleUsed'real3' then --\
?>	if (vectorField) {\
		((global real3*)dest)[dstindex] = value.vreal3;\
	} else\
<? end --\
?>	{\
		dest[dstindex] = value.vreal;\
	}

]], 	{
			solver = self,
		}))
	end

--[[
ok here's my dilemma
U group (which I considered the default display group)
uses default "pickComponent" as its name.
its group has a name, but its pickComponent fake-group doesn't.
but the einstein stuff needs to be associated only with 'U' ...
... so looks like I can't use the default for 'U' anymore?
--]]
	local alreadyAddedComponentForGroup = {}
	local function addPickComponetForGroup(group)
		local name = self:getPickComponentNameForGroup(group)
		if alreadyAddedComponentForGroup[name] then return end
		alreadyAddedComponentForGroup[name] = true
		lines:insert((self.eqn:template([[
static inline void <?=name?>(
	constant <?=solver_t?> const * const solver,
	global <?=group.bufferType?> const * const buf,
	int const component,
	int * const vectorField,
	displayValue_t * const value,
	int4 const i,
	int const index,
	global <?=cell_t?> const * const cellBuf
) {
	real3 const x = cellBuf[index].pos;
	switch (component) {
<?
for i,component in ipairs(solver.displayComponentFlatList) do
	if not component.onlyFor
	or (group.name == component.onlyFor)
	then
		if solver:isModuleUsed(component.base) then
			local label = (component.base or 'real')..' '..component.name
?>	case <?=i?>:	//<?=label?>
		{
			<?=component.code:gsub('\n', '\n\t\t\t')?>
			*vectorField = <?= solver:isVarTypeAVectorField(component.type) and '1' or '0' ?>;
			break;
		}
<?
		end
	end
end
	-- should default give some default component, like the first real?
	-- or should default give a bad value to indicate an error has happened?
?>	default:
		*vectorField = 0;
		value->vreal = 0.987654321;
		break;
	}
}
]], 	{
			name = name,
			group = group,
		})))
	end
	addPickComponetForGroup{
		bufferType = self.eqn.symbols.cons_t,
	}

--[=[
-- ok here's another idea for saving lines of code
-- add some predefined functions for getters for real, real3, real3s3, cplx, cplx3, real3x3
-- and, if your variable is just a member of a struct of one of these types, use that.
-- (this means adding in extra params for display: the offset and the struct size)

	for _,ctype in ipairs{'real', 'real3', 'real3s3', 'cplx', 'cplx3', 'real3x3'} do
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
	for _,group in ipairs(self.displayVarGroups) do
		if self.app.useGLSharing then
			group.toTexKernelName = self.app:uniqueName'calcDisplayVarToTex'
		end
		group.toBufferKernelName = self.app:uniqueName'calcDisplayVarToBuffer'
	end

	local lastDisplayVarIndex = 0
	for _,texVsBuf in ipairs(table{'Buffer'}:append(
		self.app.useGLSharing and {'Tex'} or nil
	)) do

		for _,group in ipairs(self.displayVarGroups) do

			-- TODO ... right now U is the only one that deviates, so, pick
			-- I guess types make a difference ... maybe ...
			addPickComponetForGroup(group)

			lines:insert(self.eqn:template([[
//<?=group.name?>
kernel void <?=kernelName?>(
	constant <?=solver_t?> const * const solver,
	<?=outputArg?>,
	global <?=group.bufferType?> const * const buf,
	int const displayVarIndex,
	int const component,
	global <?=cell_t?> const * const cellBuf,

	// this is an ugly ugly hack,
	// because calculating avg and stddev use a buffer that includes ghost cells
	// so unless I zero the border, it'll skew the averages
	// but in fact the only reason i'm not zeroing the border is for the display i think?
	// which sounds like a stupid reason anyways ... i should just use gl wrap
	int const zeroBorder
<?
if require 'hydro.solver.meshsolver':isa(solver) then
-- TODO is anyone even using this?  faces? or extraArgs?
?>,
	global <?=solver.coord.face_t?> const * const faces	//[numFaces]<?
end ?><?=group.extraArgs and #group.extraArgs > 0
		and ',\n\t'..table.concat(group.extraArgs, ',\n\t')
		or '' ?>
) {
	<?=SETBOUNDS?>(0,0);
<? if not require 'hydro.solver.meshsolver':isa(solver) then
?>	bool const oob = <?=OOB?>(solver->numGhost, solver->numGhost);
	int4 dsti = i;
	int dstindex = index;
	real3 x = cellBuf[index].pos;
<? for j=0,solver.dim-1 do
?>	i.s<?=j?> = clamp(i.s<?=j?>, solver->numGhost, solver->gridSize.s<?=j?> - solver->numGhost - 1);
<? end
?>	index = INDEXV(i);
<? else	-- mesh
?>	bool const oob = false;
	int dstindex = index;
	real3 x = cellBuf[index].pos;
<? end 		-- mesh vs grid
?>	displayValue_t value = {.ptr={0}};

<?=addTab(group.codePrefix or '')
?>
	global <?=cell_t?> const * const cell = cellBuf + index;

	int vectorField = 0;
	if (!(zeroBorder && oob)) {
		switch (displayVarIndex) {
]], 			table({
					group = group,
					kernelName = ({
						Tex = group.toTexKernelName,
						Buffer = group.toBufferKernelName,
					})[texVsBuf] or error'here',
					outputArg = ({
						Tex =
							-- nvidia needed 'write_only', but I don't want to write only -- I want to accumulate and do other operations
							-- TODO if I do accumulate, then I will need to ensure the buffer is initialized to zero ...
							--'write_only '..
							-- annnnd intel correctly failed with write_only set.  so it goes back out.  until I run this on nvidia again, and work around their trash quality support of OpenCL.
							(self.dim == 3 and 'image3d_t' or 'image2d_t')..' tex',
						Buffer = 'global real* dest',
					})[texVsBuf] or error'here',
					addTab = addTab,
				}, args)
			))

			for _,var in ipairs(group.vars) do
				if var.originalVar then
					var.displayVarIndex = var.originalVar.displayVarIndex
				else
					if not var.displayVarIndex then
						var.displayVarIndex = lastDisplayVarIndex
						lastDisplayVarIndex = lastDisplayVarIndex + 1
					end
					lines:insert('		//'..var.name)
					lines:insert('		case '..var.displayVarIndex..': {')

					-- hmm, this will change its success the next time through this test
					-- so this is destructive, only run it once per display var?
					-- or better yet TODO somewhere earlier, maybe before the 'prefix func' stuff,
					-- prepend 'var.codePrefix' onto 'var.code'
					--var.code = addTab(var.code)

					--[[
					if var.enabled
					or (var.vecVar and var.vecVar.enabled)
					then
					--]]do
						local env = {
							solver = self,
							var = var,
							texVsBuf = texVsBuf,
							name = ({
								Tex = group.toTexKernelName,
								Buffer = group.toBufferKernelName,
							})[texVsBuf] or error'here',
							addTab = addTab,
						}
						lines:insert(template([[
<?=addTab(addTab(addTab(addTab(var.code))))
?>				vectorField = <?=solver:isVarTypeAVectorField(var.type) and '1' or '0'?>;
]], env))
					end
					lines:insert'				break;'
					lines:insert'			}'
					lines:insert''
				end
			end	-- var

--[[
TODO ok here's an issue:
Sometimes group.bufferType is not cons_t, such as for the group of reduceBuf, in which case it is 'real'.
But in this case we are still calling the same pickComponent() as UBuf, and its bufferType is going to differ.
--]]
			lines:insert(template([[
		}
		<?=solver:getPickComponentNameForGroup(group)?>(solver, buf, component, &vectorField, &value, i, index, cellBuf);
	}

	END_DISPLAYFUNC_<?=texVsBuf:upper()?>()
}
]],				{
					solver = self,
					group = group,
					texVsBuf = texVsBuf,
				})
			)
		end	-- group
	end	-- texVsBuf

	local code = lines:concat'\n'

	-- in case any MODULE_* cmds were in there, align them back with the lhs
	code = string.split(code, '\n'):mapi(function(l)
		return (l:gsub('^%s*//// MODULE_', '//// MODULE_'))
	end):concat'\n'

	return code
end


function SolverBase:refreshInitStateProgram()
	self.eqn.initCond:refreshInitStateProgram(self)
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

	if self.eqn.resetState then
		self.eqn:resetState()
		self:boundary()
	end

	for _,op in ipairs(self.ops) do
		if op.resetState then
			op:resetState()
			self:boundary()
		end
	end

	self:constrainU()
end

function SolverBase:constrainU()
	if self.constrainUKernelObj then
		self.constrainUKernelObj(self.solverBuf, self.UBuf, self.cellBuf)
		if cmdline.printBufs == 2 then
			print()
			print('constrainU')
			self:printBuf(self.UBufObj)
		end
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

function SolverBase:convertToSIUnitsCode(units)
	units = units or '1'
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
	expr = symmath.clone(expr)()
	expr = expr:map(function(ex)
		if symmath.op.pow:isa(ex) then
			local power = ex[2].value
			assert(type(power) == 'number')
			if power == math.floor(power) and power > 0 then
				if power == 1 then return ex[1] end
				return symmath.op.mul(table{ex[1]}:rep(power):unpack())
			end
		end
	end)

	local Ccode = symmath.export.C:toCode{
		output = {expr},
		input = {
			{__solver_meter = m},
			{__solver_second = s},
			{__solver_kilogram = kg},
			{__solver_coulomb = C},
			{__solver_kelvin = K},
		},
	}
	Ccode = Ccode:gsub('__solver_', 'solver->')

	local luaFunc, luaCode = symmath.export.Lua:toFunc{
		output = {expr},
		input = {m, s, kg, C, K},
	}

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
	self.solver = assert(args.solver)

	self.name = assert(args.name)

	self.codePrefix = args.codePrefix

	-- the type of the input buffer we're deriving display vars from
	self.bufferType = assert(args.bufferType)

	-- where to look in the solver for the buffer
	self.bufferField = args.bufferField

	-- maybe this should be in args too?
	-- or - instead of buffer - how about all the kernel's args?
	-- but the reason I have to store the field here is that the buffer isn't made yet
	-- TODO? make display vars after buffers so I can store the buffer here?
	self.getBuffer = args.getBuffer
	if not self.getBuffer then
		assert(self.bufferField, "expected bufferField or getBuffer")
		self.getBuffer = function()
			return self.solver[self.bufferField]
		end
	end

	self.extraArgs = args.extraArgs

	self.vars = table(args.vars)
end


local DisplayVar = class()

-- this is the default DisplayVar
SolverBase.DisplayVar = DisplayVar

--[[
component is going to be one of several pathways to modify the data
I'm doing this in hopes to reduce the number of display kernels
The downside is that I now need to consider all permutations up front, rather than recursively create display kernels

TODO buf (dest) shouldn't have ghost cells
and dstIndex should be based on the size without ghost cells

why would I bother write to the ghost cells?
the only reason I can think of is for good subtexel lookup when rendering
--]]
function DisplayVar:init(args)
	self.code = assert(args.code)
	self.name = assert(args.name)
	self.solver = assert(args.solver)
	self.bufferType = args.bufferType	-- or self.bufferType
assert(not args.codePrefix, "move this to DisplayVarGroup")

	assert(not args.displayCode, "I was pretty sure no one was overriding DisplayVar displayCode anymore")
	--self.displayCode = args.displayCode 	-- or self.displayCode

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
	self.useLog = args.useLog or cmdline.display_useLog or false
	self.color = vec3d(math.random(), math.random(), math.random()):normalize()
	self.heatMapFixedRange = cmdline.display_fixedRange or false	-- args.name ~= 'error'
	self.heatMapValueMin = 0
	self.heatMapValueMax = 1

	self.extraArgs = args.extraArgs
end

local intptr = ffi.new('int[1]', 0)
local function int(x)
	intptr[0] = x
	return intptr
end
function DisplayVar:setArgs(kernel)
	local buffer = assert(self.group.getBuffer(), "failed to find buffer for var "..tostring(self.name))
	kernel:setArg(0, self.solver.solverBuf)
	kernel:setArg(2, buffer)
	kernel:setArg(3, int(self.displayVarIndex))
	kernel:setArg(4, int(self.component))
	kernel:setArg(5, self.solver.cellBuf)
end

function DisplayVar:setToTexArgs()
	self:setArgs(self.group.calcDisplayVarToTexKernelObj.obj)
end

function DisplayVar:setToBufferArgs()
	self:setArgs(self.group.calcDisplayVarToBufferKernelObj.obj)
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

	onlyFor = the component is only applied to the specified buffer.  hack for einstein solvers which uses norms based on metrics based on buffer fields.
--]]
function SolverBase:createDisplayComponents()
	self.displayComponents = table()
	self:addDisplayComponents('real', {
		{name = 'default'},
	})
	self:addDisplayComponents('real3', {
		{name = 'default', type = 'real3', magn='mag'},
		{name = 'mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3),0,0);'},
		{name = 'mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(value->vreal3, x),0,0);
]]},
		{name = 'x', code = 'value->vreal3 = _real3(value->vreal3.x,0,0);'},
		{name = 'y', code = 'value->vreal3 = _real3(value->vreal3.y,0,0);'},
		{name = 'z', code = 'value->vreal3 = _real3(value->vreal3.z,0,0);'},
		{name = 'xy arg', code = 'value->vreal3 = _real3(atan2(value->vreal3.y, value->vreal3.x),0,0);'},
		{name = 'yz arg', code = 'value->vreal3 = _real3(atan2(value->vreal3.z, value->vreal3.y),0,0);'},
		{name = 'zx arg', code = 'value->vreal3 = _real3(atan2(value->vreal3.x, value->vreal3.z),0,0);'},
	})
	self:addDisplayComponents('real3s3', {
		{name = 'xx', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.xx,0,0,0,0,0);'},
		{name = 'xy', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.xy,0,0,0,0,0);'},
		{name = 'xz', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.xz,0,0,0,0,0);'},
		{name = 'yy', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.yy,0,0,0,0,0);'},
		{name = 'yz', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.yz,0,0,0,0,0);'},
		{name = 'zz', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.zz,0,0,0,0,0);'},
		{name = 'norm', code = 'value->vreal3s3 = _real3s3(sqrt(real3s3_dot(value->vreal3s3, value->vreal3s3)), 0,0,0,0,0);'},
		{name = 'tr', code = 'value->vreal3s3 = _real3s3(real3s3_trace(value->vreal3s3), 0,0,0,0,0);'},
		{name = 'det', code = 'value->vreal3s3 = _real3s3(real3s3_det(value->vreal3s3), 0,0,0,0,0);'},

		{name = 'x', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.xx, value->vreal3s3.xy, value->vreal3s3.xz, 0,0,0);', type = 'real3', magn='x mag'},
		{name = 'y', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.xy, value->vreal3s3.yy, value->vreal3s3.yz, 0,0,0);', type = 'real3', magn='y mag'},
		{name = 'z', code = 'value->vreal3s3 = _real3s3(value->vreal3s3.xz, value->vreal3s3.yz, value->vreal3s3.zz, 0,0,0);', type = 'real3', magn='z mag'},
		{name = 'x mag', code = 'value->vreal3s3 = _real3s3(real3_len(real3s3_x(value->vreal3s3)), 0,0,0,0,0);'},
		{name = 'y mag', code = 'value->vreal3s3 = _real3s3(real3_len(real3s3_y(value->vreal3s3)), 0,0,0,0,0);'},
		{name = 'z mag', code = 'value->vreal3s3 = _real3s3(real3_len(real3s3_z(value->vreal3s3)), 0,0,0,0,0);'},
		{name = 'x mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3s3 = _real3s3(coordLen(real3s3_x(value->vreal3s3), x), 0,0,0,0,0);
]]},
		{name = 'y mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3s3 = _real3s3(coordLen(real3s3_y(value->vreal3s3), x), 0,0,0,0,0);
]]},
		{name = 'z mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3s3 = _real3s3(coordLen(real3s3_z(value->vreal3s3), x), 0,0,0,0,0);
]]},
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
		{name = 'mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coord_g_ll?>
value->vcplx3 = _cplx3(cplx_from_real(cplx3_weightedLenSq(value->vcplx3, coord_g_ll(x))), cplx_zero, cplx_zero);
]]},

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
		{name = 're mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vcplx3 = _cplx3(cplx_from_real(coordLen(cplx3_re(value->vcplx3), x)), cplx_zero, cplx_zero);
]]},
		{name = 'im', code = 'value->vreal3 = cplx3_im(value->vcplx3); *(real3*)(value+3) = real3_zero;', type = 'real3', magn='im mag'},
		{name = 'im mag', code = 'value->vcplx3 = _cplx3(cplx_from_real(real3_len(cplx3_im(value->vcplx3))), cplx_zero, cplx_zero);'},
		{name = 'im mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vcplx3 = _cplx3(cplx_from_real(coordLen(cplx3_im(value->vcplx3), x)), cplx_zero, cplx_zero);
]]},
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
		{name = 'tr metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coord_g_ll?>
value->vreal3x3 = _real3x3(real3x3_real3s3_dot(value->vreal3x3, coord_g_ll(x)), 0,0,0,0,0,0,0,0);
]]},

		{name = 'x', code = 'value->vreal3 = value->vreal3x3.x; value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='x mag'},
		{name = 'y', code = 'value->vreal3 = value->vreal3x3.y; value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='y mag'},
		{name = 'z', code = 'value->vreal3 = value->vreal3x3.z; value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='z mag'},
		{name = 'x mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3x3.x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'y mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3x3.y), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'z mag', code = 'value->vreal3 = _real3(real3_len(value->vreal3x3.z), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'x mag metrc', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(value->vreal3x3.x, x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;
]]},
		{name = 'y mag metrc', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(value->vreal3x3.y, x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;
]]},
		{name = 'z mag metrc', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(value->vreal3x3.z, x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;
]]},

		{name = 'T x', code = 'value->vreal3 = _real3(value->vreal3x3.x.x, value->vreal3x3.y.x, value->vreal3x3.z.x); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='T x mag'},
		{name = 'T y', code = 'value->vreal3 = _real3(value->vreal3x3.x.y, value->vreal3x3.y.y, value->vreal3x3.z.y); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='T y mag'},
		{name = 'T z', code = 'value->vreal3 = _real3(value->vreal3x3.x.z, value->vreal3x3.y.z, value->vreal3x3.z.z); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;', type = 'real3', magn='T z mag'},
		{name = 'T x mag', code = 'value->vreal3 = _real3(real3_len(_real3(value->vreal3x3.x.x, value->vreal3x3.y.x, value->vreal3x3.z.x)), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T y mag', code = 'value->vreal3 = _real3(real3_len(_real3(value->vreal3x3.x.y, value->vreal3x3.y.y, value->vreal3x3.z.y)), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T z mag', code = 'value->vreal3 = _real3(real3_len(_real3(value->vreal3x3.x.z, value->vreal3x3.y.z, value->vreal3x3.z.z)), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;'},
		{name = 'T x mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(_real3(value->vreal3x3.x.x, value->vreal3x3.y.x, value->vreal3x3.z.x), x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;
]]},
		{name = 'T y mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(_real3(value->vreal3x3.x.y, value->vreal3x3.y.y, value->vreal3x3.z.y), x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;
]]},
		{name = 'T z mag metric', code = self.eqn:template[[
//// MODULE_DEPENDS: <?=coordLen?>
value->vreal3 = _real3(coordLen(_real3(value->vreal3x3.x.z, value->vreal3x3.y.z, value->vreal3x3.z.z), x), 0,0); value->vreal3x3.y = real3_zero; value->vreal3x3.z = real3_zero;
]]},
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
	component.code = component.code or ''
	component.type = component.type or 'real'
	self.displayComponents[basetype]:insert(component)
end
function SolverBase:finalizeDisplayComponents()
	-- build a 1-based enum of all components
	self.displayComponentFlatList = table()
	self.displayComponentNames = table()
--[[ I could do this, but key order is arbitrary, and cache doesn't like that ...
	for basetype,components in pairs(self.displayComponents) do
--]]
-- [[ ... so I do this instead
	for _,basetype in ipairs(table.keys(self.displayComponents):sort()) do
		local components = self.displayComponents[basetype]
--]]
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
function SolverBase:getPickComponentNameForGroup(group)
	local name = 'pickComponent'
	if group
	and group.name	-- exclude the fake group created for default components
	then
		name = name..'_'
			-- TODO further sanitization?
			..group.name:gsub(' ', '_')
	end
	return name
end


function SolverBase:isVarTypeAVectorField(vartype)
	return vartype == 'real3' or vartype == 'cplx'
	-- real3s3 and cplx3 are too complex to merely be vector fields.
	--  maybe I'll add another display for them later.
end

-- still used by gr-hd-separate to add 'extraArgs'
function SolverBase:getUBufDisplayVarsArgs()
	return {
		codePrefix = self.eqn:getDisplayVarCodePrefix(),
		bufferType = self.eqn.symbols.cons_t,
		bufferField = 'UBuf',
	}
end


function SolverBase:addUBufDisplayVars()
	-- TODO make this getUBufDisplayVarGroupArgs, and make the vars getter separate
	local args = self:getUBufDisplayVarsArgs()

	local group = self:newDisplayVarGroup{
		name = 'U',
		bufferField = args.bufferField,
		bufferType = args.bufferType,
		codePrefix = args.codePrefix,
		extraArgs = args.extraArgs,
	}

	args.bufferField = nil
	args.bufferType = nil
	args.codePrefix =  nil
	args.extraArgs = nil

	args.group = group
	args.vars = self.eqn:getDisplayVars()

	self:addDisplayVarGroup(args, self.DisplayVar_U)
end


-- TODO this is the only function that calls DisplayVarGroup ctor
function SolverBase:newDisplayVarGroup(args)
	local displayVarGroup = DisplayVarGroup(table(args,  {solver=self}))
	self.displayVarGroups:insert(displayVarGroup)
	return displayVarGroup
end

-- TODO proper module with depends and all
function SolverBase:addDisplayVarGroup(args, cl)
	cl = cl or self.DisplayVar

	if not args.group then
		args.group = self:newDisplayVarGroup{
			name = args.name,
			bufferType = args.bufferType,
			bufferField = args.bufferField,
			getBuffer = args.getBuffer,
			codePrefix = args.codePrefix,
			extraArgs = args.extraArgs,
		}
	end

	args.codePrefix = nil
	args.bufferType = nil
	args.bufferField = nil
	args.getBuffer = nil
	args.extraArgs = nil

	local group = args.group

	local enableScalar = true
	local enableVector = true
enableVector = false

	for i,var in ipairs(args.vars) do
		xpcall(function()
			local units = var.units
			var = table(var)
			var.units = nil

			local name = assert(var.name, "expected to find name in "..tolua(var))
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
		end, function(err)
			io.stderr:write(err..'\n'..debug.traceback()..'\n')
			io.stderr:flush()
		end)
	end

	return args.group
end


function SolverBase:addDisplayVars()
	self:addUBufDisplayVars()

	-- might contain nonsense :-p
	-- TODO it also might contain vector components
	self:addDisplayVarGroup{
		name = 'reduce',
		bufferType = 'real',
		getBuffer = function() return self.reduceBuf end,
		vars = {{name='0', code='value.vreal = buf[index];'}},
	}

	local cellStructVars = self.coord.cellStruct.fields:filter(function(var)
		return var.type ~= 'int'
	end)
	for _,var in ipairs(cellStructVars) do
		self.solverModulesEnabled[var.type] = true
	end
	self:addDisplayVarGroup{
		name = 'cell',
		bufferField = 'cellBuf',
		bufferType = self.coord.cell_t,
		vars = self:createDisplayVarArgsForStructVars(
			cellStructVars,
			'cell'
		)
	}

-- [[ use for debugging only for the time being
	-- TODO make this flexible for our integrator
	-- if I put this here then integrator isn't created yet
	-- but if I put this after integrator then display variable init has already happened
	-- also TODO - either only use the UBuf's state variables,
	-- or don't regen all the display var code somehow and just bind derivbuf to the ubuf functions
	do	--if self.integrator.derivBufObj then
		local args = self:getUBufDisplayVarsArgs()

		local group = self:newDisplayVarGroup{
			name = 'deriv',
			codePrefix = args.codePrefix,
			bufferType = args.bufferType,
			getBuffer = function()
				local int = self.integrator
				-- Euler's deriv buffer
				if int.derivBufObj then return int.derivBufObj.obj end
				-- RK4's first deriv buffer
				if int.derivBufObjs and int.derivBufObjs[1] then return int.derivBufObjs[1].obj end
				-- BE's deriv buffer
				if int.krylov_dUdtObj then return int.krylov_dUdtObj.obj end
				print"HERE"
			end,
			extraArgs = args.extraArgs,
		}

		args.codePrefix = nil
		args.bufferType = nil
		args.bufferField = nil
		args.extraArgs = nil

		args.group = group
		args.vars = self:createDisplayVarArgsForStructVars(self.eqn.consStruct.fields[1].type.fields)

		-- why in addUBufDisplayVars() do I make a new group and assign args.group to it?
		self:addDisplayVarGroup(args, self.DisplayVar_U)
	end
--]]
end

--[[
accepts a list of struct var info {name=..., [type=..., units=...]}
returns a list of display var construction info
--]]
function SolverBase:createDisplayVarArgsForStructVars(structVars, ptrName, namePrefix)
	-- initialize structForType
	-- TODO put the real3x3s3Struct initialization somewhere else
	-- TODO is the structForType table the same as typeInfoForCode table within hydro/code/struct.lua?
	if not self.structForType.real3x3s3 then
		local real3x3s3Struct = Struct{
			name = 'real3x3s3',
			fields = {
				{name='x', type='real3s3'},
				{name='y', type='real3s3'},
				{name='z', type='real3s3'},
			},
			cdef = false,
		}.class
		self.structForType.real3x3s3 = real3x3s3Struct
	end

	-- should I always force it to be a ptr, hence always using -> ?
	ptrName = ptrName or 'U'
	ptrName = ptrName .. '->'

	local results = table()	-- array of ctor args for DisplayVars
	for _,var in ipairs(structVars) do
		local substruct = self.structForType[var.type]
		if substruct then
			assert(not var.units)	-- struct which are being recursively called shouldn't have units if their fields have units
			results:append(
				self:createDisplayVarArgsForStructVars(substruct.fields, '(&'..ptrName..var.name..')', var.name)
			)
		else
			results:insert{
				name = (namePrefix and (namePrefix..' ') or '')..var.name,
				code = 'value.v' .. var.type .. ' = ' .. ptrName .. var.name .. ';',
				type = var.type,
				units = var.units,

				-- if a display var has 'field' set then use a predefined calcDisplayVar function to just read the field directly (without any computations required)
				-- ... unless it has units too ... in which case ... I'll be scaling the units
				-- ... of course I could do the scaling after reading the value ...
				field = var.name,
			}
		end
	end
	return results
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
				local dupvar = table(var):setmetatable(getmetatable(var))
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
	self.displayVarForName = self.displayVars:mapi(function(var)
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
	if self.eqn.createDisplayComponents then
		self.eqn:createDisplayComponents()
	end
	self:finalizeDisplayComponents()

	self.displayVarGroups = table()
	self:addDisplayVars()
	self:finalizeDisplayVars()
end

--[[
used by the display code to dynamically adjust ranges
this returns raw values, not scaled by units

TODO
this uses calcDisplayVarToBuffer
which is also used by the display code
but the display code makes room for ghost cells for some reason,
i forget why, something to do with textures and border wrapping arguments

however when it comes to min/max, we don't need ghost cells ...
and when it comes to average and stddev, the ghost cells skew the results

TODO this caches by 't', but what about if dt=0?  or if we do something to the buffer without changing 't' ?
--]]
function SolverBase:calcDisplayVarRange(var, componentIndex)
	componentIndex = componentIndex or var.component
	if var.lastTime == self.t then
		return var.lastMin, var.lastMax
	end
	var.lastTime = self.t
	-- invalidate any values associated with the lastTime
	var.lastMin = nil
	var.lastMax = nil
	var.lastAvg = nil
	var.lastStdDev = nil

	var:setToBufferArgs(componentIndex)

	local channels = 1
	-- this size stuff is very GridSolver-based
	local volume = self.numCells
	local sizevec = var.group.getBuffer().sizevec
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

	return min, max
end

end

-- used by the output to print out avg, min, max
function SolverBase:calcDisplayVarRangeAndAvg(var, componentIndex)
	componentIndex = componentIndex or var.component

	if var.lastTime ~= self.t then
		-- this will update lastTime if necessary
		self:calcDisplayVarRange(var, componentIndex)
		-- displayVarGroup has already set up the appropriate args
	end

	if not var.lastAvg then
		-- TODO the display var is being recalc'd a few times for min/max/avg/stddev?
		-- would it be better to cache and save it?
		-- I'm really not sure ... which is more expensive?
		-- memory or the few calcs each displayvar requires?
		self:calcDisplayVarToBuffer(var, componentIndex, true)

		-- duplicated in calcDisplayVarRange
		local size = self.numCells
		-- size is including the border for reduce, but the contents are all zeroes, so they won't affect the average
		-- very ugly I know.  I should change this so the reduce buf doesn't cover ghost cells
		-- or at least doesn't during the calc display var stuff
		-- (maybe other ops like grav potential do use reduceBuf and need border?)
		local sizeForAvg = self.volumeWithoutBorder
		local sizevec = var.group.getBuffer().sizevec
		if sizevec then
			size = tonumber(sizevec:volume())
			sizeForAvg = size
		end
		-- TODO this is also averagging in the ghost cells ... hmm ....
		var.lastAvg = fromreal(self.reduceSum(nil, size)) / sizeForAvg
	end

	return var.lastMin, var.lastMax, var.lastAvg
end

function SolverBase:calcDisplayVarRangeAndAvgAndStdDev(var, componentIndex)
	componentIndex = componentIndex or var.component

	-- this will update lastTime if necessary
	self:calcDisplayVarRangeAndAvg(var, componentIndex)

	if not var.lastStdDev then
		assert(var.lastAvg)
		self:calcDisplayVarToBuffer(var, componentIndex, true)
		self.subtractAndSquareKernelObj(self.solverBuf, self.reduceBuf, real(var.lastAvg))

		-- duplicated in calcDisplayVarRange
		local size = self.numCells
		local sizeForAvg = self.volumeWithoutBorder
		local sizevec = var.group.getBuffer().sizevec
		if sizevec then
			size = tonumber(sizevec:volume())
			sizeForAvg = size
		end
		var.lastStdDev = math.sqrt(fromreal(self.reduceSum(nil, size)) / sizeForAvg)
	end

	return var.lastMin, var.lastMax, var.lastAvg, var.lastStdDev
end

-------------------------------------------------------------------------------
--                              gui                                          --
-------------------------------------------------------------------------------


SolverBase.fpsNumSamples = 30


function SolverBase:calcDT()
	local dt
	-- calc cell wavespeeds -> dts
	if self.useFixedDT then
		dt = self.fixedDT
	else
		assert(self.calcDTKernelObj)
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
	local thisTime = tonumber(getTime())
	if not self.fpsSampleCount then self.fpsSampleCount = 0 end
	if not self.lastFrameTime then self.lastFrameTime = thisTime end

	local tick = cmdline.tick or 0

	self.fpsSampleCount = self.fpsSampleCount + 1

	if (self.showFPS or cmdline.trackvars)
	and thisTime - self.lastFrameTime >= tick
	then
		local plotsOnExit = self.plotsOnExit
		if cmdline.plotOnExit then
			if not plotsOnExit then
				plotsOnExit = table()
				self.plotsOnExit = plotsOnExit
			end
		end

		local deltaTime = thisTime - self.lastFrameTime
		self.fps = self.fpsSampleCount / deltaTime

		--io.write(tostring(self), ' ')
		local sep = ''
		if self.showFPS then
			io.write(sep, 'fps=', self.fps)
			sep = '\t'
		end
		io.write(sep, 't=', self.t)
		sep = '\t'
		if cmdline.trackvars then
			if cmdline.plotOnExit then
				plotsOnExit.t = plotsOnExit.t or table()
				plotsOnExit.t:insert(self.t)
			end
			local varnames = string.split(cmdline.trackvars, ','):mapi(string.trim)
			if varnames:find'dt' then
				io.write(sep, 'dt=', self.dt or 0)	-- the first frame it won't be there ... unless I move this ...
			end
			for _,varname in ipairs(varnames) do
				if varname ~= 'dt' then
					local var = assert(self.displayVarForName[varname], "couldn't find "..varname)
					local ymin, ymax, yavg, ystddev = self:calcDisplayVarRangeAndAvgAndStdDev(var)
					if var.showInUnits and var.units then
						local unitScale = self:convertToSIUnitsCode(var.units).func()
						ymin = ymin * unitScale
						yavg = yavg * unitScale
						ymax = ymax * unitScale
						ystddev = ystddev * unitScale
					end
					io.write(sep, varname, '=[', ymin, ' ', yavg, ' ±', ystddev, ' ', ymax, ']')

					if cmdline.plotOnExit then
						local key = varname..' min'
						plotsOnExit[key] = plotsOnExit[key] or table()
						plotsOnExit[key]:insert(ymin)
						local key = varname..' avg'
						plotsOnExit[key] = plotsOnExit[key] or table()
						plotsOnExit[key]:insert(yavg)
						local key = varname..' stddev'
						plotsOnExit[key] = plotsOnExit[key] or table()
						plotsOnExit[key]:insert(ystddev)
						local key = varname..' max'
						plotsOnExit[key] = plotsOnExit[key] or table()
						plotsOnExit[key]:insert(ymax)
					end
				end
			end
		end
		print()

		self.lastFrameTime = thisTime
		self.fpsSampleCount = 0
	end


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

	self.solverPtr.t = self.t
	self.solverPtr.dt = self.dt
	self:refreshSolverBuf()

	if cmdline.testAccuracy then
		local err = self:calcExactError()
		if #self.app.solvers > 1 then
			io.write(self.name,'\t')
		end
		print('t='..self.t..' L1-error='..err)
	end

	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	-- TODO make this the responsibility of every operation that acts on UBuf to call
	-- ... to prevent redundant calls (like here)
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

		if self.addSourceKernelObj then
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
			self:constrainU()
			op:step(dt)
			if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
		end
	end
end


-- check for nans
-- expects buf to be of type cons_t, made up of numStates real variables
function SolverBase:checkFinite(buf)
	local ptr0size = tonumber(ffi.sizeof(buf.type))
	local realSize = tonumber(ffi.sizeof'real')
	local ptrsPerReal = ptr0size / realSize
	assert(ptrsPerReal == math.floor(ptrsPerReal))
	local size = buf.count * ptrsPerReal

	-- TODO calculate flags once and store?
	local checkNaNFlags = string.split(tostring(self.checkNaNs), ','):mapi(function(v) return true,v end):setmetatable(nil)

	if checkNaNFlags.gpu then
		if size ~= self.numCells * self.eqn.numStates then
			error("expected size="..size.." to be "..(self.numCells * self.eqn.numStates))
		end
		-- right now findNaNs is hardcoded to numCells
		-- and if you change that, make sure that it doesn't write past 'reduceBufObj' (which is numCells in size)
		-- or ... allocate a bigger buffer
		self.reduceBufObj:fill()
		self.findNaNsKernelObj.obj:setArg(2, buf.obj)
		self.findNaNsKernelObj()	-- convert buf into reduceBuf 0 or 1 if it is finite or not
		local max = fromreal(self.reduceMax(nil, self.numCells))

		if max == 0 then return true end

--[[ return false ... or ... only now, print out all entries
		return false, 'found non-finite offsets and numbers'
			..' at t='..self.t
--]]
	end

	-- don't free the original ptr too soon
	local ptrorig = buf:toCPU()
	local ptr = ffi.cast('real*', ptrorig)

	local found
	local function callback(i)
		assert(i >= 0 and i < size)
		local x = fromreal(ptr[i])
		if not math.isfinite(x) then
			found = found or table()
			local ins = {i, x}
			-- for certain bufs show the field
			-- TODO associate each type with the array of fields creating the struct, then reverse lookup on arbitrary types to find the field
			--if buf == self.UBufObj then
			if buf.type == self.UBufObj.type	-- self.eqn.symbols.cons_t
			-- then find the factor that the count is from the UBufObj (fluxBufObj will be 'dim' factor, UBufObj will be 1)
			and buf.count == self.UBufObj.count
			then
				local vars = self.eqn.consStruct.fields[1].type.fields
				local numScalars = HydroStruct.countScalars{vars=self.eqn.consStruct.fields[1].type.fields}
--assert(numScalars == ptrsPerReal)
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

	if checkNaNFlags.noghost	-- set to string via cmdline
	and size == self.numCells * ptrsPerReal
	and not require 'hydro.solver.meshsolver':isa(self)
	then
		if self.dim == 1 then
			for i=tonumber(self.numGhost),tonumber(self.gridSize.x-self.numGhost-1) do
				for e=0,ptrsPerReal-1 do
					callback(e + ptrsPerReal * i)
				end
			end
		elseif self.dim == 2 then
			for i=tonumber(self.numGhost),tonumber(self.gridSize.x-self.numGhost-1) do
				for j=tonumber(self.numGhost),tonumber(self.gridSize.y-self.numGhost-1) do
					for e=0,ptrsPerReal-1 do
						callback(e + ptrsPerReal * (i + self.gridSize.x * j))
					end
				end
			end
		elseif self.dim == 3 then
			for i=tonumber(self.numGhost),tonumber(self.gridSize.x-self.numGhost-1) do
				for j=tonumber(self.numGhost),tonumber(self.gridSize.y-self.numGhost-1) do
					for k=tonumber(self.numGhost),tonumber(self.gridSize.z-self.numGhost-1) do
						for e=0,ptrsPerReal-1 do
							callback(e + ptrsPerReal * (i + self.gridSize.x * (j + self.gridSize.y * k)))
						end
					end
				end
			end
		else
			error("here")
		end
	else
		for i=0,size-1 do
			callback(i)
		end
	end

	if not found then return true end

--	self:printBuf(nil, ptr)
	return false, 'found non-finite offsets and numbers'
		..(checkNaNFlags.all
			and (': '..tolua(found))
			or (', first entry: '..tolua(found[1]))
		)
		..' at t='..self.t
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
	if buf.type == self.eqn.symbols.cons_t then
		local maxdigitlen = #tostring(self.numCells-1)
--print('maxdigitlen', maxdigitlen)
		local realsPerCell = math.floor(size / self.numCells)
--print('realsPerCell', realsPerCell)
		colmax = colmax or realsPerCell
--print('colmax', colmax)
		if colmax > realsPerCell then
			error("got too many realsPerCell\n"..tolua{
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
--print('colsize', colsize)
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
		var.group.calcDisplayVarToTexKernelObj()

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

		local sizevec = var.group.getBuffer().sizevec or self.texSize
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
		glreport'here'	-- this intermittantly reports on my AMD immediately after I re-randomize the palette
		-- use texSize because meshsolver dim isn't reliable
		if self.texSize.z == 1 then
			gl.glTexSubImage2D(tex.target, 0, 0, 0, sizevec.x, sizevec.y, format, gltype, destPtr)
		else
			for z=0,tex.depth-1 do
				gl.glTexSubImage3D(tex.target, 0, 0, 0, z, sizevec.x, sizevec.y, 1, format, gltype, destPtr + channels * sizevec.x * sizevec.y * z)
			end
		end
		glreport'here'
		tex:unbind()
		glreport'here'
	end
end

-- this is abstracted because accumBuf might want to be used ...
function SolverBase:calcDisplayVarToBuffer(var, componentIndex, zeroBorder)
	componentIndex = componentIndex or var.component
	local component = self.displayComponentFlatList[componentIndex]
	local vectorField = self:isVarTypeAVectorField(component.type)
	local channels = vectorField and 3 or 1

	-- duplicated in calcDisplayVarRange
	local volume = self.numCells
	local sizevec = var.group.getBuffer().sizevec
	if sizevec then
		volume = tonumber(sizevec:volume())
	end

	if self.displayVarAccumFunc	then
		self.cmds:enqueueCopyBuffer{src=self.accumBuf, dst=self.reduceBuf, size=ffi.sizeof(self.app.real) * volume * channels}
	end
	var:setToBufferArgs()
	var.group.calcDisplayVarToBufferKernelObj.obj:setArg(1, self.reduceBuf)
	var.group.calcDisplayVarToBufferKernelObj.obj:setArg(3, int(var.displayVarIndex))
	var.group.calcDisplayVarToBufferKernelObj.obj:setArg(4, int(componentIndex))
	var.group.calcDisplayVarToBufferKernelObj.obj:setArg(5, self.cellBuf)
	var.group.calcDisplayVarToBufferKernelObj.obj:setArg(6, int(zeroBorder and 1 or 0))
	var.group.calcDisplayVarToBufferKernelObj()
	if self.displayVarAccumFunc then
		self.cmds:enqueueCopyBuffer{src=self.reduceBuf, dst=self.accumBuf, size=ffi.sizeof(self.app.real) * volume * channels}
	end
end


function SolverBase:updateGUIParams()
	ig.igText('t: '..self.t)

	-- hmm put fps somewhere else, or put ms update here
	ig.igText('fps: '..(self.fps and tostring(self.fps) or ''))

	ig.luatableTooltipCheckbox('check NaNs', self, 'checkNaNs')

	ig.luatableTooltipCheckbox('use fixed dt', self, 'useFixedDT')
	ig.igSameLine()

	ig.luatableTooltipInputFloatAsText('fixed dt', self, 'fixedDT')
	ig.luatableTooltipInputFloatAsText('CFL', self, 'cfl')


	if self.allowAccum then
		if ig.luatableTooltipCheckbox('accum', self, 'displayVarAccumFunc') then
			self:refreshSolverProgram()
			self.cmds:enqueueFillBuffer{buffer=self.accumBuf, size=ffi.sizeof(self.app.real) * self.numCells * 3}
		end
	end

	if ig.luatableTooltipCombo('integrator', self, 'integratorIndex', integratorNames) then
		self:refreshIntegrator()
	end

	-- I think I'll display my GMRES # steps to converge / epsilon error ...
	if self.integrator.updateGUI then
		ig.igSameLine()
		ig.igPushID_Str'integrator'
		if ig.igCollapsingHeader':' then
			self.integrator:updateGUI()
		end
		ig.igPopID()
	end

	for i,op in ipairs(self.ops) do
		if op.updateGUI then
			ig.igPushID_Int(i)
			op:updateGUI()
			ig.igPopID()
		end
	end

	if ig.luatableTooltipCombo('flux limiter', self, 'fluxLimiter', self.app.limiterNames) then
		self:refreshSolverProgram()
	end
end

function SolverBase:updateGUIEqnSpecific()
--[[ TODO why is this crashing
	if ig.luatableTooltipCombo('init state', self, 'initCondIndex', self.eqn.initCondNames) then
		-- TODO hmm ... the whole point of making a separate initCondProgram was to be able to refresh it without rebuilding all of the solver ...
		-- TODO try again once initCond_t is separated from solver_t
		self:refreshEqnInitState()
	end
--]]
	for _,var in ipairs(self.eqn.initCond.guiVars) do
		var:updateGUI(self)
	end
	for _,var in ipairs(self.eqn.guiVars) do
		var:updateGUI(self)
	end
	-- this is a one-off
	if self.eqn.updateGUI then
		self.eqn:updateGUI()
	end
end

do
	local function handle(self, var, title)
		local anyChanged = false
		ig.igPushID_Str(title)

		var.enabled = not not var.enabled
		local enableChanged = ig.luatableTooltipCheckbox('enabled', var, 'enabled')
		anyChanged = anyChanged or enableChanged
		ig.igSameLine()

		anyChanged = ig.luatableTooltipCheckbox('log', var, 'useLog') or anyChanged
		ig.igSameLine()

		anyChanged = ig.luatableTooltipCheckbox('units', var, 'showInUnits') or anyChanged
		ig.igSameLine()

		anyChanged = ig.luatableTooltipCheckbox('fixed range', var, 'heatMapFixedRange') or anyChanged
		ig.igSameLine()

		--ig.luatableTooltipCombo('component', var, 'component', self.displayComponentNames)

		local name = var.name
		if var.units then name = name..' '..var.units end
		if ig.igCollapsingHeader(name) then
			local unitScale = 1
			if var.units and var.showInUnits then -- convert our ranges from raw to units
				unitScale = self:convertToSIUnitsCode(var.units).func()
				var.heatMapValueMin = var.heatMapValueMin * unitScale
				var.heatMapValueMax = var.heatMapValueMax * unitScale
			end
			anyChanged = ig.luatableTooltipInputFloatAsText('value min', var, 'heatMapValueMin') or anyChanged
			anyChanged = ig.luatableTooltipInputFloatAsText('value max', var, 'heatMapValueMax') or anyChanged
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

		do
			local dim = self.app.displayDim
			if dim == 2 then
				ig.igPushID_Str'2D'
				if self.app.display2DMethodsEnabled.Graph then
					local draw2DGraph = self.draw2DGraph
					if draw2DGraph then
						ig.luatableTooltipInputInt('graph step', draw2DGraph, 'step')
					end
				end
				ig.igPopID()
			elseif dim == 3 then
				ig.igPushID_Str'3D'

				if self.app.display3DMethodsEnabled.Slices then
--[[ currently in app, currently disabled
					if useClipPlanes then
						ig.luatableRadioButton("rotate camera", self.app, 'rotateClip', 0)
						for i,clipInfo in ipairs(clipInfos) do
							ig.igPushID_Str('clip '..i)
							ig.luatableTooltipCheckbox('clip', clipInfo, 'enabled')
							ig.igSameLine()
							ig.luatableRadioButton('rotate', self.app, 'rotateClip', i)
							ig.igSameLine()
							if ig.igButton('reset') then
								clipInfo.plane = makeDefaultPlane(i)
							end
							ig.igPopID()
						end
					end
--]]
					local draw3DSlice = self.draw3DSlice
					if draw3DSlice then
						ig.luatableTooltipSliderFloat('alpha', draw3DSlice, 'alpha', 0, 1)
						ig.luatableTooltipSliderFloat('gamma', draw3DSlice, 'alphaGamma', 0, 1)
						ig.luatableTooltipCheckbox('isobars', draw3DSlice, 'useIsos')
						if draw3DSlice.useIsos then
							ig.luatableTooltipInputInt('num isobars', draw3DSlice, 'numIsobars')
						end
						ig.luatableTooltipCheckbox('lighting', draw3DSlice, 'useLighting')
						ig.luatableTooltipCheckbox('pointcloud', draw3DSlice, 'usePoints')
						if not self.draw3DSlice.usePoints then
							ig.luatableTooltipInputInt('num slices', draw3DSlice, 'numSlices')
						end
					end
				end

				ig.igPopID()
			end

			do
				ig.igPushID_Str'Vector'

				--ig.luatableCheckbox('vector field', self, 'enableVectorField')
				if self.drawVectorArrows then
					ig.luatableTooltipInputFloatAsText('vector field scale', self.drawVectorArrows, 'scale')
					--ig.luatableTooltipSliderFloat('vector field scale', self.drawVectorArrows, 'scale', 0, 100, nil, 10)

					ig.luatableTooltipInputInt('vector field step', self.drawVectorArrows, 'step')
					self.drawVectorArrows.step = math.max(self.drawVectorArrows.step, 1)
				end
				if self.drawVectorLIC then
					ig.luatableTooltipInputInt('LIC steps', self.drawVectorLIC, 'integralMaxIter')
				end

				ig.igPopID()
			end
		end

		if self.guiDisplayFilterStr == nil then self.guiDisplayFilterStr = '' end
		if self.guiDisplayFilterEnabledVars == nil then self.guiDisplayFilterEnabledVars = false end
		local refresh

		ig.luatableTooltipInputText('filter', self, 'guiDisplayFilterStr')
		ig.igSameLine()
		ig.luatableTooltipCheckbox('filter enabled', self, 'guiDisplayFilterEnabledVars')

		for i,displayVarGroup in ipairs(self.displayVarGroups) do
			ig.igPushID_Str('display '..i)
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
function SolverBase:checkStructSizes()
	local typeinfos = table{
		'real',
		'real2',
		'real3',
		'real4',
	}
-- [=[ automatically use types in the code modules
	local moduleNames = table(
		self.sharedModulesEnabled,
		self.solverModulesEnabled,
		self.initModulesEnabled
	):keys():sort()
	for _,module in ipairs(self.modules:getDependentModules(moduleNames:unpack())) do
		for _,struct in ipairs(module.structs) do
			--print('checking for struct '..struct.name..' from module '..module.name)
			if not typeinfos:find(struct)
			and not typeinfos:find(struct.name)
			then
if self.app.verbose then print('adding struct '..struct.name..' from module '..module.name) end
				typeinfos:insert(struct)
			end
		end
	end
--]=]


	local varcount = 0
	for _,typeinfo in ipairs(typeinfos) do
		varcount = varcount + 1
		if type(typeinfo) == 'string' then
			typeinfo = ffi.typeof(typeinfo)
		end
		if Struct:isa(typeinfo) then
			for _,field in typeinfo:fielditer() do
				varcount = varcount + 1
			end
		end
	end

	local cmd = self.cmds
	local _1x1_domain = self.app.env:domain{size={1}, dim=1}
	local resultPtr = ffi.new('size_t[?]', varcount)
	local resultBuf = self.app.env:buffer{name='result', type='size_t', count=varcount, data=resultPtr}

print('shared modules: '..moduleNames:concat', ')
	local codePrefix =
		vec3d.code..'\n'
		..self.modules:getTypeHeader(moduleNames:unpack())
--[=[
<?
local Struct = require 'struct'
local HydroStruct = require 'hydro.code.struct'
?>
	local testStructProgramObj = self.Program{
		name = 'checkStructSizes',
		code = table{
			codePrefix,
			template([[
kernel void checkStructSizes(
	global int* resultBuf
) {

#define offsetof __builtin_offsetof

<?
local index = 0
for i,typeinfo in ipairs(typeinfos) do
	if type(typeinfo) == 'string' then
?>	result[<?=index?>] = sizeof(<?=typeinfo?>);
<?
		index = index + 1
	else
?>	result[<?=index?>] = sizeof(<?=typeinfo.name?>);
<?
		index = index + 1
		if Struct:isa(typeinfo) then
			for _,field in ipairs(typeinfo.fields) do
?>	result[<?=index?>] = offsetof(<?=typeinfo.name?>, <?=field.name?>);
<?
				index = index + 1
			end
		end
	end
end
?>

}
]], 	{
			typeinfos = typeinfos,
		})
		}:concat'\n',
	}
	testStructProgramObj:compile()
--]=]
	local body = template([[
<?
local Struct = require 'struct'
local HydroStruct = require 'hydro.code.struct'
?>
#define offsetof __builtin_offsetof

<?
local index = 0
for i,typeinfo in ipairs(typeinfos) do
	if type(typeinfo) == 'string' then
?>	result[<?=index?>] = sizeof(<?=typeinfo?>);
<?
		index = index + 1
	else
		if Struct:isa(typeinfo) then
?>	result[<?=index?>] = sizeof(<?=typeinfo.name?>);
<?
			index = index + 1
			for _,field in ipairs(typeinfo.fields) do
?>	result[<?=index?>] = offsetof(<?=typeinfo.name?>, <?=field.name?>);
<?
				index = index + 1
			end
		end
	end
end
?>
]], {
		typeinfos = typeinfos,
	})
print'codePrefix:'
print(require 'template.showcode'(codePrefix))
print'body:'
print(require 'template.showcode'(body))
	require 'cl.obj.kernel'{
		env = self.app.env,
		domain = _1x1_domain,
		argsOut = {resultBuf},
		header = codePrefix,
		showCodeOnError = true,
		body = body,
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
			if Struct:isa(typeinfo) then
				local ffisize = tostring(ffi.sizeof(typeinfo.name))
				print('sizeof('..typeinfo.name..'): OpenCL='..clsize..', ffi='..ffisize..(clsize == ffisize and '' or ' -- !!!DANGER!!!'))

				for _,field in ipairs(typeinfo.fields) do
					local cloffset = tostring(resultPtr[index]):match'%d+'
					index = index + 1
					local ffioffset = tostring(ffi.offsetof(typeinfo.name, field.name))
					print('offsetof('..typeinfo.name..', '..field.name..'): OpenCL='..cloffset..', ffi='..ffioffset..(cloffset == ffioffset and '' or ' -- !!!DANGER!!!'))
				end
			end
		end
	end
	print('done')
	os.exit()
end
--]]

return SolverBase
