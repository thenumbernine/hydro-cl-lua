--[[
foundation class for solvers that use unstructured meshes from file

NACA 0012 from...
https://turbmodels.larc.nasa.gov/naca0012_grids.html
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local math = require 'ext.math'
local path = require 'ext.path'
local vec2i = require 'vec-ffi.vec2i'
local vec3sz = require 'vec-ffi.vec3sz'
local vec3f = require 'vec-ffi.vec3f'
local vec3d = require 'vec-ffi.vec3d'
local ig = require 'imgui'
local gl = require 'gl'
local glreport = require 'gl.report'
local template = require 'template'
local clnumber = require 'cl.obj.number'
local SolverBase = require 'hydro.solver.solverbase'
local time, getTime = table.unpack(require 'hydro.util.time')
local real = require 'hydro.real'
local vector = require 'ffi.cpp.vector'

local half = require 'cl.obj.half'
local toreal, fromreal = half.toreal, half.fromreal


local MeshSolver = SolverBase:subclass()

MeshSolver.name = 'meshsolver'

-- for compat in some display stuff
MeshSolver.numGhost = 0 

function MeshSolver:getSymbolFields()
	return MeshSolver.super.getSymbolFields(self):append{
		'boundaryCons',
		'getEdgeStates',
		-- also in gridsolver:
		'OOB',
		'SETBOUNDS',
		'SETBOUNDS_NOGHOST',
		-- also in fvsolver:
		'calcFlux',
		'calcFluxForInterface',
		'calcDerivFromFlux',
	}
end

--[[
args:
	meshfile = name of mesh file to use

NOTICE initCond is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:initMeshVars(args)
	self.showVertexes = cmdline.showVertexes or false
	self.showFaces = cmdline.showFaces or false
	self.showNormals = cmdline.showNormals or false
	self.showValues = cmdline.showValues == nil and true or cmdline.showValues	-- default to true
	self.drawCellScale = 1

	-- TODO make this a param of gridsolver as well
	local fluxName = assert(args.flux, "expected flux")
	local fluxClass = require('hydro.flux.'..fluxName)
	local fluxArgs = table(args.fluxArgs, {solver=self}):setmetatable(nil)
	self.flux = fluxClass(fluxArgs)

	MeshSolver.super.initMeshVars(self, args)

	-- alright, by this point 'gridSize' has become synonymous with global_size() ...
	-- I'm going to have to determine a vec4 globalSize for my kernels once they exceed the max 1D size
	-- (just like I'm already doing with the Tex2D size)
	-- TODO same with mins and maxs, since initCond uses them
	-- move them from mesh to this
	self.solverStructFields:append{
		{name='gridSize', type='int4'},
		{name='stepsize', type='int4'},
	}

	local meshType = assert(args.mesh.type, "expected mesh type")
	local meshFactoryClass = require('hydro.mesh.'..meshType)
	local meshFactory = meshFactoryClass(args.mesh)

	self.dim = meshFactoryClass.dim

--[[ script for testing mesh build times:
require 'ext.range'(1,6,.25):mapi(function(x)
	return true, math.floor(2^x)
end):keys():sort():mapi(function(i)
	local cmd = './run.lua sys=console meshsize='..i
	print('>'..cmd)
	print(os.execute(cmd))
end)
--]]
	-- make sure to create mesh after self.coord, since mesh modifies the coord.cell_t structure
	time('creating mesh', function()
		self.mesh = meshFactory:createMesh(self)
	end)
for k,v in pairs(self.mesh.times) do
	print(k,v)
end
--os.exit()

	-- here's our boundary ctor info hack:
	self.boundaryArgs = args.mesh.boundary

	-- ok now convert from lua tables to arrays
	-- or do this in mesh:calcAux() ?

	self.mins = vec3d(-1, -1, -1)
	self.maxs = vec3d(1, 1, 1)
	for i=0,self.dim-1 do
		self.mins.s[i] = math.huge
		self.maxs.s[i] = -math.huge
	end
	for _,v in ipairs(self.mesh.vtxs) do
		for i=0,self.dim-1 do
			self.mins.s[i] = math.min(self.mins.s[i], v.s[i])
			self.maxs.s[i] = math.max(self.maxs.s[i], v.s[i])
		end
	end

	-- save the smallest dx as well - used for vector arrow display
	self.mindx = math.huge
	for _,f in ipairs(self.mesh.faces) do
		self.mindx = math.min(self.mindx, f.cellDist)
	end

--[[
big problem
in gridsolver, these vars are dependent on the grid size
but in mesh solver they are dependent on the mesh
however the mesh is dependent on the face_t etc types
and those types are defined in initCodeModules
so what do I do?

is it safe to postpone declaring numCells until just before postInit is called, when buffers are allocated?
no, because numCells determines maxWorkGroupSize, which is needed by getSizePropsForWorkGroupSize,
	and I think eventually by the initCodeModule stuff (esp the init_dx stuff in gridsolver)
--]]

	-- TODO put these in solver_t
	self.numCells = assert(self.mesh.numCells, "did your MeshFactory remember to call :calcAux()?")
	self.numFaces = assert(self.mesh.numFaces, "did your MeshFactory remember to call :calcAux()?")
	self.numVtxs = assert(self.mesh.numVtxs, "did your MeshFactory remember to call :calcAux()?")
	self.numCellFaceIndexes = self.mesh.numCellFaceIndexes
	self.numCellVtxIndexes = self.mesh.numCellVtxIndexes
	self.numFaceVtxIndexes = self.mesh.numFaceVtxIndexes

	-- no longer is dim * numCells the number of interfaces -- it is now dependent on the mesh
	-- maybe I should rename this to numFaces?
end

function MeshSolver:postInit()
	MeshSolver.super.postInit(self)
end

function MeshSolver:initObjs(args)
	MeshSolver.super.initObjs(self, args)

	if self.fluxLimiter ~= 1 then
		print("overriding and removing fluxLimiter=="..self.fluxLimiter.." for meshsolver")
		self.fluxLimiter = 1
	end

	self.boundaryMethods = {}
	for key,boundary in ipairs(self.boundaryArgs or {}) do
		local name, boundaryObjArgs
		if type(boundary) == 'string' then
			name = boundary
		elseif type(boundary) == 'table' then
			name = boundary.name
			boundaryObjArgs = boundary.args
		else
			error("unknown boundary ctor: "..tostring(boundary))
		end
		local boundaryClass = assert(self.boundaryOptionForName[name], "failed to find boundary method named "..name)
		
		-- 'key' should be a 1-based integer ...
		self.boundaryMethods[key] = boundaryClass(boundaryObjArgs)
	end
	-- and now we're done with self.boundaryArgs
	self.boundaryArgs = nil
end

function MeshSolver:createCellStruct()
	-- here's the mesh-specific stuff
	self.coord.cellStructFields:append{
--[[ moved to coord/coord to be in gridsolver and meshsolver
		{type='real', name='volume'},	--volume of the cell
--]]			
		{type='int', name='faceOffset'},
		{type='int', name='faceCount'},
		{type='int', name='vtxOffset'},
		{type='int', name='vtxCount'},
	}

	-- here is the mesh-specific face_t fields
	self.coord.faceStructFields:append{
		{type='vec2i_t', name='cells'},	--indexes of cells
		{type='int', name='vtxOffset'},
		{type='int', name='vtxCount'},
		{type='int', name='boundaryMethodIndex'},
	}
end


-- MeshSolver Boundary
-- similar to GridSolver Boundary except GridSolver uses the grid for things like periodic/linear/quadratic cell lookup
-- this is all dependenton the mesh / cell / face structs


local Boundary = class()
MeshSolver.Boundary = Boundary
function Boundary:assignDstSrc(dst, src, args)
	local lines = table()
	if args.fields then
		for _,field in ipairs(args.fields) do
			lines:insert('('..dst..')->'..field..' = ('..src..')->'..field..';')
		end
	else
		lines:insert('*('..dst..') = *('..src..');')
	end
	return lines:concat'\n'
end

local BoundaryNone = Boundary:subclass()
BoundaryNone.name = 'none'
function BoundaryNone:getCode(args)
	return ''
end

-- here is where BoundaryPeriodic would normally go
-- but 'periodic' is not so clear for unstructured grids
-- so we'd need some kind of "read from arbitrary neighbor cell"

local BoundaryMirror = Boundary:subclass()
BoundaryMirror.name = 'mirror'
BoundaryMirror.restitution = 1
function BoundaryMirror:init(args)
	if args then
		self.restitution = args.restitution
	end
end
function BoundaryMirror:getCode(args)
	local dst = assert(args.dst)
	local src = assert(args.src)

	local lines = table()
	
	lines:insert(self:assignDstSrc(dst, src, args))
	
	local solver = args.solver
	local eqn = solver.eqn
	
	lines:insert(template([[
{
	real3 const x = (<?=face?>)->pos;
	real3 const n = (<?=face?>)->normal;	// TODO pos vs neg?
]], {
		face = assert(args.face),
	}))
	for _,var in ipairs(eqn.consStruct.vars) do
		if var.type == 'real' 
		or var.type == 'cplx'
		then
			-- do nothing
		elseif var.type == 'real3' 
		or var.type == 'cplx3'
		then
			-- TODO looks very similar to the reflect code in meshsolver
			lines:insert(template([[
	(<?=dst?>)-><?=field?> = <?=vec3?>_sub(
		(<?=dst?>)-><?=field?>,
		<?=vec3?>_<?=scalar?>_mul(
			<?=vec3?>_from_real3(n),
			<?=scalar?>_real_mul(
				<?=vec3?>_real3_dot(
					(<?=dst?>)-><?=field?>,
					n
				), 
				<?=restitutionPlusOne?>
			)
		)
	);
]], 		{
				restitutionPlusOne = clnumber(self.restitution + 1),
				vec3 = var.type,
				scalar = var.type == 'cplx3' and 'cplx' or 'real',
				field = var.name,
				dst = dst,
			}))
		else
			error("need to support reflect() for type "..var.type)
		end
	end
	lines:insert[[
}
]]
	return lines:concat'\n'
end

local BoundaryFixed = Boundary:subclass()
BoundaryFixed.name = 'fixed'
function BoundaryFixed:init(args)
	-- fixed values to use
	self.fixedCode = assert(args.fixedCode)
end
function BoundaryFixed:getCode(args)
	local dst = assert(args.dst)
	return self:fixedCode(args, dst)
end
MeshSolver.BoundaryFixed = BoundaryFixed 

local BoundaryFreeFlow = Boundary:subclass()
BoundaryFreeFlow.name = 'freeflow'
function BoundaryFreeFlow:getCode(args)
	local dst = assert(args.dst)
	local src = assert(args.src)
	return self:assignDstSrc(dst, src, args)
end

function MeshSolver:createBoundaryOptions()
	self:addBoundaryOptions{
		BoundaryNone,
		BoundaryMirror,
		BoundaryFixed,
		BoundaryFreeFlow,
	}
end

function MeshSolver:createSolverBuf()
	MeshSolver.super.createSolverBuf(self)

	-- TODO most of this matches with gridsolver

	self.solverPtr.gridSize.x = assert(self.numCells)
	self.solverPtr.gridSize.y = 1
	self.solverPtr.gridSize.z = 1
	self.solverPtr.gridSize.w = 1
	self.solverPtr.stepsize.x = 1
	self.solverPtr.stepsize.y = self.solverPtr.gridSize.x
	self.solverPtr.stepsize.z = self.solverPtr.gridSize.x * self.solverPtr.gridSize.y
	self.solverPtr.stepsize.w = self.solverPtr.gridSize.x * self.solverPtr.gridSize.y * self.solverPtr.gridSize.z

	-- while we're here, write all gui vars to the solver_t
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
end

function MeshSolver:initDraw()
	assert(self.app.targetSystem ~= 'console')

	local mesh = self.mesh

	-- TODO move this to hydro/view and set it once for all shaders
	
	-- make a GPU version of the mesh
	-- TODO triangulate ... then include a mapping from triangle to source cell,
	-- and then just copy from display buf to a texture,
	-- then store per vertex the lookup to the texture.
	-- Also, if you are going to scale cells, then you must store a unique vertex per cell here.
	-- this means no dynamic mesh (without more coding).
	local glvtxs = vector'vec3f_t'			-- vertex position
	local glvtxcenters = vector'vec3f_t'	-- center of cell for this vertex
	local glcellindex = vector'float'		-- 0-based index of cell for this vertex
	time('creating display mesh', function()
		local function addTri(va,vb,vc, ci,c)
			glvtxs:push_back(vec3f(va:unpack()))
			glvtxs:push_back(vec3f(vb:unpack()))
			glvtxs:push_back(vec3f(vc:unpack()))
			for j=0,2 do
				glvtxcenters:push_back(vec3f(c.pos:unpack()))
				glcellindex:push_back(ci)
			end
		end
		if self.dim == 2 then
			for ci,c in ipairs(mesh.cells) do
				local va = mesh.vtxs.v[mesh.cellVtxIndexes.v[0 + c.vtxOffset]]
				local vb = mesh.vtxs.v[mesh.cellVtxIndexes.v[1 + c.vtxOffset]]
				for vi=2,c.vtxCount-1 do
					local vc = mesh.vtxs.v[mesh.cellVtxIndexes.v[vi + c.vtxOffset]]
					addTri(vc,vb,va, ci, c)
					vb = vc
				end
			end
		elseif self.dim == 3 then
			for ci,c in ipairs(mesh.cells) do
				for fi=0,c.faceCount-1 do
					local f = mesh.faces.v[mesh.cellFaceIndexes.v[fi + c.faceOffset]]
					local va = mesh.vtxs.v[mesh.faceVtxIndexes.v[0 + f.vtxOffset]]
					local vb = mesh.vtxs.v[mesh.faceVtxIndexes.v[1 + f.vtxOffset]]
					for vi=2,f.vtxCount-1 do
						local vc = mesh.vtxs.v[mesh.faceVtxIndexes.v[vi + f.vtxOffset]]
						addTri(va,vb,vc, ci, c)
						vb = vc
					end
				end
			end
		end
	end)
	
	self.glvtxs = glvtxs
	self.glvtxcenters = glvtxcenters
	self.glcellindex = glcellindex
	self.numGlVtxs = #glvtxs

	local GLArrayBuffer = require 'gl.arraybuffer'
	
	self.glvtxArrayBuffer = GLArrayBuffer{
		data = glvtxs.v,
		size = #glvtxs * ffi.sizeof(glvtxs.type), 
	}:unbind()
	self.glvtxcenterArrayBuffer = GLArrayBuffer{
		data = glvtxcenters.v,
		size = #glvtxcenters * ffi.sizeof(glvtxcenters.type), 
	}:unbind()
	self.glcellindexArrayBuffer = GLArrayBuffer{
		data = glcellindex.v,
		size = #glcellindex * ffi.sizeof(glcellindex.type),
	}:unbind()

	self.drawPointsShader = self.GLProgram{
		name = 'draw_points',
		vertexCode = [[
#version 460

uniform float drawCellScale;
uniform mat4 modelViewProjectionMatrix;

in vec3 vtx;
in vec3 vtxcenter;

void main() {
	vec3 v = (vtx - vtxcenter) * drawCellScale + vtxcenter;
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
}
]],
		fragmentCode = [[
#version 460

out vec4 fragColor;
void main() {
	fragColor = vec4(1., 1., 1., 1.);
}
]],
		attrs = {
			vtx = self.glvtxArrayBuffer,
			vtxcenter = self.glvtxcenterArrayBuffer,
		}
	}:useNone()
end

-- TODO organize this between SolverBase and MeshSolver
function MeshSolver:createBuffers()
	local app = self.app

	-- set texSize before calling super
	if app.targetSystem ~= 'console' then
--print('numCells', self.numCells)		
		local maxTex2DSize = vec3sz(
			self.device:getInfo'CL_DEVICE_IMAGE2D_MAX_WIDTH',
			self.device:getInfo'CL_DEVICE_IMAGE2D_MAX_HEIGHT',
			1)
--print('maxTex2DSize', maxTex2DSize)	
		local maxTex3DSize = vec3sz(
			self.device:getInfo'CL_DEVICE_IMAGE3D_MAX_WIDTH',
			self.device:getInfo'CL_DEVICE_IMAGE3D_MAX_HEIGHT',
			self.device:getInfo'CL_DEVICE_IMAGE3D_MAX_DEPTH')
--print('maxTex3DSize', maxTex3DSize)		
		self.texSize = vec3sz()
		-- TODO if texSize >= max gl size then overflow into the next dim
		if self.numCells <= maxTex2DSize.x then
--print('using 1D texture')			
			self.texSize = vec3sz(self.numCells, 1, 1)
		else
			local sx = math.min(math.ceil(math.sqrt(self.numCells)), tonumber(maxTex2DSize.x))
			local sy = math.ceil(self.numCells / tonumber(self.texSize.x))
			if sx <= maxTex2DSize.x and sy <= maxTex2DSize.y then
--print('using 2D texture')			
				self.texSize = vec3sz(sx, sy, 1)
			else
				local sz = math.min(math.ceil(math.cbrt(self.numCells)), tonumber(maxTex3DSize.z))
				local sxy = math.ceil(self.numCells / sz)
				local sy = math.min(math.ceil(math.sqrt(sxy)), tonumber(maxTex3DSize.y))
				local sx = math.ceil(sxy / sy)
				if sx <= maxTex3DSize.x then
--print('using 3D texture')			
					self.texSize = vec3sz(sx, sy, sz)
				else
					error("couldn't fit cell buffer into texture.  max 2d size " .. maxTex2DSize .. ", max 3d size " .. maxTex3DSize)
				end
			end
		end
	end

	MeshSolver.super.createBuffers(self)

	self:clalloc('vtxBuf', 'real3', self.numVtxs)
	self:clalloc('cellBuf', self.coord.cell_t, self.numCells)
	self:clalloc('faceBuf', self.coord.face_t, self.numFaces)
	self:clalloc('cellFaceIndexesBuf', 'int', self.numCellFaceIndexes)
	
	-- specific to FiniteVolumeSolver
	self:clalloc('fluxBuf', self.eqn.symbols.cons_t, self.numFaces)

-- [[
	-- here or somewhere before 
	-- convert the mesh faces and cells from meshface_t, meshcell_t
	-- to face_t and cell_t
	-- hmm, but if any custom fields are included, I'll bet the subclass will want to populate those ...
	local mesh = self.mesh
	local faces = vector(self.coord.face_t, self.numFaces)
	for i=0,self.numFaces-1 do
		faces.v[i].pos = mesh.faces.v[i].pos
		faces.v[i].normal = mesh.faces.v[i].normal
		faces.v[i].normal2 = mesh.faces.v[i].normal2
		faces.v[i].normal3 = mesh.faces.v[i].normal3
		faces.v[i].area = mesh.faces.v[i].area
		faces.v[i].cellDist = mesh.faces.v[i].cellDist
		faces.v[i].cells = mesh.faces.v[i].cells
		faces.v[i].vtxOffset = mesh.faces.v[i].vtxOffset
		faces.v[i].vtxCount = mesh.faces.v[i].vtxCount
		faces.v[i].boundaryMethodIndex = mesh.faces.v[i].boundaryMethodIndex
	end
	mesh.faces = faces
	
	local cells = vector(self.coord.cell_t, self.numCells)
	for i=0,self.numCells-1 do
		cells.v[i].pos = mesh.cells.v[i].pos
		cells.v[i].volume = mesh.cells.v[i].volume
		cells.v[i].faceOffset = mesh.cells.v[i].faceOffset
		cells.v[i].faceCount = mesh.cells.v[i].faceCount
		cells.v[i].vtxOffset = mesh.cells.v[i].vtxOffset
		cells.v[i].vtxCount = mesh.cells.v[i].vtxCount
	end
	mesh.cells = cells
--]]
end

function MeshSolver:finalizeCLAllocs()
	MeshSolver.super.finalizeCLAllocs(self)

	self.vtxBufObj:fromCPU(self.mesh.vtxs.v)
	self.cellBufObj:fromCPU(self.mesh.cells.v)
	self.faceBufObj:fromCPU(self.mesh.faces.v)
	self.cellFaceIndexesBufObj:fromCPU(self.mesh.cellFaceIndexes.v)
	self.cmds:finish()
end

function MeshSolver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	-- numCells is the number of cells
	-- maybe I should rename it to numCells
	local localSize1d = math.min(maxWorkGroupSize, self.numCells)

	-- cell domain is the default
	self.domain = self.app.env:domain{
		size = {self.numCells},
		dim = 1,
	}

	-- face domain for flux computations
	-- ... and dt ?
	self.faceDomain = self.app.env:domain{
		size = {self.numFaces},
		dim = 1,
	}

	return {
		localSize1d = localSize1d, 
	}
end

function MeshSolver:refreshSolverProgram()
	MeshSolver.super.refreshSolverProgram(self)

	-- solverbase:
	if self.addSourceKernelObj then
		self.addSourceKernelObj.obj:setArg(3, self.cellBuf)
	end
	if self.constrainUKernelObj then
		self.constrainUKernelObj.obj:setArg(2, self.cellBuf)
	end

	-- fvsolver:
	self.calcFluxKernelObj = self.solverProgramObj:kernel{name=self.symbols.calcFlux, domain=self.faceDomain}
	self.calcFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcFluxKernelObj.obj:setArg(1, self.fluxBuf)
	self.calcFluxKernelObj.obj:setArg(2, self.UBuf)
	self.calcFluxKernelObj.obj:setArg(4, self.cellBuf)
	self.calcFluxKernelObj.obj:setArg(5, self.faceBuf)
	self.calcFluxKernelObj.obj:setArg(6, self.cellFaceIndexesBuf)

	self.calcDerivFromFluxKernelObj = self.solverProgramObj:kernel(self.symbols.calcDerivFromFlux)
	self.calcDerivFromFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(3, self.cellBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(4, self.faceBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
end

function MeshSolver:initCodeModules()
	MeshSolver.super.initCodeModules(self)
	
	self.flux:initCodeModules()

-- [[ TODO this used to be in Mesh:getMeshTypeCode
	-- TODO real3 vs vec3f/vec3d ...
	-- TODO what if real3 isn't defined yet?
	local vec2i = require 'vec-ffi.vec2i'
	-- module dependencies are built by struct fields
	-- and mesh adds vec2i to the face_t fields
	-- so make sure there is a vec2i modules
	-- TODO move this to hydro.code.math?
	self.modules:add{
		name = 'vec2i_t',
		typecode = [[

//// BEGIN EXCLUDE FOR FFI_CDEF

]]..vec2i.code..[[

//// END EXCLUDE FOR FFI_CDEF

]],
	}
--]]

	self.modules:addFromMarkup(
		self.eqn:template(path'hydro/solver/meshsolver.cl':read())
	)

	self.solverModulesEnabled[self.symbols.calcFlux] = true
	self.solverModulesEnabled[self.symbols.calcDerivFromFlux] = true
end

function MeshSolver:refreshCalcDTKernel()
	MeshSolver.super.refreshCalcDTKernel(self)
	-- TODO combine these, and offset one into the other?
	-- because I'm going to need more than just these...
	self.calcDTKernelObj.obj:setArg(3, self.cellBuf)
	self.calcDTKernelObj.obj:setArg(4, self.faceBuf)
	self.calcDTKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
end

function MeshSolver:calcDT()
	if not self.useFixedDT then
		self.calcDTKernelObj.obj:setArg(3, self.cellBuf)
		self.calcDTKernelObj.obj:setArg(4, self.faceBuf)
		self.calcDTKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
	end
	return MeshSolver.super.calcDT(self)
end

-- same as hydro/solver/fvsolver without PLM or CTU
-- TODO should MeshSolver only be a finite volume, or should we support finite difference as well?
-- or should we only support finite difference via flux calculations (like I do in fdsolver)?
function MeshSolver:calcDeriv(derivBufObj, dt)
	local dtArg = real(dt)

	self.calcFluxKernelObj.obj:setArg(3, dtArg)
	self.calcFluxKernelObj()

	self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFromFluxKernelObj()
end

function MeshSolver:setBoundaryMethods()
	print'TODO MeshSolver:setBoundaryMethods'
end

function MeshSolver:boundary()
end

-- hmm, if something else require's solverbase and uses SolverBase.DisplayVar before this does 
-- then won't it get the wrong class? (compared to if this require's first before it does?)
local DisplayVar = MeshSolver.DisplayVar
local MeshSolverDisplayVar = DisplayVar:subclass()
MeshSolver.DisplayVar = MeshSolverDisplayVar

function MeshSolverDisplayVar:setArgs(kernel)
	MeshSolverDisplayVar.super.setArgs(self, kernel)
	kernel:setArg(7, self.solver.faceBuf)
end

function MeshSolver:updateGUIParams()
	MeshSolver.super.updateGUIParams(self)

	ig.luatableTooltipCheckbox('show vertexes', self, 'showVertexes')
	ig.igSameLine()
	ig.luatableTooltipCheckbox('show faces', self, 'showFaces')
	ig.igSameLine()
	ig.luatableTooltipCheckbox('show normals', self, 'showNormals')
	ig.igSameLine()
	ig.luatableTooltipCheckbox('show cell values', self, 'showValues')
	
	ig.luatableTooltipInputFloatAsText('cell scale', self, 'drawCellScale')
end

return MeshSolver 
