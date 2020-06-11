--[[
foundation class for solvers that use unstructured meshes from files

NACA 0012 from...
https://turbmodels.larc.nasa.gov/naca0012_grids.html
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local vec3sz = require 'vec-ffi.vec3sz'
local vec3f = require 'vec-ffi.vec3f'
local vec3d = require 'vec-ffi.vec3d'
local matrix_ffi = require 'matrix.ffi'
local ig = require 'ffi.imgui'
local gl = require 'gl'
local glreport = require 'gl.report'
local GLProgram = require 'gl.program'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLVertexArray = require 'gl.vertexarray'
local tooltip = require 'hydro.tooltip'
local SolverBase = require 'hydro.solver.solverbase'
local time, getTime = table.unpack(require 'hydro.util.time')
local real = require 'hydro.real'
local vector = require 'hydro.util.vector'

matrix_ffi.real = 'float'	-- default matrix_ffi type


local MeshSolver = class(SolverBase)

--[[
args:
	meshfile = name of mesh file to use

NOTICE initState is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:initL1(args)
	self.showVertexes = false
	self.showFaces = false
	self.showCells = true
	self.drawCellScale = 1

	MeshSolver.super.initL1(self, args)


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
	time('creating mesh', function()
		self.mesh = meshFactory:createMesh(self)
	end)
for k,v in pairs(self.mesh.times) do
	print(k,v)
end
--os.exit()

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

	-- TODO put these in solver_t
	self.numCells = assert(self.mesh.numCells, "did your MeshFactory remember to call :calcAux()?")
	self.numFaces = assert(self.mesh.numFaces, "did your MeshFactory remember to call :calcAux()?")
	self.numVtxs = assert(self.mesh.numVtxs, "did your MeshFactory remember to call :calcAux()?")
	self.numCellFaceIndexes = self.mesh.numCellFaceIndexes
	self.numCellVtxIndexes = self.mesh.numCellVtxIndexes
	self.numFaceVtxIndexes = self.mesh.numFaceVtxIndexes

	-- no longer is dim * numCells the number of interfaces -- it is now dependent on the mesh
	-- maybe I should rename this to numFaces?
	
	local solver = self
	local Program = class(require 'cl.obj.program')
	if ffi.os == 'Windows' then
		os.execute'mkdir cache-cl 2> nul'
	else
		os.execute'mkdir cache-cl 2> /dev/null'
	end
	function Program:init(args)
		args.env = solver.app.env
		args.domain = solver.domain
		args.cacheFile = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
		Program.super.init(self, args)
	end
	self.Program = Program
end

function MeshSolver:initDraw()
	--[[
	next issue ... how to render this?
	how to store it?
	where do textures come into play?
	
	UBuf will be size of # of elements
	
	how about texture, and how about rendering?
	texsize will have to be > #elems ... using size is fine as long as we need a gridsize
	--]]

	if self.app.targetSystem == 'console' then return end

	local drawShaderCode = assert(file['hydro/draw/mesh_heatmap.shader'])
	self.drawShader = GLProgram{
		vertexCode = template(drawShaderCode, {
			vertexShader = true,
			solver = self,
		}),
		fragmentCode = template(drawShaderCode, {
			fragmentShader = true,
			solver = self,
		}),
		uniforms = {
			tex = 0,
			gradientTex = 1,
			drawCellScale = 1,
		},
	}
	
	self.modelViewMatrix = matrix_ffi.zeros(4,4)
	self.projectionMatrix = matrix_ffi.zeros(4,4)
	self.modelViewProjectionMatrix = matrix_ffi.zeros(4,4)

	-- make a GPU version of the mesh
	-- TODO triangulate ... then include a mapping from triangle to source cell,
	-- and then just copy from display buf to a texture,
	-- then store per vertex the lookup to the texture.
	-- Also, if you are going to scale cells, then you must store a unique vertex per cell here.
	-- this means no dynamic mesh (without more coding).
	local glvtxs = vector'vec3f_t'			-- vertex position
	local glvtxcenters = vector'vec3f_t'	-- center of cell for this vertex
	local glcellindex = vector'float'		-- 0-based index of cell for this vertex
	local function addTri(va,vb,vc, ci,c)
		glvtxs:push_back(vec3f(va:unpack()))
		glvtxs:push_back(vec3f(vb:unpack()))
		glvtxs:push_back(vec3f(vc:unpack()))
		for j=0,2 do
			glvtxcenters:push_back(vec3f(c.pos:unpack()))
			glcellindex:push_back(ci)
		end
	end
	local mesh = self.mesh
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

	self.glvtxs = glvtxs
	self.glvtxcenters = glvtxcenters
	self.glcellindex = glcellindex
	self.numGlVtxs = #glvtxs


	self.glVtxBuffer = GLArrayBuffer{
		data = glvtxs.v,
		size = #glvtxs * ffi.sizeof(glvtxs.type), 
	}
	self.glVtxCenterBuffer = GLArrayBuffer{
		data = glvtxcenters.v,
		size = #glvtxcenters * ffi.sizeof(glvtxcenters.type), 
	}
	self.glCellIndexBuffer = GLArrayBuffer{
		data = glcellindex.v,
		size = #glcellindex * ffi.sizeof(glcellindex.type),
	}


-- [[ uses VertexArray
	self.drawShaderVtxVertexArray = GLVertexArray{
		buffer = self.glVtxBuffer,
		size = 3,
		type = gl.GL_FLOAT,
	}
	self.drawShaderVtxVertexArray:bind()
	self.drawShader:setAttr('vtx', self.drawShaderVtxVertexArray)
	self.drawShaderVtxVertexArray:unbind()

	self.drawShaderVtxCenterVertexArray = GLVertexArray{
		buffer = self.glVtxCenterBuffer,
		size = 3,
		type = gl.GL_FLOAT,
	}
	self.drawShaderVtxCenterVertexArray:bind()
	self.drawShader:setAttr('vtxcenter', self.drawShaderVtxCenterVertexArray)
	self.drawShaderVtxCenterVertexArray:unbind()

	self.drawShaderCellIndexVertexArray = GLVertexArray{
		buffer = self.glCellIndexBuffer,
		size = 1,
		type = gl.GL_FLOAT,
	}
	self.drawShaderCellIndexVertexArray:bind()
	self.drawShader:setAttr('cellindex', self.drawShaderCellIndexVertexArray)
	self.drawShaderCellIndexVertexArray:unbind()
--]]

-- [[ doesn't use VertexArrays, but uses EnableVertexAttribArray instead
	self.glVtxBuffer:bind()
	self.drawShaderVtxVertexArray:setAttr(self.drawShader.attrs.vtx.loc)
	self.glVtxBuffer:unbind()
	self.glVtxCenterBuffer:bind()
	self.drawShaderVtxCenterVertexArray:setAttr(self.drawShader.attrs.vtxcenter.loc)
	self.glVtxCenterBuffer:unbind()
	self.glCellIndexBuffer:bind()
	self.drawShaderCellIndexVertexArray:setAttr(self.drawShader.attrs.cellindex.loc)
	self.glCellIndexBuffer:unbind()
--]]
end

-- TODO organize this between SolverBase and MeshSolver
function MeshSolver:createBuffers()
	local app = self.app
	
	MeshSolver.super.createBuffers(self)

	self:clalloc('vtxBuf', 'real3', self.numVtxs)
	self:clalloc('cellsBuf', 'cell_t', self.numCells)
	self:clalloc('facesBuf', 'face_t', self.numFaces)
	self:clalloc('cellFaceIndexesBuf', 'int', self.numCellFaceIndexes)
	
	-- specific to FiniteVolumeSolver
	self:clalloc('fluxBuf', self.eqn.cons_t, self.numFaces)

	-- NOTICE this all looks a lot like GridSolver's createBuffers
	if app.targetSystem ~= 'console' then

		local maxTex2DSize = vec3sz(
			self.device:getInfo'CL_DEVICE_IMAGE2D_MAX_WIDTH',
			self.device:getInfo'CL_DEVICE_IMAGE2D_MAX_HEIGHT',
			1)
		local maxTex3DSize = vec3sz(
			self.device:getInfo'CL_DEVICE_IMAGE3D_MAX_WIDTH',
			self.device:getInfo'CL_DEVICE_IMAGE3D_MAX_HEIGHT',
			self.device:getInfo'CL_DEVICE_IMAGE3D_MAX_DEPTH')
		self.texSize = vec3sz()
		-- TODO if texSize >= max gl size then overflow into the next dim
		if self.numCells <= maxTex2DSize.x then
			self.texSize = vec3sz(self.numCells, 1, 1)
		else
			local sx = math.min(math.ceil(math.sqrt(self.numCells)), maxTexSize2D.x)
			local sy = math.ceil(self.numCells / tonumber(self.texSize.x))
			if sx <= maxTex2DSize.x and sy <= maxTe2DSize.y then
				self.texSize = vec3sz(sx, sy, 1)
			else
				local sz = math.min(math.ceil(math.cbrt(self.numCells)), maxTexSize3D.z)
				local sxy = math.ceil(self.numCells / sz)
				local sy = math.min(math.ceil(math.sqrt(sxy)), maxTexSize3D.y)
				local sx = math.ceil(sxy / sy)
				if sx >= maxTexSize3D.x then
					error("couldn't fit cell buffer into texture.  max 2d size " .. maxTex2DSize .. ", max 3d size " .. maxTex3DSize)
				end
				self.texSize = vec3sz(sx, sy, sz)
			end
		end


		local GLTex2D = require 'gl.tex2d'
		local GLTex3D = require 'gl.tex3d'
		local cl = self.texSize.z == 1 and GLTex2D or GLTex3D
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
			self.calcDisplayVarToTexPtr = ffi.new(app.real..'[?]', self.texSize:volume() * 3)
		end
	end
end

function MeshSolver:finalizeCLAllocs()
	MeshSolver.super.finalizeCLAllocs(self)

	self.vtxBufObj:fromCPU(self.mesh.vtxs.v)
	self.cellsBufObj:fromCPU(self.mesh.cells.v)
	self.facesBufObj:fromCPU(self.mesh.faces.v)
	self.cellFaceIndexesBufObj:fromCPU(self.mesh.cellFaceIndexes.v)
	self.cmds:finish()
end

function MeshSolver:checkStructSizes_getTypes()
	local typeinfos = MeshSolver.super.checkStructSizes_getTypes(self)
	typeinfos:append{
		face_t,
		cell_t,
	}
	return typeinfos
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

function MeshSolver:resetState()
	MeshSolver.super.resetState(self)
--self:printBuf(self.UBufObj)
end

function MeshSolver:refreshSolverProgram()
	MeshSolver.super.refreshSolverProgram(self)

	self.calcFluxKernelObj = self.solverProgramObj:kernel{name='calcFlux', domain=self.faceDomain}
	self.calcFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcFluxKernelObj.obj:setArg(1, self.fluxBuf)
	self.calcFluxKernelObj.obj:setArg(2, self.UBuf)
	self.calcFluxKernelObj.obj:setArg(4, self.cellsBuf)
	self.calcFluxKernelObj.obj:setArg(5, self.facesBuf)
	self.calcFluxKernelObj.obj:setArg(6, self.cellFaceIndexesBuf)

	self.calcDerivFromFluxKernelObj = self.solverProgramObj:kernel'calcDerivFromFlux'
	self.calcDerivFromFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(3, self.cellsBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(4, self.facesBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
end

function MeshSolver:getSolverCode()
	return table{
		MeshSolver.super.getSolverCode(self),
	
		-- TODO flux scheme HERE
		-- this is the 'calcFluxKernel' code part of all the fvsolver subclasses
		-- so maybe make that its own module that plugs into both fvsolver and this
		template(file['hydro/solver/roe.cl'], {
			solver = self,
			eqn = self.eqn,
		}),

		-- finite volume integration
		template(file['hydro/solver/calcDerivFV.cl'], {
			solver = self,
			eqn = self.eqn,
		}),
	}:concat'\n'
end



function MeshSolver:createCodePrefix()
	MeshSolver.super.createCodePrefix(self)

	local lines = table{
		self.codePrefix,
		self.mesh.meshTypeCode,
	}

	lines:insert[[
#define OOB(lhs,rhs) (index >= get_global_size(0))
#define SETBOUNDS(lhs,rhs)	\
	int index = get_global_id(0); \
	int4 i = (int4)(index,0,0,0);	\
	if (OOB(0,0)) return;
#define SETBOUNDS_NOGHOST()	\
	int index = get_global_id(0); \
	int4 i = (int4)(index,0,0,0);
#define cell_x(i) real3_zero
]]

	self.codePrefix = lines:concat'\n'
end

function MeshSolver:refreshCalcDTKernel()
	MeshSolver.super.refreshCalcDTKernel(self)
	-- TODO combine these, and offset one into the other?
	-- because I'm going to need more than just these...
	self.calcDTKernelObj.obj:setArg(3, self.cellsBuf)
	self.calcDTKernelObj.obj:setArg(4, self.facesBuf)
	self.calcDTKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
end

function MeshSolver:calcDT()
	if not self.useFixedDT then
		self.calcDTKernelObj.obj:setArg(3, self.cellsBuf)
		self.calcDTKernelObj.obj:setArg(4, self.facesBuf)
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

function MeshSolver:boundary()
end

-- hmm, if something else require's solverbase and uses SolverBase.DisplayVar before this does 
-- then won't it get the wrong class? (compared to if this require's first before it does?)
local DisplayVar = MeshSolver.DisplayVar
local MeshSolverDisplayVar = class(DisplayVar)
MeshSolver.DisplayVar = MeshSolverDisplayVar

function MeshSolverDisplayVar:setArgs(kernel)
	MeshSolverDisplayVar.super.setArgs(self, kernel)
	kernel:setArg(4, self.solver.cellsBuf)
	kernel:setArg(5, self.solver.facesBuf)
end

function MeshSolver:calcDisplayVarToTex(var, componentIndex)
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
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex
		
		local channels = vectorField and 3 or 1
		local format = vectorField and gl.GL_RGB or gl.GL_RED
	
		var:setToBufferArgs()
		
		self:calcDisplayVarToBuffer(var)

		local sizevec = self.texSize
		local volume = tonumber(sizevec:volume())
		
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

-- TODO move this into draw/2d_heatmap as a meshsolver pathway?
-- FPS of a 64x64 quad mesh:
-- sys=console: ~900
-- sys=imgui: ~170
function MeshSolver:display(varName, ar)
	local app = self.app
	if app.targetSystem == 'console' then return end	

	local var = self.displayVarForName[varName]
	if not var then return end

-- [[ this is from the draw/*.lua
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = self:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end
--]]	

	self:calcDisplayVarToTex(var)

-- TODO use app:displayVector, but that needs meshsolver displayTex first
--if vectorField then return end
	
	local view = app.view
	view:projection(ar)
	view:modelview()

	gl.glEnable(gl.GL_DEPTH_TEST)
--	gl.glEnable(gl.GL_CULL_FACE)
	
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projectionMatrix.ptr)
	self.modelViewProjectionMatrix:mul(self.projectionMatrix, self.modelViewMatrix)

	-- draw the mesh 
	local mesh = self.mesh

	if self.showVertexes then
		gl.glPointSize(3)
		gl.glColor3f(1,1,1)
		gl.glBegin(gl.GL_POINTS)
		for _,v in ipairs(mesh.vtxs) do
			gl.glVertex3d(v:unpack())
		end
		gl.glEnd()
		gl.glPointSize(1)
	end

	if self.showFaces then
		for fi,f in ipairs(mesh.faces) do
			gl.glBegin(gl.GL_LINE_LOOP)
			for vi=0,f.vtxCount-1 do
				local v = mesh.vtxs.v[mesh.faceVtxIndexes.v[vi + f.vtxOffset]]
				gl.glVertex3d(v:unpack())
			end
			gl.glEnd()
		end
	end
	
	-- if show cells ...
	-- TODO in my mesh solver, here I pick the color according to the display value
	-- but in my grid solvers I just overlay the whole thing with a giant texture
	-- so how to combine the two? use gradient tex.  
	-- in gridshader I have a kernel for updating tex/buffer to display values based on the cell state values
	-- i can do the same, and just update the 'displayValue' inside each cell for this.
	-- in fact, before drawing this, calcDisplayVarToBuffer will fill reduceBuf with the values (associated 1-1 with each cell?)
	if self.showCells then
		
		local gradientTex = app.gradientTex
		self.drawShader:use()
		gl.glUniform1i(self.drawShader.uniforms.useLog.loc, 0)
		gl.glUniform1f(self.drawShader.uniforms.valueMin.loc, valueMin)
		gl.glUniform1f(self.drawShader.uniforms.valueMax.loc, valueMax)
		gl.glUniform1f(self.drawShader.uniforms.drawCellScale.loc, self.drawCellScale)
		gl.glUniformMatrix4fv(self.drawShader.uniforms.modelViewProjectionMatrix.loc, 1, 0, self.modelViewProjectionMatrix.ptr)
		self.tex:bind(0)
		gradientTex:bind(1)
--[[ doesn't work yet
		self.drawShaderVtxVertexArray:bind()
		self.drawShaderVtxCenterVertexArray:bind()
		self.drawShaderValueVertexArray:bind()
		gl.glDrawArrays(gl.GL_TRIANGLES, 0, self.numGlVtxs * 3)
		self.drawShaderVtxVertexArray:unbind()
		self.drawShaderVtxCenterVertexArray:unbind()
		self.drawShaderValueVertexArray:unbind()
--]]
-- [[ 180 fps @ 64x64 quad grid .. TODO has bugs now and doesn't render correctly.
		gl.glEnableVertexAttribArray(self.drawShader.attrs.cellindex.loc)
		gl.glEnableVertexAttribArray(self.drawShader.attrs.vtxcenter.loc)
		gl.glEnableVertexAttribArray(self.drawShader.attrs.vtx.loc)
		gl.glDrawArrays(gl.GL_TRIANGLES, 0, self.numGlVtxs * 3)
		gl.glDisableVertexAttribArray(self.drawShader.attrs.cellindex.loc)
		gl.glDisableVertexAttribArray(self.drawShader.attrs.vtxcenter.loc)
		gl.glDisableVertexAttribArray(self.drawShader.attrs.vtx.loc)
--]]
--[[ 115 fps @ 64x64 quad grid
		gl.glBegin(gl.GL_TRIANGLES)
		for i=0,self.numGlVtxs-1 do
			gl.glVertexAttrib1f(self.drawShader.attrs.cellindex.loc, self.glcellindex.v[i])
			gl.glVertexAttrib3f(self.drawShader.attrs.vtxcenter.loc, self.glvtxcenters.v[i]:unpack())
			gl.glVertexAttrib3f(self.drawShader.attrs.vtx.loc, self.glvtxs.v[i]:unpack())
		end
		gl.glEnd()
--]]
		gradientTex:unbind(1)
		self.tex:unbind(0)
		self.drawShader:useNone()
	end

	gl.glDisable(gl.GL_DEPTH_TEST)
--	gl.glDisable(gl.GL_CULL_FACE)

-- [[ also in draw/*.lua
	local gradientValueMin = valueMin
	local gradientValueMax = valueMax
	local showName = varName
	if var.showInUnits and var.units then
		local unitScale = self:convertToSIUnitsCode(var.units).func()
		gradientValueMin = gradientValueMin * unitScale
		gradientValueMax = gradientValueMax * unitScale
		showName = showName..' ('..var.units..')'
	end
	app:drawGradientLegend(ar, showName, gradientValueMin, gradientValueMax)
--]]

	glreport'here'
end

function MeshSolver:updateGUIParams()
	MeshSolver.super.updateGUIParams(self)

	tooltip.checkboxTable('show vertexes', self, 'showVertexes')
	ig.igSameLine()
	tooltip.checkboxTable('show faces', self, 'showFaces')
	ig.igSameLine()
	tooltip.checkboxTable('show cells', self, 'showCells')
	
	tooltip.numberTable('cell scale', self, 'drawCellScale')
end

return MeshSolver 
