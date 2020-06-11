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
local GLAttribute = require 'gl.attribute'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLVertexArray = require 'gl.vertexarray'
local tooltip = require 'hydro.tooltip'
local SolverBase = require 'hydro.solver.solverbase'
local time, getTime = table.unpack(require 'hydro.util.time')
local real = require 'hydro.real'
local vector = require 'hydro.util.vector'

matrix_ffi.real = 'float'	-- default matrix_ffi type


local MeshSolver = class(SolverBase)

-- for compat in some display stuff
MeshSolver.numGhost = 0 

--[[
args:
	meshfile = name of mesh file to use

NOTICE initState is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:initL1(args)
	self.showVertexes = false
	self.showFaces = false
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

function MeshSolver:initDraw()
	if self.app.targetSystem == 'console' then return end

	local vectorArrowCode = file['hydro/draw/vector_arrow.shader']
	self.vectorArrowShader = self.GLProgram{
		name = 'vector_arrow',
		vertexCode = template(vectorArrowCode, {
			solver = self,
			vertexShader = true,
		}),
		fragmentCode = template(vectorArrowCode, {
			solver = self,
			fragmentShader = true,
		}),
		uniforms = {
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
			scale = 1,
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

	
	local heatMapCode = assert(file['hydro/draw/mesh_heatmap.shader'])
	self.heatMap2DShader = self.GLProgram{
		vertexCode = template(heatMapCode, {
			vertexShader = true,
			solver = self,
		}),
		fragmentCode = template(heatMapCode, {
			fragmentShader = true,
			solver = self,
		}),
		uniforms = {
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
			drawCellScale = 1,
		},
		
		attrs = {
			vtx = GLArrayBuffer{
				data = glvtxs.v,
				size = #glvtxs * ffi.sizeof(glvtxs.type), 
			},
			vtxcenter = GLArrayBuffer{
				data = glvtxcenters.v,
				size = #glvtxcenters * ffi.sizeof(glvtxcenters.type), 
			},
			cellindex = GLArrayBuffer{
				data = glcellindex.v,
				size = #glcellindex * ffi.sizeof(glcellindex.type),
			},
		},
		createVAO = true,	-- tells shader to make VAO upon construction
	}
end

function MeshSolver:getHeatMap2DShader(var)
	return self.heatMap2DShader
end

-- TODO organize this between SolverBase and MeshSolver
function MeshSolver:createBuffers()
	local app = self.app

	-- set texSize before calling super
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
	end

	MeshSolver.super.createBuffers(self)

	self:clalloc('vtxBuf', 'real3', self.numVtxs)
	self:clalloc('cellsBuf', 'cell_t', self.numCells)
	self:clalloc('facesBuf', 'face_t', self.numFaces)
	self:clalloc('cellFaceIndexesBuf', 'int', self.numCellFaceIndexes)
	
	-- specific to FiniteVolumeSolver
	self:clalloc('fluxBuf', self.eqn.cons_t, self.numFaces)
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

function MeshSolver:updateGUIParams()
	MeshSolver.super.updateGUIParams(self)

	tooltip.checkboxTable('show vertexes', self, 'showVertexes')
	ig.igSameLine()
	tooltip.checkboxTable('show faces', self, 'showFaces')
	
	tooltip.numberTable('cell scale', self, 'drawCellScale')
end

-- TODO move this into draw/2d_heatmap as a meshsolver pathway

-- this is a mirror of Draw2DHeatmap:drawSolverWithVar
local function drawSolverWithVar(app, solver, var, heatMap2DShader)
	solver:calcDisplayVarToTex(var)

	local tex = solver:getTex()
	tex:bind(0)
	if app.displayBilinearTextures then
		gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
	else
		gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST)
	end

--[[ 110 fps: glVertexAttrib prim calls
	gl.glBegin(gl.GL_TRIANGLES)
	for i=0,solver.numGlVtxs-1 do
		gl.glVertexAttrib1f(heatMap2DShader.attrs.cellindex.loc, solver.glcellindex.v[i])
		gl.glVertexAttrib3f(heatMap2DShader.attrs.vtxcenter.loc, solver.glvtxcenters.v[i]:unpack())
		gl.glVertexAttrib3f(heatMap2DShader.attrs.vtx.loc, solver.glvtxs.v[i]:unpack())
	end
	gl.glEnd()
--]]
--[[ 130 fps: glBindBuffer / glVertexAttribPointer / glEnableVertexAttribArray
	solver.heatMap2DShader:setAttrs(solver.heatMapShaderAttrs)
	gl.glDrawArrays(gl.GL_TRIANGLES, 0, solver.numGlVtxs * 3)
--]]
-- [[ 150fps: glVertexArray
	solver.heatMap2DShader.vao:bind()
	gl.glDrawArrays(gl.GL_TRIANGLES, 0, solver.numGlVtxs * 3)
	solver.heatMap2DShader.vao:unbind()
--]]
	tex:unbind(0)
end
	
local function showDisplayVar(app, solver, var, varName, ar)
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end

--	gl.glEnable(gl.GL_DEPTH_TEST)
--	gl.glEnable(gl.GL_CULL_FACE)

	local heatMap2DShader = solver:getHeatMap2DShader()
	heatMap2DShader:use()
	
	local gradientTex = app.gradientTex
	gradientTex:bind(1)

	-- useCoordMap is missing from MeshSolver
	gl.glUniform1i(heatMap2DShader.uniforms.useLog.loc, 0)
	gl.glUniform1f(heatMap2DShader.uniforms.valueMin.loc, valueMin)
	gl.glUniform1f(heatMap2DShader.uniforms.valueMax.loc, valueMax)
	-- drawCellScale isn't present in GridSolver
	gl.glUniform1f(heatMap2DShader.uniforms.drawCellScale.loc, solver.drawCellScale)

	-- this is only in MeshSolver...
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, solver.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, solver.projectionMatrix.ptr)
	solver.modelViewProjectionMatrix:mul(solver.projectionMatrix, solver.modelViewMatrix)
	gl.glUniformMatrix4fv(heatMap2DShader.uniforms.modelViewProjectionMatrix.loc, 1, 0, solver.modelViewProjectionMatrix.ptr)
	
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	drawSolverWithVar(app, solver, var, heatMap2DShader)

	gl.glDisable(gl.GL_BLEND)

	gradientTex:unbind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)
	heatMap2DShader:useNone()

--	gl.glDisable(gl.GL_DEPTH_TEST)
--	gl.glDisable(gl.GL_CULL_FACE)

	local gradientValueMin = valueMin
	local gradientValueMax = valueMax
	local showName = varName
	if var.showInUnits and var.units then
		local unitScale = solver:convertToSIUnitsCode(var.units).func()
		gradientValueMin = gradientValueMin * unitScale
		gradientValueMax = gradientValueMax * unitScale
		showName = showName..' ('..var.units..')'
	end
	app:drawGradientLegend(ar, showName, gradientValueMin, gradientValueMax)
end

function MeshSolver:display(varName, ar)
	local app = self.app
	if app.targetSystem == 'console' then return end	
	
	local view = app.view
	view:projection(ar)
	view:modelview()

	local var = self.displayVarForName[varName]
	if not var then return end

	-- if it's a vector field then let app handle it.
	local component = self.displayComponentFlatList[var.component]
	if self:isVarTypeAVectorField(component.type) then return end

	if self.showVertexes then
		local mesh = self.mesh
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
		local mesh = self.mesh
		for fi,f in ipairs(mesh.faces) do
			gl.glBegin(gl.GL_LINE_LOOP)
			for vi=0,f.vtxCount-1 do
				local v = mesh.vtxs.v[mesh.faceVtxIndexes.v[vi + f.vtxOffset]]
				gl.glVertex3d(v:unpack())
			end
			gl.glEnd()
		end
	end

	-- from here on it's showDisplayVar
	showDisplayVar(app, self, var, varName, ar)
	glreport'here'
end

return MeshSolver 
