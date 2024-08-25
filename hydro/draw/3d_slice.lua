local table = require 'ext.table'
local path = require 'ext.path'
local assertgt = require 'ext.assert'.gt
local assertlen = require 'ext.assert'.len
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local vector = require 'ffi.cpp.vector-lua'
local CartesianCoordinateSystem = require 'hydro.coord.cartesian'
local Draw = require 'hydro.draw.draw'


local Draw3DSlice = Draw:subclass()

-- 2D
local vertexesInQuad = {{0,0},{1,0},{1,1},{0,1}}

--[[
looks great for flat space
TODO for curved space: provide a coordMapInv function (might have to be manual to account for domains of rotations)
 and then call this as we march through volumes
 and treat out-of-bound values as fully transparent
--]]
Draw3DSlice.usePoints = cmdline.displayPointCloud
Draw3DSlice.numIsobars = 20
Draw3DSlice.useLighting = false
Draw3DSlice.alpha = .15
Draw3DSlice.alphaGamma = 1
Draw3DSlice.numSlices = 255

Draw3DSlice.useIsos = true
if cmdline.isobars ~= nil then
	Draw3DSlice.useIsos = cmdline.isobars
end

function Draw3DSlice:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = solver.app
	local shader = solver.volumeSliceShader
	local uniforms = shader.uniforms
	if require 'hydro.solver.meshsolver':isa(solver) then return end

	app.view:setup(ar)

	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end

	solver:calcDisplayVarToTex(var)

	shader:use()
	local tex = solver:getTex(var)
		:bind(0)
		:setParameter(gl.GL_TEXTURE_MAG_FILTER, app.displayBilinearTextures and gl.GL_LINEAR or gl.GL_NEAREST)

	app.gradientTex:bind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.alpha.loc, self.alpha)
	gl.glUniform1f(uniforms.alphaGamma.loc, self.alphaGamma)
	gl.glUniform1i(uniforms.useIsos.loc, self.useIsos)
	gl.glUniform1f(uniforms.numIsobars.loc, self.numIsobars)
	gl.glUniform1i(uniforms.useLighting.loc, self.useLighting)

	if self.usePoints then
		local numGhost = solver.numGhost

		assertgt(solver.gridSize.x, 2 * numGhost)
		assertgt(solver.gridSize.y, 2 * numGhost)
		assertgt(solver.gridSize.z, 2 * numGhost)
		local numVtxs = (solver.gridSize - 2 * numGhost):volume()

		solver.draw3DSlicePtVtxs = solver.draw3DSlicePtVtxs or vector'vec3f_t'

		if not solver.draw3DSlicePtVtxBuf
		or #solver.draw3DSlicePtVtxs ~= numVtxs
		then
			solver.draw3DSlicePtVtxs:resize(numVtxs)
			solver.draw3DSlicePtVtxs:resize(0)

			for i=numGhost+1,tonumber(solver.gridSize.x-numGhost) do
				for j=numGhost+1,tonumber(solver.gridSize.y-numGhost) do
					for k=numGhost+1,tonumber(solver.gridSize.z-numGhost) do
						solver.draw3DSlicePtVtxs:emplace_back():set(
							(i - numGhost - .5)/tonumber(solver.gridSize.x - 2*numGhost),
							(j - numGhost - .5)/tonumber(solver.gridSize.y - 2*numGhost),
							(k - numGhost - .5)/tonumber(solver.gridSize.z - 2*numGhost))
					end
				end
			end
			assertlen(solver.draw3DSlicePtVtxs, numVtxs)

			--[[ TODO cuz I think we can't use draw without VAOs in GL core ...
			self.draw3DSlicePtVtxBuf = GLArrayBuffer{
				data = solver.draw3DSlicePtVtxs.v,
				size = #solver.draw3DSlicePtVtxs * ffi.sizeof(solver.draw3DSlicePtVtxs.type),
			}:unbind()
			--]]
		end

		gl.glEnable(gl.GL_DEPTH_TEST)
		gl.glPointSize(2)

		gl.glEnableVertexAttribArray(shader.attrs.vertex.loc)
		gl.glVertexAttribPointer(shader.attrs.vertex.loc, 3, gl.GL_FLOAT, false, 0, solver.draw3DSlicePtVtxs.v)
		gl.glDrawArrays(gl.GL_POINTS, 0, numVtxs)
		gl.glDisableVertexAttribArray(shader.attrs.vertex.loc)

		gl.glDisable(gl.GL_DEPTH_TEST)
	else

		gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
		gl.glEnable(gl.GL_BLEND)

		local n = self.numSlices
		local fwd = -app.frustumView.angle:zAxis()

		local fwddir
		local jmin, jmax, jdir

		fwddir = select(2, table{fwd:unpack()}:mapi(math.abs):sup())

		if fwd.s[fwddir-1] < 0 then
			jmin, jmax, jdir = 0, n, 1
		else
			jmin, jmax, jdir = n, 0, -1
		end

		gl.glUniform3f(uniforms.normal.loc,
			fwddir == 1 and jdir or 0,
			fwddir == 2 and jdir or 0,
			fwddir == 3 and jdir or 0)

		gl.glBegin(gl.GL_QUADS)
		for j=jmin,jmax,jdir do
			local f = j/n
			for _,vtx in ipairs(vertexesInQuad) do
				if fwddir == 1 then
					gl.glVertex3f(f, vtx[1], vtx[2])
				elseif fwddir == 2 then
					gl.glVertex3f(vtx[1], f, vtx[2])
				elseif fwddir == 3 then
					gl.glVertex3f(vtx[1], vtx[2], f)
				end
			end
		end
		gl.glEnd()

		gl.glDisable(gl.GL_BLEND)
	end

	app.gradientTex:unbind(1)
	tex:unbind(0)
	shader:useNone()

	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

function Draw3DSlice:display(varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = solver.app
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax, useLog)
	end
end

function Draw3DSlice:prepareShader()
	local solver = self.solver
	if solver.volumeSliceShader then return end
	local app = solver.app

	local volumeSliceCode = assert(path'hydro/draw/3d_slice.glsl':read())

	solver.volumeSliceShader = solver.GLProgram{
		name = '3d_slice',
		vertexCode = solver.eqn:template(volumeSliceCode, {
			draw = self,
			vertexShader = true
		}),
		fragmentCode = solver.eqn:template(volumeSliceCode, {
			draw = self,
			fragmentShader = true,
			-- TODO move this from app, or make it a field of app?
			clipInfos = app.useClipPlanes and app.clipInfos or nil,
		}),
		uniforms = {
			volTex = 0,
			gradientTex = 1,
			valueMin = 0,
			valueMax = 0,
		},
	}:useNone()
end

return Draw3DSlice
