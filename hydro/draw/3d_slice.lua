local table = require 'ext.table'
local class = require 'ext.class'
local file = require 'ext.file'
local template = require 'template'
local gl = require 'ffi.OpenGL'
local CartesianCoordinateSystem = require 'hydro.coord.cartesian'
local Draw = require 'hydro.draw.draw'


local Draw3DSlice = class(Draw)

-- 2D
local vertexesInQuad = {{0,0},{1,0},{1,1},{0,1}}

--[[
looks great for flat space
TODO for curved space: provide a coordMapInv function (might have to be manual to account for domains of rotations)
 and then call this as we march through volumes 
 and treat out-of-bound values as fully transparent
--]]
Draw3DSlice.usePoints = false
Draw3DSlice.useIsos = true
Draw3DSlice.numIsobars = 20
Draw3DSlice.useLighting = false
Draw3DSlice.alpha = .15
Draw3DSlice.alphaGamma = 1
Draw3DSlice.numSlices = 255


function Draw3DSlice:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax, useLog)
	if require 'hydro.solver.meshsolver'.is(solver) then return end

	app.view:setup(ar)

	if useClipPlanes then
		for i,clipInfo in ipairs(clipInfos) do
			gl.glClipPlane(gl.GL_CLIP_PLANE0+i-1, clipInfo.plane.s)
		end
	end
	
	
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

	local shader = solver.volumeSliceShader
	local uniforms = shader.uniforms

	shader:use()
	local tex = solver:getTex(var)
	tex:bind(0)
	if app.displayBilinearTextures then
		gl.glTexParameteri(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
	else
		gl.glTexParameteri(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST)
	end
	
	app.gradientTex:bind(1)
	
	self:setupDisplayVarShader(shader, app, solver, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.alpha.loc, self.alpha)
	gl.glUniform1f(uniforms.alphaGamma.loc, self.alphaGamma)
	gl.glUniform3f(uniforms.solverMins.loc, solver.mins:unpack())
	gl.glUniform3f(uniforms.solverMaxs.loc, solver.maxs:unpack())
	gl.glUniform1i(uniforms.useIsos.loc, self.useIsos)
	gl.glUniform1f(uniforms.numIsobars.loc, self.numIsobars)
	gl.glUniform1i(uniforms.useLighting.loc, self.useLighting)
	gl.glUniform1f(uniforms.numGhost.loc, solver.numGhost)
	gl.glUniform3f(uniforms.texSize.loc, solver.texSize:unpack())

	if useClipPlanes then
		for i,info in ipairs(clipInfos) do
			gl.glUniform1i(uniforms['clipEnabled'..i].loc, info.enabled and 1 or 0)
		end
	end

	if self.usePoints then
		gl.glEnable(gl.GL_DEPTH_TEST)
		gl.glPointSize(2)
		gl.glBegin(gl.GL_POINTS)
		local numGhost = solver.numGhost
		for i=numGhost+1,tonumber(solver.gridSize.x-numGhost) do
			for j=numGhost+1,tonumber(solver.gridSize.y-numGhost) do
				for k=numGhost+1,tonumber(solver.gridSize.z-numGhost) do
					gl.glVertex3d(
						(i - numGhost - .5)/tonumber(solver.gridSize.x - 2*numGhost),
						(j - numGhost - .5)/tonumber(solver.gridSize.y - 2*numGhost),
						(k - numGhost - .5)/tonumber(solver.gridSize.z - 2*numGhost))
				end
			end
		end
		gl.glEnd()
		gl.glDisable(gl.GL_DEPTH_TEST)
	
	else
	
		gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
		gl.glEnable(gl.GL_BLEND)

		local n = self.numSlices
		local fwd = -app.frustumView.angle:conjugate():zAxis()
		
		local fwddir
		local jmin, jmax, jdir


		-- hack for picking order of axis for non-Cartesian
		if app.display_useCoordMap
		and (require 'hydro.coord.sphere'.is(solver.coord) 
			or require 'hydro.coord.sphere-log-radial'.is(solver.coord))
		then
			fwddir = 1
			jmin, jmax, jdir = 0, n, 1
			gl.glUniform3f(uniforms.normal.loc, (-fwd):unpack())
		else
		
			fwddir = select(2, table{fwd:unpack()}:map(math.abs):sup())

			if fwd.s[fwddir-1] < 0 then
				jmin, jmax, jdir = 0, n, 1
			else
				jmin, jmax, jdir = n, 0, -1
			end
			
			gl.glUniform3f(uniforms.normal.loc, 
				fwddir == 1 and jdir or 0, 
				fwddir == 2 and jdir or 0, 
				fwddir == 3 and jdir or 0)
		end

		if CartesianCoordinateSystem.is(solver.coord) then
			-- [[	single quad
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
			--]]
		else
			-- [[	use a grid, so curved coordinates can be seen
			for j=jmin,jmax,jdir do
				local f = j/n
				local xres = 20
				local yres = 20
				for ybase=1,yres-1 do
					gl.glBegin(gl.GL_TRIANGLE_STRIP)
					for x=1,xres do
						for y=ybase,ybase+1 do
							if fwddir == 1 then
								gl.glVertex3f(f, x/xres, y/yres)
							elseif fwddir == 2 then
								gl.glVertex3f(x/xres, f, y/yres)
							elseif fwddir == 3 then
								gl.glVertex3f(x/xres, y/yres, f)
							end
						end
					end
					gl.glEnd()
				end
			end
			--]]
		end
	
		gl.glDisable(gl.GL_BLEND)
	end

	app.gradientTex:unbind(1)
	tex:unbind(0)
	shader:useNone()

	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

function Draw3DSlice:display(app, solvers, varName, ar, xmin, xmax, ymin, ymax, useLog)
	for _,solver in ipairs(solvers) do 
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:prepareShader(solver)
			self:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax, useLog)
		end
	end
end

function Draw3DSlice:prepareShader(solver)
	if solver.volumeSliceShader then return end 
	
	local volumeSliceCode = assert(file['hydro/draw/3d_slice.shader'])
	
	solver.volumeSliceShader = solver.GLProgram{
		name = '3d_slice',
		vertexCode = template(volumeSliceCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			vertexShader = true
		}),
		fragmentCode = template(volumeSliceCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			fragmentShader = true,
			-- TODO move this from app, or make it a field of app?
			clipInfos = useClipPlanes and clipInfos or nil,
		}),
		uniforms = {
			volTex = 0,
			gradientTex = 1,
			valueMin = 0,
			valueMax = 0,
		},
	}
end

return Draw3DSlice
