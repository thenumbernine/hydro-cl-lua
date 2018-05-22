local gl = require 'ffi.OpenGL'
local table = require 'ext.table'
local CartesianGeom = require 'geom.cartesian'

return function(HydroCLApp)

-- 2D
local vertexesInQuad = {{0,0},{1,0},{1,1},{0,1}}

--[[
looks great for flat space
TODO for curved space: provide a coordMapInv function (might have to be manual to account for domains of rotations)
 and then call this as we march through volumes 
 and treat out-of-bound values as fully transparent
--]]
HydroCLApp.display3D_Slice_usePoints = false
HydroCLApp.display3D_Slice_useIsos = true
HydroCLApp.display3D_Slice_numIsobars = 20
HydroCLApp.display3D_Slice_useLighting = false
HydroCLApp.display3D_Slice_alpha = .15
HydroCLApp.display3D_Slice_alphaGamma = 1
HydroCLApp.display3D_Slice_numSlices = 255
function HydroCLApp:display3D_Slice(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	for _,solver in ipairs(solvers) do 
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			
			self.view:projection(ar)
			self.view:modelview()

if useClipPlanes then
	for i,clipInfo in ipairs(clipInfos) do
		gl.glClipPlane(gl.GL_CLIP_PLANE0+i-1, clipInfo.plane:ptr())
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

			solver.volumeSliceShader:use()
			solver:getTex(var):bind(0)
			self.gradientTex:bind(1)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.alpha.loc, self.display3D_Slice_alpha)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.alphaGamma.loc, self.display3D_Slice_alphaGamma)
			gl.glUniform3f(solver.volumeSliceShader.uniforms.mins.loc, solver.mins:unpack())
			gl.glUniform3f(solver.volumeSliceShader.uniforms.maxs.loc, solver.maxs:unpack())
			gl.glUniform1i(solver.volumeSliceShader.uniforms.useLog.loc, var.useLog)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.valueMin.loc, valueMin)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.valueMax.loc, valueMax)
			gl.glUniform1i(solver.volumeSliceShader.uniforms.useIsos.loc, self.display3D_Slice_useIsos)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.numIsobars.loc, self.display3D_Slice_numIsobars)
			gl.glUniform1i(solver.volumeSliceShader.uniforms.useLighting.loc, self.display3D_Slice_useLighting)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.numGhost.loc, solver.numGhost)
			gl.glUniform3f(solver.volumeSliceShader.uniforms.texSize.loc, solver.gridSize:unpack())

if useClipPlanes then
			for i,info in ipairs(clipInfos) do
				gl.glUniform1i(solver.volumeSliceShader.uniforms['clipEnabled'..i].loc, info.enabled and 1 or 0)
			end
end

			if self.display3D_Slice_usePoints then
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

				local n = self.display3D_Slice_numSlices
				local fwd = -self.frustumView.angle:conjugate():zAxis()
				local fwddir = select(2, table(fwd):map(math.abs):sup())

				local jmin, jmax, jdir
				if fwd[fwddir] < 0 then
					jmin, jmax, jdir = 0, n, 1
				else
					jmin, jmax, jdir = n, 0, -1
				end
				gl.glUniform3f(solver.volumeSliceShader.uniforms.normal.loc, 
					fwddir == 1 and jdir or 0, 
					fwddir == 2 and jdir or 0, 
					fwddir == 3 and jdir or 0)
						
				if CartesianGeom.is(solver.geometry) then
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

			self.gradientTex:unbind(1)
			solver:getTex(var):unbind(0)
			solver.volumeSliceShader:useNone()
		
			self:drawGradientLegend(ar, varName, valueMin, valueMax)
		end
	end
end

end
