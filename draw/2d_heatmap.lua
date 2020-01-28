local gl = require 'ffi.OpenGL'
local class = require 'ext.class'


local Draw2DHeatmap = class()

function Draw2DHeatmap:drawSolverWithVar(app, solver, var, heatMap2DShader, xmin, xmax, ymin, ymax)
-- hmm ... this is needed for sub-solvers
local origSolver = var.solver
var.solver = solver
	
	solver:calcDisplayVarToTex(var)
	
	gl.glUniform2f(heatMap2DShader.uniforms.solverMins.loc, solver.mins[1], solver.mins[2])
	gl.glUniform2f(heatMap2DShader.uniforms.solverMaxs.loc, solver.maxs[1], solver.maxs[2])

	local tex = solver:getTex(var)
	local size = var.getBuffer().sizevec or solver.gridSize
	gl.glUniform2f(heatMap2DShader.uniforms.texCoordMax.loc, 
		tonumber(size.x) / tex.width,
		tonumber(size.y) / tex.height)
	tex:bind(0)
	if app.displayBilinearTextures then
		gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
	else
		gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST)
	end

	gl.glBegin(gl.GL_QUADS)
	gl.glVertex2d(xmin, ymin)
	gl.glVertex2d(xmax, ymin)
	gl.glVertex2d(xmax, ymax)
	gl.glVertex2d(xmin, ymax)
	gl.glEnd()
	
	tex:unbind(0)

var.solver = origSolver
end

function Draw2DHeatmap:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax)
	-- TODO allow a fixed, manual colormap range
	-- NOTICE with AMR this will only get from the root node
	--  which should at least have blitters of the children
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end
	
	local heatMap2DShader = solver:getHeatMap2DShader(var)
	heatMap2DShader:use()
	app.gradientTex:bind(1)

	gl.glUniform1i(heatMap2DShader.uniforms.useCoordMap.loc, app.display_useCoordMap)
	gl.glUniform1i(heatMap2DShader.uniforms.useLog.loc, var.useLog)
	gl.glUniform1f(heatMap2DShader.uniforms.valueMin.loc, valueMin)
	gl.glUniform1f(heatMap2DShader.uniforms.valueMax.loc, valueMax)
	
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)
	
	self:drawSolverWithVar(app, solver, var, heatMap2DShader, xmin, xmax, ymin, ymax)

-- [[
	if solver.amr then
		for k,subsolver in pairs(solver.amr.child) do
			self:drawSolverWithVar(app, subsolver, var, heatMap2DShader, xmin, xmax, ymin, ymax)
		end
	end
--]]

	gl.glDisable(gl.GL_BLEND)

	app.gradientTex:unbind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)
	heatMap2DShader:useNone()

--	gl.glDisable(gl.GL_DEPTH_TEST)

	-- TODO only draw the first
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

--	gl.glEnable(gl.GL_DEPTH_TEST)
end

function Draw2DHeatmap:display(app, solvers, varName, ar, graph_xmin, graph_ymin, graph_xmax, graph_ymax)
	app.view:projection(ar)
	app.view:modelview()
	if app.view.getOrthoBounds then
		xmin, xmax, ymin, ymax = app.view:getOrthoBounds(ar)
	else
		xmin, xmax, ymin, ymax = graph_xmin, graph_ymin, graph_xmax, graph_ymax
	end

--	gl.glEnable(gl.GL_DEPTH_TEST)
	
	local gridz = 0	--.1

	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex3f(x*xstep,ymin, gridz)
		gl.glVertex3f(x*xstep,ymax, gridz)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex3f(xmin,y*ystep, gridz)
		gl.glVertex3f(xmax,y*ystep, gridz)
	end
	gl.glEnd()
	
	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex3f(xmin, 0, gridz)
	gl.glVertex3f(xmax, 0, gridz)
	gl.glVertex3f(0, ymin, gridz)
	gl.glVertex3f(0, ymax, gridz)
	gl.glEnd()
			
	-- NOTICE overlays of multiple solvers won't be helpful.  It'll just draw over the last solver.
	-- I've got to rethink the visualization
	for _,solver in ipairs(solvers) do 
		if not require 'solver.meshsolver'.is(solver) then
			local var = solver.displayVarForName[varName]
			if var and var.enabled then
				self:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax)
			end
		end
	end
--	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:display2D_Heatmap(...)
		if not self.draw2DHeatmap then self.draw2DHeatmap = Draw2DHeatmap() end
		return self.draw2DHeatmap:display(self, ...)
	end
end
