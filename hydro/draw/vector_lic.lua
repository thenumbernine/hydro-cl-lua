local gl = require 'ffi.OpenGL'
local ffi = require 'ffi'
local class = require 'ext.class'
local GLTex2D = require 'gl.tex2d'


local DrawVectorLIC = class()

DrawVectorLIC.integralMaxIter = 10

function DrawVectorLIC:drawSolverWithVar(app, solver, var, vectorLICShader, xmin, xmax, ymin, ymax)
-- hmm ... this is needed for sub-solvers
local origSolver = var.solver
var.solver = solver
	
	solver:calcDisplayVarToTex(var)
	
	gl.glUniform2f(vectorLICShader.uniforms.solverMins.loc, solver.mins.x, solver.mins.y)
	gl.glUniform2f(vectorLICShader.uniforms.solverMaxs.loc, solver.maxs.x, solver.maxs.y)

	local tex = solver:getTex(var)
	local size = var.getBuffer().sizevec or solver.texSize
	gl.glUniform2f(vectorLICShader.uniforms.texCoordMax.loc, 
		tonumber(size.x) / tex.width,
		tonumber(size.y) / tex.height)
	tex:bind(0)
	self.noiseTex:bind(2)
	if app.displayBilinearTextures then
		gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
	else
		gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST)
	end

	gl.glBegin(gl.GL_QUADS)
	gl.glVertex3d(xmin, ymin, app.displayFixedZ)
	gl.glVertex3d(xmax, ymin, app.displayFixedZ)
	gl.glVertex3d(xmax, ymax, app.displayFixedZ)
	gl.glVertex3d(xmin, ymax, app.displayFixedZ)
	gl.glEnd()
	
	tex:unbind(0)
	self.noiseTex:unbind(2)

var.solver = origSolver
end

function DrawVectorLIC:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax)
	-- TODO allow a fixed, manual colormap range
	-- NOTICE with AMR this will only get from the root node
	--  which should at least have blitters of the children
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		local component = solver.displayComponentFlatList[var.component]
		local vectorField = solver:isVarTypeAVectorField(component.type)
		if vectorField then
			-- calc range of magnitude of vector variable
			valueMin, valueMax = solver:calcDisplayVarRange(var, component.magn)
		else
			valueMin, valueMax = solver:calcDisplayVarRange(var)
		end
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end

	if not self.noiseTex then
		local noiseVol = app.drawVectorLICNoiseSize * app.drawVectorLICNoiseSize * 4
		local noiseData = ffi.new('float[?]', noiseVol)
		for i=0,noiseVol-1 do
			noiseData[i] = math.random()
		end
		self.noiseTex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,
			width = app.drawVectorLICNoiseSize,
			height = app.drawVectorLICNoiseSize,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			data = noiseData,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_LINEAR,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}
	end

	local vectorLICShader = solver.vectorLICShader
	vectorLICShader:use()
	app.gradientTex:bind(1)

	gl.glUniform1i(vectorLICShader.uniforms.useCoordMap.loc, app.display_useCoordMap)
	gl.glUniform1i(vectorLICShader.uniforms.useLog.loc, var.useLog)
	gl.glUniform1f(vectorLICShader.uniforms.valueMin.loc, valueMin)
	gl.glUniform1f(vectorLICShader.uniforms.valueMax.loc, valueMax)
	gl.glUniform1i(vectorLICShader.uniforms.integralMaxIter.loc, self.integralMaxIter)
	
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)
	
	self:drawSolverWithVar(app, solver, var, vectorLICShader, xmin, xmax, ymin, ymax)

-- [[
	if solver.amr then
		for k,subsolver in pairs(solver.amr.child) do
			self:drawSolverWithVar(app, subsolver, var, vectorLICShader, xmin, xmax, ymin, ymax)
		end
	end
--]]

	gl.glDisable(gl.GL_BLEND)

	app.gradientTex:unbind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)
	vectorLICShader:useNone()

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

function DrawVectorLIC:display(app, solvers, varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
	app.view:projection(ar)
	app.view:modelview()
	if app.view.getOrthoBounds then
		xmin, xmax, ymin, ymax = app.view:getOrthoBounds(ar)
	else
		xmin, xmax, ymin, ymax = graph_xmin, graph_xmax, graph_ymin, graph_ymax
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
		if not require 'hydro.solver.meshsolver'.is(solver) then
			local var = solver.displayVarForName[varName]
			if var and var.enabled then
				self:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax)
			end
		end
	end
--	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:displayVector_LIC(...)
		if not self.drawVectorLIC then self.drawVectorLIC = DrawVectorLIC() end
		return self.drawVectorLIC:display(self, ...)
	end
end
