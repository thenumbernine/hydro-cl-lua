local class = require 'ext.class'
local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'
local GLPingPong = require 'gl.pingpong'
local GLProgram = require 'gl.program'
local template = require 'template'
local clnumber = require 'cl.obj.number'

local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
}

local DrawVectorField = class()

function DrawVectorField:displayVectorField(app, solver, varName, xmin, xmax, ymin, ymax, useLog)
	--if app.displayDim ~= 2 then return end
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		--[[
		for each solver, for each var, I'll need a vector field state to keep track of the line offset within the cell and length
		--]]
		solver.vectorFieldStates = solver.vectorFieldStates or {}
		solver.vectorFieldStates[var.name] = solver.vectorFieldStates[var.name] or {}
		local state = solver.vectorFieldStates[var.name]

		local w = tonumber(solver.gridSize.x)
		local h = tonumber(solver.gridSize.x)
		if not state.pingpong
		or (state.pingpong.width ~= w or state.pingpong.height ~= h)
		then
			state.pingpong = GLPingPong{
				width = w,
				height = h,
				internalFormat = gl.GL_RGBA,
				format = gl.GL_RGBA,
				minFilter = gl.GL_NEAREST,
				magFilter = gl.GL_NEAREST,
				type = gl.GL_FLOAT,
			}
		end

		--[[
		now do a FBO render of the vector field state
		determine offsets and lengths of lines
		channels: ofsx, ofsy, 0, magn
		--]]
		if not solver.updateVectorFieldShader then
			solver.updateVectorFieldShader = GLProgram{
				vertexCode = [[
varying vec2 pos;

void main() {
	gl_Position = gl_Vertex;
	pos = gl_Vertex.xy;
}
]],
				fragmentCode = template([[
varying vec2 pos;

uniform sampler2D fieldTex;
uniform sampler2D stateTex;

#define lenSq(x)	dot(x,x)

void main() {
	vec2 tc = pos.xy * .5 + .5;

	vec4 field_c = texture2D(fieldTex, tc);
	vec4 field_xR = texture2D(fieldTex, tc + vec2(<?=dy?>, 0.));
	vec4 field_xL = texture2D(fieldTex, tc - vec2(<?=dy?>, 0.));
	vec4 field_yR = texture2D(fieldTex, tc + vec2(0., <?=dy?>));
	vec4 field_yL = texture2D(fieldTex, tc - vec2(0., <?=dy?>));
	float lenSq_c = lenSq(field_c.xy);
	float lenSq_xR = lenSq(field_xR.xy);
	float lenSq_xL = lenSq(field_xL.xy);
	float lenSq_yR = lenSq(field_yR.xy);
	float lenSq_yL = lenSq(field_yL.xy);
	
	vec4 state_c = texture2D(stateTex, tc);
	vec4 state_xR = texture2D(stateTex, tc + vec2(<?=dy?>, 0.));
	vec4 state_xL = texture2D(stateTex, tc - vec2(<?=dy?>, 0.));
	vec4 state_yR = texture2D(stateTex, tc + vec2(0., <?=dy?>));
	vec4 state_yL = texture2D(stateTex, tc - vec2(0., <?=dy?>));

	float value = sqrt(lenSq_c);
	gl_FragColor.a = value;

	float xRy = state_xR.y - field_xR.y * state_xR.x / field_xR.x;
	float xLy = state_xL.y + field_xL.y * (1. - state_xL.x) / field_xL.x;
	float yRx = state_yR.x - field_yR.x * state_yR.y / field_yR.y;
	float yLx = state_yL.x + field_yL.x * (1. - state_yL.y) / field_yL.y;

	vec2 x;
	//TODO search from best to worst and using those inside the [0,1] bounds
	//instead of finding the best and hoping it is within the [0,1] bounds
	//Should I bother with the center cell?
	if (lenSq_xR > lenSq_xL && lenSq_xR > lenSq_yL && lenSq_xR > lenSq_yR && xRy >= 0. && xRy <= 1. && lenSq_xR > lenSq_c) {
		//reposition the cell center so it points at the intersection
		x = vec2(0., xRy);
		
	} else if (lenSq_xL > lenSq_xR && lenSq_xL > lenSq_yL && lenSq_xL > lenSq_yR && xLy >= 0. && xLy <= 1. && lenSq_xL > lenSq_c) {
		x = vec2(1., xLy);
	} else if (lenSq_yR > lenSq_yL && lenSq_yR > lenSq_xL && lenSq_yR > lenSq_xR && yRx >= 0. && yRx <= 1. && lenSq_yR > lenSq_c) {
		x = vec2(yRx, 0.);
	} else if (lenSq_yL > lenSq_yR && lenSq_yL > lenSq_xL && lenSq_yL > lenSq_xR && yLx >= 0. && yLx <= 1. && lenSq_yL > lenSq_c) {
		x = vec2(yLx, 1.);
	}
	
}
]], {
	dx = clnumber(1 / solver.gridSize.x),
	dy = clnumber(1 / solver.gridSize.y),
}),
				uniforms = {
					fieldTex = 0,
					stateTex = 1,
				},
			}
		end

		
		local valueMin, valueMax
		if var.heatMapFixedRange then
			valueMin = var.heatMapValueMin
			valueMax = var.heatMapValueMax
		else
			local component = solver.displayComponentFlatList[var.component]
			local vectorField = solver:isVarTypeAVectorField(var.type)
			if vectorField then
				valueMin, valueMax = solver:calcDisplayVarRange(var, component.magn)
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
			end
			var.heatMapValueMin = valueMin
			var.heatMapValueMax = valueMax
		end
		
		solver:calcDisplayVarToTex(var)	

		-- solver:getTex(var) now has the vector field
		-- now we have to render-to-tex with input as solver:getTex(var), and the last iteration's state tex
		-- ... using a fbo that updates the state tex
		
		state.pingpong:draw{
			shader = solver.updateVectorFieldShader,
			texs = {
				solver:getTex(var),
				state.pingpong:prev(),
			},
			callback = function()
				gl.glBegin(gl.GL_TRIANGLE_STRIP)
				for _,v in ipairs{{0,0},{1,0},{0,1},{1,1}} do
					gl.glTexCoord2d(v[1], v[2])
					gl.glVertex2d(v[1], v[2])
				end
				gl.glEnd()
			end,	
		}


		solver.vectorFieldShader:use()
		gl.glUniform1i(solver.vectorFieldShader.uniforms.useLog.loc, var.useLog)
		gl.glUniform1i(solver.vectorFieldShader.uniforms.displayDim.loc, app.displayDim)
		-- [[ this gives the l1 bounds of the vector field
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMin.loc, valueMin)
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMax.loc, valueMax)
		--]]
		--[[ it'd be nice instead to get the norm bounds ... 
		-- but looking at the reduce calculations, the easiest way to do that is
		-- to associate each vector display shader with a norm display shader
		-- and then just reduce that
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMin.loc, 0)
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMax.loc, math.max(math.abs(valueMin), math.abs(valueMax)))
		--]]
		solver:getTex(var):bind(0)
		app.gradientTex:bind(1)
		
		gl.glUniform3f(solver.vectorFieldShader.uniforms.mins.loc, solver.mins:unpack())
		gl.glUniform3f(solver.vectorFieldShader.uniforms.maxs.loc, solver.maxs:unpack())
		-- how to determine scale?
		--local scale = app.displayVectorField_scale * (valueMax - valueMin)
		--local scale = app.displayVectorField_scale / (valueMax - valueMin)
		local scale = app.displayVectorField_scale 
			* app.displayVectorField_step 
			* math.min(
				(solver.maxs.x - solver.mins.x) / tonumber(solver.gridSize.x),
				(solver.maxs.y - solver.mins.y) / tonumber(solver.gridSize.y),
				(solver.maxs.z - solver.mins.z) / tonumber(solver.gridSize.z))
		gl.glUniform1f(solver.vectorFieldShader.uniforms.scale.loc, scale) 

		local step = app.displayVectorField_step

		--[[ goes just slightly faster.  24 vs 23 fps.
		app.vectorField_displayList = app.vectorField_displayList or {}
		local glCallOrDraw = require 'gl.call'
		glCallOrDraw(app.vectorField_displayList, function()
		--]]	
			gl.glBegin(gl.GL_LINES)
			for k=0,tonumber(solver.sizeWithoutBorder.z-1),step do
				for j=0,tonumber(solver.sizeWithoutBorder.y-1),step do
					for i=0,tonumber(solver.sizeWithoutBorder.x-1),step do
						local tx = (i + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
						local ty = (j + .5 + (solver.dim > 1 and solver.numGhost or app.displayFixedY * tonumber(solver.gridSize.z))) / tonumber(solver.gridSize.y)
						local tz = (k + .5 + (solver.dim > 2 and solver.numGhost or app.displayFixedZ * tonumber(solver.gridSize.z))) / tonumber(solver.gridSize.z)
						gl.glMultiTexCoord3f(gl.GL_TEXTURE0, tx, ty, tz)	
						local x = (i + .5) / tonumber(solver.sizeWithoutBorder.x)
						local y = (j + .5) / tonumber(solver.sizeWithoutBorder.y)
						local z = (k + .5) / tonumber(solver.sizeWithoutBorder.z)
						gl.glMultiTexCoord3f(gl.GL_TEXTURE1, x, y, z)
						for _,q in ipairs(arrow) do
							gl.glVertex2f(q[1], q[2])
						end
					end
				end
			end
			gl.glEnd()
		--[[
		end)
		--]]

		app.gradientTex:unbind(1)
		solver:getTex(var):unbind(0)
		solver.vectorFieldShader:useNone()
	end
end

function DrawVectorField:display(app, solvers, ar, ...)
	app.view:projection(ar)
	app.view:modelview()

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)

	for _,solver in ipairs(solvers) do
		self:displayVectorField(app, solver, ...)
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:displayVectorField(...)
		if not self.drawVectorField2 then self.drawVectorField2 = DrawVectorField() end
		return self.drawVectorField2:display(self, ...)
	end
end
