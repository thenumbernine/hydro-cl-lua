local ffi = require 'ffi'
local class = require 'ext.class'
local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'
local GLPingPong = require 'gl.pingpong'
local GLProgram = require 'gl.program'
local template = require 'template'

local Draw2DLIC = class()

function Draw2DLIC:showDisplayVar(app, solver, var, xmin, xmax, ymin, ymax, useLog)

	--[[
	for each solver, for each var, I'll need a vector field state to keep track of the line offset within the cell and length
	--]]
	solver.display2D_LIC_State = solver.display2D_LIC_State or {}
	solver.display2D_LIC_State[var.name] = solver.display2D_LIC_State[var.name] or {}
	local state = solver.display2D_LIC_State[var.name]

	local w = tonumber(solver.gridSize.x)
	local h = tonumber(solver.gridSize.x)
	if not state.pingpong
	or (state.pingpong.width ~= w or state.pingpong.height ~= h)
	then
		local data
		data = ffi.new('unsigned char[?]', w * h * 4)
		for j=0,h-1 do
			for i=0,w-1 do
				for k=0,3 do
					data[k+4*(i+w*j)] = math.random(0,255)
				end
			end
		end
		state.pingpong = GLPingPong{
			width = w,
			height = h,
			internalFormat = gl.GL_RGBA,
			format = gl.GL_RGBA,
			type = gl.GL_UNSIGNED_BYTE,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_NEAREST,
			data = data,
		}
	end

	--[[
	now do a FBO render of the vector field state
	determine offsets and lengths of lines
	channels: ofsx, ofsy, 0, magn
	--]]
	if not state.updateVectorField2Shader then
		state.updateVectorField2Shader = GLProgram{
			vertexCode = [[
varying vec2 pos;

void main() {
	gl_Position = gl_Vertex;
	pos = gl_Vertex.xy;
}
]],
			fragmentCode = template([[
<?
local clnumber = require 'cl.obj.number'
local dx = 1 / tonumber(solver.gridSize.x)
local dy = 1 / tonumber(solver.gridSize.y)
?>
varying vec2 pos;
uniform sampler2D varTex;
uniform sampler2D stateTex;
uniform float alpha;

void main() {
	vec2 dx = vec2(<?=clnumber(dx)?>, <?=clnumber(dy)?>);
	vec2 ofs = pos - texture2D(varTex, pos).xy * dx;

	gl_FragColor = mix(
		texture2D(stateTex, pos),
		texture2D(stateTex, ofs),
		alpha);
}
]], {
	solver = solver,
}),
			uniforms = {
				varTex = 0,
				stateTex = 1,
			},
		}
	end

	if not state.vectorLICShader then
		local env = {
			clnumber = require 'cl.obj.number',
			app = app,
			solver = solver,
			coord = solver.coord,
		}
		state.vectorLICShader = GLProgram{
			vertexCode = template([[
varying vec3 viewCoord;

<?=coord:getCoordMapInvGLSLCode()?>

void main() {
	viewCoord = gl_Vertex.xyz;
	// disregard 'z' for rendering
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xy, 0., 1.);
}
]], env),
			fragmentCode = template([[
varying vec3 viewCoord;

<?=coord:getCoordMapInvGLSLCode()?>

//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useCoordMap;
uniform bool useLog;
uniform float valueMin, valueMax;

uniform vec2 texCoordMax;
uniform <?=solver.dim == 3 and 'sampler3D' or 'sampler2D'?> tex;
uniform sampler1D gradientTex;

uniform vec2 solverMins, solverMaxs;

uniform sampler2D stateTex;

void main() {

	vec2 gridCoord = viewCoord.xy;
	if (useCoordMap) {
		gridCoord = coordMapInv(vec3(gridCoord, 0.)).xy;
	} else {
		gridCoord = .5 * (gridCoord + 1.) * (solverMaxs.xy - solverMins.xy) + solverMins.xy; 
	}

	if (gridCoord.x < solverMins.x || gridCoord.x > solverMaxs.x ||
		gridCoord.y < solverMins.y || gridCoord.y > solverMaxs.y
	) {
		discard;
	}

	vec3 texCoord;
	texCoord.xy = vec2(
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * <?=clnumber(solver.sizeWithoutBorder.x)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.x)?>,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * <?=clnumber(solver.sizeWithoutBorder.y)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.y)?>
	) * texCoordMax;
	texCoord.z = viewCoord.z;

<? if solver.dim == 3 then
?>	float value = texture3D(tex, texCoord).r;
<? else
?>	float value = texture2D(tex, texCoord.xy).r;
<? end
?>	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
	value = (value * <?=clnumber(app.gradientTex.width-1)?> + .5) / <?=clnumber(app.gradientTex.width)?>;
	gl_FragColor = texture1D(gradientTex, value) * texture2D(stateTex, texCoord.xy);
}
]],	env),
			uniforms = {
				tex = 0,
				gradientTex = 1,
				stateTex = 2,
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

-- [[
	state.pingpong:swap()
	state.pingpong:draw{
		viewport = {0, 0, state.pingpong.width, state.pingpong.height},
		resetProjection = true,
		shader = state.updateVectorField2Shader,
		texs = {
			solver:getTex(var),
			state.pingpong:prev(),
		},
		callback = function()
			gl.glBegin(gl.GL_TRIANGLE_STRIP)
			for _,v in ipairs{{0,0},{1,0},{0,1},{1,1}} do
				gl.glTexCoord2d(v[1], v[2])
				gl.glVertex2d(2*v[1]-1, 2*v[2]-1)
			end
			gl.glEnd()
		end,
	}
--]]

	state.vectorLICShader:use()
	gl.glUniform1i(state.vectorLICShader.uniforms.useLog.loc, var.useLog)
	-- [[ this gives the l1 bounds of the vector field
	gl.glUniform1f(state.vectorLICShader.uniforms.valueMin.loc, valueMin)
	gl.glUniform1f(state.vectorLICShader.uniforms.valueMax.loc, valueMax)
	--]]
	--[[ it'd be nice instead to get the norm bounds ... 
	-- but looking at the reduce calculations, the easiest way to do that is
	-- to associate each vector display shader with a norm display shader
	-- and then just reduce that
	gl.glUniform1f(state.vectorLICShader.uniforms.valueMin.loc, 0)
	gl.glUniform1f(state.vectorLICShader.uniforms.valueMax.loc, math.max(math.abs(valueMin), math.abs(valueMax)))
	--]]
	solver:getTex(var):bind(0)
	app.gradientTex:bind(1)
	state.pingpong:cur():bind(2)
	
	gl.glUniform3f(state.vectorLICShader.uniforms.solverMins.loc, solver.mins:unpack())
	gl.glUniform3f(state.vectorLICShader.uniforms.solverMaxs.loc, solver.maxs:unpack())

	gl.glBegin(gl.GL_QUADS)
	gl.glVertex3d(xmin, ymin, app.displayFixedZ)
	gl.glVertex3d(xmax, ymin, app.displayFixedZ)
	gl.glVertex3d(xmax, ymax, app.displayFixedZ)
	gl.glVertex3d(xmin, ymax, app.displayFixedZ)
	gl.glEnd()

	state.pingpong:cur():unbind(2)
	app.gradientTex:unbind(1)
	solver:getTex(var):unbind(0)
	state.vectorLICShader:useNone()
end

function Draw2DLIC:display(app, solvers, ar, varName, ...)
	app.view:projection(ar)
	app.view:modelview()

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)

	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then	
			self:showDisplayVar(app, solver, var, ...)
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:displayVector_LIC(...)
		if not self.drawVectorLIC then self.drawVectorLIC = Draw2DLIC() end
		return self.drawVectorLIC:display(self, ...)
	end
end
