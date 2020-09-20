#version 460

<? 
local clnumber = require 'cl.obj.number'
local inout = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set inout to")
?>

<?=inout?> float cellindexv;

<? if vertexShader then ?>

uniform float drawCellScale;
uniform mat4 modelViewProjectionMatrix;

attribute vec3 vtx;
attribute vec3 vtxcenter;
attribute float cellindex;

void main() {
	vec3 v = (vtx - vtxcenter) * drawCellScale + vtxcenter;
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
	cellindexv = cellindex;
}

<? end
if fragmentShader then ?>

<?=solver:getGradientGLSLCode()?>
<?=draw:getCommonGLSLFragCode(solver)?>

<? 
if require 'gl.tex2d'.is(solver.tex) then 
	if solver.texSize.y == 1 then
?>
float getValue() {
	float tc = (cellindexv + .5) * <?=clnumber(1 / tonumber(solver.texSize.x))?>;
	return texture2D(tex, vec2(tc, .5)).r;
}
<?
	else
?>
float getValue() {
	vec2 texSize = vec2(<?=clnumber(solver.texSize.x)?>, <?=clnumber(solver.texSize.y)?>);
	float i = cellindexv;
	vec2 tc;
	tc.x = mod(i, texSize.x);
	i -= tc.x;
	tc.x += .5;
	tc.x /= texSize.x;
	tc.y = (i + .5) / texSize.y;
	return texture2D(tex, tc).r;
}
<? 
	end
elseif require 'gl.tex3d'.is(solver.tex) then 
?>
function getValue() {
	vec3 texSize = vec3(<?=clnumber(solver.texSize.x)?>, <?=clnumber(solver.texSize.y)?>, <?=clnumber(solver.texSize.z)?>);
	float i = cellindexv;
	vec3 tc;
	
	tc.x = mod(i, texSize.x);
	i -= tc.x;
	tc.x += .5;
	tc.x /= texSize.x;
	
	tc.y = fmod(i, texSize.y);
	i -= tc.y;
	tc.y += .5;
	tc.y /= texSize.y;
	
	tc.z = (i + .5) / texSize.z;
	return texture3D(tex, tc).r;
}
<? 
else
	error("don't know how to handle your tex lua obj class")
end 
?>

out vec4 fragColor;

void main() {
	float value = getValue();
	fragColor = getGradientColor(value);
}
<? end ?>
