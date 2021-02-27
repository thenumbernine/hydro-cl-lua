#version 460

<? 
local varying = vertexShader and "out"
		or fragmentShader and "in"
		or error("don't know what to set varying to")
?>

<?=varying?> float cellindexv;

<?=draw:getModuleCodeGLSL(coordMapGLSL, coordMapInvGLSL)?>
<?=draw:getCommonGLSLFragCode()?>

<?
if vertexShader then
?>

uniform float drawCellScale;

attribute vec3 vtx;
attribute vec3 vtxcenter;
attribute float cellindex;

void main() {
	vec3 v = (vtx - vtxcenter) * drawCellScale + vtxcenter;
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
	cellindexv = cellindex;
}

<?
end
if fragmentShader then 
	if require "gl.tex2d":isa(solver.tex) then 
		if solver.texSize.y == 1 then
?>
float getValue() {
	float tc = (cellindexv + .5) / texSize.x;
	return texture2D(tex, vec2(tc, .5)).r;
}
<?
		else
?>
float getValue() {
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
	elseif require "gl.tex3d":isa(solver.tex) then 
?>
float getValue() {
	float i = cellindexv;
	vec3 tc;
	
	tc.x = mod(i, texSize.x);
	i -= tc.x;
	tc.x += .5;
	tc.x /= texSize.x;
	
	tc.y = mod(i, texSize.y);
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
<?
end
?>
