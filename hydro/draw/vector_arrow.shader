#version 460

<?
local coord = solver.coord
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>

<?=draw:getModuleCodeGLSL(eqn.symbols.coordMapGLSL, eqn.symbols.coordMapInvGLSL, eqn.symbols.cartesianFromCoord)?>
<?=draw:getCommonGLSLFragCode()?>

<? if vertexShader then ?>

<?=varying?> vec4 color;

attribute vec2 vtx;
attribute vec3 gridCoord;	//in hydro/draw/draw.lua I said gridCoord was already half-off, but this is integers starting at 0

uniform float scale;

void main() {
	//hmm, to rotate the sampled slice in coordinate space, I would have to 
	// 1) map from grid coord to world coords
	// 2) pad with fixed coords.  
	// 3) rotate
	// 4) map back from world coords to grid coords
	// 5) map back from grid coords to texel coords
	// ... seems redundant.
	vec3 texCoord = (gridCoord + .5) / sizeWithoutBorder;
	vec3 chartCoord = texToChartCoord(texCoord);
	vec3 center = chartCoord;
	if (displayDim <= 1) chartCoord.y = displayFixed.x;
	if (displayDim <= 2) chartCoord.z = displayFixed.y;
	chartCoord = chartToWorldCoord(chartCoord);
	chartCoord = quatRotate(displaySliceAngle, chartCoord);
	chartCoord = worldToChartCoord(chartCoord);
	texCoord = chartToTexCoord(chartCoord);
#if 0
	if (texCoord.x < 0. || texCoord.x > 1.) discard;
	if (displayDim > 1 && (texCoord.y < 0. || texCoord.y > 1.)) discard;
	if (displayDim > 2 && (texCoord.z < 0. || texCoord.z > 1.)) discard;
#endif
	vec3 dir = getTex(texCoord).rgb;

	dir = cartesianFromCoord(dir, center);	
	float value = length(dir);
	dir /= value;

	vec3 tv;
	if (displayDim < 3) {
		tv = vec3(-dir.y, dir.x, 0.);
	} else {	//displayDim == 3
		vec3 vx = vec3(0., -dir.z, dir.y);
		vec3 vy = vec3(dir.z, 0., -dir.x);
		vec3 vz = vec3(-dir.y, dir.x, 0.);
		float lxsq = dot(vx,vx);
		float lysq = dot(vy,vy);
		float lzsq = dot(vz,vz);
		if (lxsq > lysq) {		//x > y
			if (lxsq > lzsq) {	//x > z, x > y
				tv = vx;
			} else {			//z > x > y
				tv = vz;
			}
		} else {				//y >= x
			if (lysq > lzsq) {	//y >= x, y > z
				tv = vy;
			} else {			// z > y >= x
				tv = vz;
			}
		}
	}

	//TODO don't just normalize dir, but scale down by the view size ... ? or the dx?
	
	value = getGradientFrac(value);
//TODO make this a flag, whether to normalize vectors or not? 
	float valuescale = scale;// * clamp(value, 0., 1.);
	value = getGradientTexCoord(value);
	color = texture1D(gradientTex, value);

	//cartesian coords
	vec3 v = chartToWorldCoord(center) + valuescale * (vtx.x * dir + vtx.y * tv);
	
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
}

<?
end
if fragmentShader then
?>

<?=varying?> vec4 color;

out vec4 fragColor;

void main() {
	fragColor = color;
}

<? end ?>
