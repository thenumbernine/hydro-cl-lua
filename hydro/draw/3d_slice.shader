#version 460

<?
local coord = solver.coord
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>

<?=varying?> vec3 texCoord;	//[0,1]^n

<? local useClipPlanes = false ?>
<? if useClipPlanes then ?>
<?=varying?> vec3 pos;		//positive after coordinate mapping, before view transform
<? end ?>

<?=coord:getModuleCodeGLSL("coordMapGLSL", "coordMapInvGLSL")?>
<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>

//in tex coords
attribute vec4 vertex;

void main() {
	texCoord = vertex.xyz;
	
	vec4 x = vertex;
	x.xyz = texToWorldCoord(x.xyz);
<? if useClipPlanes then ?>
	pos = x.xyz;
<? end ?>
	gl_Position = modelViewProjectionMatrix * x;
}

<? end
if fragmentShader then ?>

out vec4 fragColor;

uniform bool useLighting;
uniform float alpha;
uniform float alphaGamma;
uniform bool useIsos;
uniform float numIsobars;

uniform vec3 normal;
<? if useClipPlanes then ?>
<? for i,clipInfo in ipairs(clipInfos) do
?>uniform bool clipEnabled<?=i?>;
<? end
?>
<? end ?>

float getVoxelValue(vec3 tc) {
	return getTex(tc).r;
}

void main() {
<? if useClipPlanes then ?>
	vec4 worldPos = gl_ModelViewMatrix * vec4(pos, 1.);
<? for i,clipInfo in ipairs(clipInfos) do
?>	if (clipEnabled<?=i?> && dot(worldPos, gl_ClipPlane[<?=i-1?>]) < 0.) discard;
<? end
end
?>	

	vec3 tc = texCoord;
	if (displayDim <= 1) tc.y = displayFixed.x;
	if (displayDim <= 2) tc.z = displayFixed.y;
	tc = texToWorldCoord(tc);
	tc = quatRotate(displaySliceAngle, tc);
	tc = worldToTexCoord(tc);
	tc = texToNoGhostCoord(tc);	//getting rid of the ghost cells

	float value = getVoxelValue(tc);
	
	float frac = getGradientFrac(value);
	float gradTC = getGradientTexCoord(frac);
	vec4 voxelColor = texture1D(gradientTex, gradTC);

	//don't bother with the gamma factor if we're using isobars
	if (useIsos) {
		float ffrac = fract(frac * numIsobars + .5);
//percentage of each isobar that is solid
#define stepThickness .1
		float isoAlpha = smoothstep(0., stepThickness / 3., ffrac) 
			- smoothstep(2./3.*stepThickness, stepThickness, ffrac);
		
		voxelColor.a = isoAlpha;
	} else {
		voxelColor.a = pow(alpha, alphaGamma);
	}

	if (useLighting) {
		const float dx = <?=1/32?>;
		const float ambient = .3;
		vec3 lightVec = vec3(0., 0., 1.);
		vec3 texGrad = vec3(
			getVoxelValue(tc + vec3(dx,0.,0.)) - getVoxelValue(tc - vec3(dx,0.,0.)),
			getVoxelValue(tc + vec3(0.,dx,0.)) - getVoxelValue(tc - vec3(0.,dx,0.)),
			getVoxelValue(tc + vec3(0.,0.,dx)) - getVoxelValue(tc - vec3(0.,0.,dx)));
		texGrad = normalMatrix * texGrad;
		texGrad = normalize(texGrad);
		float lum = max(ambient, abs(dot(lightVec, texGrad)));	
		voxelColor.rgb *= lum;
	}

	//calculate normal in screen coordinates
	vec4 n = modelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	fragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}

<? end ?>
