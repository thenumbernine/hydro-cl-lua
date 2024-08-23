<?=draw.glslVersion?>

<?
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>

<?=varying?> vec3 texCoord;	//[0,1]^n

<? if app.useClipPlanes then ?>
<?=varying?> vec3 pos;		//positive after coordinate mapping, before view transform
<? end ?>

<?=draw:getModuleCodeGLSL(coordMapGLSL, coordMapInvGLSL)?>
<?=draw:getCommonGLSLFragCode()?>

<? if vertexShader then ?>

//in tex coords
in vec4 vertex;

void main() {
	vec4 x = vertex;

	//map from [0,1]^3 to [cartesianMin, cartesianMax]
	x.xyz = x.xyz * (cartesianMax - cartesianMin) + cartesianMin;
	texCoord = x.xyz;

<? if app.useClipPlanes then ?>
	pos = x.xyz;
<? end ?>
	gl_Position = mvProjMat * x;
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
<? if app.useClipPlanes then ?>
<? for i,clipInfo in ipairs(app.clipInfos) do
?>uniform bool clipEnabled<?=i?>;
uniform vec4 clipPlane<?=i?>;
<? end
?>
<? end ?>

float getVoxelValue(vec3 tc) {
	return getTex(tc).r;
}

void main() {
<? if app.useClipPlanes then ?>
	vec4 pos4 = vec4(pos, 1.);
<? for i,clipInfo in ipairs(app.clipInfos) do
?>	if (clipEnabled<?=i?> && (dot(pos4, clipPlane<?=i?>) < 0.)) discard;
<? end
end
?>	

	vec3 tc = texCoord;
	if (displayDim <= 1) tc.y = displayFixed.x;
	if (displayDim <= 2) tc.z = displayFixed.y;
	tc = quatRotate(displaySliceAngle, tc);
	tc = worldToTexCoord(tc);	//inverse map from world coord to chart coord, then to tex coords
	tc = texToNoGhostCoord(tc);	//getting rid of the ghost cells
	
	if (tc.x < 0 || tc.x > 1 || 
		tc.y < 0 || tc.y > 1 || 
		tc.z < 0 || tc.z > 1) discard;

	float value = getVoxelValue(tc);
	
	float frac = getGradientFrac(value);
	float gradTC = getGradientTexCoord(frac);
	vec4 voxelColor = texture(gradientTex, vec2(gradTC, .5));

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
	vec4 n = mvProjMat * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	fragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}

<? end ?>
