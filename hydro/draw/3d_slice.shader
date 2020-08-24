<?
local clnumber = require 'cl.obj.number'
local coord = solver.coord
local app = solver.app
?>
varying vec3 texCoord;	//[0,1]^n

<? local useClipPlanes = false ?>
<? if useClipPlanes then ?>
varying vec3 pos;		//positive after coordinate mapping, before view transform
<? end ?>

<? if vertexShader then ?>
<?=coord:getModuleCodeGLSL'coordMapGLSL'?>

uniform bool useCoordMap;

uniform vec3 solverMins, solverMaxs;
void main() {
	texCoord = gl_Vertex.xyz;
	
	vec4 x = gl_Vertex;
	if (useCoordMap) {
		x.xyz *= solverMaxs - solverMins;
		x.xyz += solverMins;
		x = vec4(coordMap(x.xyz), x.w);
	} else {
		x.xyz *= 2.;
		x.xyz -= 1.;
	}
<? if useClipPlanes then ?>
	pos = x.xyz;
<? end ?>
	gl_Position = gl_ModelViewProjectionMatrix * x;
}

<? end
if fragmentShader then ?>

<?=solver:getGradientGLSLCode()?>

<? if solver.dim == 3 then
?>uniform sampler3D tex;
<? else
?>uniform sampler2D tex;
<? end
?>

uniform float numGhost;
uniform vec3 texSize;

uniform vec3 normal;
uniform float alpha;
uniform float alphaGamma;

uniform bool useIsos;
uniform float numIsobars;
<? if useClipPlanes then ?>
<? for i,clipInfo in ipairs(clipInfos) do
?>uniform bool clipEnabled<?=i?>;
<? end
?>
<? end ?>
uniform bool useLighting;

float getTex(vec3 tc) {
	//using texel coordinates as-is
	//return texture3D(tex, tc).r;

	//getting rid of the ghost cells
	tc = tc * ((texSize - 2. * numGhost) / texSize) + (numGhost / texSize);
<? if solver.dim == 3 then
?>	return texture3D(tex, tc).r;
<? else
?>	return texture2D(tex, tc.xy).r;
<? end
?>}

void main() {
<? if useClipPlanes then ?>
	vec4 worldPos = gl_ModelViewMatrix * vec4(pos, 1.);
<? for i,clipInfo in ipairs(clipInfos) do
?>	if (clipEnabled<?=i?> && dot(worldPos, gl_ClipPlane[<?=i-1?>]) < 0.) discard;
<? end
?>
<? end
?>	float value = getTex(texCoord);
	
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
			getTex(texCoord + vec3(dx,0.,0.)) - getTex(texCoord - vec3(dx,0.,0.)),
			getTex(texCoord + vec3(0.,dx,0.)) - getTex(texCoord - vec3(0.,dx,0.)),
			getTex(texCoord + vec3(0.,0.,dx)) - getTex(texCoord - vec3(0.,0.,dx)));
		texGrad = gl_NormalMatrix * texGrad;
		texGrad = normalize(texGrad);
		float lum = max(ambient, abs(dot(lightVec, texGrad)));	
		voxelColor.rgb *= lum;
	}

	//calculate normal in screen coordinates
	vec4 n = gl_ModelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	gl_FragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}

<? end ?>
