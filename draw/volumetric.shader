<?
local clnumber = require 'cl.obj.number'
?>

varying vec3 texCoordStart;
varying vec3 vertexStart;
varying vec3 eye;

<? if vertexShader then ?>
<?=solver.coord:getCoordMapGLSLCode()?>

uniform vec3 mins, maxs;
void main() {
	//this only needs to be done once per render
	//i.e. uniform?
	eye = (gl_ModelViewMatrixInverse * vec4(0., 0., 0., 1.)).xyz;
	texCoordStart = gl_MultiTexCoord0.xyz;
	
	vec4 x = gl_Vertex;
	x.xyz *= maxs - mins;
	x.xyz += mins;
	x = vec4(coordMap(x.xyz), x.w);	
	
	vertexStart = x.xyz;	//this means we have to invert coordMap as we travel through the cartesian space
	gl_Position = gl_ModelViewProjectionMatrix * x;
}

<? end
if fragmentShader then ?>

uniform sampler3D tex;
uniform sampler1D gradientTex;
uniform int maxiter;
uniform vec3 oneOverDx;

uniform float numGhost;
uniform vec3 texSize;

uniform bool useLog;
uniform float valueMin, valueMax;

uniform float alpha;
uniform float alphaGamma;
uniform bool useIsos;
uniform bool useLighting;
uniform float numIsobars;

#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>

float getValue(vec3 p) {
	float value = texture3D(tex, p).r;
	if (useLog) value = log(1. + abs(value)) * _1_LN_10;
	return value;
}

float logmap(float x) { 
	return log(1. + abs(x)) * _1_LN_10;
}

float calcFrac(float value) {
	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		return (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		return (value - valueMin) / (valueMax - valueMin);
	}
}

float getTex(vec3 tc) {
	//using texel coordinates as-is
	//return texture3D(tex, tc).r;

	//getting rid of the ghost cells
	tc = tc * ((texSize - 2. * numGhost) / texSize) + (numGhost / texSize);
	return texture3D(tex, tc).r;
}

void main() {
	vec3 p = texCoordStart;
	vec4 result = vec4(0., 0., 0., 1.);
		
	float value = getTex(p);
	float frac = calcFrac(value);
	float gradTC = frac * <?=clnumber((app.gradientTex.width-1) / app.gradientTex.width)?>;
	vec4 voxelColor = texture1D(gradientTex, gradTC);
	
	voxelColor.a = alpha;// * min(1., mod(value * 4., 3.));
	voxelColor.a *= min(1., length(voxelColor.rgb));
	result.rgb += result.a * voxelColor.a * voxelColor.rgb;
	result.a *= 1. - voxelColor.a;
	
	vec3 step = vertexStart - eye;
	step = normalize(step);
	step /= float(maxiter);
	//step /= oneOverDx;
	for (int i = 2; i <= maxiter; i++) {
		p += step;
		if (p.x < 0. || p.y < 0. || p.z < 0. ||
			p.x > 1. || p.y > 1. || p.z > 1.) break;

		//notice if you store 1-alpha in the color you trace through the volume
		//then you can forward-trace your blending
		//instead of having to backward-trace
		//(as you would when rendering transparent stuff on top of each other)
		//this will allow you to bailout early if your transparency ever hits fully opaque
		value = getTex(p);
		frac = calcFrac(value);
		gradTC = frac * <?=clnumber((app.gradientTex.width-1) / app.gradientTex.width)?>;
		voxelColor = texture1D(gradientTex, gradTC);
		//voxelColor.a /= min(1., length(voxelColor.rgb));

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
				getTex(p + vec3(dx,0.,0.)) - getTex(p - vec3(dx,0.,0.)),
				getTex(p + vec3(0.,dx,0.)) - getTex(p - vec3(0.,dx,0.)),
				getTex(p + vec3(0.,0.,dx)) - getTex(p - vec3(0.,0.,dx)));
			texGrad = gl_NormalMatrix * texGrad;
			texGrad = normalize(texGrad);
			float lum = max(ambient, abs(dot(lightVec, texGrad)));	
			voxelColor.rgb *= lum;
		}

		result.rgb += result.a * voxelColor.a * voxelColor.rgb;
		result.a *= 1. - voxelColor.a;

		//if (result.a < .01) break;
	}
	gl_FragColor = vec4(result.rgb, 1. - result.a);
}

<? end ?>
