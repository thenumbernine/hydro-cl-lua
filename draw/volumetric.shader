varying vec3 texCoordStart;
varying vec3 vertexStart;
varying vec3 eye;

<? if vertexShader then ?>

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
uniform sampler1D gradient;
uniform int maxiter;
uniform vec3 oneOverDx;
uniform float scale;
uniform bool useLog;
uniform float alpha;
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>

float getValue(vec3 p) {
	float value = texture3D(tex, p).r;
	if (useLog) value = log(1. + abs(value)) * _1_LN_10;
	value *= scale;
	return value;
}

void main() {
	vec3 p = texCoordStart;
	vec4 result = vec4(0., 0., 0., 1.);
	float value = getValue(p); 
	
	float voxelAlpha = alpha;// * min(1., mod(value * 4., 3.));
	vec3 voxelColor = texture1D(gradient, value).rgb;
	voxelAlpha *= min(1., length(voxelColor));
	result.rgb += result.a * voxelAlpha * voxelColor;
	result.a *= 1. - voxelAlpha;
	
	vec3 step = vertexStart - eye;
	step = normalize(step) / float(maxiter);
	step /= oneOverDx;
	for (int i = 2; i <= maxiter; i++) {
		p += step;
		if (p.x < 0. || p.y < 0. || p.z < 0. ||
			p.x > 1. || p.y > 1. || p.z > 1.) break;

		//notice if you store 1-alpha in the color you trace through the volume
		//then you can forward-trace your blending
		//instead of having to backward-trace
		//(as you would when rendering transparent stuff on top of each other)
		//this will allow you to bailout early if your transparency ever hits fully opaque
		value = getValue(p);
		
		voxelAlpha = alpha;// * min(1., mod(value * 4., 3.));
		voxelColor = texture1D(gradient, value).rgb;
		voxelAlpha *= min(1., length(voxelColor));
		result.rgb += result.a * voxelAlpha * voxelColor;
		result.a *= 1. - voxelAlpha;

		if (result.a < .01) break;
	}
result.a = 0.;
	gl_FragColor = vec4(result.rgb, 1. - result.a);
}

<? end ?>
