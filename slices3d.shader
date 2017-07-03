varying vec3 texCoord;

<?
if vertexShader then
?>

uniform vec3 mins, maxs;
void main() {
	texCoord = gl_Vertex.xyz;
	
	vec4 x = gl_Vertex;
	x.xyz *= maxs - mins;
	x.xyz += mins;
	x = vec4(coordMap(x.xyz), x.w);
	gl_Position = gl_ModelViewProjectionMatrix * x;
}

<?
end
if fragmentShader then
?>

//1/log(10)
#define _1_LN_10 0.4342944819032517611567811854911269620060920715332
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useLog;
uniform float valueMin, valueMax;
uniform sampler3D tex;
uniform sampler1D gradientTex;

uniform vec3 normal;
uniform float alpha;
uniform float alphaGamma;

void main() {
	//vec4 worldPos = gl_ModelViewMatrix * vec4(pos,1.);

	float value = texture3D(tex, texCoord).r;
	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
	//TODO insert the gradient tex size
	value = value * <?=clnumber(gradTexWidth-1)?> / <?=clnumber(gradTexWidth)?>;
	vec4 voxelColor = texture1D(gradientTex, value);
	voxelColor.a = pow(alpha, alphaGamma);
	
	//calculate normal in screen coordinates
	vec4 n = gl_ModelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	gl_FragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}

<?
end
?>
