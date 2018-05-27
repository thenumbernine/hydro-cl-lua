varying vec2 texCoord;

<?  if vertexShader then ?>

void main() {
	texCoord = gl_MultiTexCoord0.xy;

	//with coordinate mapping
	vec4 x = vec4(coordMap(gl_Vertex.xyz), gl_Vertex.w);
	//without (rectangular display, even with curvilinear coordinates)	
	//vec4 x = gl_Vertex;
	
	gl_Position = gl_ModelViewProjectionMatrix * x;
}

<? end
if fragmentShader then ?>

//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useLog;
uniform float valueMin, valueMax;
uniform sampler2D tex;
uniform sampler1D gradientTex;

void main() {
	float value = texture2D(tex, texCoord).r;
	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
	value = (value * <?=clnumber(gradTexWidth-1)?> + .5) / <?=clnumber(gradTexWidth)?>;
	gl_FragColor = texture1D(gradientTex, value);
}

<? end ?>
