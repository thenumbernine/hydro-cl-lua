varying vec2 texCoord;

#ifdef VERTEX_SHADER

void main() {
	texCoord = gl_MultiTexCoord0.xy;
	vec4 x = vec4(coordMap(gl_Vertex.xyz), gl_Vertex.w);
	gl_Position = gl_ModelViewProjectionMatrix * x;
}

#endif	//VERTEX_SHADER

#ifdef FRAGMENT_SHADER

//1/log(10)
#define _1_LN_10 0.4342944819032517611567811854911269620060920715332
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useLog;
uniform float valueMin, valueMax;
uniform sampler2D tex;
uniform sampler1D gradient;
uniform float alpha;

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
	gl_FragColor = vec4(texture1D(gradient, value).rgb, alpha);
}

#endif	//FRAGMENT_SHADER
