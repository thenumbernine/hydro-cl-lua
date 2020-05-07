<?
local clnumber = require 'cl.obj.number'
local app = solver.app
local coord = solver.coord
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
varying vec3 viewCoord;

<?=coord:getCoordMapInvGLSLCode()?>

<? if vertexShader then ?>

void main() {
	viewCoord = gl_Vertex.xyz;
	// disregard 'z' for rendering
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xy, 0., 1.);
}

<? end
if fragmentShader then ?>

//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useCoordMap;
uniform bool useLog;
uniform float valueMin, valueMax;

uniform sampler1D gradientTex;
uniform sampler2D noiseTex;
uniform vec2 solverMins, solverMaxs;
uniform vec2 texCoordMax;
uniform int integralMaxIter;

<? if solver.dim == 3 then ?>
uniform sampler3D tex;
vec4 getTex(vec3 texCoord) {
	return texture3D(tex, texCoord);
}
<? else ?>
uniform sampler2D tex;
vec4 getTex(vec3 texCoord) {
	return texture2D(tex, texCoord.xy);
}
<? end ?>

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	//TODO make this a flag:
	vec2 gridCoord = viewCoord.xy;
	if (useCoordMap) {
		gridCoord = coordMapInv(vec3(gridCoord, 0.)).xy;
	} else {
		gridCoord = .5 * (gridCoord + 1.) * (solverMaxs.xy - solverMins.xy) + solverMins.xy; 
	}

	if (gridCoord.x < solverMins.x || gridCoord.x > solverMaxs.x ||
		gridCoord.y < solverMins.y || gridCoord.y > solverMaxs.y
	) {
		discard;
	}

	vec3 texCoord;
	texCoord.xy = vec2(
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * <?=clnumber(solver.sizeWithoutBorder.x)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.x)?>,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * <?=clnumber(solver.sizeWithoutBorder.y)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.y)?>
	) * texCoordMax;
	texCoord.z = viewCoord.z;

<? local ds = 1/solver.app.drawVectorLICNoiseSize ?>
	float licMag = texture2D(noiseTex, texCoord.xy).r;
	float totalWeight = 1.;
	<? for dir=-1,1,2 do ?>{
		vec3 pos = texCoord;
		for (int iter = 0; iter < integralMaxIter; ++iter) {
			float f = float(iter + 1) / float(integralMaxIter+1);
			float weight = smoothstep(1., 0., f);
			vec3 dPos_ds = normalize(getTex(pos).rgb);
			pos += dPos_ds * <?=ds * dir?>;
			licMag += texture2D(noiseTex, pos.xy).r * weight;
			totalWeight += weight;
		}
	}<? end ?>

	licMag /= totalWeight;

	//add some contrast
	//licMag = smoothstep(0., 1., licMag);


	//color by magnitude
	float fieldMagn = length(getTex(texCoord));
	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		fieldMagn = (logmap(fieldMagn) - logValueMin) / (logValueMax - logValueMin);
	} else {
		fieldMagn = (fieldMagn - valueMin) / (valueMax - valueMin);
	}
	fieldMagn = (fieldMagn * <?=clnumber(app.gradientTex.width-1)?> + .5) / <?=clnumber(app.gradientTex.width)?>;
	vec4 licColorTimesGradColor = texture1D(gradientTex, fieldMagn);
	licColorTimesGradColor.rgb *= licMag;

	gl_FragColor = licColorTimesGradColor;
}

<? end ?>
