<?
local clnumber = require 'cl.obj.number'
local app = solver.app
local coord = solver.coord
?>

varying vec2 viewCoord;

<?=coord:getCoordMapGLSLCode()?>
<?=coord:getCoordMapInvGLSLCode()?>

<? if vertexShader then ?>

void main() {
	viewCoord = gl_Vertex.xy;
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}

<? end
if fragmentShader then ?>

//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useLog;
uniform float valueMin, valueMax;
uniform vec2 texCoordMax;
uniform sampler2D tex;
uniform sampler1D gradientTex;

uniform vec2 solverMins, solverMaxs;

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	vec2 gridCoord = coordMapInv(vec3(viewCoord.xy, 0.)).xy;
	
	if (gridCoord.x < solverMins.x || gridCoord.x > solverMaxs.x ||
		gridCoord.y < solverMins.y || gridCoord.y > solverMaxs.y
	) {
		discard;
	}

const float epsilon = 1e-2;
if (abs(gridCoord.x - solverMins.x) < epsilon ||
	abs(gridCoord.x - solverMaxs.x) < epsilon ||
	abs(gridCoord.y - solverMins.y) < epsilon ||
	abs(gridCoord.y - solverMaxs.y) < epsilon)
{
	gl_FragColor = vec4(1.,1.,1.,1.);
	return;
}

	vec2 texCoord = vec2(
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * <?=clnumber(solver.sizeWithoutBorder.x)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.x)?>,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * <?=clnumber(solver.sizeWithoutBorder.y)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.y)?>
	) * texCoordMax;

	float value = texture2D(tex, texCoord).r;
	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
	value = (value * <?=clnumber(app.gradientTex.width-1)?> + .5) / <?=clnumber(app.gradientTex.width)?>;
	gl_FragColor = texture1D(gradientTex, value);
}

<? end ?>
