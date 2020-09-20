<?
local clnumber = require 'cl.obj.number'
local coord = solver.coord
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
varying vec3 viewCoord;

<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>

void main() {
	viewCoord = gl_Vertex.xyz;
	// disregard 'z' for rendering
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xy, 0., 1.);
}

<? end
if fragmentShader then ?>

<?=
solver.coord:getModuleCodeGLSL(
	'coordMapInvGLSL',
	solver.coord.vectorComponent == 'cartesian'
	and 'cartesianToCoord'	-- coord_cartesianToCoord
	or 'coord_conn_apply23'
)
?>

uniform sampler2D noiseTex;
uniform int integralMaxIter;

vec3 getTexCoordForGridCoord(vec3 gridCoord) {
	vec3 texCoord = (
		(gridCoord - solverMins) / (solverMaxs - solverMins) 
		* vec3(<?=
			clnumber(tonumber(solver.sizeWithoutBorder.x))?>, <?=
			clnumber(tonumber(solver.sizeWithoutBorder.y))?>, <?=
			clnumber(tonumber(solver.sizeWithoutBorder.z))
		?>) 
		+ vec3(<?=
			clnumber(tonumber(solver.numGhost))?>, <?=
			clnumber(tonumber(solver.numGhost))?>, <?=
			clnumber(tonumber(solver.numGhost))
		?>)
	) * vec3(<?=
		clnumber(1 / tonumber(solver.gridSize.x))?>, <?=
		clnumber(1 / tonumber(solver.gridSize.y))?>, <?=
		clnumber(1 / tonumber(solver.gridSize.z))
	?>);
	texCoord.z = viewCoord.z;
	return texCoord;
}

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	//TODO make this a flag:
	vec3 gridCoord = vec3(viewCoord.xy, 0.);
	if (useCoordMap) {
		gridCoord = coordMapInv(gridCoord);
	} else {
		gridCoord = .5 * (gridCoord + 1.) * (solverMaxs - solverMins) + solverMins; 
	}
	
	if (gridCoord.x < solverMins.x || gridCoord.x > solverMaxs.x ||
		gridCoord.y < solverMins.y || gridCoord.y > solverMaxs.y
	) {
		discard;
	}

<? local ds = 1/solver.app.drawVectorLICNoiseSize ?>
	float licMag = texture2D(noiseTex, getTexCoordForGridCoord(gridCoord).xy).r;
	float totalWeight = 1.;
	<? for dir=-1,1,2 do ?>{
		vec3 pos = gridCoord;
		vec3 vel = vec3(0., 0., 0.);
		vec3 last_dPos_ds = vec3(0., 0., 0.);
		for (int iter = 0; iter < integralMaxIter; ++iter) {
			
			vec3 dPos_ds = getTex(getTexCoordForGridCoord(pos)).xyz;
			
<? if solver.coord.vectorComponent == 'cartesian' then ?>
			//If the vector components are cartesian but our coords are chart-based
			//then convert the vector back to the chart based coordinates.
			vec3 dVel_ds = vec3(0., 0., 0.);
			dPos_ds = coord_cartesianToCoord(dPos_ds, pos);
<? else ?>
			// u'^i = -Gamma^i_jk u^j u^k
			vec3 dVel_ds = -coord_conn_apply23(
				vel,	//u^j
				vel,	//u^k
				pos		//x
			);
<? end ?>

			if (iter > 0 && dot(dPos_ds, last_dPos_ds) < 0.) break;

			float uMag = length(dPos_ds);
			if (uMag == 0.) break;

			pos += dPos_ds * (<?=ds * dir?> / uMag);
			vel += dVel_ds * (<?=ds * dir?> / uMag);
			
			float f = float(iter + 1) / float(integralMaxIter+1);
			float weight = smoothstep(1., 0., f);
			licMag += texture2D(noiseTex, getTexCoordForGridCoord(pos).xy).r * weight;
			totalWeight += weight;
			last_dPos_ds = dPos_ds;
		}
	}<? end ?>
	
	licMag /= totalWeight;
	
	//add some contrast
	licMag = smoothstep(.2, .8, licMag);
	
	//color by magnitude
	float fieldMagn = length(getTex(getTexCoordForGridCoord(gridCoord)).xyz);
// where did these flickers come from?  I'm getting an unnderrun of values...
fieldMagn = clamp(fieldMagn, valueMin, valueMax);
	vec4 licColorTimesGradColor = getGradientColor(fieldMagn);
	licColorTimesGradColor.rgb *= licMag;
	
	gl_FragColor = licColorTimesGradColor;
}

<? end ?>
