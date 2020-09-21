#version 460

<?
local coord = solver.coord
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
<?=varying?> vec3 viewCoord;

<?=coord:getModuleCodeGLSL("coordMapGLSL", "coordMapInvGLSL")?>
<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>

attribute vec4 vertex;

void main() {
	vec4 v = vertex;
	viewCoord = v.xyz;
	gl_Position = modelViewProjectionMatrix * v;
}

<? end
if fragmentShader then ?>

out vec4 fragColor;

uniform sampler2D noiseTex;
uniform int integralMaxIter;

<?=
solver.coord:getModuleCodeGLSL(
	'coordMapInvGLSL',
	solver.coord.vectorComponent == 'cartesian'
	and 'cartesianToCoord'	-- coord_cartesianToCoord
	or 'coord_conn_apply23'
)
?>

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	//TODO make this a flag:
	vec3 chartCoord = viewCoord;
	chartCoord.z = displayFixed.y;
	chartCoord = quatRotate(displaySliceAngle, chartCoord);
	chartCoord = worldToChartCoord(chartCoord);

	vec3 texCoord = chartToTexCoord(chartCoord);
#if 1
	if (texCoord.x < 0. || texCoord.x > 1.) discard;
	if (displayDim > 1 && (texCoord.y < 0. || texCoord.y > 1.)) discard;
	if (displayDim > 2 && (texCoord.z < 0. || texCoord.z > 1.)) discard;
#endif

<? local ds = 1/solver.app.drawVectorLICNoiseSize ?>
	float licMag = texture2D(noiseTex, chartToNoGhostCoord(chartCoord).xy).r;
	float totalWeight = 1.;
	<? for dir=-1,1,2 do ?>{
		vec3 pos = chartCoord;
		vec3 vel = vec3(0., 0., 0.);
		vec3 last_dPos_ds = vec3(0., 0., 0.);
		for (int iter = 0; iter < integralMaxIter; ++iter) {
			
			vec3 dPos_ds = getTex(chartToNoGhostCoord(pos)).xyz;
			
<? if solver.coord.name == 'cartesian' then ?>
			vec3 dVel_ds = vec3(0., 0., 0.);
<? elseif solver.coord.vectorComponent == 'cartesian' then ?>
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
			licMag += texture2D(noiseTex, chartToNoGhostCoord(pos).xy).r * weight;
			totalWeight += weight;
			last_dPos_ds = dPos_ds;
		}
	}<? end ?>
	
	licMag /= totalWeight;
	
	//add some contrast
	licMag = smoothstep(.2, .8, licMag);
	
	//color by magnitude
	float fieldMagn = length(getTex(chartToNoGhostCoord(chartCoord)).xyz);
// where did these flickers come from?  I'm getting an unnderrun of values...
fieldMagn = clamp(fieldMagn, valueMin, valueMax);
	vec4 licColorTimesGradColor = getGradientColor(fieldMagn);
	licColorTimesGradColor.rgb *= licMag;
	
	fragColor = licColorTimesGradColor;
}

<? end ?>
