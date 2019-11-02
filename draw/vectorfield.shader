<?
local clnumber = require 'cl.obj.number'
local coord = solver.coord
local app = solver.app
?>
varying vec4 color;

<? -- coordMapInv isn't used by volumetric.shader, but in coord/sphere-log-radial.lua it does have important predefined functions (sinh, etc) which are needed by coordMap ?>
<?=coord:getCoordMapInvGLSLCode()?>

<? if vertexShader then ?>
<?=coord:getCoordMapGLSLCode()?>

uniform vec3 mins, maxs;
uniform float scale;

<? if solver.dim < 3 then ?>
uniform sampler2D tex;
<? else ?>
uniform sampler3D tex;
<? end ?>


//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useLog;
uniform float valueMin, valueMax;
uniform sampler1D gradientTex;

void main() {
	
	//manifold coords
	vec3 pt = gl_MultiTexCoord1.xyz * (maxs - mins) + mins;

<? if solver.dim < 3 then ?> 
	
	vec3 dir = texture2D(tex, gl_MultiTexCoord0.xy).rgb;
	dir = cartesianFromCoord(dir, pt);	
	float value = length(dir);
	dir /= value;

	vec3 tv = vec3(-dir.y, dir.x, 0.);

<? else ?>
	
	vec3 dir = texture3D(tex, gl_MultiTexCoord0.xyz).rgb;
	dir = cartesianFromCoord(dir, pt);
	float value = length(dir);
	dir /= value;
	
	vec3 vx = vec3(0., -dir.z, dir.y);
	vec3 vy = vec3(dir.z, 0., -dir.x);
	vec3 vz = vec3(-dir.y, dir.x, 0.);
	float lxsq = dot(vx,vx);
	float lysq = dot(vy,vy);
	float lzsq = dot(vz,vz);
	vec3 tv;
	if (lxsq > lysq) {		//x > y
		if (lxsq > lzsq) {	//x > z, x > y
			tv = vx;
		} else {			//z > x > y
			tv = vz;
		}
	} else {				//y >= x
		if (lysq > lzsq) {	//y >= x, y > z
			tv = vy;
		} else {			// z > y >= x
			tv = vz;
		}
	}

<? end ?>

	//TODO don't just normalize dir, but scale down by the view size ... ? or the dx?
	//TODO put this code somewhere that heatmap2d, vectorfield, and whatever else can get to it
	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
//TODO make this a flag, whether to normalize vectors or not? 
	float valuescale = scale;// * clamp(value, 0., 1.);
	value = (value * <?=clnumber(app.gradientTex.width-1)?> + .5) / <?=clnumber(app.gradientTex.width)?>;
	color = texture1D(gradientTex, value);


	vec2 offset = gl_Vertex.xy;

	//cartesian coords
	vec3 v = coordMap(pt) + valuescale * (offset.x * dir + offset.y * tv);
	
	gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * vec4(v, 1.);
}

<?
end
if fragmentShader then
?>

void main() {
	gl_FragColor = color;
}

<? end ?>
