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

uniform vec3 solverMins, solverMaxs;

<? if not vectorField2 then ?>
uniform float scale;
<? end ?>

<? if solver.dim < 3 then ?>
uniform sampler2D tex;
<? else ?>
uniform sampler3D tex;
<? end ?>

uniform int displayDim;

<?=solver:getGradientGLSLCode()?>

<? if vectorField2 then ?>
uniform sampler2D offsetTex;
<? end ?>

void main() {
	
	//manifold coords
	vec3 pt = gl_MultiTexCoord1.xyz;
	
<? 
local clnumber = require 'cl.obj.number'
if vectorField2 then 
	local dx = 1 / tonumber(solver.gridSize.x)
	local dy = 1 / tonumber(solver.gridSize.y) 
?>
	vec4 offsetField = texture2D(offsetTex, gl_MultiTexCoord1.xy);
	vec2 delta = vec2(<?=clnumber(dx)?>, <?=clnumber(dy)?>);
	pt.xy += (offsetField.xy - vec2(.5, .5)) * delta;
<? 
end -- vectorField2
?>

	pt *= (solverMaxs - solverMins);
	pt += solverMins;

<? if solver.dim < 3 then ?>
	vec3 dir = texture2D(tex, gl_MultiTexCoord0.xy).rgb;
<? else ?>
	vec3 dir = texture3D(tex, gl_MultiTexCoord0.xyz).rgb;
<? end ?>

	dir = cartesianFromCoord(dir, pt);	
	float value = length(dir);
	dir /= value;

	vec3 tv;
	if (displayDim < 3) {
		tv = vec3(-dir.y, dir.x, 0.);
	} else {	//displayDim == 3
		vec3 vx = vec3(0., -dir.z, dir.y);
		vec3 vy = vec3(dir.z, 0., -dir.x);
		vec3 vz = vec3(-dir.y, dir.x, 0.);
		float lxsq = dot(vx,vx);
		float lysq = dot(vy,vy);
		float lzsq = dot(vz,vz);
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
	}

	//TODO don't just normalize dir, but scale down by the view size ... ? or the dx?
	
	value = getGradientFrac(value);
//TODO make this a flag, whether to normalize vectors or not? 
<? if vectorField2 then ?>
	float valuescale = offsetField.w;
<? else ?>
	float valuescale = scale;// * clamp(value, 0., 1.);
<? end ?>
	value = getGradientTexCoord(value);
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
