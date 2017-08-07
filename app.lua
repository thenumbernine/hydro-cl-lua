--[[
predefined vars:
	dim=
	gridSize=
	fluxLimiter=
	boundary=
	integrator=
	eqn=
	mins=
	maxs=
	float= set to false to override double precision
--]]
local cmdline = {}
--[[
for _,k in ipairs{
	'dim', 'gridSize', 'fluxLimiter', 'boundary', 
	'integrator', 'eqn', 'mins', 'maxs', 'float'
} do
	cmdline[k] = _G[k]
	_G[k] = nil
	if cmdline[k] == nil then 
		cmdline[k] = _G[k:lower()] 
		_G[k:lower()] = nil
	end
end
--]]
for _,w in ipairs(arg or {}) do
	local k,v = w:match'^(.-)=(.*)$'
	if k then
		cmdline[k] = assert(loadstring('return '..v))()
	else
		cmdline[w] = true
	end
end

local bit = require 'bit'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local cl = require 'ffi.OpenCL'
local gl = require 'gl'
local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local math = require 'ext.math'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local ImGuiApp = require 'imguiapp'
local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'
local GLProgram = require 'gl.program'
local GLGradientTex = require 'gl.gradienttex'
local GLTex2D = require 'gl.tex2d'
local glreport = require 'gl.report'
local Font = require 'gui.font'
local vec4d = require 'ffi.vec.vec4d'
local vec3d = require 'ffi.vec.vec3d'
local tooltip = require 'tooltip'

local CartesianGeom = require 'geom.cartesian'
local CylinderGeom = require 'geom.cylinder'

local HydroCLApp = class(ImGuiApp)

HydroCLApp.title = 'Hydrodynamics in OpenCL'

-- list from https://en.wikipedia.org/wiki/Flux_limiter
HydroCLApp.limiters = table{
	{name='donor cell', code='return 0.;'},
	{name='Lax-Wendroff', code='return 1.;'},
	{name='Beam-Warming', code='return r;'},
	{name='Fromm', code='return .5 * (1. + r);'},
	{name='CHARM', code='return max(0., r) * (3. * r + 1.) / ((r + 1.) * (r + 1.));'},
	{name='HCUS', code='return max(0., 1.5 * (r + fabs(r)) / (r + 2.));'},
	{name='HQUICK', code='return max(0., 2. * (r + fabs(r)) / (r + 3.));'},
	{name='Koren', code='return max(0., min(2. * r, min((1. + 2. * r) / 3., 2.)));'},
	{name='minmod', code='return max(0., min(r, 1.));'},
	{name='Oshker', code='return max(0., min(r, 1.5));	//replace 1.5 with 1 <= beta <= 2'},
	{name='ospre', code='return .5 * (r * r + r) / (r * r + r + 1.);'},
	{name='smart', code='return max(0., min(2. * r, min(.25 + .75 * r, 4.)));'},
	{name='Sweby', code='return max(0., max(min(1.5 * r, 1.), min(r, 1.5)));	//replace 1.5 with 1 <= beta <= 2'},
	{name='UMIST', code='return max(0., min(min(2. * r, .75 + .25 * r), min(.25 + .75 * r, 2.)));'},
	{name='van Albada 1', code='return (r * r + r) / (r * r + 1.);'},
	{name='van Albada 2', code='return 2. * r / (r * r + 1.);'},
	{name='van Leer', code='return max(0., r) * 2. / (1. + r);'},	--return (r + fabs(r)) / (1. + fabs(r));
	{name='monotized central', code='return max(0., min(2., min(.5 * (1. + r), 2. * r)));'},
	{name='superbee', code='return max((real)0., (real)max((real)min((real)1., (real)2. * r), (real)min((real)2., r)));'},
	{name='Barth-Jespersen', code='return .5 * (r + 1.) * min(1., min(4. * r / (r + 1.), 4. / (r + 1.)));'},
}
HydroCLApp.limiterNames = HydroCLApp.limiters:map(function(limiter) return limiter.name end)

-- setup for the solver
-- override this for specific experiments
-- this can't override float vs double precision yet
function HydroCLApp:setup()
	
	--[[ self-gravitation simulation of planet earth?
	cmdline = {
		dim = 1,
		mins = {-12e+6},
		maxs = {12e+6},
		gridSize = 256,
		initState = 'self-gravitation test 1',
	}
	--]]
	
	
	-- create this after 'real' is defined
	--  specifically the call to 'refreshGridSize' within it
	local dim = 2
	local args = {
		app = self, 
		eqn = cmdline.eqn,
		dim = cmdline.dim or dim,
		
		integrator = cmdline.integrator or 'forward Euler',	
		--integrator = 'Runge-Kutta 2',
		--integrator = 'Runge-Kutta 2 Heun',
		--integrator = 'Runge-Kutta 2 Ralston',
		--integrator = 'Runge-Kutta 3',
		--integrator = 'Runge-Kutta 4',
		--integrator = 'Runge-Kutta 4, 3/8ths rule',
		--integrator = 'Runge-Kutta 2, TVD',
		--integrator = 'Runge-Kutta 2, non-TVD',
		--integrator = 'Runge-Kutta 3, TVD',
		--integrator = 'Runge-Kutta 4, TVD',
		--integrator = 'Runge-Kutta 4, non-TVD',
		--integrator = 'backward Euler',

		--fixedDT = .0001,
		cfl = .25/dim,

		fluxLimiter = cmdline.fluxLimiter or 'superbee',

		--usePLM = true,	-- piecewise-linear slope limiter
		--fluxLimiter = 'donor cell',
		--slopeLimiter = 'minmod',
		
		-- [[ Cartesian
		geometry = 'cartesian',
		mins = cmdline.mins or {-1, -1, -1},
		maxs = cmdline.maxs or {1, 1, 1},
		-- 256^2 = 2^16 = 2 * 32^3
		gridSize = ({
			{256,1,1},
			{256,256,1},
			{32,32,32},
		})[dim],
		boundary = {
			xmin=cmdline.boundary or 'freeflow',
			xmax=cmdline.boundary or 'freeflow',
			ymin=cmdline.boundary or 'freeflow',
			ymax=cmdline.boundary or 'freeflow',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]
		--[[ cylinder
		geometry = 'cylinder',
		mins = cmdline.mins or {.1, 0, -1},
		maxs = cmdline.maxs or {1, 2*math.pi, 1},
		gridSize = ({
			{128, 1, 1}, -- 1D
			{32, 128, 1}, -- 2D
			{16, 64, 16}, -- 3D
		})[dim],
		boundary = {
			xmin=cmdline.boundary or 'freeflow',		-- hmm, how to treat the r=0 boundary ...
			xmax=cmdline.boundary or 'freeflow',
			ymin=cmdline.boundary or 'periodic',
			ymax=cmdline.boundary or 'periodic',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]
		--[[ sphere
		geometry = 'sphere',
		mins = cmdline.mins or {0, -math.pi, .5},
		maxs = cmdline.maxs or {math.pi, math.pi, 1},
		gridSize = {
			cmdline.gridSize or 64,
			cmdline.gridSize or 128,
			cmdline.gridSize or 64,
		},
		boundary = {
			xmin=cmdline.boundary or 'freeflow',
			xmax=cmdline.boundary or 'freeflow',
			ymin=cmdline.boundary or 'freeflow',
			ymax=cmdline.boundary or 'freeflow',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]
		--[[ sphere1d
		geometry = 'sphere1d',
		mins = cmdline.mins or {1, -math.pi, .5},
		maxs = cmdline.maxs or {100, math.pi, 1},
		gridSize = {
			cmdline.gridSize or 256,
			cmdline.gridSize or 128,
			cmdline.gridSize or 64,
		},
		boundary = {
			xmin=cmdline.boundary or 'freeflow',
			xmax=cmdline.boundary or 'freeflow',
			ymin=cmdline.boundary or 'freeflow',
			ymax=cmdline.boundary or 'freeflow',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]

		--useGravity = true,

		-- no initial state means use the first
		--initState = cmdline.initState,
		
		-- Euler / SRHD / MHD initial states:
		--initState = 'constant',
		--initState = 'constant with motion',
		--initState = 'linear',
		--initState = 'gaussian',
		--initState = 'advect wave',
		--initState = 'sphere',
		--initState = 'rarefaction wave',
		
		initState = 'Sod',
		--initState = 'Sedov',
		--initState = 'Kelvin-Hemholtz',
		--initState = 'Rayleigh-Taylor',
		--initState = 'Colella-Woodward',
		--initState = 'double mach reflection',
		--initState = 'square cavity',
		--initState = 'shock bubble interaction',

		--initState = 'configuration 1',
		--initState = 'configuration 2',
		--initState = 'configuration 3',
		--initState = 'configuration 4',
		--initState = 'configuration 5',
		--initState = 'configuration 6',

		-- self-gravitation tests:
		--initState = 'self-gravitation test 1',
		--initState = 'self-gravitation test 1 spinning',
		--initState = 'self-gravitation test 2',
		--initState = 'self-gravitation test 2 orbiting',
		--initState = 'self-gravitation test 4',
		--initState = 'self-gravitation soup',
		
		-- those designed for SRHD / GRHD:
		--initState = 'relativistic shock reflection',			-- not working.  these initial conditions are constant =P
		--initState = 'relativistic blast wave test problem 1',
		--initState = 'relativistic blast wave test problem 2',
		--initState = 'relativistic blast wave interaction',
	
		-- MHD-only init states: (that use 'b')
		--initState = 'Brio-Wu',
		--initState = 'Orszag-Tang',
		
		-- EM:
		--initState = 'Maxwell default',
		--initState = 'Maxwell scattering around cylinder',
		--initState = 'Maxwell wire',
		--initState = 'Maxwell FDTD test',
		
		--initState = 'two-fluid EMHD soliton ion',
		--initState = 'two-fluid EMHD soliton electron',
		--initState = 'two-fluid EMHD soliton maxwell',
	
		-- GR
		--initState = 'gauge shock wave',
		--initState = 'Alcubierre warp bubble',
		--initState = 'Schwarzschild black hole',
		--initState = 'black hole - isotropic',
		--initState = 'binary black holes - isotropic',
		--initState = 'stellar model',
		--initState = 'stellar model 2',
		--initState = 'stellar model 3',
	}
	
	-- HD - Roe
	--self.solvers:insert(require 'solver.euler-roe'(args))

	-- the same as solver.euler-roe:
	-- TODO specify behavior operations (selfgrav, nodiv, etc) in eqn, and apply them to the solver
	--self.solvers:insert(require 'solver.selfgrav'(require 'solver.roe')(table(args, {eqn='euler'})))

	-- HD - Burgers
	-- f.e. and b.e. are working, but none of the r.k. integrators 
	-- PLM isn't implemented yet
	-- neither is source term / poisson stuff
	--self.solvers:insert(require 'solver.euler-burgers'(args))

	-- SRHD.  
	-- rel blast wave 1 & 2 works in 1D at 256 with superbee flux lim
	-- rel blast wave interaction works with superbee flux lim in 1D works at 256, fails at 1024 with float (works with double)
	-- 	256x256 double fails with F.E., RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths,
	-- rel blast wave 1 doesn't work in 64x64. with superbee flux lim
	-- rel blast wave 2 with superbee flux lim, Roe solver, works at 64x64 with forward euler
	-- 	at 256x256 fails with F.E, RK2, RK2-non-TVD., RK3-TVD, RK4, RK4-TVD, RK4-non-TVD 
	--    but works with RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths
	-- Kelvin-Hemholtz works for all borderes freeflow, float precision, 256x256, superbee flux limiter
	--self.solvers:insert(require 'solver.srhd-roe'(args))
	
	-- GRHD
	-- right now this is just like srhd except extended by Font's eqns
	-- this has plug-ins for ADM metric alpha, beta, gammas, but I need to make a composite solver to combine it with GR equations. 
	self.solvers:insert(require 'solver.grhd-roe'(args))

	-- GRHD+GR
	-- here's the GRHD solver with the BSSNOK plugged into it
	--self.solvers:insert(require 'solver.gr-hd-separate-behavior'(args))

	-- M+HD. 
	-- with superbee flux lim:  
	-- Brio-Wu works in 1D at 256, works in 2D at 64x64 in a 1D profile in the x and y directions.
	-- Orszag-Tang with forward Euler integrator fails at 64x64 around .7 or .8
	-- 		but works with 'Runge-Kutta 4, TVD' integrator at 64x64
	-- 		RK4-TVD fails at 256x256 at just after t=.5
	--		and works fine with backwards Euler 
	-- when run alongside HD Roe solver, curves don't match (different heat capacity ratios?)
	--		but that could be because of issues with simultaneous solvers.
	--self.solvers:insert(require 'solver.mhd-roe'(args))
	
	-- EM
	--self.solvers:insert(require 'solver.maxwell-roe'(args))
	
	-- EM+HD
	-- I'm having some memory issues with two solvers running simultanously .. 
	--self.solvers:insert(require 'solver.twofluid-emhd-separate-roe'(args))
	-- so to try and get around that, here the two are combined into one solver:
	--self.solvers:insert(require 'solver.twofluid-emhd-roe'(args))

	-- GR
	--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm1d_v1'})))
	--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm1d_v2'})))
	--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d'})))
	--
	-- the BSSNOK solver works similar to the adm3d for the warp bubble simulation
	--  but something gets caught up in the freeflow boundary conditions, and it explodes
	-- so I have set constant Minkowski boundary conditions?
	-- the BSSNOK solver sometimes explodes / gets errors / nonzero Hamiltonian constraint for forward euler
	-- however they tend to not explode with backward euler ... though these numerical perturbations still appear, but at least they don't explode
	--self.solvers:insert(require 'solver.bssnok-fd'(args))
	
	-- TODO GR+HD by combining the SR+HD 's alphas and gammas with the GR's alphas and gammas
end

local useClipPlanes
if useClipPlanes then
-- TODO put all the display3D_slices stuff in its own file
local rotateClip = ffi.new('int[1]', 0)
local function makeDefaultPlane(i)
	assert(i >= 1 and i <= 4)
	local plane = vec4d(0,0,0,0)
	plane:ptr()[math.min(i,3)-1] = -1
	return plane
end
local clipInfos = range(4):map(function(i)
	return {
		enabled = i == 3,
		plane = makeDefaultPlane(i),
	}
end)
end

function HydroCLApp:initGL(...)
	HydroCLApp.super.initGL(self, ...)

	-- TODO favor cl_khr_gl_sharing, cl_khr_fp64, cl_khr_3d_image_writes
	self.env = CLEnv{
		verbose = true,
		precision = cmdline.float and 'float' or nil,
		cpu = cmdline.cpu,
	}
	print(self.env.platform:getName())
	print(self.env.device:getName())

	self.is64bit = self.env.real == 'double'
	self.useGLSharing = self.env.useGLSharing
	self.device = self.env.device
	self.ctx = self.env.ctx
	self.cmds = self.env.cmds
	self.real = self.env.real

	ffi.cdef('typedef '..self.real..' real;')

	ffi.cdef(file['math.h'])

	self.solvers = table()

	self:setup()

	-- This only looks good when overlaying vector fields on top of other graphs.
	-- When it comes to separate variables, they usually look better apart.
	self.displayAllTogether = self.solvers[1] and self.solvers[1].dim > 1 or false


	local gradTexWidth = 1024
	self.gradientTex = GLGradientTex(gradTexWidth, {
	-- [[ white, rainbow, black
		{0,0,0,.5},	-- black ... ?
		{0,0,1,1},	-- blue
		{0,1,1,1},	-- cyan
		{0,1,0,1},	-- green
		{1,1,0,1},	-- yellow
		{1,.5,0,1},	-- orange
		{1,0,0,1},	-- red
		{1,1,1,1},	-- white
	--]]
	--[[ stripes 
		range(32):map(function(i)
			return ({
				{0,0,0,0},
				{1,1,1,1},
			})[i%2+1]
		end):unpack()
	--]]
	}, false)
	-- don't wrap the colors, but do use GL_REPEAT
	self.gradientTex:setWrap{s = gl.GL_REPEAT}

	local graphShaderCode = file['graph.shader']
	self.graphShader = GLProgram{
		vertexCode = '#define VERTEX_SHADER\n'..graphShaderCode,
		fragmentCode = '#define FRAGMENT_SHADER\n'..graphShaderCode,
		uniforms = {
			tex = 0,
			scale = 1,
			ambient = 1,
		},
	}
	
	local code = file['heatmap2d.shader']
	for _,solver in ipairs(self.solvers) do
		solver.heatMap2DShader = GLProgram{
			vertexCode = template(
				table{
					solver.geometry:getCoordMapGLSLCode(),
					code, 
				}:concat'\n',
				{
					vertexShader = true,
				}
			),
			fragmentCode = template(code, {
				fragmentShader = true,
				clnumber = clnumber,
				gradTexWidth = gradTexWidth,
			}),
			uniforms = {
				valueMin = 0,
				valueMax = 0,
				tex = 0,
				gradientTex = 1,
			},
		}

		if solver.dim == 3 then
			-- raytracing (stalling)
			
			local maxiter = math.max(tonumber(solver.gridSize.x), tonumber(solver.gridSize.y), tonumber(solver.gridSize.z))
			local code = file['volumetric.shader']
			local volumeRayShader = GLProgram{
				vertexCode = template(
					table{
						solver.geometry:getCoordMapGLSLCode(),
						code,
					}:concat'\n',
					{
						vertexShader = true,
					}
				),
				fragmentCode = template(code, {
					fragmentShader = true,
				}),
				uniforms = {
					tex = 0,
					gradient = 1,
					maxiter = maxiter,
					oneOverDx = {(solver.maxs - solver.mins):unpack()},
				},
			}
			solver.volumeRayShader = volumeRayShader

			-- volume slices

			local code = file['slices3d.shader']
			solver.volumeSliceShader = GLProgram{
				vertexCode = template(
					table{
						solver.geometry:getCoordMapGLSLCode(),
						code,
					}:concat'\n',
					{
						vertexShader = true,
					}
				),
				fragmentCode = template(code, {
					fragmentShader = true,
					clnumber = clnumber,
					gradTexWidth = gradTexWidth,
					clipInfos = useClipPlanes and clipInfos or nil,
				}),
				uniforms = {
					volTex = 0,
					gradientTex = 1,
					valueMin = 0,
					valueMax = 0,
				},
			}
		end
	
		solver.vectorFieldShader = GLProgram{
			vertexCode = template([[
uniform vec3 mins, maxs;
uniform float scale;

<? if dim < 3 then ?>
uniform sampler2D tex;
<? else ?>
uniform sampler3D tex;
<? end ?>

void main() {
<? if dim < 3 then ?> 
	vec3 dir = texture2D(tex, gl_MultiTexCoord0.xy).rgb;
	vec3 tv = vec3(-dir.y, dir.x, 0.);
<? else ?>
	vec3 dir = texture3D(tex, gl_MultiTexCoord0.xyz).rgb;
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

	vec2 offset = gl_Vertex.xy;
	vec3 v = gl_MultiTexCoord1.xyz * (maxs - mins) + mins + scale * (offset.x * dir + offset.y * tv);
	gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * vec4(v, 1.);
}
]], {dim=solver.dim}),
			fragmentCode = [[
void main() {
	//TODO uniform & use graph color
	gl_FragColor = vec4(1.,1.,1.,1.);
}
]],
			uniforms = {
				scale = 1,
				tex = 0,
			},
		}

	end

	self.isobarShader = GLProgram{
		vertexCode = [[
varying vec4 color;
void main() {
	color = gl_Color;
	gl_Position = ftransform();
}
]],
		fragmentCode = [[
varying vec4 color;
void main() {
	float dx_dz = dFdx(gl_FragCoord.z);
	float dy_dz = dFdy(gl_FragCoord.z);

	vec3 n = normalize(vec3(dx_dz, dy_dz, 10.));

	gl_FragColor = vec4(color.rgb * n.z, color.a);
}
]],
	}

	-- [[ need to get image loading working
	local fonttex = GLTex2D{
		filename = 'font.png',
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		magFilter = gl.GL_LINEAR,
	}
	if not pcall(function()
		gl.glGenerateMipmap(gl.GL_TEXTURE_2D) 
	end) then
		gl.glTexParameteri(fonttex.target, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST)
		gl.glTexParameteri(fonttex.target, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR) 
	end
	self.font = Font{tex = fonttex}
	--]]


	self.orthoView = require 'view.ortho'()
	self.frustumView = require 'view.frustum'()
	self.view = (#self.solvers > 0 and self.solvers[1].dim == 3) and self.frustumView or self.orthoView
end


--[[
rendering:
graph variables - many - 1D and 2D ... and hypercoordinates in 3D?
heatmap variable - 2D only
volumetric variable - 3D only 
point cloud? i.e. cheap volumetric?
slices? another cheap volumetric?
--]]

--[[

display options
	ortho vs frustum for 2D, (ortho only for 1D, frustum-only for 3D?)
	graph variable(s?) for 1D and 2D
	1D graph: whether the y axis is auto-scaling or fixed
	heatmap variable(s?) for 2D
	2D heatmap: whether the color range is auto-scaling or fixed
	2D graph: whether the z-axis is auto-scaling or fixed
	volumetric variable(s?) for 3D
	3D volumetric: whether the color range is auto-scaling or fixed
	3D volumetric alpha / gamma values

gui:
	checkbox next to each variable for quickly adding and removing views (based on the first / primary var in the view)
	list of views, with their details next to them, for adding and removing more variables

--]]

local View = class()

local Display = class()

function Display:init(args)
	self.ortho = true

	self.mins = vec3(-1,-1,-1)
	self.maxs = vec3(1,1,1)
	
	--[[
	var name, does it contribute to the graph space autoscale range? 
	--]]
	self.graphVars = table()
	
	--[[
	var name, heatmap scale min/max, heatmap autoscale? texture, alpha
	--]]
	self.heatMapVars = table()
end

-- whether to output a text file that holds all variable ranges
-- TODO how about some way to track any variable in the gui?
local dumpFile = {
	enabled = false,
}

function dumpFile:update(app, t)
	if not self.enabled then return end

	local f = self.file
	if not f then
		f = io.open('var-ranges.txt', 'w')
		self.file = f
		
		-- don't change any vars while outputting or else your rows won't match your header
		f:write'#t'
		
		for _,solver in ipairs(app.solvers) do 
			--[[ display variables:
			for _,var in ipairs(solver.displayVars) do
				if var.enabled then
					f:write('\t',var.name..'_min')
					f:write('\t',var.name..'_max')
				end
			end
			--]]
			-- [[ gmres error
			f:write'\tgmres_err'
			f:write'\tgmres_iter'
			--]]
		end
		f:write'\n'
		f:flush()
	end
	
	f:write(t)
	for _,solver in ipairs(app.solvers) do 
		--[[ display variables:
		for _,var in ipairs(solver.displayVars) do
			if var.enabled then
				local ymin, ymax = solver:calcDisplayVarRange(var)
				f:write(('\t%.16e'):format(ymin))
				f:write(('\t%.16e'):format(ymax))
			end
		end
		--]]
		-- [[ gmres error
		f:write('\t'..solver.integrator.last_err)
		f:write('\t'..solver.integrator.last_iter)
		--]]
	end
	f:write'\n'
	f:flush()
end


HydroCLApp.running = nil

function HydroCLApp:update(...)
	if self.running then
		if self.running == 'step' then 
			print('performing single step...')
			self.running = nil
		end

		-- update the one furthest behind
		local oldestSolver = self.solvers:inf(function(a,b)
			return a.t < b.t
		end)
		if oldestSolver then 
			oldestSolver:update() 
	
			-- TODO should the time be oldestSolver.t after oldestSolver just updated?
			-- or - if dumpFile is enabled - should we re-search-out the oldest solver and use its time?
			dumpFile:update(self, oldestSolver.t)
		
			if self.exitTime and oldestSolver.t > self.exitTime then
				self:requestExit()
			end
		end
	else	
		-- clear all 'lastFrameTime's of solvers so the rough fps calcs don't get messed with
		for _,solver in ipairs(self.solvers) do
			solver.lastFrameTime = nil
		end
	end

	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	
	local w, h = self:size()

	local varNamesEnabled = table()
	local varNamesEnabledByName = {}
	for _,solver in ipairs(self.solvers) do
		for i,var in ipairs(solver.displayVars) do
			if var.enabled then
				if not varNamesEnabledByName[var.name] then
					varNamesEnabled:insert(var.name)
					varNamesEnabledByName[var.name] = true
				end
			end
		end
	end

	local graphsWide = math.ceil(math.sqrt(#varNamesEnabled))
	local graphsHigh = math.ceil(#varNamesEnabled/graphsWide)
	local graphCol = 0
	local graphRow = 0

	if self.displayAllTogether then
		graphsWide = 1
		graphsHigh = 1
	end

	local ar = (w / graphsWide) / (h / graphsHigh)

	local useLog
	local vectorField
	for _,varName in ipairs(varNamesEnabled) do
		local xmin, xmax, ymin, ymax
		for _,solver in ipairs(self.solvers) do
			local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
			
			if varIndex
			and var.enabled
			--and solver.visiblePtr and solver.visiblePtr[0] 
			then
				useLog = var.useLog
				vectorField = var.vectorField
				
				local solverxmin, solverxmax = solver.mins[1], solver.maxs[1]
				solverxmin, solverxmax = 1.1 * solverxmin - .1 * solverxmax, 1.1 * solverxmax - .1 * solverxmin
				if solver.dim > 1 then
					local center = .5 * (solverxmin + solverxmax)
					solverxmin = (solverxmin - center) * ar + center
					solverxmax = (solverxmax - center) * ar + center
				end

				local solverymin, solverymax
				if solver.dim > 1 then
					solverymin, solverymax = solver.mins[2], solver.maxs[2]
					solverymin, solverymax = 1.1 * solverymin - .1 * solverymax, 1.1 * solverymax - .1 * solverymin
				else
					solverymin, solverymax = solver:calcDisplayVarRange(var)
		
					if useLog then
						solverymin = math.log(solverymin, 10)
						solverymax = math.log(solverymax, 10)
						if not math.isfinite(solverymin) then solverymin = -math.huge end
						if not math.isfinite(solverymax) then solverymax = -math.huge end
						solverymin = math.max(-30, solverymin)
						solverymax = math.max(-30, solverymax)
						solverymax = math.max(solverymax, solverymin + 1e-5)
					end			
					
					if solverymin
					and solverymax 
					and solverymin == solverymin 
					and solverymax == solverymax 
					then
						local base = 10	-- round to nearest base-10
						local scale = 10 -- ...with increments of 10
						solverymin, solverymax = 1.1 * solverymin - .1 * solverymax, 1.1 * solverymax - .1 * solverymin
						local newymin = (solverymin < 0 and -1 or 1) * (math.abs(solverymin) == math.huge and 1e+100 or base ^ math.log(math.abs(solverymin), base))
						local newymax = (solverymax < 0 and -1 or 1) * (math.abs(solverymax) == math.huge and 1e+100 or base ^ math.log(math.abs(solverymax), base))
						solverymin, solverymax = newymin, newymax
						do
							local minDeltaY = 1e-5
							local deltaY = solverymax - solverymin
							if deltaY < minDeltaY then
								solverymax = solverymax + .5 * minDeltaY
								solverymin = solverymin - .5 * minDeltaY
							end
						end
					end
				end
				
				xmin = xmin or solverxmin
				xmax = xmax or solverxmax
				ymin = ymin or solverymin
				ymax = ymax or solverymax
				
				if xmin and solverxmin then xmin = math.min(xmin, solverxmin) end
				if xmax and solverxmax then xmax = math.max(xmax, solverxmax) end
				if ymin and solverymin then ymin = math.min(ymin, solverymin) end
				if ymax and solverymax then ymax = math.max(ymax, solverymax) end
			end
		end
		
		if not xmin or not xmax or xmin ~= xmin or xmax ~= xmax then
			xmin = -5
			xmax = 5
		end
		if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
			ymin = -5
			ymax = 5
		end

		if not self.displayAllTogether then
			gl.glViewport(
				graphCol / graphsWide * w,
				(1 - (graphRow + 1) / graphsHigh) * h,
				w / graphsWide,
				h / graphsHigh)
		end

		-- TODO maybe find the first solver for this var and use it to choose 1D,2D,3D
		local dim = self.solvers[1].dim
		
		if not vectorField then
			if dim == 1 then
				self:display1D(self.solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
			elseif dim == 2 then
				self:display2D(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
			elseif dim == 3 then
				self:display3D(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
			end
		else
			--if self.enableVectorField then
			self:displayVectorField(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
			--end
		end
	
		if not self.displayAllTogether then
			graphCol = graphCol + 1
			if graphCol == graphsWide then
				graphCol = 0
				graphRow = graphRow + 1
			end
		end
	end

	gl.glViewport(0,0,w,h)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(0, w/h, 0, 1, -1, 1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	if self.font then
		local solverNames = self.solvers:map(function(solver)
			return {
				text = ('(%.3f) '):format(solver.t)..solver.name,
				color = solver.color,
			}
		end)
		local fontSizeX = .02
		local fontSizeY = .02
		local maxlen = solverNames:map(function(solverName)
			return self.font:draw{
				text = solverName.text,
				fontSize = {fontSizeX, -fontSizeY},
				dontRender = true,
				multiLine = false,
			}
		end):inf()
		for i,solverName in ipairs(solverNames) do
			self.font:draw{
				pos = {w/h-maxlen,fontSizeY*(i+1)},
				text = solverName.text,
				color = {solverName.color[1], solverName.color[2], solverName.color[3],1},
				fontSize = {fontSizeX, -fontSizeY},
				multiLine = false,
			}
		end
	end

	HydroCLApp.super.update(self, ...)
end

function HydroCLApp:showDisplayVar1D(solver, varIndex, var)
	solver:calcDisplayVarToTex(var)	
	-- display

	self.graphShader:use()
	solver:getTex(var):bind()

	gl.glUniform1f(self.graphShader.uniforms.scale.loc, 1)
	gl.glUniform1f(self.graphShader.uniforms.offset.loc, 0)
	gl.glUniform1f(self.graphShader.uniforms.ambient.loc, 1)
	gl.glUniform1i(self.graphShader.uniforms.useLog.loc, var.useLog)
	gl.glUniform2f(self.graphShader.uniforms.xmin.loc, solver.mins[1], 0)
	gl.glUniform2f(self.graphShader.uniforms.xmax.loc, solver.maxs[1], 0)
	gl.glUniform1i(self.graphShader.uniforms.axis.loc, solver.dim)
	gl.glUniform2f(self.graphShader.uniforms.size.loc, solver.gridSize.x, solver.gridSize.y)

	gl.glColor3f(table.unpack((#self.solvers > 1 and solver or var).color))
	gl.glBegin(gl.GL_LINE_STRIP)
	local step = 1
	for i=2,tonumber(solver.gridSize.x)-2,step do
		local x = (i+.5)/tonumber(solver.gridSize.x)
		gl.glVertex2f(x, 0)
	end
	gl.glEnd()
	
	solver:getTex(var):unbind()
	self.graphShader:useNone()
end

function HydroCLApp:display1D(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex2f(x*xstep,ymin)
		gl.glVertex2f(x*xstep,ymax)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex2f(xmin,y*ystep)
		gl.glVertex2f(xmax,y*ystep)
	end
	gl.glEnd()

	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex2f(xmin, 0)
	gl.glVertex2f(xmax, 0)
	gl.glVertex2f(0, ymin)
	gl.glVertex2f(0, ymax)
	gl.glEnd()

	-- display here
	for _,solver in ipairs(solvers) do
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			self:showDisplayVar1D(solver, varIndex, var)
		end

		if self.font then
			local fontSizeX = (xmax - xmin) * .05
			local fontSizeY = (ymax - ymin) * .05
			local ystep = ystep * 2
			for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
				local realY = useLog and 10^y or y
				self.font:draw{
					pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
					text=tostring(realY),
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
			self.font:draw{
				pos={xmin, ymax},
				text=varName,
				color = {1,1,1,1},
				fontSize={fontSizeX, -fontSizeY},
				multiLine=false,
			}
		end
	end
end

function HydroCLApp:display2D_Heatmap(solvers, varName, ar, graph_xmin, graph_ymin, graph_xmax, graph_ymax)
	self.view:projection(ar)
	self.view:modelview()
	if self.view.getOrthoBounds then
		xmin, xmax, ymin, ymax = self.view:getOrthoBounds(ar)
	else
		xmin, xmax, ymin, ymax = graph_xmin, graph_ymin, graph_xmax, graph_ymax
	end

--	gl.glEnable(gl.GL_DEPTH_TEST)
	
	local gridz = 0	--.1

	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex3f(x*xstep,ymin, gridz)
		gl.glVertex3f(x*xstep,ymax, gridz)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex3f(xmin,y*ystep, gridz)
		gl.glVertex3f(xmax,y*ystep, gridz)
	end
	gl.glEnd()
	
	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex3f(xmin, 0, gridz)
	gl.glVertex3f(xmax, 0, gridz)
	gl.glVertex3f(0, ymin, gridz)
	gl.glVertex3f(0, ymax, gridz)
	gl.glEnd()
			
	-- NOTICE overlays of multiple solvers won't be helpful.  It'll just draw over the last solver.
	-- I've got to rethink the visualization
	for _,solver in ipairs(solvers) do 
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			-- TODO allow a fixed, manual colormap range
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end

			solver:calcDisplayVarToTex(var)
	
			solver.heatMap2DShader:use()
			gl.glUniform1i(solver.heatMap2DShader.uniforms.useLog.loc, var.useLog)
			gl.glUniform1f(solver.heatMap2DShader.uniforms.valueMin.loc, valueMin)
			gl.glUniform1f(solver.heatMap2DShader.uniforms.valueMax.loc, valueMax)
			solver:getTex(var):bind(0)
			self.gradientTex:bind(1)
	
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
			gl.glEnable(gl.GL_BLEND)
	
			local gridScale = 4
			local udivs = math.ceil(tonumber(solver.gridSize.x)/gridScale)
			local vdivs = math.ceil(tonumber(solver.gridSize.y)/gridScale)
			if CartesianGeom.is(solver.geometry) then
				udivs, vdivs = 1, 1
			elseif CylinderGeom.is(solver.geometry) then
				udivs = 1
			end	
			for vbase=0,vdivs-1 do
				gl.glBegin(gl.GL_TRIANGLE_STRIP)
				for ui=0,udivs do
					for vofs=0,1 do
						local vi = vbase+vofs
						local u = ui/udivs
						local v = vi/vdivs
						local tx = (u * tonumber(solver.sizeWithoutBorder.x) + solver.numGhost) / tonumber(solver.gridSize.x)
						local ty = (v * tonumber(solver.sizeWithoutBorder.y) + solver.numGhost) / tonumber(solver.gridSize.y)
						gl.glTexCoord2d(tx, ty)
						gl.glVertex2d(
							u * solver.maxs[1] + (1 - u) * solver.mins[1],
							v * solver.maxs[2] + (1 - v) * solver.mins[2])
					end
				end
				gl.glEnd()
			end
	
			gl.glDisable(gl.GL_BLEND)
			
			self.gradientTex:unbind(1)
			solver:getTex(var):unbind(0)
			solver.heatMap2DShader:useNone()
			
--			gl.glDisable(gl.GL_DEPTH_TEST)

			local palwidth = (xmax - xmin) * .15 
			self.gradientTex:enable()
			self.gradientTex:bind()
			gl.glBegin(gl.GL_QUADS)
			gl.glColor3f(1,1,1)
			gl.glTexCoord1f(0) gl.glVertex2f(xmin, ymin)
			gl.glTexCoord1f(0) gl.glVertex2f(xmin + palwidth, ymin)
			gl.glTexCoord1f(1) gl.glVertex2f(xmin + palwidth, ymax)
			gl.glTexCoord1f(1) gl.glVertex2f(xmin, ymax)
			gl.glEnd()
			self.gradientTex:unbind()
			self.gradientTex:disable()
	
			if self.font then
				local fontSizeX = (xmax - xmin) * .025
				local fontSizeY = (ymax - ymin) * .025
				local ystep = 10^(math.log(ymax - ymin, 10) - 1.5)
				for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
					local value = (y - ymin) * (valueMax - valueMin) / (ymax - ymin) + valueMin
					local absvalue = math.abs(value)
					self.font:draw{
						pos={xmin * .99 + xmax * .01, y + fontSizeY * .5},
						text=(
							(absvalue > 1e+5 or absvalue < 1e-5)
							and ('%.5e'):format(value) or ('%.5f'):format(value)),
						color = {1,1,1,1},
						fontSize={fontSizeX, -fontSizeY},
						multiLine=false,
					}
				end
				self.font:draw{
					pos = {xmin, ymax, gridz},
					text=varName,
					color = {1,1,1,1},
					fontSize = {fontSizeX, -fontSizeY},
					multiLine = false,
				}
			end
			
--			gl.glEnable(gl.GL_DEPTH_TEST)
		end
	end
	
--	gl.glDisable(gl.GL_DEPTH_TEST)
end

function HydroCLApp:display2D_Graph(solvers, varName, ar, graph_xmin, graph_ymin, graph_xmax, graph_ymax)
	self.view:projection(ar)
	self.view:modelview()
	gl.glColor3f(1,1,1)
	gl.glEnable(gl.GL_DEPTH_TEST)

	for _,solver in ipairs(solvers) do 
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			-- TODO allow a fixed, manual colormap range
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end

			-- TODO gui this somewhere
			local step = 1
			local ambient = .3
			
			-- TODO where to specify using the heatmap gradient vs using the variable/solver color
			gl.glColor3f(table.unpack((#self.solvers > 1 and solver or var).color))

			solver:calcDisplayVarToTex(var)
	
			self.graphShader:use()
			solver:getTex(var):bind()

			local scale = 1 / (valueMax - valueMin)
			local offset = valueMin
			gl.glUniform1f(self.graphShader.uniforms.scale.loc, scale)
			gl.glUniform1f(self.graphShader.uniforms.offset.loc, offset)
			
			gl.glUniform1f(self.graphShader.uniforms.ambient.loc, ambient)
			gl.glUniform1i(self.graphShader.uniforms.axis.loc, solver.dim)
			gl.glUniform1i(self.graphShader.uniforms.useLog.loc, var.useLog)
			gl.glUniform2f(self.graphShader.uniforms.size.loc, solver.gridSize.x, solver.gridSize.y)
			gl.glUniform2f(self.graphShader.uniforms.xmin.loc, solver.mins[1], solver.mins[2])
			gl.glUniform2f(self.graphShader.uniforms.xmax.loc, solver.maxs[1], solver.maxs[2])

			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

			for ybase=2,tonumber(solver.gridSize.y)-2-step,step do
				gl.glBegin(gl.GL_TRIANGLE_STRIP)
				for x=2,tonumber(solver.gridSize.x)-2,step do
					for yofs=0,step,step do
						local y = ybase + yofs
						gl.glVertex2d(
							(x + .5) / tonumber(solver.gridSize.x),
							(y + .5) / tonumber(solver.gridSize.y))
					end
				end
				gl.glEnd()
			end
			
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
			
			self.graphShader:useNone()
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
end

local display2DMethods = table{
	{Heatmap = HydroCLApp.display2D_Heatmap},
	{Graph = HydroCLApp.display2D_Graph},
}
local display2DMethodNames = display2DMethods:map(function(kv)
	return (next(kv))
end)
function HydroCLApp:display2D(...)
	self.display2DMethod = self.display2DMethod or ffi.new('int[1]', 0)
	select(2, next(display2DMethods[ self.display2DMethod[0]+1 ]))(self, ...)
end

-- 2D
local vertexesInQuad = {{0,0},{1,0},{1,1},{0,1}}

-- 3D
local vertexesInCube = {
	0,0,0,
	1,0,0,
	0,1,0,
	1,1,0,
	0,0,1,
	1,0,1,
	0,1,1,
	1,1,1,
}

local quadsInCube = {
	0,1,3,2,
	4,6,7,5,
	1,5,7,3,
	0,2,6,4,
	0,4,5,1,
	2,3,7,6,
}

--[[
looks great for flat space
TODO for curved space: provide a coordMapInv function (might have to be manual to account for domains of rotations)
 and then call this as we march through volumes 
 and treat out-of-bound values as fully transparent
--]]
HydroCLApp.display3D_Slice_usePoints = false
HydroCLApp.display3D_Slice_useIsos = true
HydroCLApp.display3D_Slice_numIsobars = 20
HydroCLApp.display3D_Slice_useLighting = false
HydroCLApp.display3D_Slice_alpha = .15
HydroCLApp.display3D_Slice_alphaGamma = 1
HydroCLApp.display3D_Slice_numSlices = 255
function HydroCLApp:display3D_Slice(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

if useClipPlanes then
	for i,clipInfo in ipairs(clipInfos) do
		gl.glClipPlane(gl.GL_CLIP_PLANE0+i-1, clipInfo.plane:ptr())
	end
end
	
	for _,solver in ipairs(solvers) do 
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end

			solver:calcDisplayVarToTex(var)	

			solver.volumeSliceShader:use()
			solver:getTex(var):bind(0)
			self.gradientTex:bind(1)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.alpha.loc, self.display3D_Slice_alpha)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.alphaGamma.loc, self.display3D_Slice_alphaGamma)
			gl.glUniform3f(solver.volumeSliceShader.uniforms.mins.loc, solver.mins:unpack())
			gl.glUniform3f(solver.volumeSliceShader.uniforms.maxs.loc, solver.maxs:unpack())
			gl.glUniform1i(solver.volumeSliceShader.uniforms.useLog.loc, var.useLog)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.valueMin.loc, valueMin)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.valueMax.loc, valueMax)
			gl.glUniform1i(solver.volumeSliceShader.uniforms.useIsos.loc, self.display3D_Slice_useIsos)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.numIsobars.loc, self.display3D_Slice_numIsobars)
			gl.glUniform1i(solver.volumeSliceShader.uniforms.useLighting.loc, self.display3D_Slice_useLighting)
			gl.glUniform1f(solver.volumeSliceShader.uniforms.numGhost.loc, solver.numGhost)
			gl.glUniform3f(solver.volumeSliceShader.uniforms.texSize.loc, solver.gridSize:unpack())

if useClipPlanes then
			for i,info in ipairs(clipInfos) do
				gl.glUniform1i(solver.volumeSliceShader.uniforms['clipEnabled'..i].loc, info.enabled and 1 or 0)
			end
end

			if self.display3D_Slice_usePoints then
				gl.glEnable(gl.GL_DEPTH_TEST)
				gl.glPointSize(2)
				gl.glBegin(gl.GL_POINTS)
				local numGhost = solver.numGhost
				for i=numGhost+1,tonumber(solver.gridSize.x-numGhost) do
					for j=numGhost+1,tonumber(solver.gridSize.y-numGhost) do
						for k=numGhost+1,tonumber(solver.gridSize.z-numGhost) do
							gl.glVertex3d(
								(i - numGhost - .5)/tonumber(solver.gridSize.x - 2*numGhost),
								(j - numGhost - .5)/tonumber(solver.gridSize.y - 2*numGhost),
								(k - numGhost - .5)/tonumber(solver.gridSize.z - 2*numGhost))
						end
					end
				end
				gl.glEnd()
				gl.glDisable(gl.GL_DEPTH_TEST)
			
			else
			
				gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
				gl.glEnable(gl.GL_BLEND)

				local n = self.display3D_Slice_numSlices
				local fwd = -self.frustumView.angle:conjugate():zAxis()
				local fwddir = select(2, table(fwd):map(math.abs):sup())

				local jmin, jmax, jdir
				if fwd[fwddir] < 0 then
					jmin, jmax, jdir = 0, n, 1
				else
					jmin, jmax, jdir = n, 0, -1
				end
				gl.glUniform3f(solver.volumeSliceShader.uniforms.normal.loc, 
					fwddir == 1 and jdir or 0, 
					fwddir == 2 and jdir or 0, 
					fwddir == 3 and jdir or 0)
						
				if CartesianGeom.is(solver.geometry) then
					-- [[	single quad
					gl.glBegin(gl.GL_QUADS)
					for j=jmin,jmax,jdir do
						local f = j/n
						for _,vtx in ipairs(vertexesInQuad) do
							if fwddir == 1 then
								gl.glVertex3f(f, vtx[1], vtx[2])
							elseif fwddir == 2 then
								gl.glVertex3f(vtx[1], f, vtx[2])
							elseif fwddir == 3 then
								gl.glVertex3f(vtx[1], vtx[2], f)
							end
						end
					end
					gl.glEnd()
					--]]
				else
					-- [[	use a grid, so curved coordinates can be seen
					for j=jmin,jmax,jdir do
						local f = j/n
						local xres = 20
						local yres = 20
						for ybase=1,yres-1 do
							gl.glBegin(gl.GL_TRIANGLE_STRIP)
							for x=1,xres do
								for y=ybase,ybase+1 do
									if fwddir == 1 then
										gl.glVertex3f(f, x/xres, y/yres)
									elseif fwddir == 2 then
										gl.glVertex3f(x/xres, f, y/yres)
									elseif fwddir == 3 then
										gl.glVertex3f(x/xres, y/yres, f)
									end
								end
							end
							gl.glEnd()
						end
					end
					--]]
				end
			
				gl.glDisable(gl.GL_BLEND)
			end

			self.gradientTex:unbind(1)
			solver:getTex(var):unbind(0)
			solver.volumeSliceShader:useNone()
		end
	end
end

function HydroCLApp:display3D_Ray(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

	for _,solver in ipairs(solvers) do
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end
			
			solver:calcDisplayVarToTex(var)	

			gl.glColor3f(1,1,1)
			for pass=0,1 do
				if pass == 0 then
					gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
				else
					gl.glEnable(gl.GL_CULL_FACE)
					gl.glCullFace(gl.GL_FRONT)
					gl.glEnable(gl.GL_DEPTH_TEST)
					gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
					gl.glEnable(gl.GL_BLEND)
					solver.volumeRayShader:use()
					gl.glUniform1f(solver.volumeRayShader.uniforms.scale.loc, 1)--scale)
					gl.glUniform1i(solver.volumeRayShader.uniforms.useLog.loc, 0)--useLog and 1 or 0)
					gl.glUniform1f(solver.volumeRayShader.uniforms.alpha.loc, 1)--alpha)
					solver:getTex(var):bind(0)
					self.gradientTex:bind(1)
				end
				gl.glBegin(gl.GL_QUADS)
				for i=1,24 do
					local x = vertexesInCube[quadsInCube[i] * 3 + 0 + 1]
					local y = vertexesInCube[quadsInCube[i] * 3 + 1 + 1]
					local z = vertexesInCube[quadsInCube[i] * 3 + 2 + 1]
					gl.glTexCoord3f(x, y, z)
					x = x * (solver.maxs[1] - solver.mins[1]) + solver.mins[1]
					y = y * (solver.maxs[2] - solver.mins[2]) + solver.mins[2]
					z = z * (solver.maxs[3] - solver.mins[3]) + solver.mins[3]
					gl.glVertex3f(x, y, z)
				end
				gl.glEnd()
				if pass == 0 then
					gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
				else
					self.gradientTex:unbind(1)
					solver:getTex(var):unbind(0)
					solver.volumeRayShader:useNone()
					gl.glDisable(gl.GL_BLEND)
					gl.glDisable(gl.GL_DEPTH_TEST)
					gl.glCullFace(gl.GL_BACK)
					gl.glDisable(gl.GL_CULL_FACE)
				end
			end
		end
	end
glreport'here'
end

--[[

6---12--7
|\      |\
| 10    8 11
7  \    |  \
|   4--9----5
|   |   |   |
2---|6--3   |
 \  3    \  5
  2 |     4 |
   \|      \|
    0---1---1


76543210
00010110

--]]

-- maps from the sum of two vertex bits to the edge index 
local edgeForTwoCorners = table.map({
	2^0 + 2^1,
	2^0 + 2^2,
	2^0 + 2^4,
	2^1 + 2^3,
	2^1 + 2^5,
	2^2 + 2^3,
	2^2 + 2^6,
	2^3 + 2^7,
	2^4 + 2^5,
	2^4 + 2^6,
	2^5 + 2^7,
	2^6 + 2^7,
}, function(v,k)
	return k,v
end)

function reverse(t)
	local nt = {}
	for i=#t,1,-1 do
		table.insert(nt, t[i])
	end
	return nt
end

local singles = {
	[2^0] = {1, 2, 3},
	[2^1] = reverse{1, 4, 5},
	[2^2] = {2, 6, 7},
	[2^3] = reverse{4, 6, 8},
	[2^4] = reverse{3, 9, 10},
	[2^5] = {5, 9, 11},
	[2^6] = reverse{7, 10, 12},
	[2^7] = {8, 11, 12},
}
-- [[ combine all singles that do not share an edge in common
local function neighbor(corner, corner2)
	return corner == bit.bxor(corner2, 1)
		or corner == bit.bxor(corner2, 2)
		or corner == bit.bxor(corner2, 4)
end
for corner=0,6 do
	for corner2=corner+1,7 do
		if not neighbor(corner, corner2) then
			singles[2^corner + 2^corner2] = table():append(singles[2^corner], singles[2^corner2])
		end
		for corner3=corner2+1,7 do
			if not neighbor(corner, corner2) 
			and not neighbor(corner2, corner3)
			and not neighbor(corner, corner3)
			then
				singles[2^corner + 2^corner2 + 2^corner3] = table():append(singles[2^corner], singles[2^corner2], singles[2^corner3])
			end
			for corner4=corner3+1,7 do
				if not neighbor(corner, corner2) 
				and not neighbor(corner, corner3)
				and not neighbor(corner, corner4)
				and not neighbor(corner2, corner3)
				and not neighbor(corner2, corner4)
				and not neighbor(corner3, corner4)
				then
					singles[2^corner + 2^corner2 + 2^corner3 + 2^corner4] = table():append(singles[2^corner], singles[2^corner2], singles[2^corner3], singles[2^corner4])
				end		
			end
		end
	end
end
--]]
edgesForInside = table(single)

local doubles = {
	[2^0 + 2^1] = {2, 3, 5, 5, 4, 2},
	[2^2 + 2^3] = {2, 4, 8, 8, 7, 2},
	[2^4 + 2^5] = reverse{3, 5, 11, 11, 10, 3},
	[2^6 + 2^7] = {7, 8, 11, 11, 10, 7},
	
	[2^0 + 2^2] = reverse{1, 3, 7, 7, 6, 1},
	[2^1 + 2^3] = {1, 5, 8, 8, 6, 1},
	[2^4 + 2^6] = {3, 7, 12, 12, 9, 3},
	[2^5 + 2^7] = reverse{5, 8, 12, 12, 9, 5},
	
	[2^0 + 2^4] = {1, 2, 10, 10, 9, 1},
	[2^1 + 2^5] = reverse{1, 4, 11, 11, 9, 1},
	[2^2 + 2^6] = {2, 6, 12, 12, 10, 2},
	[2^3 + 2^7] = reverse{4, 6, 12, 12, 11, 4},
}
-- [[
for _,pair in ipairs{
	{ {0,1}, {6,7} },
	{ {4,5}, {2,3} },
	
	{ {0,2}, {5,7} },
	{ {1,3}, {4,6} },
	
	{ {0,4}, {3,7} },
	{ {1,5}, {2,6} },
} do
	local pa, pb = table.unpack(pair)
	local ka = 2^pa[1] + 2^pa[2]
	local kb = 2^pb[1] + 2^pb[2]
	doubles[ka + kb] = table():append(doubles[ka], doubles[kb])
end
--]]
edgesForInside = table(edgesForInside, doubles)

local triples = {
--[[ something in here is off
	--bad
	[2^0 + 2^1 + 2^2] = {3, 5, 7, 4, 6, 7, 7, 5, 4},	
	[2^0 + 2^1 + 2^3] = {3, 5, 9, 3, 9, 6, 6, 2, 3},
	[2^0 + 2^2 + 2^3] = {3, 8, 7, 1, 4, 8, 8, 3, 1},
	[2^1 + 2^2 + 2^3] = {5, 8, 7, 1, 5, 7, 7, 2, 1},

	--bad
	[2^4 + 2^5 + 2^6] = {3, 7, 5, 5, 7, 12, 12, 11, 5},
	[2^4 + 2^5 + 2^7] = {3, 9, 5, 3, 10, 12, 12, 8, 3},
	[2^4 + 2^6 + 2^7] = {3, 7, 8, 3, 8, 11, 11, 9, 3},
	[2^5 + 2^6 + 2^7] = {5, 7, 8, 5, 9, 10, 10, 7, 5},
--]]
	
	--good
	[2^0 + 2^2 + 2^6] = {1, 6, 12, 1, 12, 10, 10, 3, 1},
	[2^0 + 2^2 + 2^4] = {1, 6, 9, 6, 7, 10, 10, 9, 6},
	[2^0 + 2^4 + 2^6] = {1, 12, 9, 1, 2, 7, 7, 12, 1},
	[2^2 + 2^4 + 2^6] = {6, 12, 9, 2, 6, 9, 9, 3, 2},

	--good
	[2^1 + 2^3 + 2^5] = {1, 9, 6, 9, 11, 8, 8, 6, 9},
	[2^1 + 2^3 + 2^7] = {1, 12, 6, 1, 5, 11, 11, 12, 1},
	[2^1 + 2^5 + 2^7] = {1, 9, 12, 1, 12, 8, 8, 4, 1},
	[2^3 + 2^5 + 2^7] = {6, 9, 12, 4, 5, 9, 9, 6, 4},

	--good
	[2^0 + 2^1 + 2^4] = {2, 10, 4, 4, 10, 9, 9, 5, 4},
	[2^0 + 2^1 + 2^5] = {2, 11, 4, 2, 3, 9, 9, 11, 2},
	[2^0 + 2^4 + 2^5] = {2, 10, 11, 1, 2, 11, 11, 5, 1},
	[2^1 + 2^4 + 2^5] = {4, 10, 11, 1, 3, 10, 10, 4, 1},

	--good
	[2^2 + 2^3 + 2^6] = {2, 4, 10, 4, 8, 12, 12, 10, 4},
	[2^2 + 2^3 + 2^7] = {2, 4, 11, 2, 11, 12, 12, 7, 2},
	[2^2 + 2^6 + 2^7] = {2, 11, 10, 2, 6, 8, 8, 11, 2},
	[2^3 + 2^6 + 2^7] = {4, 11, 10, 4, 10, 7, 7, 6, 4},
}
edgesForInside = table(edgesForInside, triples)

local quads = {
	[2^0 + 2^1 + 2^2 + 2^3] = {3, 5, 8, 8, 7, 3},
	[2^0 + 2^2 + 2^4 + 2^6] = {1, 6, 12, 12, 9, 1},
	[2^0 + 2^1 + 2^4 + 2^5] = reverse{2, 4, 11, 11, 10, 2},
}
edgesForInside = table(edgesForInside, quads)

for _,k in ipairs(table.keys(edgesForInside)) do
	edgesForInside[bit.band(0xff, bit.bnot(k))] = reverse(edgesForInside[k])
end

--[[
for i=0,255 do
	if not edgesForInside[i] then print(('0x%x'):format(i)) end
end
--]]

function HydroCLApp:display3D_Isosurface(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

	-- draw wireframe
	for _,solver in ipairs(solvers) do
		local volumeRayShader = solver.volumeRayShader
		gl.glColor3f(1,1,1)
		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
		gl.glBegin(gl.GL_QUADS)
		for i=1,24 do
			local x = vertexesInCube[quadsInCube[i] * 3 + 0 + 1]
			local y = vertexesInCube[quadsInCube[i] * 3 + 1 + 1]
			local z = vertexesInCube[quadsInCube[i] * 3 + 2 + 1]
			gl.glTexCoord3f(x, y, z)
			x = x * (solver.maxs[1] - solver.mins[1]) + solver.mins[1]
			y = y * (solver.maxs[2] - solver.mins[2]) + solver.mins[2]
			z = z * (solver.maxs[3] - solver.mins[3]) + solver.mins[3]
			gl.glVertex3f(x, y, z)
		end
		gl.glEnd()
		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
	end
			
	gl.glColor3f(1,1,1)
	gl.glEnable(gl.GL_CULL_FACE)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	gl.glEnable(gl.GL_BLEND)
	
	self.isobarShader:use()
	gl.glBegin(gl.GL_TRIANGLES)

	for _,solver in ipairs(solvers) do 
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end
			
			solver:calcDisplayVarToTex(var)	
			
			assert(not self.useGLSharing, "I still need to code in the GL sharing version")
			local dest = ffi.cast('float*', solver.calcDisplayVarToTexPtr)
			local cornerValues = {}
			local cornerBars = {}
			local edgeVtxs = {}

			local numBars = 1

			for k=solver.numGhost,tonumber(solver.gridSize.z)-solver.numGhost-1 do
				for j=solver.numGhost,tonumber(solver.gridSize.y)-solver.numGhost-1 do
					for i=solver.numGhost,tonumber(solver.gridSize.x)-solver.numGhost-1 do
						local ofs = i + solver.gridSize.x * (j + solver.gridSize.y * k)
					
						local cellMinBar = math.huge
						local cellMaxBar = -math.huge

						for corner=0,7 do
							local cornerOfs = ofs
							for n=0,2 do
								if bit.band(bit.rshift(corner, n), 1) == 1 then
									cornerOfs = cornerOfs + solver.stepSize:ptr()[n]
								end
							end
							
							local cornerValue = dest[cornerOfs]
							cornerValues[corner] = cornerValue
							
							local cornerBar = math.floor( (cornerValue - valueMin) / (valueMax - valueMin) * numBars - .5 )
							cornerBars[corner] = cornerBar
							cellMinBar = math.min(cellMinBar, cornerBar)
							cellMaxBar = math.max(cellMaxBar, cornerBar)
						end

						if not (cellMaxBar < 0 or cellMinBar > numBars) then 
							for bar = math.max(cellMinBar, 0), math.min(cellMaxBar, numBars-1) do
								local isoValue = (bar + .5) / numBars * (valueMax - valueMin) + valueMin
						
								local inside = 0
								for corner=0,7 do
									if cornerValues[corner] > isoValue then
										inside = bit.bor(inside, bit.lshift(1, corner))
									end
								end

								-- now for each 0/1 pair along each edge, calculate the fractions based on the fractions of values across corners
								for corner=0,6 do
									for n=0,2 do
										local sign = bit.band(bit.rshift(corner, n), 1) == 0 and 1 or -1
										local nextCorner = bit.bxor(corner, bit.lshift(1, n))
										if nextCorner > corner then
											--v.x, v.y, v.z = i, j, k
											local v = vec3d(i,j,k)
											for m=0,2 do
												if bit.band(corner, bit.lshift(1,m)) ~= 0 then
													v:ptr()[m] = v:ptr()[m] + 1
												end
											end
											
											local frac = (isoValue - cornerValues[corner]) / (cornerValues[nextCorner] - cornerValues[corner])
											v:ptr()[n] = v:ptr()[n] + sign * frac
										
											local twoCornerIndex = bit.bor(bit.lshift(1, corner), bit.lshift(1, nextCorner))
											local edge = edgeForTwoCorners[twoCornerIndex]
											edgeVtxs[edge] = v 
										end
									end
								end
								
								local edges = edgesForInside[inside]
								if edges then
									for _,edge in ipairs(edges) do
										local x,y,z = edgeVtxs[edge]:unpack()
										
										x = (x - solver.numGhost) / tonumber(solver.gridSize.x - 2 * solver.numGhost) * (solver.maxs[1] - solver.mins[1]) + solver.mins[1]
										y = (y - solver.numGhost) / tonumber(solver.gridSize.y - 2 * solver.numGhost) * (solver.maxs[2] - solver.mins[2]) + solver.mins[2]
										z = (z - solver.numGhost) / tonumber(solver.gridSize.z - 2 * solver.numGhost) * (solver.maxs[3] - solver.mins[3]) + solver.mins[3]
									
										gl.glColor3f( table.unpack(({
											{1,0,0},
											{0,1,0},
											[0] = {0,0,1},
										})[_%3]) )
										gl.glVertex3d(x,y,z)
									end
								end
							end
						end
					end
				end
			end	
		end
	end
	
	gl.glEnd()
	self.isobarShader:useNone()
	gl.glDisable(gl.GL_DEPTH_TEST)
	gl.glDisable(gl.GL_CULL_FACE)
	gl.glDisable(gl.GL_BLEND)
end

local display3DMethods = table{
	{Slices = HydroCLApp.display3D_Slice},
	{Raytrace = HydroCLApp.display3D_Ray},
	{Isosurfaces = HydroCLApp.display3D_Isosurface},
}
local display3DMethodNames =  display3DMethods:map(function(kv)
	return (next(kv))
end)
function HydroCLApp:display3D(...)
	self.display3DMethod = self.display3DMethod or ffi.new('int[1]', 0)
	select(2, next(display3DMethods[ self.display3DMethod[0]+1 ]))(self, ...)
end

local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
}
HydroCLApp.displayVectorField_scale = .1
HydroCLApp.displayVectorField_step = 4
function HydroCLApp:displayVectorField(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)

	for _,solver in ipairs(solvers) do
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex and var.enabled then
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end
			
			solver:calcDisplayVarToTex(var)	
					
			solver.vectorFieldShader:use()
			solver:getTex(var):bind(0)
			
			gl.glUniform3f(solver.vectorFieldShader.uniforms.mins.loc, solver.mins:unpack())
			gl.glUniform3f(solver.vectorFieldShader.uniforms.maxs.loc, solver.maxs:unpack())
			-- how to determine scale?
			local scale = self.displayVectorField_scale * (valueMax - valueMin)
			gl.glUniform1f(solver.vectorFieldShader.uniforms.scale.loc, scale) 

			local step = self.displayVectorField_step
			gl.glBegin(gl.GL_LINES)
			for k=0,tonumber(solver.sizeWithoutBorder.z-1),step do
				for j=0,tonumber(solver.sizeWithoutBorder.y-1),step do
					for i=0,tonumber(solver.sizeWithoutBorder.x-1),step do
						local tx = (i + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
						local ty = (j + .5 + (solver.dim > 1 and solver.numGhost or 0)) / tonumber(solver.gridSize.y)
						local tz = (k + .5 + (solver.dim > 2 and solver.numGhost or 0)) / tonumber(solver.gridSize.z)
						gl.glMultiTexCoord3f(gl.GL_TEXTURE0, tx, ty, tz)	
						local x = (i + .5) / tonumber(solver.sizeWithoutBorder.x)
						local y = (j + .5) / tonumber(solver.sizeWithoutBorder.y)
						local z = (k + .5) / tonumber(solver.sizeWithoutBorder.z)
						gl.glMultiTexCoord3f(gl.GL_TEXTURE1, x, y, z)
						for _,q in ipairs(arrow) do
							gl.glVertex2f(q[1], q[2])
						end
					end
				end
			end
			gl.glEnd()
	
			solver:getTex(var):unbind(0)
			solver.vectorFieldShader:useNone()
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
	
glreport'here'
end

function HydroCLApp:updateGUI()
	if ig.igCollapsingHeader'simulation' then
		if ig.igButton(self.running and 'Stop' or 'Start') then
			self.running = not self.running
		end
		ig.igSameLine()
		if ig.igButton'Step' then
			self.running = 'step'
		end
		ig.igSameLine()
		if ig.igButton'Reset' then
			print'resetting...'
			for _,solver in ipairs(self.solvers) do
				solver:resetState()
			end
			self.running = nil
		end

		if ig.igButton'Save' then
			-- save as cfits 
			for i,solver in ipairs(self.solvers) do
				solver:save(tostring(i))
			end
		end

		tooltip.checkboxTable('dump to file', dumpFile, 'enabled')
		ig.igSameLine()

		tooltip.checkboxTable('stack graphs', self, 'displayAllTogether')

		if ig.igRadioButtonBool('ortho', self.view == self.orthoView) then
			self.view = self.orthoView
		end
		ig.igSameLine()
		if ig.igRadioButtonBool('frustum', self.view == self.frustumView) then
			self.view = self.frustumView
		end
	
		-- TODO flag for separate/combined displays (esp for ortho view)

		-- TODO flag to toggle slice vs volume display
		-- or maybe checkboxes for each kind?
		
		if self.solvers[1] then
			local dim = self.solvers[1].dim
			if dim == 2 then
				ig.igPushIdStr'2D'
				ig.igCombo('Display Method', self.display2DMethod, display2DMethodNames)
				ig.igPopId()
			elseif dim == 3 then
				ig.igPushIdStr'3D'
				ig.igCombo('Display Method', self.display3DMethod, display3DMethodNames)
				
				-- if we're doing 3D slice display 
				if HydroCLApp.display3D_Slice == select(2, next(display3DMethods[self.display3DMethod[0]+1])) then

if useClipPlanes then
					ig.igRadioButton("rotate camera", rotateClip, 0)
					for i,clipInfo in ipairs(clipInfos) do
						ig.igPushIdStr('clip '..i)
						tooltip.checkbox('clip', clipInfo.enabled)
						ig.igSameLine()
						ig.igRadioButton('rotate', rotateClip, i)
						ig.igSameLine()
						if ig.igButton('reset') then
							clipInfo.plane = makeDefaultPlane(i)
						end
						ig.igPopId()
					end				
end					
					tooltip.sliderTable('alpha', self, 'display3D_Slice_alpha', 0, 1)
					tooltip.sliderTable('gamma', self, 'display3D_Slice_alphaGamma', 0, 1)
					tooltip.checkboxTable('isobars', self, 'display3D_Slice_useIsos')
					if self.display3D_Slice_useIsos then
						tooltip.intTable('num isobars', self, 'display3D_Slice_numIsobars')
					end
					tooltip.checkboxTable('lighting', self, 'display3D_Slice_useLighting')
					tooltip.checkboxTable('pointcloud', self, 'display3D_Slice_usePoints')
					if not self.display3D_Slice_usePoints then
						tooltip.intTable('num slices', self, 'display3D_Slice_numSlices')
					end
				end
				ig.igPopId()
			end
		
			--ig.igCheckbox('vector field', self.enableVectorField)
		
			tooltip.numberTable('vector field scale', self, 'displayVectorField_scale')
			--tooltip.sliderTable('vector field scale', self, 'displayVectorField_scale', 0, 100, nil, 10)
			
			tooltip.intTable('vector field step', self, 'displayVectorField_step')
		end
	end
	
	for i,solver in ipairs(self.solvers) do
		ig.igPushIdStr('solver '..i)
		if ig.igCollapsingHeader(solver.name) then
			-- TODO new window for each
			solver:updateGUI()
		end
		ig.igPopId()
	end
end

local leftButtonDown
local rightButtonDown
local leftShiftDown
local rightShiftDown
local leftGuiDown
local rightGuiDown
function HydroCLApp:event(event, ...)
	HydroCLApp.super.event(self, event, ...)
	local canHandleMouse = not ig.igGetIO()[0].WantCaptureMouse
	local canHandleKeyboard = not ig.igGetIO()[0].WantCaptureKeyboard
	local shiftDown = leftShiftDown or rightShiftDown
	local guiDown = leftGuiDown or rightGuiDown
	if event.type == sdl.SDL_MOUSEMOTION then
		if canHandleMouse then
			local dx = event.motion.xrel
			local dy = event.motion.yrel
			if leftButtonDown and not guiDown then
				if shiftDown then
					if dx ~= 0 or dy ~= 0 then
						self.view:mouseZoom(-dy, dy)
					end
				else
					if dx ~= 0 or dy ~= 0 then
						self.view:mousePan(dx, dy, self:size())
					end
				end
			end
		end
	elseif event.type == sdl.SDL_MOUSEBUTTONDOWN then
		if event.button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = true
		elseif event.button.button == sdl.SDL_BUTTON_RIGHT then
			rightButtonDown = true
		end
	elseif event.type == sdl.SDL_MOUSEBUTTONUP then
		if event.button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = false
		elseif event.button.button == sdl.SDL_BUTTON_RIGHT then
			rightButtonDown = false
		end
	elseif event.type == sdl.SDL_KEYDOWN then
		if event.key.keysym.sym == sdl.SDLK_LSHIFT then
			leftShiftDown = true
		elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
			rightShiftDown = true
		elseif event.key.keysym.sym == sdl.SDLK_LGUI then
			leftGuiDown = true
		elseif event.key.keysym.sym == sdl.SDLK_RGUI then
			rightGuiDown = true
		end
	elseif event.type == sdl.SDL_KEYUP then
		if event.key.keysym.sym == sdl.SDLK_LSHIFT then
			leftShiftDown = false
		elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
			rightShiftDown = false
		elseif event.key.keysym.sym == sdl.SDLK_LGUI then
			leftGuiDown = false
		elseif event.key.keysym.sym == sdl.SDLK_RGUI then
			rightGuiDown = false
		elseif canHandleKeyboard then
			if event.key.keysym.sym == sdl.SDLK_SPACE then
				self.running = not self.running
			elseif event.key.keysym.sym == ('u'):byte() then
				self.running = 'step'
			elseif event.key.keysym.sym == ('r'):byte() then
				print'resetting...'
				for _,solver in ipairs(self.solvers) do
					solver:resetState()
				end
				self.running = nil
			end
		end
	end
end

return HydroCLApp
