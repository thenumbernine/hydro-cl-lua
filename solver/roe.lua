local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'ffi.OpenGL'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local vec3sz = require 'solver.vec3sz'
local vec3 = require 'vec.vec3'
local clnumber = require 'clnumber'
local showcode = require 'showcode'
local processcl = require 'processcl'
local CLImageGL = require 'cl.imagegl'
local CLProgram = require 'cl.program'
local GLTex2D = require 'gl.tex2d'
local GLTex3D = require 'gl.tex2d'
local glreport = require 'gl.report'


local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}

local function xs_to_rs(code)
	return (code:gsub('{x(%d)}', function(i)
		return 'r.'..xs[i+0]
	end))
end

local function getCode_real3_to_real3(name, exprs)
	return 'inline real3 '..name..'(real3 r) { return _real3('
		..range(3):map(function(i)
			return xs_to_rs(exprs[i])
		end):concat', '..'); }'
end

local function getCode_define_i3_to_real3(name, exprs)
	return '#define '..name..'(i) _real3('
		..exprs:map(function(expr,i)
			for j=1,3 do
				expr = expr:gsub('{x'..j..'}', 'cell_x'..(j-1)..'(i.'..xs[j]..')')
			end
			return expr
		end):concat', '..')'
end

local Solver = class()

Solver.name = 'Roe'
Solver.numGhost = 2

-- enable these to verify accuracy
-- disable these to save on allocation / speed
Solver.checkFluxError = true 
Solver.checkOrthoError = true 

Solver.integrators = require 'int.all'
Solver.integratorNames = Solver.integrators:map(function(integrator) return integrator.name end)

--[[
args:
	app
	gridSize
	dim (optional, required if gridSize is a vec3sz)
--]]
function Solver:init(args)
	assert(args)
	local gridSize = assert(args.gridSize)
	if type(gridSize) == 'number' then 
		self.gridSize = vec3sz(gridSize,1,1)
		assert(not args.dim or args.dim == 1)
		self.dim = 1
	elseif type(gridSize) == 'cdata' 
	and ffi.typeof(gridSize) == vec3sz 
	then
		self.gridSize = vec3sz(gridSize)
		self.dim = assert(args.dim)
	elseif type(gridSize) == 'table' then
		self.gridSize = vec3sz(table.unpack(gridSize))
		assert(not args.dim or args.dim <= #gridSize)
		self.dim = args.dim or #gridSize
	else
		error("can't understand args.gridSize type "..type(args.gridSize).." value "..tostring(args.gridSize))
	end

	for i=self.dim,2 do self.gridSize:ptr()[i] = 1 end

	self.color = vec3(math.random(), math.random(), math.random()):normalize()

	self.mins = vec3(table.unpack(args.mins or {-1, -1, -1}))
	self.maxs = vec3(table.unpack(args.maxs or {1, 1, 1}))

	self.t = 0

	self:createEqn(args.eqn)
	
	self.name = self.eqn.name..' '..self.name

	self.app = assert(args.app)
	
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	self.maxWorkGroupSize = tonumber(self.app.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')

	self.offset = vec3sz(0,0,0)
	self.localSize1d = math.min(self.maxWorkGroupSize, tonumber(self.gridSize:volume()))
	
--	self.localSize = self.dim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
	-- TODO better than constraining by math.min(gridSize),
	-- look at which gridSizes have the most room, and double them accordingly, until all of maxWorkGroupSize is taken up
	self.localSize = vec3sz(1,1,1)
	local rest = self.maxWorkGroupSize
	local localSizeX = math.min(gridSize[1], 2^math.ceil(math.log(rest^(1/self.dim),2)))
	self.localSize.x = localSizeX
	if self.dim > 1 then
		rest = rest / localSizeX
		if self.dim == 2 then
			self.localSize.y = math.min(gridSize[2], rest)
		elseif self.dim == 3 then
			local localSizeY = math.min(gridSize[2], 2^math.ceil(math.log(math.sqrt(rest),2)))
			self.localSize.y = localSizeY
			self.localSize.z = math.min(gridSize[3], rest / localSizeY)
		end
	end

	self.useFixedDT = ffi.new('bool[1]', false)
	self.fixedDT = ffi.new('float[1]', 0)
	self.cfl = ffi.new('float[1]', .5)
	
	self.initStatePtr = ffi.new('int[1]', (table.find(self.eqn.initStateNames, args.initState) or 1)-1)
	
	self.integratorPtr = ffi.new('int[1]', (self.integratorNames:find(args.integrator) or 1)-1)
	
	self.slopeLimiterPtr = ffi.new('int[1]', (self.app.slopeLimiterNames:find(args.slopeLimiter) or 1)-1)

	self.boundaryMethods = {}
	for i=1,3 do
		for _,minmax in ipairs(minmaxs) do
			local var = xs[i]..minmax
			self.boundaryMethods[var] = ffi.new('int[1]', self.app.boundaryMethods:find(
				(args.boundary or {})[var] or 'freeflow'
			)-1)
		end
	end

	self.geometry = require('geom.'..args.geometry){solver=self}

	self:refreshGridSize()
end

function Solver:getCoordMapCode()
	return table{
		getCode_real3_to_real3('coordMap', range(3):map(function(i)
			return self.geometry.uCode[i] or '{x'..i..'}'
		end)),
	}:concat'\n'
end

function Solver:getCoordMapGLSLCode()
	return (self:getCoordMapCode()
		:gsub('inline%S*', '')
		:gsub('_real3', 'vec3')
		:gsub('real3', 'vec3')
	)
end

-- this is the general function - which just assigns the eqn provided by the arg
-- but it can be overridden for specific equations
function Solver:createEqn(eqn)
	self.eqn = require('eqn.'..assert(eqn))(self)
end

function Solver:refreshGridSize()

	self.volume = tonumber(self.gridSize:volume())
	self.dxs = vec3(range(3):map(function(i)
		return (self.maxs[i] - self.mins[i]) / tonumber(self.gridSize:ptr()[i-1])	
	end):unpack())

	self:refreshIntegrator()	-- depends on eqn & gridSize

	self:createDisplayVars()	-- depends on eqn
	
	-- depends on eqn & gridSize
	self.buffers = table()
	self:createBuffers()
	self:finalizeCLAllocs()
	
	self:createCodePrefix()		-- depends on eqn, gridSize, displayVars

	self:refreshInitStateProgram()
	self:refreshCommonProgram()
	self:refreshSolverProgram()
	self:refreshDisplayProgram()
	self:refreshBoundaryProgram()

	self:resetState()
end

function Solver:refreshIntegrator()
	self.integrator = self.integrators[self.integratorPtr[0]+1](self)
end

local ConvertToTex = class()
Solver.ConvertToTex = ConvertToTex

ConvertToTex.type = 'real'	-- default

ConvertToTex.displayCode = [[
__kernel void <?=name?>(
	<?=input?>,
	int displayVar,
	const __global <?=type?>* buf
) {
	SETBOUNDS(0,0);
	int dstindex = index;
	int4 dsti = i;
	
	//now constrain
	if (i.x < 2) i.x = 2;
	if (i.x > gridSize_x - 2) i.x = gridSize_x - 2;
#if dim >= 2
	if (i.y < 2) i.y = 2;
	if (i.y > gridSize_y - 2) i.y = gridSize_y - 2;
#endif
#if dim >= 3
	if (i.z < 2) i.z = 2;
	if (i.z > gridSize_z - 2) i.z = gridSize_z - 2;
#endif
	//and recalculate read index
	index = INDEXV(i);
	
	int side = 0;
	int intindex = side + dim * index;
	real value = 0;
<?=body?>
<?=output?>
}
]]

function ConvertToTex:init(args)
	local solver = assert(args.solver)	
	self.name = assert(args.name)
	self.solver = solver
	self.type = args.type	-- or self.type
	self.displayCode = args.displayCode	-- or self.displayCode
	self.displayBodyCode = args.displayBodyCode
	self.vars = table()
	for i,name in ipairs(args.vars) do
		self.vars:insert{
			convertToTex = self,
			name = self.name..'_'..name,
			enabled = ffi.new('bool[1]', 
				self.name == 'U' and (solver.dim==1 or i==1)
				or (self.name == 'error' and solver.dim==1)
			),
			useLogPtr = ffi.new('bool[1]', args.useLog or false),
			color = vec3(math.random(), math.random(), math.random()):normalize(),
			--heatMapTexPtr = ffi.new('int[1]', 0),	-- hsv, isobar, etc ...
			heatMapFixedRangePtr = ffi.new('bool[1]', self.name ~= 'error'),
			heatMapValueMinPtr = ffi.new('float[1]', 0),
			heatMapValueMaxPtr = ffi.new('float[1]', 1),
		}
	end
end

function ConvertToTex:getCode()
	return self.displayCode 
end

function ConvertToTex:setArgs(kernel, var)
	kernel:setArg(1, ffi.new('int[1]', var.globalIndex))
	kernel:setArg(2, self.solver[var.convertToTex.name..'Buf'])
end

function ConvertToTex:setToTexArgs(var)
	self:setArgs(self.calcDisplayVarToTexKernel, var)
end

function ConvertToTex:setToBufferArgs(var)
	self:setArgs(self.calcDisplayVarToBufferKernel, var)
end

function Solver:addConvertToTex(args, class)
	class = class or ConvertToTex
	self.convertToTexs:insert(class(
		table(args, {
			solver = self,
		})
	))
end

function Solver:createDisplayVars()
	self.convertToTexs = table()
	self:addConvertToTexs()
	self.displayVars = table()
	for _,convertToTex in ipairs(self.convertToTexs) do
		self.displayVars:append(convertToTex.vars)
	end
	-- set the index in the master list of all display vars
	for i,var in ipairs(self.displayVars) do
		var.globalIndex = i
	end
end

function Solver:addConvertToTexUBuf()
	self:addConvertToTex{
		name = 'U',
		type = 'cons_t',
		displayBodyCode = [[
	const __global cons_t* U = buf + index;
]]..self.eqn:getCalcDisplayVarCode(),
		vars = assert(self.eqn.displayVars),
	}
end

function Solver:addConvertToTexs()
	self:addConvertToTexUBuf()
	self:addConvertToTex{
		name = 'wave',
		vars = range(0, self.eqn.numWaves-1),
		displayBodyCode = [[
	const __global real* wave = buf + intindex * numWaves;
	value = wave[displayVar - displayFirst_wave];
]],
	}
	
	local eigenDisplayVars = self.eqn:getEigenInfo().displayVars
	if eigenDisplayVars and #eigenDisplayVars > 0 then
		self:addConvertToTex{
			name = 'eigen',
			type = 'eigen_t',
			vars = eigenDisplayVars,
			displayBodyCode = [[
	const __global eigen_t* eigen = buf + intindex;
]]..self.eqn:getCalcDisplayVarEigenCode(),
		}
	end

	self:addConvertToTex{
		name = 'deltaUTilde', 
		vars = range(0,self.eqn.numWaves-1),
		displayBodyCode = [[
	const __global real* deltaUTilde = buf + intindex * numWaves;
	value = deltaUTilde[displayVar - displayFirst_deltaUTilde];
]],
	}
	self:addConvertToTex{
		name = 'rTilde',
		vars = range(0,self.eqn.numWaves-1),
		displayBodyCode = [[
	const __global real* rTilde = buf + intindex * numWaves;
	value = rTilde[displayVar - displayFirst_rTilde];
]],
	}
	self:addConvertToTex{
		name = 'flux', 
		vars = range(0,self.eqn.numStates-1),
		displayBodyCode = [[
	const __global real* flux = buf + intindex * numStates;
	value = flux[displayVar - displayFirst_flux];
]],
	}
	-- might contain nonsense :-p
	self:addConvertToTex{
		name = 'reduce', 
		vars = {'0'},
		displayBodyCode = [[
	value = buf[index];
]],
	}
	if self.checkFluxError or self.checkOrthoError then	
		self:addConvertToTex{
			name = 'error', 
			vars = {'ortho', 'flux'},
			type = 'error_t',
			useLog = true,
			displayBodyCode = [[
	if (displayVar == display_error_ortho) {
		value = buf[intindex].ortho;
	} else if (displayVar == display_error_flux) {
		value = buf[intindex].flux;
	}
]],
		}
	end
end

local errorType = 'error_t'
local errorTypeCode = 'typedef struct { real ortho, flux; } '..errorType..';'

-- my best idea to work around the stupid 8-arg max kernel restriction
Solver.allocateOneBigStructure = false

function Solver:clalloc(name, size)
	self.buffers:insert{name=name, size=size}	
end

function Solver:finalizeCLAllocs()
	local total = 0
	for _,buffer in ipairs(self.buffers) do
		buffer.offset = total
		local name = buffer.name
		local size = buffer.size
		total = total + size
		print('allocating '..name..' size '..size..' total '..total)
		if not self.allocateOneBigStructure then
			self[name] = self.app.ctx:buffer{rw=true, size=size}
		end
	end
	if self.allocateOneBigStructure then
		self.oneBigBuf = self.app.ctx:buffer{rw=true, size=total}
	end
end

function Solver:createBuffers()
	local realSize = ffi.sizeof(self.app.real)

	-- to get sizeof
	ffi.cdef(self.eqn:getEigenInfo().typeCode)
	ffi.cdef(errorTypeCode)

	-- should I put these all in one AoS?
	self:clalloc('UBuf', self.volume * self.eqn.numStates * realSize)
	self:clalloc('waveBuf', self.volume * self.dim * self.eqn.numWaves * realSize)
	self:clalloc('eigenBuf', self.volume * self.dim * ffi.sizeof'eigen_t')
	self:clalloc('deltaUTildeBuf', self.volume * self.dim * self.eqn.numWaves * realSize)
	self:clalloc('rTildeBuf', self.volume * self.dim * self.eqn.numWaves * realSize)
	self:clalloc('fluxBuf', self.volume * self.dim * self.eqn.numStates * realSize)
	
	-- debug only
	if self.checkFluxError then
		self:clalloc('fluxXformBuf', self.volume * self.dim * ffi.sizeof'fluxXform_t')
	end	
	if self.checkFluxError or self.checkOrthoError then
		local errorTypeSize = ffi.sizeof(errorType)
		self:clalloc('errorBuf', self.volume * self.dim * errorTypeSize)
	end
	
	self:clalloc('reduceBuf', self.volume * realSize)
	self:clalloc('reduceSwapBuf', self.volume * realSize / self.localSize1d)
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- CL/GL interop

	self.tex = (self.dim < 3 and GLTex2D or GLTex3D){
		width = self.gridSize.x,
		height = self.gridSize.y,
		depth = self.gridSize.z,
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
	}

	if self.app.useGLSharing then
		self.texCLMem = CLImageGL{context=self.app.ctx, tex=self.tex, write=true}
	else
		self.calcDisplayVarToTexPtr = ffi.new(self.app.real..'[?]', self.volume)
		
		--[[ PBOs?
		self.calcDisplayVarToTexPBO = ffi.new('gl_int[1]', 0)
		gl.glGenBuffers(1, self.calcDisplayVarToTexPBO)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, self.calcDisplayVarToTexPBO[0])
		gl.glBufferData(gl.GL_PIXEL_UNPACK_BUFFER, self.tex.width * self.tex.height * ffi.sizeof(self.app.real) * 4, nil, gl.GL_STREAM_READ)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
		--]]
	end
end

function Solver:createCodePrefix()
	local lines = table()
	if self.dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end
	if self.app.real == 'double' then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_fp64 : enable'
	end

	-- real types in CL natives: 1,2,4,8
	lines:append(table{'',2,4,8}:map(function(n)
		return 'typedef '..self.app.real..n..' real'..n..';'
	end))

	lines:insert(self.app.real3TypeCode)
	lines:insert[[
#define _real3(a,b,c) (real3){.s={a,b,c}}

static inline real real3_dot(real3 a, real3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline real3 real3_scale(real3 a, real s) {
	return _real3(a.x * s, a.y * s, a.z * s);
}

static inline real3 real3_add(real3 a, real3 b) {
	return _real3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline real3 real3_sub(real3 a, real3 b) {
	return _real3(a.x - b.x, a.y - b.y, a.z - b.z);
}

]]

	lines:append{
		'#define dim '..self.dim,
		'#define numGhost '..self.numGhost,
		'#define numStates '..self.eqn.numStates,
		'#define numWaves '..self.eqn.numWaves,
	}:append(xs:map(function(x,i)
		return '#define mins_'..x..' '..clnumber(self.mins[i])..'\n'
			.. '#define maxs_'..x..' '..clnumber(self.maxs[i])..'\n'
	end)):append{
		'constant real3 mins = _real3(mins_x, '..(self.dim<2 and '0' or 'mins_y')..', '..(self.dim<3 and '0' or 'mins_z')..');', 
		'constant real3 maxs = _real3(maxs_x, '..(self.dim<2 and '0' or 'maxs_y')..', '..(self.dim<3 and '0' or 'maxs_z')..');', 
	}:append{
		'#define dx_min '..clnumber(math.min(table.unpack(self.dxs, 1, self.dim))),
	}:append(xs:map(function(name,i)
		return '#define gridSize_'..name..' '..tonumber(self.gridSize[name])
	end)):append{
		'constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);',
		'constant int4 stepsize = (int4)(1, gridSize_x, gridSize_x * gridSize_y, gridSize_x * gridSize_y * gridSize_z);',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		INDEX((i).x, (i).y, (i).z)',
	}:append(range(3):map(function(i)
		return (('#define grid_dx{i} ((maxs_{x} - mins_{x}) / (real)(gridSize_{x} - '..(2*self.numGhost)..'))')
			:gsub('{i}', i-1)
			:gsub('{x}', xs[i]))
	end)):append(range(3):map(function(i)
		return (('#define cell_x{i}(i) ((real)(i + '..clnumber(.5-self.numGhost)..') * grid_dx{i} + mins_'..xs[i]..')')
			:gsub('{i}', i-1))
	end)):append{
		'#define CELL_X(i) _real3(cell_x0(i.x), cell_x1(i.y), cell_x2(i.z));',

		self.eqn.getTypeCode and self.eqn:getTypeCode() or '',

		-- run here for the code, and in buffer for the sizeof()
		self.eqn:getEigenInfo().typeCode or '',
		
		-- bounds-check macro
		'#define OOB(lhs,rhs) (i.x < lhs || i.x >= gridSize_x - rhs'
			.. (self.dim < 2 and '' or ' || i.y < lhs || i.y >= gridSize_y - rhs')
			.. (self.dim < 3 and '' or ' || i.z < lhs || i.z >= gridSize_z - rhs')
			.. ')',
		
		-- define i, index, and bounds-check
		'#define SETBOUNDS(lhs,rhs)	\\',
		'\tint4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0); \\',
		'\tif (OOB(lhs,rhs)) return; \\',
		'\tint index = INDEXV(i);',
	}

	-- TODO replace this with geom code
	lines:append(range(self.dim):map(function(i)
		local code = self.geometry.dxCodes[i]
		for j=1,3 do
			code = code:gsub(
				'{x'..j..'}',
				'cell_x'..(j-1)..'(i.'..xs[j]..')')
		end
		return '#define dx_at'..(i-1)..'(i) (grid_dx'..(i-1)..' * ('..code..'))'
	end))

	lines:append(range(3):map(function(i)
		return getCode_define_i3_to_real3(
			'e'..(i-1)..'_at', 
			range(3):map(function(j)
				local code = (self.geometry.eCode[i] or {})[j] or '0'
				return code
			end))
	end)):append(range(3):map(function(i)
		return getCode_define_i3_to_real3(
			'e'..(i-1)..'unit_at', 
			range(3):map(function(j)
				return (self.geometry.eUnitCode[i] or {})[j] or '0'
			end))
	end))

	lines:append(self.displayVars:map(function(var,i)
		return '#define display_'..var.name..' '..i
	end))

	-- output the first and last indexes of display vars associated with each buffer
	for _,convertToTex in ipairs(self.convertToTexs) do
		assert(convertToTex.vars[1], "failed to find vars for convertToTex "..convertToTex.name)
		lines:insert('#define displayFirst_'..convertToTex.name..' display_'..convertToTex.vars[1].name)
		lines:insert('#define displayLast_'..convertToTex.name..' display_'..convertToTex.vars:last().name)
	end

	lines:append{
		self.checkFluxError and '#define checkFluxError' or '',
		self.checkOrthoError and '#define checkOrthoError' or '',
		self.allocateOneBigStructure and '#define allocateOneBigStructure' or '',
		errorTypeCode or '',
		self.eqn:getEigenInfo().code or '',
		self:getCoordMapCode() or '',
		
		-- this is dependent on coord map / length code
		self.eqn:getCodePrefix(self) or '',
	}
	
	self.codePrefix = lines:concat'\n'

	print('codePrefix:\n',self.codePrefix)
end

function Solver:refreshInitStateProgram()
	local initStateCode = table{
		self.codePrefix,
		self.eqn:getInitStateCode(self),
	}:concat'\n'
	self.initStateProgram = CLProgram{context=self.app.ctx, devices={self.app.device}, code=initStateCode}
	self.initStateKernel = self.initStateProgram:kernel('initState', self.UBuf)
end

function Solver:resetState()
	self.app.cmds:finish()
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initStateKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:finish()
end

function Solver:getCalcDTCode()
	return processcl(assert(file['solver/calcDT.cl']), {solver=self})
end

function Solver:refreshCommonProgram()
	-- code that depend on real and nothing else
	-- TODO move to app, along with reduceBuf

	local commonCode = table():append
		{self.app.is64bit and '#pragma OPENCL EXTENSION cl_khr_fp64 : enable' or nil
	}:append{
		'typedef '..self.app.real..' real;',
		-- used to find the min/max of a buffer
		processcl(assert(file['solver/reduce.cl']), {
			name = 'reduceMin',
			initValue = 'INFINITY',
			op = function(x,y) return 'min('..x..', '..y..')' end,
		}),
		processcl(assert(file['solver/reduce.cl']), {
			name = 'reduceMax',
			initValue = '-INFINITY',
			op = function(x,y) return 'max('..x..', '..y..')' end,
		}),	
		-- used by the integrators
		[[
__kernel void multAdd(
	__global real* a,
	const __global real* b,
	const __global real* c,
	real d
) {
	size_t i = get_global_id(0);
	if (i >= get_global_size(0)) return;
	a[i] = b[i] + c[i] * d;
}	
]],
	}:concat'\n'

	self.commonProgram = CLProgram{context=self.app.ctx, devices={self.app.device}, code=commonCode}

	for _,name in ipairs{'Min', 'Max'} do
		self['reduce'..name..'Kernel'] = self.commonProgram:kernel(
			'reduce'..name,
			self.reduceBuf,
			{ptr=nil, size=self.localSize1d * ffi.sizeof(self.app.real)},
			ffi.new('int[1]', self.volume),
			self.reduceSwapBuf)
	end

	-- used by the integrators
	self.multAddKernel = self.commonProgram:kernel'multAdd'
end

function Solver:getSolverCode()
	local slopeLimiterCode = 'real slopeLimiter(real r) {'
		.. self.app.slopeLimiters[1+self.slopeLimiterPtr[0]].code 
		.. '}'
		
	return table{
		self.codePrefix,
		slopeLimiterCode,
		
		'typedef struct { real min, max; } range_t;',
		self.eqn:getSolverCode(self) or '',

		self:getCalcDTCode() or '',
		
		processcl(assert(file['solver/solver.cl'], {solver=self}))
	}:concat'\n'
end

-- depends on buffers
function Solver:refreshSolverProgram()
	local code = self:getSolverCode()
	self.solverProgram = CLProgram{context=self.app.ctx, code=code}
	local success, message = self.solverProgram:build{self.app.device}
	if not success then
		-- show code
		print(showcode(code))
		-- show errors
--		message = string.split(string.trim(message),'\n'):filter(function(line) return line:find'error' end):concat'\n'
		error(message)	
	end

	self.calcDTKernel = self.solverProgram:kernel('calcDT', self.reduceBuf, self.UBuf)
	
	self.calcEigenBasisKernel = self.solverProgram:kernel('calcEigenBasis', self.waveBuf, self.eigenBuf, self.UBuf, self.fluxXformBuf)
	self.calcDeltaUTildeKernel = self.solverProgram:kernel('calcDeltaUTilde', self.deltaUTildeBuf, self.UBuf, self.eigenBuf)
	self.calcRTildeKernel = self.solverProgram:kernel('calcRTilde', self.rTildeBuf, self.deltaUTildeBuf, self.waveBuf)
	self.calcFluxKernel = self.solverProgram:kernel('calcFlux', self.fluxBuf, self.UBuf, self.waveBuf, self.eigenBuf, self.deltaUTildeBuf, self.rTildeBuf)
	
	self.calcDerivFromFluxKernel = self.solverProgram:kernel'calcDerivFromFlux'
	self.calcDerivFromFluxKernel:setArg(1, self.fluxBuf)
	if self.eqn.useSourceTerm then
		self.addSourceKernel = self.solverProgram:kernel'addSource'
		self.addSourceKernel:setArg(1, self.UBuf)
	end

	if self.checkFluxError or self.checkOrthoError then
		self.calcErrorsKernel = self.solverProgram:kernel('calcErrors', self.errorBuf, self.waveBuf, self.eigenBuf, self.fluxXformBuf)	
	end
end

function Solver:refreshDisplayProgram()

	local lines = table{
		self.codePrefix
	}

	if self.app.useGLSharing then
		for _,convertToTex in ipairs(self.convertToTexs) do
			lines:append{
				processcl(convertToTex:getCode(), {
					name = 'calcDisplayVarToTex_'..convertToTex.name,
					input = '__write_only '..(self.dim == 3 and 'image3d_t' or 'image2d_t')..' tex',
					output = '	write_imagef(tex, '
						..(self.dim == 3 and '(int4)(dsti.x, dsti.y, dsti.z, 0)' or '(int2)(dsti.x, dsti.y)')
						..', (float4)(value, 0., 0., 0.));',
					body = convertToTex.displayBodyCode or '',
					type = convertToTex.type,
					eigenBody = self.eqn:getCalcDisplayVarEigenCode(),
				})
			}
		end
	end

	for _,convertToTex in ipairs(self.convertToTexs) do
		lines:append{
			processcl(convertToTex:getCode(), {
				name = 'calcDisplayVarToBuffer_'..convertToTex.name,
				input = '__global real* dest',
				output = '	dest[dstindex] = value;',
				body = convertToTex.displayBodyCode or '',
				type = convertToTex.type,
				eigenBody = self.eqn:getCalcDisplayVarEigenCode(),
			})
		-- end display code
		}
	end

	local code = lines:concat'\n'
	self.displayProgram = CLProgram{context=self.app.ctx, code=code}
	local success, message = self.displayProgram:build{self.app.device}
	if not success then
		print(showcode(code))
		error(message)	
	end

	if self.app.useGLSharing then
		for _,convertToTex in ipairs(self.convertToTexs) do
			convertToTex.calcDisplayVarToTexKernel = self.displayProgram:kernel('calcDisplayVarToTex_'..convertToTex.name, self.texCLMem)
		end
	end

	for _,convertToTex in ipairs(self.convertToTexs) do
		convertToTex.calcDisplayVarToBufferKernel = self.displayProgram:kernel('calcDisplayVarToBuffer_'..convertToTex.name, self.reduceBuf)
	end
end

function Solver:createBoundaryProgramAndKernel(args)
	local bufType = assert(args.type)
	
	local lines = table()
	lines:insert(self.codePrefix)
	lines:insert((([[
__kernel void boundary(
	__global {type}* buf
) {
]]):gsub('{type}', bufType)))
	if self.dim == 1 then 
		lines:insert[[
	if (get_global_id(0) != 0) return;
]]
	elseif self.dim == 2 then
		lines:insert[[
	int i = get_global_id(0);
]]
	elseif self.dim == 3 then
		lines:insert[[
	int2 i = (int2)(get_global_id(0), get_global_id(1));
]]
	end
	
	-- 1D: use a small 1D kernel and just run once 
	-- 2D: use a 1D kernel the size of the max dim 
		
	lines:insert[[
	for (int j = 0; j < numGhost; ++j) {
]]

	for side=1,self.dim do

		if self.dim == 2 then
			if side == 1 then
				lines:insert[[
	if (i < gridSize_y) {
]]
			elseif side == 2 then
				lines:insert[[
	if (i < gridSize_x) {
]]
			end
		elseif self.dim == 3 then
			if side == 1 then
				lines:insert[[
	if (i.x < gridSize_y && i.y < gridSize_z) {
]]
			elseif side == 2 then
				lines:insert[[
	if (i.x < gridSize_x && i.y < gridSize_z) {
]]
			elseif side == 3 then
				lines:insert[[
	if (i.x < gridSize_x && i.y < gridSize_y) {
]]
			end
		end

		local function index(j)
			if self.dim == 1 then
				return j
			elseif self.dim == 2 then
				if side == 1 then
					return 'INDEX('..j..',i,0)'
				elseif side == 2 then
					return 'INDEX(i,'..j..',0)'
				end
			elseif self.dim == 3 then
				if side == 1 then
					return 'INDEX('..j..',i.x,i.y)'
				elseif side == 2 then
					return 'INDEX(i.x,'..j..',i.y)'
				elseif side == 3 then
					return 'INDEX(i.x,i.y,'..j..')'
				end
			end
			error'TODO'
		end

		local x = xs[side]

		local gridSizeSide = 'gridSize_'..xs[side]
		local rhs = gridSizeSide..'-numGhost+j'
		lines:insert(({
			periodic = '\t\tbuf['..index'j'..'] = buf['..index(gridSizeSide..'-2*numGhost+j')..'];',
			mirror = table{
				'\t\tbuf['..index'j'..'] = buf['..index'2*numGhost-1-j'..'];',
			}:append(table.map((args.mirrorVars or {})[side] or {}, function(var)
				return '\t\tbuf['..index'j'..'].'..var..' = -buf['..index'j'..'].'..var..';'
			end)):concat'\n',
			freeflow = '\t\tbuf['..index'j'..'] = buf['..index'numGhost'..'];',
		})[args.methods[x..'min']])

		lines:insert(({
			periodic = '\t\tbuf['..index(rhs)..'] = buf['..index'numGhost+j'..'];',
			mirror = table{
				'\t\tbuf['..index(rhs)..'] = buf['..index(gridSizeSide..'-numGhost-1-j')..'];',
			}:append(table.map((args.mirrorVars or {})[side] or {}, function(var)
				return '\t\tbuf['..index(rhs)..'].'..var..' = -buf['..index(rhs)..'].'..var..';'
			end)):concat'\n',
			freeflow = '\t\tbuf['..index(rhs)..'] = buf['..index(gridSizeSide..'-numGhost-1')..'];',
		})[args.methods[x..'max']])
	
		if self.dim > 1 then
			lines:insert'\t}'
		end
	end
	
	lines:insert'\t}'
	lines:insert'}'

	local code = lines:concat'\n'

	local boundaryProgram = CLProgram{context=self.app.ctx, devices={self.app.device}, code=code} 
	local boundaryKernel = boundaryProgram:kernel'boundary'
	return boundaryProgram, boundaryKernel
end

function Solver:refreshBoundaryProgram()
	self.boundaryProgram, self.boundaryKernel = 
		self:createBoundaryProgramAndKernel{
			type = 'cons_t',
			-- remap from enum/combobox int values to names
			methods = table.map(self.boundaryMethods, function(v,k)
				return self.app.boundaryMethods[1+v[0]], k
			end),
			mirrorVars = self.eqn.mirrorVars,
		}
	self.boundaryKernel:setArg(0, self.UBuf)
end

-- assumes the buffer is already in the kernel's arg
function Solver:applyBoundaryToBuffer(kernel)
	-- 1D:
	if self.dim == 1 then
		self.app.cmds:enqueueNDRangeKernel{kernel=kernel, globalSize=self.localSize1d, localSize=self.localSize1d}
	elseif self.dim == 2 then
		local maxSize = math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y))
		self.app.cmds:enqueueNDRangeKernel{kernel=kernel, globalSize=maxSize, localSize=math.min(self.localSize1d, maxSize)}
	elseif self.dim == 3 then
		print'TODO'
	else
		error("can't run boundary for dim "..tonumber(self.dim))
	end
end

function Solver:boundary()
	self:applyBoundaryToBuffer(self.boundaryKernel)
end

function Solver:reduceMin()
	return self:reduce(self.reduceMinKernel)
end

function Solver:reduceMax()
	return self:reduce(self.reduceMaxKernel)
end

function Solver:reduce(kernel)
	local reduceSize = self.volume
	local dst = self.reduceSwapBuf
	local src = self.reduceBuf
	local iter = 0

	-- TODO should be the min of the rounded-up/down? power-of-2 of reduceSize
	local reduceLocalSize1D = math.min(reduceSize, self.maxWorkGroupSize)
	
	while reduceSize > 1 do
		iter = iter + 1
		local nextSize = math.floor(reduceSize / reduceLocalSize1D)
		if 0 ~= bit.band(reduceSize, reduceLocalSize1D - 1) then 
			nextSize = nextSize + 1 
		end
		local reduceGlobalSize = math.max(reduceSize, reduceLocalSize1D)
		kernel:setArg(0, src)
		kernel:setArg(2, ffi.new('int[1]', reduceSize))
		kernel:setArg(3, dst)
		self.app.cmds:enqueueNDRangeKernel{kernel=kernel, dim=1, globalSize=reduceGlobalSize, localSize=math.min(reduceGlobalSize, reduceLocalSize1D)}
		--self.app.cmds:finish()
		dst, src = src, dst
		reduceSize = nextSize
	end
	self.app.cmds:enqueueReadBuffer{buffer=src, block=true, size=ffi.sizeof(self.app.real), ptr=self.reduceResultPtr}
	return self.reduceResultPtr[0]
end

function Solver:calcDeriv(derivBuf, dt)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcEigenBasisKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	if self.checkFluxError or self.checkOrthoError then
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcErrorsKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	end
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDeltaUTildeKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcRTildeKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.calcFluxKernel:setArg(6, ffi.new('real[1]', dt))
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}

	-- calcDerivFromFlux zeroes the derivative buffer
	self.calcDerivFromFluxKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}

	-- addSource adds to the derivative buffer
	if self.eqn.useSourceTerm then
		self.addSourceKernel:setArg(0, derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.addSourceKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	end
end

function Solver:update()
	self:boundary()

	-- calc cell wavespeeds -> dts
	if self.useFixedDT[0] then
		dt = tonumber(self.fixedDT[0])
	else
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDTKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		dt = tonumber(self.cfl[0]) * self:reduceMin()
	end
	self.fixedDT[0] = dt

	self:step(dt)
end

function Solver:step(dt)
	self.integrator:integrate(dt, function(derivBuf)
		self:calcDeriv(derivBuf, dt)
	end)

	self.t = self.t + dt
end

function Solver:calcDisplayVarToTex(varIndex)
	local var = self.displayVars[varIndex]
	local convertToTex = var.convertToTex
	if self.app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		self.app.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	
		convertToTex:setToTexArgs(var)
		self.app.cmds:enqueueNDRangeKernel{kernel=convertToTex.calcDisplayVarToTexKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		self.app.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		self.app.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex
		
		convertToTex:setToBufferArgs(var)
		self.app.cmds:enqueueNDRangeKernel{kernel=convertToTex.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		self.app.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(self.app.real) * self.volume, ptr=ptr}
		if self.app.is64bit then
			for i=0,self.volume-1 do
				ffi.cast('float*',ptr)[i] = ffi.cast('double*',ptr)[i]
			end
		end
		tex:bind()
		if self.dim < 3 then
			gl.glTexSubImage2D(gl.GL_TEXTURE_2D, 0, 0, 0, tex.width, tex.height, gl.GL_RED, gl.GL_FLOAT, ptr)
		else
			for z=0,tex.depth-1 do
				gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, tex.width, tex.height, 1, gl.GL_RED, gl.GL_FLOAT, ptr + tex.width * tex.height * z)
			end
		end
		tex:unbind()
		glreport'here'
	end
end

function Solver:calcDisplayVarRange(varIndex)
	local var = self.displayVars[varIndex]
	local convertToTex = var.convertToTex
	
	convertToTex:setToBufferArgs(var)
	
	self.app.cmds:enqueueNDRangeKernel{kernel=convertToTex.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	local ymin = self:reduceMin()
	self.app.cmds:enqueueNDRangeKernel{kernel=convertToTex.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	local ymax = self:reduceMax()
	return ymin, ymax
end

return Solver
