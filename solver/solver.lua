local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local math = require 'ext.math'
local vec3 = require 'vec.vec3'
local CLImageGL = require 'cl.imagegl'
local CLBuffer = require 'cl.obj.buffer'
local GLTex2D = require 'gl.tex2d'
local GLTex3D = require 'gl.tex3d'
local glreport = require 'gl.report'
local clnumber = require 'clnumber'
local template = require 'template'
local vec3sz = require 'ffi.vec.vec3sz'

local function getn(...)
	local t = {...}
	t.n = select('#', ...)
	return t
end
local function time(name, cb)
	print(name..'...')
	local startTime = os.clock()
	local result = getn(cb())
	local endTime = os.clock()
	print('...done '..name..' ('..(endTime - startTime)..'s)')
	return table.unpack(result, 1, result.n)
end

local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}

local function convertParams(code)
	code = code:gsub('{x^(%d)}', function(i)
		return 'x.'..xs[i+0]
	end)
	code = code:gsub('{v^(%d)}', function(i)
		return 'v.'..xs[i+0]
	end)
	return code
end

local function getCode_real3_to_real(name, code)
	return template([[
inline real <?=name?>(real3 x) {
	return <?=code?>;
}]], {
		name = name,
		code = convertParams(code),
	})
end

-- f(x) where x is a point in the coordinate chart
local function getCode_real3_to_real3(name, exprs)
	return template([[
inline real3 <?=name?>(real3 x) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] and convertParams(exprs[i]) or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end

-- f(v,x) where x is a point on the coordinate chart and v is most likely a tensor
local function getCode_real3_real3_to_real(name, expr)
	return template([[
inline real <?=name?>(real3 v, real3 x) {
	return <?=convertParams(expr)?>;
}]], {
		name = name,
		expr = expr,
		convertParams = convertParams,
	})
end

local function getCode_real3_real3_to_real3(name, exprs)
	return template([[
inline real3 <?=name?>(real3 v, real3 x) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] and convertParams(exprs[i]) or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end

local function getCode_real3_to_sym3(name, exprs)
	return template([[
inline sym3 <?=name?>(real3 x) {
	return (sym3){
<? for i=1,3 do
	for j=i,3 do
?>		.<?=xs[i]..xs[j]?> = <?=exprs[i] and exprs[i][j] 
			and convertParams(exprs[i][j]) or '0.'?>,
<?	end
end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end

local Solver = class()

Solver.name = 'Solver'
Solver.numGhost = 2

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

	self:createEqn(args.eqn)
	
	self.name = self.eqn.name..' '..self.name

	self.app = assert(args.app)
	
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	self.maxWorkGroupSize = tonumber(self.app.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')

	self.offset = vec3sz(0,0,0)
	self.localSize1d = math.min(self.maxWorkGroupSize, tonumber(self.gridSize:volume()))

	if self.dim == 3 then
		local localSizeX = math.min(gridSize[1], 2^math.ceil(math.log(self.maxWorkGroupSize,2)/2))
		local localSizeY = self.maxWorkGroupSize / localSizeX
		self.localSize2d = {localSizeX, localSizeY}
	end

--	self.localSize = self.dim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
	-- TODO better than constraining by math.min(gridSize),
	-- look at which gridSizes have the most room, and double them accordingly, until all of maxWorkGroupSize is taken up
	self.localSize = vec3sz(1,1,1)
	local rest = self.maxWorkGroupSize
	local localSizeX = math.min(gridSize[1], 2^math.ceil(math.log(rest,2)/self.dim))
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

	self.useFixedDT = false
	self.fixedDT = .001
	self.cfl = .5	--/self.dim
	self.initStatePtr = ffi.new('int[1]', (table.find(self.eqn.initStateNames, args.initState) or 1)-1)
	self.integratorPtr = ffi.new('int[1]', (self.integratorNames:find(args.integrator) or 1)-1)
	self.fluxLimiter = ffi.new('int[1]', (self.app.limiterNames:find(args.fluxLimiter) or 1)-1)
	
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

	self.usePLM = args.usePLM
	assert(not self.usePLM or self.fluxLimiter[0] == 0, "are you sure you want to use flux and slope limiters at the same time?")
	self.slopeLimiter = ffi.new('int[1]', (self.app.limiterNames:find(args.slopeLimiter) or 1)-1)

	self:refreshGridSize()
end

function Solver:getCoordMapCode()
	return table{
		getCode_real3_to_real3('coordMap', range(3):map(function(i)
			return self.geometry.uCode[i] or '{x^'..i..'}'
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

	self.domain = self.app.env:domain{
		size = {self.gridSize:unpack()},
		dim = dim,
	}

	-- don't include the ghost cells as a part of the grid coordinate space
	self.sizeWithoutBorder = vec3sz(self.gridSize:unpack())
	for i=0,self.dim-1 do
		self.sizeWithoutBorder:ptr()[i] = self.sizeWithoutBorder:ptr()[i] - 2 * self.numGhost
	end

	self.volume = tonumber(self.gridSize:volume())
	self.dxs = vec3(range(3):map(function(i)
		return (self.maxs[i] - self.mins[i]) / tonumber(self.sizeWithoutBorder:ptr()[i-1])
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
kernel void <?=name?>(
	<?=input?>,
	const global <?= convertToTex.type ?>* buf
	<?= #convertToTex.extraArgs > 0 
		and ','..table.concat(convertToTex.extraArgs, ',\n\t')
		or '' ?>
) {
	SETBOUNDS(0,0);
	int dstindex = index;
	int4 dsti = i;
	
	real3 x = cell_x(i);
	real3 xInt[<?=solver.dim?>];
<? for i=0,solver.dim-1 do
?>	xInt[<?=i?>] = x;
	xInt[<?=i?>].s<?=i?> -= .5 * grid_dx<?=i?>;
<? end
?>
	//now constrain
	if (i.x < 2) i.x = 2;
	if (i.x > gridSize_x - 2) i.x = gridSize_x - 2;
<? if solver.dim >= 2 then ?>
	if (i.y < 2) i.y = 2;
	if (i.y > gridSize_y - 2) i.y = gridSize_y - 2;
<? end
if solver.dim >= 3 then ?>
	if (i.z < 2) i.z = 2;
	if (i.z > gridSize_z - 2) i.z = gridSize_z - 2;
<? end ?>
	//and recalculate read index
	index = INDEXV(i);
	
	int side = 0;
	int indexInt = side + dim * index;
	real value = 0;

<?= convertToTex.varCodePrefix or '' ?>
<?= var.code ?>

<?= output ?>
}
]]

ConvertToTex.extraArgs = {}

function ConvertToTex:init(args)
	local solver = assert(args.solver)	
	self.name = assert(args.name)
	self.solver = solver
	self.type = args.type	-- or self.type
	self.extraArgs = args.extraArgs
	self.displayCode = args.displayCode	-- or self.displayCode
	self.displayBodyCode = args.displayBodyCode
	self.varCodePrefix = args.varCodePrefix
	self.vars = table()
	for i,var in ipairs(args.vars) do
		assert(type(var) == 'table', "failed on var "..self.name)
		local name, code = next(var)
		self.vars:insert{
			convertToTex = self,
			code = code,
			name = self.name..'_'..name,
			enabled = ffi.new('bool[1]', 
				self.name == 'U' and (solver.dim==1 or i==1)
				or (self.name == 'error' and solver.dim==1)
			),
			useLogPtr = ffi.new('bool[1]', 
				args.useLog or false
			),
			color = vec3(math.random(), math.random(), math.random()):normalize(),
			--heatMapTexPtr = ffi.new('int[1]', 0),	-- hsv, isobar, etc ...
			heatMapFixedRangePtr = ffi.new('bool[1]', false),	-- self.name ~= 'error'
			heatMapValueMinPtr = ffi.new('float[1]', 0),
			heatMapValueMaxPtr = ffi.new('float[1]', 1),
		}
	end
end

function ConvertToTex:setArgs(kernel, var)
	kernel:setArg(1, self.solver[var.convertToTex.name..'Buf'])
end

function ConvertToTex:setToTexArgs(var)
	self:setArgs(var.calcDisplayVarToTexKernel, var)
end

function ConvertToTex:setToBufferArgs(var)
	self:setArgs(var.calcDisplayVarToBufferKernel, var)
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
end

function Solver:addConvertToTexUBuf()
	self:addConvertToTex{
		name = 'U',
		type = self.eqn.cons_t,
		varCodePrefix = self.eqn:getDisplayVarCodePrefix(),
		vars = assert(self.eqn:getDisplayVars()),
	}
end

function Solver:addConvertToTexs()
	self:addConvertToTexUBuf()
	
	-- might contain nonsense :-p
	self:addConvertToTex{
		name = 'reduce', 
		vars = {{['0'] = 'value = buf[index];'}},
	}
end

-- my best idea to work around the stupid 8-arg max kernel restriction
-- this is almost as bad of an idea as using OpenCL was to begin with
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
		if not self.allocateOneBigStructure then
			if size % ffi.sizeof(self.app.env.real) ~= 0 then
				print()
				print'!!!!!!!!!!! WARNING !!!!!!!!!!!'
				print(' unaligned buffer: '..name)
				print()
				size = size + ffi.sizeof(self.app.env.real)
			end
		end
		total = total + size
		if not self.allocateOneBigStructure then
			local bufObj = CLBuffer{
				env = self.app.env,
				name = name,
				type = 'real',
				size = size / ffi.sizeof(self.app.env.real),
			}
			self[name..'Obj'] = bufObj
			self[name] = bufObj.obj
		end
	end
	if self.allocateOneBigStructure then
		self.oneBigBuf = self.app.ctx:buffer{rw=true, size=total}
	end
end

function Solver:getConsLRTypeCode()
	return template([[
typedef union {
	<?=eqn.cons_t?> LR[2];
	struct {
		<?=eqn.cons_t?> L, R;
	};
} <?=eqn.consLR_t?>;
]], {
	eqn = self.eqn,
})
end

function Solver:createBuffers()
	local realSize = ffi.sizeof(self.app.real)

	-- to get sizeof
	ffi.cdef(self.eqn:getTypeCode())
	ffi.cdef(self:getConsLRTypeCode())

	-- for twofluid, cons_t has been renamed to euler_maxwell_t and maxwell_cons_t
	if ffi.sizeof(self.eqn.cons_t) ~= self.eqn.numStates * ffi.sizeof'real' then
	   error('expected sizeof('..self.eqn.cons_t..') to be '
		   ..self.eqn.numStates..' * sizeof(real) = '..(self.eqn.numStates * ffi.sizeof'real')
		   ..' but found '..ffi.sizeof(self.eqn.cons_t))
	end

	-- should I put these all in one AoS?
	
	-- two fluid uses multiple solvers, each with their own cons_t ... 
	-- but in ffi, typedefs can't be undefined / redefined
	self:clalloc('UBuf', self.volume * ffi.sizeof(self.eqn.cons_t))
	
	if self.usePLM then
		self:clalloc('ULRBuf', self.volume * self.dim * ffi.sizeof(self.eqn.consLR_t))
	end

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	self:clalloc('reduceBuf', self.volume * realSize)
	self:clalloc('reduceSwapBuf', self.volume * realSize / self.localSize1d)
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- CL/GL interop

	local cl = self.dim < 3 and GLTex2D or GLTex3D
	self.tex = cl{
		width = tonumber(self.gridSize.x),
		height = tonumber(self.gridSize.y),
		depth = tonumber(self.gridSize.z),
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

	-- real3
	lines:insert(self.app.real3TypeCode)
	lines:insert(self.app.real3Code)
	lines:insert(self.app.sym3TypeCode)
	lines:insert(self.app.sym3Code)

	lines:append{
		'#define geometry_'..self.geometry.name..' 1',
		'#define dim '..self.dim,
		'#define numGhost '..self.numGhost,
		'#define numStates '..self.eqn.numStates,
		'#define numWaves '..self.eqn.numWaves,
	}:append(xs:map(function(x,i)
	-- coordinate space = u,v,w
	-- cartesian space = x,y,z
	-- min and max in coordinate space
		return '#define mins_'..x..' '..clnumber(self.mins[i])..'\n'
			.. '#define maxs_'..x..' '..clnumber(self.maxs[i])..'\n'
	end)):append{
		'constant real3 mins = _real3(mins_x, '..(self.dim<2 and '0' or 'mins_y')..', '..(self.dim<3 and '0' or 'mins_z')..');', 
		'constant real3 maxs = _real3(maxs_x, '..(self.dim<2 and '0' or 'maxs_y')..', '..(self.dim<3 and '0' or 'maxs_z')..');', 
	}:append(xs:map(function(name,i)
	-- grid size
		return '#define gridSize_'..name..' '..tonumber(self.gridSize[name])
	end)):append{
		'constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);',
		'constant int4 stepsize = (int4)(1, gridSize_x, gridSize_x * gridSize_y, gridSize_x * gridSize_y * gridSize_z);',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		indexForInt4ForSize(i, gridSize_x, gridSize_y, gridSize_z)',
	
	--[[
	naming conventions ...
	* the grid indexes i_1..i_n that span 1 through gridSize_1..gridSize_n
	(between the index and the coordinate space:)
		- grid_dx? is the change in coordinate space wrt the change in grid space
		- cell_x(i) calculates the coordinates at index i
	* the coordinates space x_1..x_n that spans mins.s1..maxs.sn
	(between the coordinate and the embedded space:)
		- vectors can be calculated from Cartesian by cartesianToCoord
		- the length of the basis vectors wrt the change in indexes is given by dx?_at(x)
		- the Cartesian length of the holonomic basis vectors is given by coordHolBasisLen.  
			This is like dx?_at except not scaled by grid_dx?
			This is just the change in embedded wrt the change in coordinate, not wrt the change in grid
		- volume_at(x) gives the volume between indexes at the coordinate x
		- the Cartesian length of a vector in coordinate space is given by coordLen and coordLenSq
	* the embedded Cartesian space ... idk what letters I should use for this.  
		Some literature uses x^I vs coordinate space x^a, but that's still an 'x', no good for programming.
		Maybe I'll use 'xc' vs 'x', like I've already started to do in the initial conditions.
	--]]
	
	}:append(range(3):map(function(i)
	-- this is the change in coordinate wrt the change in code
	-- delta in coordinate space along one grid cell
		return (('#define grid_dx{i} ((maxs_{x} - mins_{x}) / (real)(gridSize_{x} - '..(2*self.numGhost)..'))')
			:gsub('{i}', i-1)
			:gsub('{x}', xs[i]))
	end)):append(range(3):map(function(i)
	-- mapping from index to coordinate 
		return (('#define cell_x{i}(i) ((real)(i + '..clnumber(.5-self.numGhost)..') * grid_dx{i} + mins_'..xs[i]..')')
			:gsub('{i}', i-1))
	end)):append{
		'#define cell_x(i) _real3(cell_x0(i.x), cell_x1(i.y), cell_x2(i.z))',

		self.eqn.getTypeCode and self.eqn:getTypeCode() or '',

		-- run here for the code, and in buffer for the sizeof()
		self.eqn:getEigenTypeCode() or '',
		
		-- bounds-check macro
		'#define OOB(lhs,rhs) (i.x < lhs || i.x >= gridSize_x - rhs'
			.. (self.dim < 2 and '' or ' || i.y < lhs || i.y >= gridSize_y - rhs')
			.. (self.dim < 3 and '' or ' || i.z < lhs || i.z >= gridSize_z - rhs')
			.. ')',
		
		-- define i, index, and bounds-check
		'#define SETBOUNDS(lhs,rhs)	\\',
		'\tint4 i = globalInt4(); \\',
		'\tif (OOB(lhs,rhs)) return; \\',
		'\tint index = INDEXV(i);',
	}

	-- dx0, ...
	-- this is the change in cartesian wrt the change in grid
	lines:append(range(self.dim):map(function(i)
		local code = self.geometry.dxCodes[i]
		for j=1,3 do
			code = code:gsub(
				'{x^'..j..'}',
				'cell_x'..(j-1)..'(i.'..xs[j]..')')
		end
		return '#define dx'..(i-1)..'_at(i) (grid_dx'..(i-1)..' * ('..code..'))'
	end))
	
	-- volume
	local volumeCode = '(' .. self.geometry.volumeCode .. ')'
	for i=1,self.dim do
		volumeCode = volumeCode .. ' * grid_dx'..(i-1)
	end
	lines:insert(getCode_real3_to_real('volume_at', volumeCode))
	
	-- coord len code: l(v) = v^i v^j g_ij
	lines:append{
		getCode_real3_real3_to_real('coordLenSq', self.geometry.uLenSqCode),
		[[
inline real coordLen(real3 r, real3 x) {
	return sqrt(coordLenSq(r, x));
}]],
	}

	lines:insert(getCode_real3_real3_to_real3('coord_conn', self.geometry.connCodes))

	--[[
	for i=0,self.dim-1 do
		lines:insert(getCode_real3_to_real('coordHolBasisLen'..i, self.geometry.eHolLenCode[i+1]))
	end
	--]]

	for i,eiCode in ipairs(self.geometry.eCode) do
		lines:insert(getCode_real3_to_real3('coordBasis'..(i-1), eiCode))
	end

	lines:insert(getCode_real3_real3_to_real3('coord_lower', self.geometry.lowerCodes))

	do
		local function addSym3Components(name, codes)
			for i=1,3 do
				for j=i,3 do
					local code = (codes[i] and codes[i][j] and convertParams(codes[i][j]) or clnumber(i==j and 1 or 0))
					lines:insert('#define '..name..(i-1)..(j-1)..'(r) '..code)
					if i ~= j then
						lines:insert('#define '..name..(j-1)..(i-1)..'(r) '..code)
					end
				end
			end
		end
		
		addSym3Components('coord_g', self.geometry.gCode)
		addSym3Components('coord_gU', self.geometry.gUCode)
		addSym3Components('coord_sqrt_gU', self.geometry.sqrt_gUCode)
		lines:insert(getCode_real3_to_sym3('coord_g', self.geometry.gCode))
		lines:insert(getCode_real3_to_sym3('coord_gU', self.geometry.gUCode))
	end

	lines:insert(template([[

//converts a vector from cartesian coordinates to grid coordinates
//by projecting the vector into the grid basis vectors 
//at x, which is in grid coordinates
real3 cartesianToCoord(real3 v, real3 x) {
	real3 vCoord;
	<? for i=0,solver.dim-1 do ?>{
		real3 e = coordBasis<?=i?>(x);
		//anholonomic normalized
		vCoord.s<?=i?> = real3_dot(e, v) / real3_len(e);
		//holonomic
		//vCoord.s<?=i?> = real3_dot(e, v) / real3_lenSq(e);
	}<? end
	for i=solver.dim,2 do ?>
	vCoord.s<?=i?> = 0.;
	<? end ?>
	return vCoord;
}
]], {
		solver = self,
	}))

	lines:append{
		-- not messing with this one yet
		self.allocateOneBigStructure and '#define allocateOneBigStructure' or '',
		
		self:getCoordMapCode() or '',
		
		-- this is dependent on coord map / length code
		self.eqn:getCodePrefix() or '',
		self:getConsLRTypeCode(),
	}

	self.codePrefix = lines:concat'\n'

print'done building solver.codePrefix'
--print(self.codePrefix)
end

function Solver:refreshInitStateProgram()
	local initStateCode = table{
		self.app.env.code,
		self.codePrefix,
		self.eqn:getInitStateCode(),
	}:concat'\n'
	time('compiling init state program', function()
		self.initStateProgram = self.app.ctx:program{devices={self.app.device}, code=initStateCode}
	end)
	self.initStateKernel = self.initStateProgram:kernel('initState', self.UBuf)
end

function Solver:resetState()
	self.app.cmds:finish()
	self.app.cmds:enqueueNDRangeKernel{
		kernel=self.initStateKernel, 
		dim=self.dim, 
		globalSize=self.gridSize:ptr(), 
		localSize=self.localSize:ptr()}
	self:boundary()
	self.app.cmds:finish()
	self.t = 0
end

function Solver:getCalcDTCode()
	if self.eqn.hasCalcDT then return end
	return template(file['solver/calcDT.cl'], {solver=self, eqn=self.eqn})
end

function Solver:refreshCommonProgram()
	-- code that depend on real and nothing else
	-- TODO move to app, along with reduceBuf

	local commonCode = table():append
		{self.app.is64bit and '#pragma OPENCL EXTENSION cl_khr_fp64 : enable' or nil
	}:append{
		'typedef '..self.app.real..' real;',
		-- used by the integrators
		[[
kernel void multAdd(
	global real* a,
	const global real* b,
	const global real* c,
	real d
) {
	size_t i = get_global_id(0);
	if (i >= get_global_size(0)) return;
	a[i] = b[i] + c[i] * d;
}	
]],
	}:concat'\n'

	time('compiling common program', function()
		self.commonProgram = self.app.ctx:program{devices={self.app.device}, code=commonCode}
	end)

	-- used by the integrators
	self.multAddKernel = self.commonProgram:kernel'multAdd'

	self.reduceMin = self.app.env:reduce{
		size = self.volume,
		op = function(x,y) return 'min('..x..', '..y..')' end,
		initValue = 'INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
	self.reduceMax = self.app.env:reduce{
		size = self.volume,
		op = function(x,y) return 'max('..x..', '..y..')' end,
		initValue = '-INFINITY',
		buffer = self.reduceBuf,
		swapBuffer = self.reduceSwapBuf,
		result = self.reduceResultPtr,
	}
end

function Solver:getSolverCode()
	local fluxLimiterCode = 'real fluxLimiter(real r) {'
		.. self.app.limiters[1+self.fluxLimiter[0]].code 
		.. '}'

	local slopeLimiterCode = 'real slopeLimiter(real r) {'
		.. self.app.limiters[1+self.slopeLimiter[0]].code 
		.. '}'
	
	return table{
		self.app.env.code,
		self.codePrefix,
		
		-- TODO move to Roe, or FiniteVolumeSolver as a parent of Roe and HLL?
		self.eqn:getEigenCode() or '',
		
		fluxLimiterCode,
		slopeLimiterCode,
		
		'typedef struct { real min, max; } range_t;',
		self.eqn:getSolverCode() or '',

		self:getCalcDTCode() or '',
	
		-- messing with this ...
		self.usePLM and template(file['solver/plm.cl'], {solver=self, eqn=self.eqn}) or '',
	}:concat'\n'
end

-- depends on buffers
function Solver:refreshSolverProgram()
	-- set pointer to the buffer holding the LR state information
	-- for piecewise-constant that is the original UBuf
	-- for piecewise-linear that is the ULRBuf
	self.getULRBuf = self.usePLM and self.ULRBuf or self.UBuf

	self.getULRArg = self.usePLM 
		and ('const global '..self.eqn.consLR_t..'* ULRBuf')
		or ('const global '..self.eqn.cons_t..'* UBuf')

	-- this code creates the const global cons_t* UL, UR variables
	-- it assumes that indexL, indexR, and side are already defined
	self.getULRCode = self.usePLM and ([[
	const global ]]..self.eqn.cons_t..[[* UL = &ULRBuf[side + dim * indexL].R;
	const global ]]..self.eqn.cons_t..[[* UR = &ULRBuf[side + dim * indexR].L;
]]) or ([[
	const global ]]..self.eqn.cons_t..[[* UL = UBuf + indexL;
	const global ]]..self.eqn.cons_t..[[* UR = UBuf + indexR;
]])

	local code = self:getSolverCode()

	time('compiling solver program', function()
		self.solverProgram = self.app.ctx:program{devices={self.app.device}, code=code}
	end)

--[[ trying to find out what causes stalls on certain kernels
do
	local cl = require 'ffi.OpenCL'
	
	local size = ffi.new('size_t[1]', 0)
	local err = cl.clGetProgramInfo(self.solverProgram.id, cl.CL_PROGRAM_BINARY_SIZES, ffi.sizeof(size), size, nil)
	if err ~= cl.CL_SUCCESS then error"failed to get binary size" end
print('size',size[0])

	local binary = ffi.new('char[?]', size[0])
	local err = cl.clGetProgramInfo(self.solverProgram.id, cl.CL_PROGRAM_BINARIES, size[0], binary, nil)
	if err ~= cl.CL_SUCCESS then error("failed to get binary: "..err) end
file['solverProgram.bin'] = ffi.string(binary, size[0])
end
--]]

	self:refreshCalcDTKernel()

	if self.usePLM then
		self.calcLRKernel = self.solverProgram:kernel(
			'calcLR',
			self.ULRBuf,
			self.UBuf)
	end
end

-- for solvers who don't rely on calcDT
function Solver:refreshCalcDTKernel()
	self.calcDTKernel = self.solverProgram:kernel(
		'calcDT',
		self.reduceBuf,
		self.UBuf)
end

function Solver:refreshDisplayProgram()

	local lines = table{
		self.app.env.code,
		self.codePrefix,
	}
	
	for _,convertToTex in ipairs(self.convertToTexs) do
		for _,var in ipairs(convertToTex.vars) do
			var.id = tostring(var):sub(10)
		end
	end

	if self.app.useGLSharing then
		for _,convertToTex in ipairs(self.convertToTexs) do
			for _,var in ipairs(convertToTex.vars) do
				if var.enabled[0] then
					lines:append{
						template(convertToTex.displayCode, {
							solver = self,
							var = var,
							convertToTex = convertToTex,
							name = 'calcDisplayVarToTex_'..var.id,
							input = '__write_only '
								..(self.dim == 3 
									and 'image3d_t' 
									or 'image2d_t'
								)..' tex',
							output = '	write_imagef(tex, '
								..(self.dim == 3 
									and '(int4)(dsti.x, dsti.y, dsti.z, 0)' 
									or '(int2)(dsti.x, dsti.y)'
								)..', (float4)(value, 0., 0., 0.));',
						})
					}
				end
			end
		end
	end

	for _,convertToTex in ipairs(self.convertToTexs) do
		for _,var in ipairs(convertToTex.vars) do
			if var.enabled[0] then
				lines:append{
					template(convertToTex.displayCode, {
						solver = self,
						var = var,
						convertToTex = convertToTex,
						name = 'calcDisplayVarToBuffer_'..var.id,
						input = 'global real* dest',
						output = '	dest[dstindex] = value;',
					})
				}
			end
		end
	end

	local code = lines:concat'\n'
	time('compiling display program', function()
		self.displayProgram = self.app.ctx:program{devices={self.app.device}, code=code}
	end)

	if self.app.useGLSharing then
		for _,convertToTex in ipairs(self.convertToTexs) do
			for _,var in ipairs(convertToTex.vars) do
				if var.enabled[0] then
					var.calcDisplayVarToTexKernel = self.displayProgram:kernel('calcDisplayVarToTex_'..var.id, self.texCLMem)
				end
			end
		end
	end

	for _,convertToTex in ipairs(self.convertToTexs) do
		for _,var in ipairs(convertToTex.vars) do
			if var.enabled[0] then
				var.calcDisplayVarToBufferKernel = self.displayProgram:kernel('calcDisplayVarToBuffer_'..var.id, self.reduceBuf)
			end
		end
	end
end

function Solver:createBoundaryProgramAndKernel(args)
	local assign = args.assign or function(a, b) return a .. ' = ' .. b end
	
	local lines = table()
	lines:insert(self.app.env.code)
	lines:insert(self.codePrefix)
	lines:insert(template([[
kernel void boundary(
	global <?=bufType?>* buf
) {
]], {
	bufType = args.type,
}))
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
		end

		local x = xs[side]

		local gridSizeSide = 'gridSize_'..xs[side]
		local rhs = gridSizeSide..'-numGhost+j'
		lines:insert(({
			periodic = '\t\t'..assign('buf['..index'j'..']', 'buf['..index(gridSizeSide..'-2*numGhost+j')..']')..';',
			mirror = table{
				'\t\t'..assign('buf['..index'j'..']', ' buf['..index'2*numGhost-1-j'..']')..';',
			}:append(table.map((args.mirrorVars or {})[side] or {}, function(var)
				return '\t\t'..'buf['..index'j'..'].'..var..' = -buf['..index'j'..'].'..var..';'
			end)):concat'\n',
			freeflow = '\t\t'..assign('buf['..index'j'..']', 'buf['..index'numGhost'..']')..';',
		})[args.methods[x..'min']])

		lines:insert(({
			periodic = '\t\t'..assign('buf['..index(rhs)..']', 'buf['..index'numGhost+j'..']')..';',
			mirror = table{
				'\t\t'..assign('buf['..index(rhs)..']', 'buf['..index(gridSizeSide..'-numGhost-1-j')..']')..';',
			}:append(table.map((args.mirrorVars or {})[side] or {}, function(var)
				return '\t\t'..'buf['..index(rhs)..'].'..var..' = -buf['..index(rhs)..'].'..var..';'
			end)):concat'\n',
			freeflow = '\t\t'..assign('buf['..index(rhs)..']', 'buf['..index(gridSizeSide..'-numGhost-1')..']')..';',
		})[args.methods[x..'max']])
	
		if self.dim > 1 then
			lines:insert'\t}'
		end
	end
	
	lines:insert'\t}'
	lines:insert'}'

	local code = lines:concat'\n'

	local boundaryProgram
	time('compiling boundary program', function()
		boundaryProgram = self.app.ctx:program{devices={self.app.device}, code=code} 
	end)
	local boundaryKernel = boundaryProgram:kernel'boundary'
	return boundaryProgram, boundaryKernel
end

function Solver:refreshBoundaryProgram()
	self.boundaryProgram, self.boundaryKernel = 
		self:createBoundaryProgramAndKernel{
			type = self.eqn.cons_t,
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
		-- xy xz yz
		local maxSizeX = math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y))
		local maxSizeY = math.max(tonumber(self.gridSize.y), tonumber(self.gridSize.z))
		self.app.cmds:enqueueNDRangeKernel{
			kernel = kernel,
			globalSize = {maxSizeX, maxSizeY},
			localSize = {
				math.min(self.localSize2d[1], maxSizeX),
				math.min(self.localSize2d[2], maxSizeY),
			},
		}
	else
		error("can't run boundary for dim "..tonumber(self.dim))
	end
end

function Solver:boundary()
	self:applyBoundaryToBuffer(self.boundaryKernel)
end

function Solver:calcDT()
	local dt
	-- calc cell wavespeeds -> dts
	if self.useFixedDT then
		dt = self.fixedDT
	else
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDTKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		dt = self.cfl * self.reduceMin()
		if not math.isfinite(dt) then
			print("got a bad dt!") -- TODO dump all buffers
		end
		self.fixedDT = dt
	end
	return dt
end

function Solver:update()
--local before = self.UBufObj:toCPU()
--print'\nself.UBufObj before boundary:' self:printBuf(self.UBufObj) print(debug.traceback(),'\n\n')	
	self:boundary()
--[[
local after = self.UBufObj:toCPU()
print'\nself.UBufObj after boundary:' self:printBuf(self.UBufObj) print(debug.traceback(),'\n\n')	
if self.app.dim == 2 then
	for i=0,tonumber(self.gridSize.x)-1 do
		for j=0,self.numGhost-1 do
			for _,s in ipairs{0,4} do
				for _,offset in ipairs{
					s + self.eqn.numStates * (i + self.gridSize.x * j),
					s + self.eqn.numStates * (j + self.gridSize.x * i),
					s + self.eqn.numStates * (i + self.gridSize.x * (self.gridSize.y - 1 - j)),
					s + self.eqn.numStates * ((self.gridSize.x - 1 - j) + self.gridSize.x * i),
				} do
					if after[offset] < 0 then
						local msg = "negative found in boundary after :boundary() was called"
						msg = msg .. "\nat position " .. tostring(offset / self.eqn.numStates)
						msg = msg .. '\nit was '..(before[offset] < 0 and 'negative' or 'positive')..' before the boundary update'
						error(msg)
					end
				end
			end
		end
	end
end
--]]	
	local dt = self:calcDT()
--print('dt',dt)	
	self:step(dt)
	self.t = self.t + dt
end

function Solver:step(dt)
	self.integrator:integrate(dt, function(derivBuf)
		self:calcDeriv(derivBuf, dt)
	end)
end

function Solver:printBuf(buf, ptr)
	ptr = ptr or buf:toCPU()
	local max = #tostring(self.volume-1)
	for i=0,self.volume-1 do
		io.write((' '):rep(max-#tostring(i)), i,':')
		for j=0,self.eqn.numStates-1 do
			io.write(' ', ptr[j + self.eqn.numStates * i])
		end 
		print()
	end
end

-- check for nans
-- expects buf to be of type cons_t, made up of numStates real variables
function Solver:checkFinite(buf)
	local ptr = buf:toCPU()
	local found
	for i=0,self.volume * self.eqn.numStates - 1 do
		if not math.isfinite(ptr[i]) then
			found = found or table()
			found:insert(i)
		end
	end
	if not found then return end
	self:printBuf(nil, ptr)
	print(found:map(tostring):concat', ')
	error'found non-finite numbers'
end

function Solver:getTex(var) 
	return self.tex
end

function Solver:calcDisplayVarToTex(var)
	local app = self.app
	local convertToTex = var.convertToTex
	if app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		app.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	
		convertToTex:setToTexArgs(var)
		app.cmds:enqueueNDRangeKernel{kernel=var.calcDisplayVarToTexKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		app.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		app.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex
		
		convertToTex:setToBufferArgs(var)
		app.cmds:enqueueNDRangeKernel{kernel=var.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		app.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * self.volume, ptr=ptr}
		if app.is64bit then
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

function Solver:calcDisplayVarRange(var)
	local convertToTex = var.convertToTex
	
	convertToTex:setToBufferArgs(var)
	
	self.app.cmds:enqueueNDRangeKernel{kernel=var.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	local ymin = self.reduceMin()
	self.app.cmds:enqueueNDRangeKernel{kernel=var.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	local ymax = self.reduceMax()
	return ymin, ymax
end

local float = ffi.new'float[1]'
local bool = ffi.new'bool[1]'

function Solver:updateGUIParams()
	if ig.igCollapsingHeader'parameters:' then
		bool[0] = self.useFixedDT
		if ig.igCheckbox('use fixed dt', bool) then
			self.useFixedDT = bool[0]
		end
		float[0] = self.fixedDT
		if ig.igInputFloat('fixed dt', float) then 
			self.fixedDT = float[0] 
		end
		float[0] = self.cfl
		if ig.igInputFloat('CFL', float) then 
			self.cfl = float[0] 
		end

		if ig.igCombo('integrator', self.integratorPtr, self.integratorNames) then
			self:refreshIntegrator()
		end

		if ig.igCombo('slope limiter', self.fluxLimiter, self.app.limiterNames) then
			self:refreshSolverProgram()
		end

		for i=1,self.dim do
			for _,minmax in ipairs(minmaxs) do
				local var = xs[i]..minmax
				if ig.igCombo(var, self.boundaryMethods[var], self.app.boundaryMethods) then
					self:refreshBoundaryProgram()
				end
			end
		end
	end
end

function Solver:updateGUIEqnSpecific()
	if ig.igCollapsingHeader'equation-specific:' then
		-- equation-specific:

		if ig.igCombo('init state', self.initStatePtr, self.eqn.initStateNames) then
			
			self:refreshInitStateProgram()
			
			-- TODO changing the init state program might also change the boundary methods
			-- ... but I don't want it to change the settings for the running scheme (or do I?)
			-- ... but I don't want it to not change the settings ...
			-- so maybe refreshing the init state program should just refresh everything?
			-- or maybe just the boundaries too?
			-- hack for now:
			self:refreshBoundaryProgram()
		end	
		
		local f = ffi.new'float[1]'
		local i = ffi.new'int[1]'
		for _,var in ipairs(self.eqn.guiVars) do
			var:updateGUI(self)
		end
	end
end

do
	-- display vars: TODO graph vars
	local function handle(var, title)
		ig.igPushIdStr(title)
		local enableChanged = ig.igCheckbox(var.name, var.enabled) 
		ig.igSameLine()
		if ig.igCollapsingHeader'' then	
			ig.igCheckbox('log', var.useLogPtr)
			ig.igCheckbox('fixed range', var.heatMapFixedRangePtr)
			ig.igInputFloat('value min', var.heatMapValueMinPtr)
			ig.igInputFloat('value max', var.heatMapValueMaxPtr)
		end
		ig.igPopId()
		return enableChanged 
	end

	-- do one for 'all'
	local function _and(a,b) return a and b end
	local fields = {'enabled', 'useLogPtr', 'heatMapFixedRangePtr', 'heatMapValueMinPtr', 'heatMapValueMaxPtr'}
	local types = {'bool', 'bool', 'bool', 'float', 'float'}
	local defaults = {true, true, true, math.huge, -math.huge}
	local combines = {_and, _and, _and, math.min, math.max}
	local all = {name='all'}
	for i=1,#fields do
		all[fields[i]] = ffi.new(types[i]..'[1]')
	end
	local original = {}
	for i,field in ipairs(fields) do
		original[field] = ffi.new(types[i]..'[1]')
	end

	function Solver:updateGUIDisplay()
		if ig.igCollapsingHeader'display:' then
			for i,convertToTex in ipairs(self.convertToTexs) do
				ig.igPushIdStr('display '..i)
				if ig.igCollapsingHeader(convertToTex.name) then				
					for i=1,#fields do
						all[fields[i]][0] = defaults[i]
					end
					for _,var in ipairs(convertToTex.vars) do
						for i,field in ipairs(fields) do
							all[field][0] = combines[i](all[field][0], var[field][0])
						end
					end
					for _,field in ipairs(fields) do
						original[field][0] = all[field][0]
					end
					handle(all, 'all')
					for _,field in ipairs(fields) do
						if all[field][0] ~= original[field][0] then
							for _,var in ipairs(convertToTex.vars) do
								var[field][0] = all[field][0]
							end
							if field == 'enabled' then
								self:refreshDisplayProgram()
							end
						end
					end

					for _,var in ipairs(convertToTex.vars) do
						if handle(var, convertToTex.name..' '..var.name) then
							self:refreshDisplayProgram()
						end
					end
				end
				ig.igPopId()
			end
		end
	end
end

function Solver:updateGUI()
	self:updateGUIParams()
	self:updateGUIEqnSpecific()
	self:updateGUIDisplay()

	-- heat map var

	-- TODO volumetric var
end

return Solver
