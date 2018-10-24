--[[
foundation class for solvers that use unstructured meshes from files

NACA 0012 from...
https://turbmodels.larc.nasa.gov/naca0012_grids.html
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local string = require 'ext.string'
local file = require 'ext.file'
local vec3 = require 'vec.vec3'
local unique = require 'util.unique'
local SolverBase = require 'solver.solverbase'


local MeshSolver = class(SolverBase)


local cellKeys = {
	'x',	-- real3
	'maxDist',
	'numSides',
}

local ifaceKeys = {
	'x',	-- real3
	'n',	-- real3
	'area',
	'dist',
}

MeshSolver.meshTypeCode = [[

enum { maxNumSides = 4 };

typedef struct cell_s {
	real3 x;		//center.  technically could be a 'realN'
	real maxDist;	//max distance of all interfaces
	int numSides;	//number of interfaces
	int ifaces[maxNumSides];	//list of interfaces. maybe I should make this a pointer to another list?
} cell_t;

typedef struct iface_s {
	real3 x;		//center.  realN.
	real3 n;		//normal pointing from first to second
	int cellIndex[2];	//indexes of cells
	real area;		//edge length / surface area
	real dist;		//dist between cell centers along 'n'	
} iface_t;
]]

--[[
args:
	meshfile = name of mesh file to use

NOTICE initState is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:preInit(args)
	MeshSolver.super.preInit(self, args)

	ffi.cdef(self.meshTypeCode)

	local meshfilename = assert(args.meshfile)
	local ls = assert(file['grids/'..meshfilename..'.p2dfmt']):trim():split'\n'
	local first = ls:remove(1)
	local m, n = ls:remove(1):trim():split'%s+':map(function(l) return tonumber(l) end):unpack()
	local x = ls:concat():trim():split'%s+':map(function(l) return tonumber(l) end)
	assert(#x == 2*m*n)
	print(m, n, m*n)
	-- [[
	local us = x:sub(1,m*n)
	local vs = x:sub(m*n+1)
	assert(#us == #vs)
	local vtxs = matrix()
	for i=1,#us do
		local u,v = us[i], vs[i]
		print(u,v)
		vtxs[#vtxs+1] = matrix{u,v}
	end

	assert(#vtxs == m*n, "expected "..#vtxs.." to equal "..(m*n))
	
	for i=1,n-1 do
		for j=1,m-1 do
			local c = addCell(
				1 + j-1 + m * (i-1),
				1 + j-1 + m * i,
				1 + j + m * i,
				1 + j + m * (i-1))
			
			c.U = consFromPrim(matrix{
				1,
				.1, 0, 0,
				1
			})
		end
	end


error'fixme'


	self.mins = matrix()
	self.maxs = matrix()
	for i=1,#size do
		self.mins[i] = math.huge
		self.maxs[i] = -math.huge
	end
	for _,v in ipairs(vtxs) do
		for i=1,#size do
			self.mins[i] = math.min(self.mins[i], v[i])
			self.maxs[i] = math.max(self.maxs[i], v[i])
		end
	end
	
	local function real3(v)
		local x = ffi.new'real3'
		x.x = 0 x.y = 0 x.z = 0
		for i=1,#v do
			x.s[i-1] = v[i]
		end
		return x
	end

	local vtxStepSize = {1, size[1], size[1] * size[2]}
	local cellStepSize = {1, size[1]-1, (size[1]-1) * (size[2]-1)}

	local cells = table.append(range(0,size[1]-2):map(function(i)
		return range(0,size[2]-2):map(function(j)
			local vtxIndexes = table{
				i + size[1] * j,
				i + 1 + size[1] * j,
				i + 1 + size[1] * (j + 1),
				i + size[1] * (j + 1),
			}
			local vtxs = vtxIndexes:map(function(index) return vtxs[index+1] end)
			local x = vtxs:sum() / #vtxs
			
			local ifaces = table()
			
			return {
				vtxIndexes = vtxIndexes, 
				vtxs = vtxs,
				x = x,
				ifaces = ifaces,
				numSides = #ifaces,
				maxDist = 1,	-- TODO
			}
		end)
	end):unpack())

	local ifaces = table.append(range(0,size[1]-3):map(function(i)
		return table.append(range(0,size[2]-3):map(function(j)
			return range(0,#size-1):map(function(k)
				local vtxIndexes = table{
					i + size[1] * j,
					i + size[1] * j + vtxStepSize[k+1],
				}
				local vtxs = vtxIndexes:map(function(index) return vtxs[index+1] end)
				local x = vtxs:sum() / #vtxs

				local cellIndexes = table{
					i + (size[1]-1) * j,
					i + (size[1]-1) * j + cellStepSize[k+1],
				}

				local cellL = cells[cellIndexes[1]+1]
				local cellR = cells[cellIndexes[2]+1]
				local delta = cellR.x - cellL.x
				local n = delta / delta:norm()

				return {
					x = x,
					n = n,
					area = (vtxs[2] - vtxs[1]):norm(),
					dist = delta:norm(),
				}
			end)
		end):unpack())
	end):unpack())

	for i,cell in ipairs(cells) do
		cell.x = real3(cell.x) 
	end
	for i,iface in ipairs(ifaces) do 
		iface.x = real3(iface.x) 
		iface.n = real3(iface.n)
	end

	self.mesh = {
		vtxs = vtxs,
		-- list of vertex indexes for each cell
		-- how about cell centers too?
		cells = cells,
		ifaces = ifaces,
	}
	-- TODO put these in solver_t
	self.numCells = #self.mesh.cells
	self.numInterfaces = (size[1] - 1) * size[2] + size[1] * (size[2] - 1)

	--[[
	next issue ... how to render this?
	how to store it?
	where do textures come into play?

	UBuf will be size of # of elements
	
	how about texture, and how about rendering?
	texsize will have to be > #elems ... using size is fine as long as we need a gridsize
	--]]

	-- no longer is dim * numCells the number of interfaces -- it is now dependent on the mesh
	-- maybe I should rename this to numInterfaces?

	local solver = self
	local Program = class(require 'cl.obj.program')
	function Program:init(args)
		args.env = solver.app.env
		args.domain = solver.domain
		args.cacheFile = 'cache-cl/'..unique(assert(args.name))
		Program.super.init(self, args)
	end
	self.Program = Program
end

function MeshSolver:createBuffers()
	MeshSolver.super.createBuffers(self)
	
	self:clalloc('cellsBuf', self.numCells * ffi.sizeof'cell_t')
	self:clalloc('ifacesBuf', self.numInterfaces * ffi.sizeof'iface_t')
end

function MeshSolver:finalizeCLAllocs()
	MeshSolver.super.finalizeCLAllocs(self)

	local cellsData = ffi.new('cell_t[?]', self.numCells)
	for i,cell in ipairs(self.mesh.cells) do
		for _,k in ipairs(cellKeys) do
			cellsData[i-1][k] = cell[k]
		end
	end
	
	local ifacesData = ffi.new('iface_t[?]', self.numInterfaces)
	for i,iface in ipairs(self.mesh.ifaces) do
		for _,k in ipairs(ifaceKeys) do
			ifacesData[i-1][k] = iface[k]
		end
	end

	self.app.cmds:enqueueWriteBuffer{buffer=self.cellsBuf, block=true, size=ffi.sizeof'cell_t' * self.numCells, ptr=cellsData}
	self.app.cmds:enqueueWriteBuffer{buffer=self.ifacesBuf, block=true, size=ffi.sizeof'iface_t' * self.numInterfaces, ptr=ifacesData}
end

function MeshSolver:refreshInitStateProgram()
	-- override this until I know what to do...
	--self.eqn.initState:refreshInitStateProgram(self)
end

function MeshSolver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	-- numCells is the number of cells
	-- maybe I should rename it to numCells
	local localSize1d = math.min(maxWorkGroupSize, self.numCells)
	
	self.domain = self.app.env:domain{
		size = {localSize1d},
		dim = self.dim,
	}
	
	return {
		localSize1d = localSize1d, 
	}
end


function MeshSolver:applyInitCond()
	assert(self.UBuf)
	-- TODO read it in from a file?  and upload it?  or something?
end

function MeshSolver:createCodePrefix()
	MeshSolver.super.createCodePrefix(self)

	local lines = table{
		self.codePrefix,
		self.meshTypeCode,
	}

	lines:insert[[
#define OOB(lhs,rhs) (index >= get_global_size(0))
#define SETBOUNDS(lhs,rhs)	\
	int index = get_global_id(0); \
	int4 i = (int4)(index,0,0,0);	\
	if (OOB(0,0)) return;
#define SETBOUNDS_NOGHOST()	\
	int index = get_global_id(0); \
	int4 i = (int4)(index,0,0,0);
#define cell_x(i) real3_zero
]]

	self.codePrefix = lines:concat'\n'
end

function MeshSolver:refreshCalcDTKernel()
	MeshSolver.super.refreshCalcDTKernel(self)
	-- TODO combine these, and offset one into the other?
	self.calcDTKernelObj.obj:setArg(2, cellsBuf)
	self.calcDTKernelObj.obj:setArg(3, ifacesBuf)
end

function MeshSolver:update()
end

return MeshSolver 
