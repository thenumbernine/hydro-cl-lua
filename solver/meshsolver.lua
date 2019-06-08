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
local SolverBase = require 'solver.solverbase'


local MeshSolver = class(SolverBase)


MeshSolver.meshTypeCode = [[

enum { maxNumSides = 4 };
enum { maxNumVtxsPerCell = 4 };	// should match 'maxNumSides' for 2D

typedef struct cell_s {
	real3 x;		//center.  technically could be a 'realN'
	real volume;	//volume of the cell
	real maxDist;	//max distance between the centers of any two cells
	int numSides;	//number of interfaces
	int ifaces[maxNumSides];	//list of interfaces. maybe I should make this a pointer to another list?
	int numVtxs;	//number of vertices
	int vtxs[maxNumVtxsPerCell];
} cell_t;

typedef struct iface_s {
	real3 x;		//center.  realN.
	real3 normal;	//normal pointing from first to second
	real area;		//edge length / surface area
	real dist;		//dist between cell centers along 'normal'
	int cellIndex[2];	//indexes of cells
	
	//2D-only
	real3 delta;
} iface_t;
]]

local function real3(v)
	local x = ffi.new'real3'
	x.x = 0 x.y = 0 x.z = 0
	for i=1,#v do
		x.s[i-1] = v[i]
	end
	return x
end

local Vertex = class()
function Vertex:init(x,y,z)
	self.x = vec3(
		tonumber(x) or 0,
		tonumber(y) or 0,
		tonumber(z) or 0)
	self.edges = table()
end

local function addVtx(vtxs, x)
	if mergeUniqueVtxs then
		for _,v in ipairs(vtxs) do
			if v.x == x then return v end
		end
	end
	local v = Vertex()
	for i=1,3 do
		v.x[i] = tonumber(x[i]) or 0
	end
	vtxs:insert(v)
	return v
end

local Edge = class()
function Edge:init(...)
	self.vtxs = table{...}
	self.cells = table()
end

local function addEdge(edges,vtxs,i,j)
	local a,b = vtxs[i], vtxs[j]
	for _,e in ipairs(edges) do
		if (e.vtxs[1] == a and e.vtxs[2] == b)
		or (e.vtxs[1] == b and e.vtxs[2] == a)
		then
			return e
		end
	end
	local e = Edge(a,b)
	edges:insert(e)
	a.edges:insert(e)
	b.edges:insert(e)
	return e
end

local Cell = class()
function Cell:init()
	self.edges = table()
	self.vtxs = table()
end

local function addCell(cells, edges, vtxs, ...)
	local c = Cell()
	local n = select('#', ...)
	local indexes = table()
	for i=1,n do
		local vtxindex = select(i, ...)
		local v = assert(vtxs[vtxindex], "couldn't find vertex "..vtxindex)
		c.vtxs:insert(v)
		indexes:insert(vtxindex)
	end
	for i=1,n do
		c.edges:insert(addEdge(edges, vtxs, indexes[i], indexes[i%n+1]))
	end
	cells:insert(c)
	c.x = c.vtxs:map(function(v) return v.x end):sum() / n
	for _,e in ipairs(c.edges) do
		e.cells:insert(c)
	end
	return c
end

-- this is 2D
-- TODO n-d functions for this
local function polyVol(...)
	local n = select('#', ...)
	local v = 0
	for i=1,n do
		local pi = select(i, ...)
		local pj = select(i%n+1, ...)
		v = v + .5 * (pi[1] * pj[2] - pi[2] * pj[1])
	end
	return v
end



--[[
args:
	meshfile = name of mesh file to use

NOTICE initState is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:initL1(args)
	MeshSolver.super.initL1(self, args)

	ffi.cdef(self.meshTypeCode)

	local meshfilename = assert(args.meshfile)
	local ls = string.split(assert(string.trim(file['grids/'..meshfilename..'.p2dfmt'])), '\n')
	local first = ls:remove(1)
	local m, n = string.split(string.trim(ls:remove(1)), '%s+'):map(function(l) return tonumber(l) end):unpack()
	local x = string.split(string.trim(ls:concat()), '%s+'):map(function(l) return tonumber(l) end)
	assert(#x == 2*m*n)
	print(m, n, m*n)
	-- [[
	local us = x:sub(1,m*n)
	local vs = x:sub(m*n+1)
	assert(#us == #vs)
	local vtxs = table()
	for i=1,#us do
		local u,v = us[i], vs[i]
		addVtx(vtxs, {u,v})
	end

	assert(#vtxs == m*n, "expected "..#vtxs.." to equal "..(m*n))
	
	local cells = table()
	local edges = table()
	for i=1,n-1 do
		for j=1,m-1 do
			local c = addCell(
				cells,
				edges,
				vtxs,
				1 + j-1 + m * (i-1),
				1 + j-1 + m * i,
				1 + j + m * i,
				1 + j + m * (i-1))
		
			--[[ TODO do this later, within the initState shader
			c.U = consFromPrim({
				1,
				.1, 0, 0,
				1
			})
			--]]
		end
	end

	for _,e in ipairs(edges) do
		local a,b = e.vtxs:unpack()
		e.x = (a.x + b.x) * .5
		e.delta = a.x - b.x
		e.length = e.delta:length()
		e.normal = vec3(-e.delta[2], e.delta[1], e.delta[3])
		e.normal = e.normal / e.normal:length()
		
		local a,b = e.cells:unpack()
		if a and b then
			-- make sure the normal points to a
			if (a.x - e.x):dot(e.normal) < 0 then
				a,b = b,a
				e.cells = table{a,b}
			end
			e.cellDist = (b.x - a.x):length()
		else
			local c = a or b
			-- for ghost state's sake:
			e.cellDist = (c.x - e.x):length() * 2
		end
	end

	for i,v in ipairs(vtxs) do
		v.index = i
	end
	
	for i,e in ipairs(edges) do
		e.index = i	
	end

	for i,c in ipairs(cells) do
		c.index = i
		c.volume = polyVol(c.vtxs[1].x, c.vtxs[2].x, c.vtxs[3].x, c.vtxs[4].x)
		assert(#c.edges <= ffi.C.maxNumSides)
		assert(#c.vtxs <= ffi.C.maxNumVtxsPerCell)
	end


	self.mins = vec3(-1, -1, -1)
	self.maxs = vec3(1, 1, 1)
	for i=1,self.dim do
		self.mins[i] = math.huge
		self.maxs[i] = -math.huge
	end
	for _,v in ipairs(vtxs) do
		for i=1,self.dim do
			self.mins[i] = math.min(self.mins[i], v.x[i])
			self.maxs[i] = math.max(self.maxs[i], v.x[i])
		end
	end

	-- convert vertices
	-- only store position
	local cpu_vtxs = ffi.new('real3[?]', #vtxs)
	for i,v in ipairs(vtxs) do
		local cpu_vtx = cpu_vtxs[i-1]
		for j=1,3 do
			cpu_vtx.s[j-1] = vtxs[i].x[j]
		end
	end

	-- convert edges
	local cpu_ifaces = ffi.new('iface_t[?]', #edges)
	for i,edge in ipairs(edges) do
		local cpu_iface = cpu_ifaces[i-1]
		for j=1,3 do
			cpu_iface.x.s[j-1] = edge.x[j]
			cpu_iface.delta.s[j-1] = edge.delta[j]
			cpu_iface.normal.s[j-1] = edge.normal[j]
		end
		cpu_iface.area = edge.length
		cpu_iface.dist = edge.cellDist

		if #edge.cells > 2 then
			error("expected interfaces to have <=2 cells but found "..#edge.cells)
		end
		for j,cell in ipairs(edge.cells) do
			cpu_iface.cellIndex[j-1] = not cell and -1 or cell.index-1
		end
	end

	local cpu_cells = ffi.new('cell_t[?]', #cells)
	for i,cell in ipairs(cells) do
		local cpu_cell = cpu_cells[i-1]
		for j=1,3 do
			cpu_cell.x.s[j-1] = cell.x[j]
		end
		cpu_cell.volume = cell.volume
		cpu_cell.numSides = #cell.edges
		for i,edge in ipairs(cell.edges) do
			cpu_cell.ifaces[i-1] = edge.index-1
		end
		cpu_cell.numVtxs = #cell.vtxs
		for i,vtx in ipairs(cell.vtxs) do
			cpu_cell.vtxs[i-1] = vtx.index-1
		end
		cpu_cell.maxDist = 0
		for i=1,#cell.edges-1 do
			for j=i+1,#cell.edges do
				local dist = (cell.edges[i].x - cell.edges[j].x):length()
				cpu_cell.maxDist = math.max(cpu_cell.maxDist, dist)
			end
		end
	end

	self.mesh = {
		vtxs = cpu_vtxs,
		-- list of vertex indexes for each cell
		-- how about cell centers too?
		cells = cpu_cells,
		ifaces = cpu_ifaces,
	}
	-- TODO put these in solver_t
	self.numCells = #cells
	self.numInterfaces = #edges
	self.numVtxs = #vtxs
	
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
		args.cacheFile = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
		Program.super.init(self, args)
	end
	self.Program = Program
end

function MeshSolver:createBuffers()
	MeshSolver.super.createBuffers(self)

	self:clalloc('vtxBuf', 'real3', self.numVtxs)
	self:clalloc('cellsBuf', 'cell_t', self.numCells)
	self:clalloc('ifacesBuf', 'iface_t', self.numInterfaces)
end

function MeshSolver:finalizeCLAllocs()
	MeshSolver.super.finalizeCLAllocs(self)

	self.app.cmds:enqueueWriteBuffer{buffer=self.cellsBuf, block=true, size=ffi.sizeof'cell_t' * self.numCells, ptr=self.mesh.cells}
	self.app.cmds:enqueueWriteBuffer{buffer=self.ifacesBuf, block=true, size=ffi.sizeof'iface_t' * self.numInterfaces, ptr=self.mesh.ifaces}
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

function MeshSolver:boundary()
end

return MeshSolver 
