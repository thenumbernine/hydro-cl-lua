--[[
foundation class for solvers that use unstructured meshes from files

NACA 0012 from...
https://turbmodels.larc.nasa.gov/naca0012_grids.html
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local string = require 'ext.string'
local file = require 'ext.file'
local SolverBase = require 'solver.solverbase'
local vec3 = require 'vec.vec3'


local MeshSolver = class(SolverBase)

--[[
args:
	meshfile = name of mesh file to use

NOTICE initState is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:init(args)
	MeshSolver.super.init(self, args)

	local meshfilename = assert(args.meshfile)
	local data = string.split(string.trim(file['grids/'..meshfilename..'.p2dfmt']), '\n')
		:map(function(l)
			return string.split(string.trim(l), '%s+')
		end)
	local line1 = data:remove(1)
	local size = data:remove(1)
	
	local vtxs = table()
	for i=1,#data do
		vtxs:append(data[i])
	end

	print('size:product()', size:product())
	print('#vtxs', #vtxs)
	print('#size', #size)
	print('#vtxs / #size', #vtxs / #size)
	assert(size:product() == #vtxs / #size)

	self.mesh = {
		vtxs = vtxs,
		-- list of vertex indexes for each cell
		-- how about cell centers too?
		cells = table.append(range(0,size[1]-2):map(function(i)
			return range(0,size[2]-2):map(function(j)
				return table{
					i + size[1] * j,
					i + 1 + size[1] * j,
					i + 1 + size[1] * (j + 1),
					i + size[1] * (j + 1),
				}
			end)
		end):unpack()),
	}

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
		Program.super.init(self, args)
	end
	self.Program = Program
end

function MeshSolver:refreshInitStateProgram()
	-- override this until I know what to do...
	--self.eqn.initState:refreshInitStateProgram(self)
end

function MeshSolver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	local numCells = #self.mesh.cells
	
	-- numCells is the number of cells
	-- maybe I should rename it to numCells
	local localSize1d = math.min(maxWorkGroupSize, numCells)
	
	self.domain = self.app.env:domain{
		size = {localSize1d},
		dim = self.dim,
	}
	
	return {
		numCells = numCells,
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
	}

	lines:insert[[
#define SETBOUNDS(lhs,rhs)	\
	int4 i = globalInt4(); \
	int index = 0;
#define cell_x(i)			_real3(0,0,0)
]]

	self.codePrefix = lines:concat'\n'
end

function MeshSolver:update()
end

return MeshSolver 
