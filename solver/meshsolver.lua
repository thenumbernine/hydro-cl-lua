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
--local vec2sz = require 'vec-ffi.vec2sz'
local vec2sz = require 'vec-ffi.create_vec2'{ctype='size_t'}
local vec2i = require 'vec-ffi.vec2i'
local vec2d = require 'vec-ffi.vec2d'
local vec3d = require 'vec-ffi.vec3d'
local SolverBase = require 'solver.solverbase'


local MeshSolver = class(SolverBase)


-- TODO real3 vs vec3f/vec3d ...
MeshSolver.meshTypeCode = [[
typedef struct iface_s {
	real3 pos;		//center.  realN.
	real3 normal;	//normal pointing from first to second
	real area;		//edge length / surface area
	real cellDist;	//dist between cell centers along 'normal'
	int cells[2];	//indexes of cells
	
	int vtxOffset;
	int vtxCount;
} iface_t;

typedef struct cell_s {
	real3 pos;		//center.  technically could be a 'realN'
	real volume;	//volume of the cell
	int faceOffset;
	int faceCount;
	int vtxOffset;
	int vtxCount;
} cell_t;
]]

-- TODO what if real3 isn't defined yet?
ffi.cdef(MeshSolver.meshTypeCode)

local iface_mt = {}
local cell_mt = {}

local iface_metatype = ffi.metatype('iface_t', iface_mt)
local cell_metatype = ffi.metatype('cell_t', cell_mt)

local function iface_t()
	local iface = iface_metatype()
	iface.pos.x = 0
	iface.pos.y = 0
	iface.pos.z = 0
	iface.normal.x = 0
	iface.normal.y = 0
	iface.normal.z = 0
	iface.area = 0
	iface.cellDist = 0
	iface.cells[0] = -1
	iface.cells[1] = -1
	iface.vtxOffset = 0
	iface.vtxCount = 0
	return iface
end

local function cell_t()
	local cell = cell_metatype()
	cell.pos.x = 0
	cell.pos.y = 0
	cell.pos.z = 0
	cell.faceOffset = 0
	cell.faceCount = 0
	cell.vtxOffset = 0
	cell.vtxCount = 0
	return cell
end

--[=[
--[[
stl vector class
have I made this before?
--]]
local vector = class()
function vector:init(args)
	self.type = args.type
	self.size = 0
	self.capacity = 0
	self:setcapacity(32)
	if args.size > 0 then
		self:resize(args.size)
	end
end
function self:setcapacity(newcap)
	if newcap <= self.capacity then return end
	local newv = ffi.new(self.type..'[?]', newcap)
	for i=0,self.size-1 do	-- TODO ffi.copy
		newv[i] = self.v[i]
	end
	self.v = newv
	self.capacity = newcap
end
function self:resize(newsize)
	self:setcapacity( math.ceil(newsize / 32) * 32 )
	self.size = newsize
end
function vector:push_back(obj)
	self:resize(self.size + 1)
	self.v[self.size - 1] = obj
end
function vector:back()
	assert(self.size > 0)
	return self.v + self.size - 1
end
function vector:insertAll(...)
	if select('#', ...) == 3 then
		local where, first, last = ...
		local i = first
		while i ~= last do
			
			i = i + 1
		end
	end
end
function vector:__len()
	return rawget(self, 'size')
end
function vector:__index(k)
	assert(k >= 0 and k < #self)
	return rawget(self, 'v')[k]
end
function vector:__newindex(k, v)
	assert(k >= 0 and k < #self)
	rawget(self, 'v')[k] = v
end
--]=]


-- fields we need:
-- real3[] vtxs
-- iface_t[] faces
-- cell_t[] cells
-- int[] cellFaceIndexes
-- int[] cellVtxIndexes
-- int[] faceVtxIndexes
local Mesh = class()

function Mesh:init(solver)
	self.solver = solver
	self.vtxs = table()	--vector{type='real3'}
	self.faces = table()	--vector{type='iface_t'}
	self.cells = table()	--vector{type='cell_t'}
	self.cellFaceIndexes = table()	--vector{type='int'}
	self.cellVtxIndexes = table()	--vector{type='int'}
	self.faceVtxIndexes = table()	--vector{type='int'}
	self.real3 = function(x,y,z)
		return ffi.new('real3', {x=x, y=y, z=z})
	end
end


--	3D - polygon


function Mesh:polygon3DVolume(vs, normal)
	local n = #vs
	local volume = 0
	for i=1,n do
		local a = vs[i]
		local b = vs[i%n+1]
		volume = volume + self:parallelepipedVolume(a, b, normal)
	end
	return .5 * volume
end

function Mesh:polygon3DCOM(vs, area, normal)
	local n = #vs
	if area == 0 then
		if n == 0 then
			error"you can't get a COM without any vertexes"
		end
		return table.sum(vs) / n
	end
	local com = vec3d()
	local a = vs[1]
	for i=3,n do
		local b = vs[i-1]
		local c = vs[i]
		com = com + (a + b + c) * normal:dot(cross(c - a, c - b))
	end
	return com * (1 / (6 * area))
end


--	2D - parallelogram


function Mesh:parallelogramVolume(a, b)
	--epsilon_ij a^i b^j
	return a.x * b.y - a.y * b.x
end


--	2D - polygon


function Mesh:polygonVolume(vs)
	local n = #vs
	local volume = 0
	for i=1,n do
		local a = vs[i]
		local b = vs[i%n+1]
		volume = volume + self:parallelogramVolume(a, b)
	end
	return .5 * volume
end

function Mesh:polygonCOM(vs, volume)
	local n = #vs
	if volume == 0 then
		if n == 0 then
			error"you can't get a COM without any vertexes"
		end
		return table.sum(vs) / n
	end
	local com = vec2d()
	for i=1,n do
		local a = vs[i]
		local b = vs[i%n+1]
		local vol = self:parallelogramVolume(a, b)
		for j=0,1 do
			com.s[j] = com.s[j] + (a.s[j] + b.s[j]) * vol
		end
	end
	return com * (1 / (6 * volume))
end



--	3D - parallelepiped


function Mesh:parallelepipedVolume(a, b, c)
	--epsilon_ijk a^i b^j c^k
	return a.x * b.y * c.z
		+ a.y * b.z * c.x
		+ a.z * b.x * c.y
		- c.x * b.y * a.z
		- c.y * b.z * a.x
		- c.z * b.x * a.y
end


--	3D - polyhedron


function Mesh:polyhedronVolume(faces)
	local volume = 0
	for _,face in ipairs(faces) do
		for i=3,#faces do
			--tri from 0, i-1, 1
			local a = face[1]
			local b = face[i-1]
			local c = face[i]
			volume = volume + self:parallelepipedVolume(a, b, c)
		end
	end
	--volume of n-sided pyramid in nD is volume of parallelogram divided by 1/n!
	return volume / 6
end

function Mesh:polyhedronCOM(faces, volume)
	if volume == 0 then
		if #faces == 0 then
			error"you can't get a COM without any vertexes"
		end
		local sum = self.real3()
		local total = 0
		for _,face in ipairs(faces) do
			for _,vtx in ipairs(face) do
				sum = sum + vtx
				total = total + 1
			end
		end
		return sum * (1 / total)
	end
	local com = self.real3()
	for _,face in ipairs(faces) do
		for i=3,#face do
			--tri from 0, i-1, 1
			local a = face[1]
			local b = face[i-1]
			local c = face[i]
			com = com + ((a + b) * (a + b) + (b + c) * (b + c) + (c + a) * (c + a)) * (b - a):cross(c - a) / 48
		end
	end
	return com / volume
end



-- vs is a Lua table
-- vs values are 0-based indexes
-- returns a 0-based index
function Mesh:addFaceForVtxs(vs, n)
	for fi=0,#self.faces-1 do
		local f = assert(self.faces[1+fi], "failed to get face at 0-based "..fi)
		if self.solver.dim == 2 then
			local va = vs[1]
			local vb = vs[2]
			if f.vtxCount == 2 then
				if (self.faceVtxIndexes[1+f.vtxOffset+0] == va and self.faceVtxIndexes[1+f.vtxOffset+1] == vb) or
					(self.faceVtxIndexes[1+f.vtxOffset+0] == vb and self.faceVtxIndexes[1+f.vtxOffset+1] == va)
				then
					return fi
				end
			end
		elseif self.solver.dim == 3 then
			for j=0,n-1 do
				--check in one direction
				local matches = true
				for i=0,n-1 do
					if self.faceVtxIndexes[1+f.vtxOffset+i] ~= vs[1+(j+i)%n] then
						matches = false;
						break
					end
				end
				if matches then return fi end
			
				--check in the other direction
				matches = true
				for i=0,n-1 do
					if self.faceVtxIndexes[1+f.vtxOffset+i] ~= vs[1+(j+n-i)%n] then
						matches = false
						break
					end
				end
				if matches then return fi end
			end
		end
	end

	local fi = #self.faces
	
	self.faces:insert(iface_t())
	local f = self.faces:last()
	
	f.vtxOffset = #self.faceVtxIndexes
	for i=1,n do
		self.faceVtxIndexes:insert(vs[i])
	end
	f.vtxCount = #self.faceVtxIndexes - f.vtxOffset

	--TODO calc these for n=3
	if self.solver.dim == 2 then
		local a = self.vtxs[1+vs[1]]
		local b = self.vtxs[1+vs[2]]
		for i=0,2 do
			f.pos.s[i] = (a.s[i] + b.s[i]) * .5
		end
		local delta = ffi.new'real3'
		for i=0,2 do
			delta.s[i] = a.s[i] - b.s[i]
		end
		local deltaLen = math.sqrt(delta.x*delta.x + delta.y*delta.y)
		f.area = deltaLen
		f.normal = self.real3(delta.y / deltaLen, -delta.x / deltaLen)
	elseif self.solver.dim == 3 then
		local polyVtxs = table()	--vector{type='real3', size=n}
		for i=1,n do
			polyVtxs[i] = self.vtxs[1+vs[i]]
		end
		for i=0,n-1 do
			local i2 = (i+1)%n
			local i3 = (i2+1)%n
			f.normal = f.normal + (polyVtxs[1+i2] - polyVtxs[1+i]):cross(polyVtxs[1+i3] - polyVtxs[1+i2])
		end
		f.normal = f.normal:normalize()
		f.area = self:polygon3DVolume(polyVtxs, f.normal)
		f.pos = self:polygon3DCOM(polyVtxs, f.area, f.normal)
	else
		error'here'
	end
	
	return fi
end

-- vs is a Lua table
-- vs values are 0-based indexes
function Mesh:addFace(vs, n, ci)
	local fi = self:addFaceForVtxs(vs, n)
--	if fi ~= -1 then return fi end
	local f = self.faces[fi+1]
	if f.cells[0] == -1 then
		f.cells[0] = ci
	elseif f.cells[1] == -1 then
		f.cells[1] = ci
	else
		error"tried to add too many cells to an edge"
	end
	return fi
end

-- vis is a Lua table of vertex indexes
-- vis values are 0-based indexes
function Mesh:addCell(vis)
	local ci = #self.cells
	self.cells:insert(cell_t())
	local c = self.cells:last()
	
	local n = #vis

	--c.vtxs = vis
	c.vtxOffset = #self.cellVtxIndexes
	self.cellVtxIndexes:append(vis)
	c.vtxCount = #self.cellVtxIndexes - c.vtxOffset

	local lastFaceSize = #self.faces
	c.faceOffset = #self.cellFaceIndexes
	if self.solver.dim == 2 then
		--face is a 1-form
		for i=1,n do
			local vtxs = {vis[i], vis[i%n+1]}
			local fi = self:addFace(vtxs, #vtxs, ci)
			do -- if (fi != -1) then
				self.cellFaceIndexes:insert(fi)
			end
		end
		local polyVtxs = table.mapi(vis, function(vi)
			return self.vtxs[1+vi]
		end)
		
		c.volume = self:polygonVolume(polyVtxs)
		local com = self:polygonCOM(polyVtxs, c.volume)
		c.pos = self.real3(com:unpack())
	
	elseif self.solver.dim == 3 then
		--face is a 2-form
		assert(#vis == 8)	--only adding cubes at the moment
--[[
6----7
/|   /|
4----5 |
| |  | |  z
| 2--|-3  ^ y
|/   |/   |/
0----1    *-> x
--]]
		--identity cube sides
		local identityCubeSides = {
			{0,4,6,2},	--x-
			{1,3,7,5},	--x+
			{0,1,5,4},	--y-
			{2,6,7,3},	--y+
			{0,2,3,1},	--z-
			{4,5,7,6},	--z+
		}

		--vector of per-face vector-of-vertexes
		--to pass to the polyhedron functions
		local cubeVtxs = table()

		for _,side in ipairs(identityCubeSides) do
			
			local thisFaceVtxIndexes = table.mapi(side, function(side_i)
				return vis[1+side_i]
			end)
			
			local fi = self:addFace(self.thisFaceVtxIndexes, #self.thisFaceVtxIndexes, ci)
			do -- if fi ~= -1 then
				local f = self.faces[fi+1]
				
				--if face area is zero then don't add it to cell's faces
				-- and don't use it later for determining cell's volume
				if f.area <= 1e-7 then
				else
					self.cellFaceIndexes:insert(fi)

					cubeVtxs:insert(
						table.mapi(thisFaceVtxIndexes, function(i)
							return vtxs[i]
						end)
					)
				end
			end
		end

		if #self.cellFaceIndexes - lastFaceSize < 4 then	--not enough sides to even form a tetrahedron
			c.volume = 0
		else
			c.volume = self:polyhedronVolume(cubeVtxs)
			c.pos = self:polyhedronCOM(cubeVtxs, c.volume)
		end
	else
		error"you are here";
	end
	
	c.faceCount = #self.cellFaceIndexes - c.faceOffset
end

function Mesh:calcAux()
	for _,f in ipairs(self.faces) do
		local a = f.cells[0]
		local b = f.cells[1]
		local cella = self.cells[1+a]
		local cellb = self.cells[1+b]
		if a ~= -1 and b ~= -1 then
			local dx = cella.pos.x - cellb.pos.x
			local dy = cella.pos.y - cellb.pos.y
			local dz = cella.pos.z - cellb.pos.z
			local nDotDelta = f.normal.x * dx + f.normal.y * dy + f.normal.z * dz
			if nDotDelta < 0 then
				a,b = b,a
				f.cells[0] = a
				f.cells[1] = b
				f.normal.x = -f.normal.x
				f.normal.y = -f.normal.y
				f.normal.z = -f.normal.z
				
				cella = self.cells[1+a]
				cellb = self.cells[1+b]
				dx = -dx
				dy = -dy
				dz = -dz
				nDotDelta = -nDotDelta 
			end
			--distance between cell centers
			--f.cellDist = (self.cells[1+b].pos - self.cells[1+a].pos).length();
			--distance projected to edge normal
			f.cellDist = math.abs(nDotDelta)
		elseif a ~= -1 then
			local dx = cella.pos.x - f.pos.x
			local dy = cella.pos.y - f.pos.y
			local dz = cella.pos.z - f.pos.z
			f.cellDist = math.sqrt(dx*dx + dy*dy + dz*dz) * 2
		else
			error"looks like you created a face that isn't touching any cells..."
		end
	end
end


local MeshFactory = class()


local P2DFMTMeshFactory = class(MeshFactory)

function P2DFMTMeshFactory:createMesh(args)
	-- TODO
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
end


local Chart2DMeshFactory = class(MeshFactory)
function Chart2DMeshFactory:init(args)
	args = args or {}
	self.size = vec2sz(args.size or {20,20})
	self.mins = vec2d(args.mins or {-1, -1})
	self.maxs = vec2d(args.maxs or {1, 1})
	self.wrap = vec2i(args.wrap or {0, 0})
	self.capmin = vec2i(args.capmin or {0, 0})
end
function Chart2DMeshFactory:coordChart(x) return x end


local Tri2DMeshFactory = class(Chart2DMeshFactory)
Tri2DMeshFactory.name = 'Tri2DMesh'
function Tri2DMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)

	local n = self.size + 1
	local step = vec2i(1, n.x)
	
	local vtxsize = n:volume()
	if self.capmin.x then vtxsize = vtxsize + 1 end
	
	--mesh.vtxs = ffi.new('real3[?]', vtxsize)
	local i = vec2i(0,0)
	for iy=0,n.y-1 do
		i.y = iy
		for ix=0,n.x-1 do
			i.x = ix
			local x = vec2d(vec2d(i) / vec2d(self.size) * (self.maxs - self.mins) + self.mins)
			local u = self:coordChart(x)
			mesh.vtxs[tonumber(1+i:dot(step))] = self.real3(u:unpack())
		end
	end
	
	local capindex = n:volume()
	if self.capmin.x then
		local sum = mesh.real3()
		for j=0,n.y-1 do
			sum = sum + mesh.vtxs[tonumber(1+ 0 + n.x * j)]
		end
		mesh.vtxs[tonumber(1+capindex)].pos = sum / tonumber(n.y)
	end
	
	local imax = vec2i()
	for j=0,1 do
		imax.s[j] = self.wrap.s[j] ~= 0 and n.s[j] or n.s[j] - 1
	end
	
	local ni = vec2i()
	for iy=0,imax.y-1 do
		i.y = iy
		ni.y = (i.y + 1) % n.y
		for ix=0,imax.x-1 do
			ni.x = (i.x + 1) % n.x
			mesh:addCell{
				tonumber(i.x + n.x * i.y),
				tonumber(ni.x + n.x * i.y),
				tonumber(ni.x + n.x * ni.y),
			}
			mesh:addCell{
				tonumber(ni.x + n.x * ni.y),
				tonumber(i.x + n.x * ni.y),
				tonumber(i.x + n.x * i.y),
			}
		end
	end
	
	mesh:calcAux()
	return mesh;
end


local Quad2DMeshFactory = class(Chart2DMeshFactory)
Quad2DMeshFactory.name = 'Quad2DMesh'
function Quad2DMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)

	local n = self.size + 1
	local step = vec2i(1, n.x)
	local vtxsize = n:volume()
	if self.capmin.x ~= 0 then vtxsize = vtxsize + 1 end
	--mesh.vtxs = ffi.new('real3[?]', vtxsize)

	local coordRangeMax = vec2i(self.size:unpack())
	if self.wrap.x ~= 0 or self.capmin.x ~= 0 then coordRangeMax.x = coordRangeMax.x + 1 end
	if self.wrap.y ~= 0 or self.capmin.y ~= 0 then coordRangeMax.y = coordRangeMax.y + 1 end

	local iofs = vec2i()
	if self.capmin.x ~= 0 then iofs.x = 1 end
	if self.capmin.y ~= 0 then iofs.y = 1 end

	local i = vec2i()
	for iy=0,tonumber(n.y)-1 do
		i.y = iy
		for ix=0,tonumber(n.x)-1 do
			i.x = ix
			local x = vec2d((i + iofs):unpack()) / vec2d(coordRangeMax:unpack()) * (self.maxs - self.mins) + self.mins
			local u = self:coordChart(x)
			mesh.vtxs[tonumber(1+i:dot(step))] = mesh.real3(u:unpack())
		end
	end
	
	local capindex = n:volume()
	if self.capmin.x ~= 0 then
		local sum = mesh.real3()
		for j=0,n.y-1 do
			sum = sum + mesh.vtxs[tonumber(1+ 0 + n.x * j)]
		end
		mesh.vtxs[tonumber(1+capindex)] = sum / tonumber(n.y)
	end

	local imax = vec2i()
	for j=0,1 do
		imax.s[j] = self.wrap.s[j] ~= 0 and n.s[j] or n.s[j]-1
	end

	local ni = vec2i()
	for iy=0,tonumber(imax.y-1) do
		i.y = iy
		ni.y = (i.y + 1) % n.y
		for ix=0,tonumber(imax.x-1) do
			i.x = ix
			ni.x = (i.x + 1) % n.x
			mesh:addCell{
				tonumber(i.x + n.x * i.y),
				tonumber(ni.x + n.x * i.y),
				tonumber(ni.x + n.x * ni.y),
				tonumber(i.x + n.x * ni.y),
			}
		end
	end

	if self.capmin.x ~= 0 then
		for j=0,imax.y-1 do
			local jn = (j + 1) % n.y
			mesh.addCell{
				tonumber(0 + n.x * j), 
				tonumber(0 + n.x * jn), 
				tonumber(capindex),
			}
		end
	end

	mesh:calcAux()
	return mesh
end


local meshFactoryClassForType = {
	Tri2DMesh = Tri2DMeshFactory,
	Quad2DMesh = Quad2DMeshFactory,
}


--[[
args:
	meshfile = name of mesh file to use

NOTICE initState is tied closely to grid mins/maxs...
so how should meshfiles use init states?
--]]
function MeshSolver:initL1(args)
	MeshSolver.super.initL1(self, args)


	local meshType = assert(args.mesh.type, "expected mesh type")
	local meshFactoryClass = assert(meshFactoryClassForType[meshType], "failed to find mesh factory for type "..meshType)
	local meshFactory = meshFactoryClass()
	self.mesh = meshFactory:createMesh(self)


	-- ok now convert from lua tables to arrays
	-- or do this in mesh:calcAux() ?


	self.mins = vec3d(-1, -1, -1)
	self.maxs = vec3d(1, 1, 1)
	for i=0,self.dim-1 do
		self.mins.s[i] = math.huge
		self.maxs.s[i] = -math.huge
	end
	for _,v in ipairs(self.mesh.vtxs) do
		for i=0,self.dim-1 do
			self.mins.s[i] = math.min(self.mins.s[i], v.s[i])
			self.maxs.s[i] = math.max(self.maxs.s[i], v.s[i])
		end
	end

	-- convert vertices
	-- only store position
	local cpu_vtxs = ffi.new('real3[?]', #self.mesh.vtxs)
	for i,v in ipairs(self.mesh.vtxs) do
		cpu_vtxs[i-1] = v
	end

	local cpu_ifaces = ffi.new('iface_t[?]', #self.mesh.faces)
	for i,iface in ipairs(self.mesh.faces) do
		cpu_ifaces[i-1] = iface
	end

	local cpu_cells = ffi.new('cell_t[?]', #self.mesh.cells)
	for i,cell in ipairs(self.mesh.cells) do
		cpu_cells[i-1] = cell
	end

	-- TODO put these in solver_t
	self.numCells = #self.mesh.cells
	self.numInterfaces = #self.mesh.faces
	self.numVtxs = #self.mesh.vtxs
	
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

	self.cmds:enqueueWriteBuffer{buffer=self.cellsBuf, block=true, size=ffi.sizeof'cell_t' * self.numCells, ptr=self.mesh.cells}
	self.cmds:enqueueWriteBuffer{buffer=self.ifacesBuf, block=true, size=ffi.sizeof'iface_t' * self.numInterfaces, ptr=self.mesh.faces}
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
	self.calcDTKernelObj.obj:setArg(2, self.cellsBuf)
	self.calcDTKernelObj.obj:setArg(3, self.ifacesBuf)
end

function MeshSolver:update()
end

function MeshSolver:boundary()
end


-- hmm, if something else require's solverbase and uses SolverBase.DisplayVar before this does 
-- then won't it get the wrong class? (compared to if this require's first before it does?)
local DisplayVar = MeshSolver.DisplayVar
local MeshSolverDisplayVar = class(DisplayVar)
MeshSolver.DisplayVar = MeshSolverDisplayVar

function MeshSolverDisplayVar:setArgs(kernel)
	MeshSolverDisplayVar.super.setArgs(self, kernel)
	kernel:setArg(4, self.solver.cellsBuf)
	kernel:setArg(5, self.solver.ifacesBuf)
end


return MeshSolver 
