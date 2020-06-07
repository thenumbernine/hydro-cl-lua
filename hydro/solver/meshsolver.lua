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
local vec3f = require 'vec-ffi.vec3f'
local vec3d = require 'vec-ffi.vec3d'
local matrix_ffi = require 'matrix.ffi'
local ig = require 'ffi.imgui'
local gl = require 'gl'
local glreport = require 'gl.report'
local GLProgram = require 'gl.program'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLVertexArray = require 'gl.vertexarray'
local template = require 'template'
local tooltip = require 'hydro.tooltip'
local SolverBase = require 'hydro.solver.solverbase'
local time, getTime = table.unpack(require 'hydro.util.time')
local real = require 'hydro.real'
local struct = require 'hydro.struct.struct'

matrix_ffi.real = 'float'	-- default matrix_ffi type

-- [=[
--[[
stl vector class
have I made this before?
--]]
local vector = class()
function vector:init(ctype, arg)
	self.type = ctype
	self.size = 0
	self.capacity = 0
	self:setcapacity(32)
	if arg ~= nil then
		if type(arg) == 'number' then
			self:resize(arg)
		elseif type(arg) == 'table' then
			self:resize(#arg)
			for i,v in ipairs(arg) do
				self.v[i-1] = v
			end
		end
	end

	-- I don't want to set __index to vector or else things will go bad within class init ... so set it here:
	--[[
	local mt = table(vector)
	function mt:__index(k)
		local mv = mt[k] 
		if mv ~= nil then return mv end
		if k < 0 or k > self.size then
			error("got out of bounds index: "..tostring(k))
		end
		return self.v[k]
	end
	function mt:__newindex(k, v)
		assert(k >= 0 and k < #self)
		rawget(self, 'v')[k] = v
	end
	setmetatable(self, mt)
	--]]
end
function vector:__ipairs()
	return coroutine.wrap(function()
		for i=0,self.size-1 do
			coroutine.yield(i, self.v[i])
		end
	end)
end
function vector:setcapacity(newcap)
	if newcap <= self.capacity then return end
	local newv = ffi.new(self.type..'[?]', newcap)
	for i=0,self.size-1 do	-- TODO ffi.copy
		newv[i] = self.v[i]
	end
	self.v = newv
	self.capacity = newcap
end

function vector:resize(newsize)
	newsize = assert(tonumber(newsize))
	self:setcapacity( (math.floor(tonumber(newsize) / 32) + 1) * 32 )
	self.size = newsize
end

function vector:push_back(obj)
	self:resize(self.size + 1)
	self.v[self.size - 1] = obj
end

-- returns a ptr to the last element
function vector:back()
	assert(self.size > 0)
	return self.v + self.size - 1
end
-- returns a ptr to the first element
function vector:begin()
	return self.v
end
-- because the __index removes access from the object itself ...
vector.data = vector.begin
-- returns a ptr past the last element
function vector:iend()
	return self.v + self.size
end

--[[
insert(where, ptr first, ptr last)
insert(where, ref value)
--]]
function vector:insert(...)
	local n = select('#', ...)
	if n == 3 then
		local where, first, last = ...
		local i = first
		while i ~= last do
			self:insert(where, i[0])
			i = i + 1
		end
	elseif n == 2 then
		local where, value = ...
		local offset = where - self.v
		assert(0 <= offset and offset <= self.size)
		self:resize(self.size + 1)
		for i=self.size-1,offset+1,-1 do
			self.v[i] = self.v[i-1]
		end
		self.v[offset] = value
	end
end
function vector:__len()
	--return rawget(self, 'size')
	return self.size
end
function vector:__tostring()
	local s = '['
	local sep = ''
	for i=0,self.size-1 do
		s = s .. sep .. self.v[i]
		sep = ', '
	end
	return s .. ']'
end

--[[
local v = vector'int'
v:push_back(10)
print(v)
print(#v)
print(v[0])
v:push_back(20)
v:push_back(30)
local ip = ffi.new('int[3]', 1,2,3)
v:insert(v.v + 1, ip, ip + 3)
print(#v)
print(v)
os.exit()
--]]
--]=]


local MeshSolver = class(SolverBase)

local face_t = struct{
	name = 'face_t',
	dontUnion = true,
	vars = {
		{type='real3', name='pos'},		--center.  realN.
		{type='real3', name='normal'},	--normal pointing from first to second
		{type='real3', name='normal2'},	--orthonormal basis around normal
		{type='real3', name='normal3'},
		{type='real', name='area'},		--edge length / surface area
		{type='real', name='cellDist'},	--dist between cell centers along 'normal'
		{type='vec2i_t', name='cells'},	--indexes of cells
		{type='int', name='vtxOffset'},
		{type='int', name='vtxCount'},
	},
}

local cell_t = struct{
	name = 'cell_t',
	dontUnion = true,
	vars = {
		{type='real3', name='pos'},		--center.  technically could be a 'realN'
		{type='real', name='volume'},	--volume of the cell
		{type='int', name='faceOffset'},
		{type='int', name='faceCount'},
		{type='int', name='vtxOffset'},
		{type='int', name='vtxCount'},
	},
}

-- TODO real3 vs vec3f/vec3d ...
MeshSolver.meshTypeCode = table{
	vec2i.typeCode,
	face_t:getTypeCode(),
	cell_t:getTypeCode()
}:concat'\n'

-- TODO what if real3 isn't defined yet?
ffi.cdef(MeshSolver.meshTypeCode)

local face_mt = {}
local cell_mt = {}

local function new_face_t()
	local face = face_t.metatype()
	face.pos.x = 0
	face.pos.y = 0
	face.pos.z = 0
	face.normal.x = 0
	face.normal.y = 0
	face.normal.z = 0
	face.normal2.x = 0
	face.normal2.y = 0
	face.normal2.z = 0
	face.normal3.x = 0
	face.normal3.y = 0
	face.normal3.z = 0
	face.area = 0
	face.cellDist = 0
	face.cells.s[0] = -1
	face.cells.s[1] = -1
	face.vtxOffset = 0
	face.vtxCount = 0
	return face
end

local function new_cell_t()
	local cell = cell_t.metatype()
	cell.pos.x = 0
	cell.pos.y = 0
	cell.pos.z = 0
	cell.faceOffset = 0
	cell.faceCount = 0
	cell.vtxOffset = 0
	cell.vtxCount = 0
	return cell
end


local Mesh = class()

function Mesh:init(solver)
	self.solver = solver
	self.vtxs = vector'real3'
	self.faces = vector'face_t'
	self.cells = vector'cell_t'
	self.cellFaceIndexes = vector'int'
	self.cellVtxIndexes = vector'int'
	self.faceVtxIndexes = vector'int'
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
		local f = assert(self.faces.v[fi], "failed to get face at 0-based "..fi)
		if self.solver.dim == 2 then
			local va = vs[1]
			local vb = vs[2]
			if f.vtxCount == 2 then
				if (self.faceVtxIndexes.v[f.vtxOffset+0] == va and self.faceVtxIndexes.v[f.vtxOffset+1] == vb) or
					(self.faceVtxIndexes.v[f.vtxOffset+0] == vb and self.faceVtxIndexes.v[f.vtxOffset+1] == va)
				then
					return fi
				end
			end
		elseif self.solver.dim == 3 then
			for j=0,n-1 do
				--check in one direction
				local matches = true
				for i=0,n-1 do
					if self.faceVtxIndexes.v[f.vtxOffset+i] ~= vs[1+(j+i)%n] then
						matches = false;
						break
					end
				end
				if matches then return fi end
			
				--check in the other direction
				matches = true
				for i=0,n-1 do
					if self.faceVtxIndexes.v[f.vtxOffset+i] ~= vs[1+(j+n-i)%n] then
						matches = false
						break
					end
				end
				if matches then return fi end
			end
		end
	end

	local fi = #self.faces
	
	self.faces:push_back(new_face_t())
	local f = self.faces:back()
	
	f.vtxOffset = #self.faceVtxIndexes
	for i=1,n do
		self.faceVtxIndexes:push_back(vs[i])
	end
	f.vtxCount = #self.faceVtxIndexes - f.vtxOffset

	--TODO calc these for n=3
	if self.solver.dim == 2 then
		local a = self.vtxs.v[vs[1]]
		local b = self.vtxs.v[vs[2]]
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
		local polyVtxs = table()	--vector('real3', n)
		for i=1,n do
			polyVtxs[i] = self.vtxs.v[vs[i]]
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
	local f = self.faces.v[fi]
	if f.cells.s[0] == -1 then
		f.cells.s[0] = ci
	elseif f.cells.s[1] == -1 then
		f.cells.s[1] = ci
	else
		error"tried to add too many cells to an edge"
	end
	return fi
end

-- vis is a Lua table of vertex indexes
-- vis values are 0-based indexes
function Mesh:addCell(vis)
	local ci = #self.cells
	self.cells:push_back(new_cell_t())
	local c = self.cells:back()
	
	local n = #vis

	--c.vtxs = vis
	c.vtxOffset = #self.cellVtxIndexes
	self.cellVtxIndexes:insert(self.cellVtxIndexes:iend(), vis:begin(), vis:iend())
	c.vtxCount = #self.cellVtxIndexes - c.vtxOffset

	local lastFaceSize = #self.faces
	c.faceOffset = #self.cellFaceIndexes
	if self.solver.dim == 2 then
		--face is a 1-form
		for i=1,n do
			local vtxs = {vis.v[i-1], vis.v[i%n]}
			local fi = self:addFace(vtxs, #vtxs, ci)
			do -- if (fi != -1) then
				self.cellFaceIndexes:push_back(fi)
			end
		end
		local polyVtxs = range(0,#vis-1):mapi(function(i)
			return self.vtxs.v[vis.v[i]]
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
				return vis.v[side_i]
			end)
			
			local fi = self:addFace(self.thisFaceVtxIndexes, #self.thisFaceVtxIndexes, ci)
			do -- if fi ~= -1 then
				local f = self.faces.v[fi]
				
				--if face area is zero then don't add it to cell's faces
				-- and don't use it later for determining cell's volume
				if f.area <= 1e-7 then
				else
					self.cellFaceIndexes:push_back(fi)

					cubeVtxs:insert(
						table.mapi(thisFaceVtxIndexes, function(i)
							return vtxs.v[i]
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
		local a = f.cells.s[0]
		local b = f.cells.s[1]
		if a ~= -1 and b ~= -1 then
			local cella = self.cells.v[a]
			local cellb = self.cells.v[b]
			local dx = cella.pos.x - cellb.pos.x
			local dy = cella.pos.y - cellb.pos.y
			local dz = cella.pos.z - cellb.pos.z
			local nDotDelta = f.normal.x * dx + f.normal.y * dy + f.normal.z * dz
			if nDotDelta < 0 then
				a,b = b,a
				f.cells.s[0] = a
				f.cells.s[1] = b
				f.normal.x = -f.normal.x
				f.normal.y = -f.normal.y
				f.normal.z = -f.normal.z
				
				cella = self.cells.v[a]
				cellb = self.cells.v[b]
				dx = -dx
				dy = -dy
				dz = -dz
				nDotDelta = -nDotDelta 
			end
			--distance between cell centers
			--f.cellDist = (self.cells.v[b].pos - self.cells.v[a].pos).length();
			--distance projected to edge normal
			f.cellDist = math.abs(nDotDelta)
		elseif a ~= -1 then
			local cella = self.cells.v[a]
			local dx = cella.pos.x - f.pos.x
			local dy = cella.pos.y - f.pos.y
			local dz = cella.pos.z - f.pos.z
			f.cellDist = math.sqrt(dx*dx + dy*dy + dz*dz) * 2
		else
			error"looks like you created a face that isn't touching any cells..."
		end

		-- here form an orthonormal basis
		local xcn = self.real3(0, f.normal.z, -f.normal.y)
		local ycn = self.real3(-f.normal.z, 0, f.normal.x)
		local zcn = self.real3(f.normal.y, -f.normal.x, 0)
		local xlenSq = xcn:lenSq()
		local ylenSq = ycn:lenSq()
		local zlenSq = zcn:lenSq()
		if xlenSq > ylenSq then
			if xlenSq > zlenSq then	-- x
				f.normal2:set(xcn:normalize())
			else					-- z 
				f.normal2:set(zcn:normalize())
			end
		elseif ylenSq > zlenSq then	-- y
			f.normal2:set(ycn:normalize())
		else						-- z
			f.normal2:set(zcn:normalize())
		end
		f.normal3:set(f.normal2:cross(f.normal):normalize())
	end

	-- now convert all our buffers to C types
	-- hmm, this is where a std::vector would be helpful ... store the size and ptr in one structure
	-- also, how will I pass all this to the GPU?
	-- do we want so many different parameters?
	-- maybe make some weird data structure where we have to skip past arrays?  
	-- maybe have a initial structure that tells us where to skip to for each of these?
	-- how about shared memory (CL 2.0)?  will that be compat with NVIDIA?
	local function tableToC(t, ty)
		local ptr = ffi.new(ty..'[?]', #t)
		for i=1,#t do
			ptr[i-1] = t[i]
		end
		return ptr, #t
	end
	
	--self.numVtxs, self.vtxs = tableToC(self.vtxs, 'real3')
	self.numVtxs = #self.vtxs
	self.numFaces = #self.faces
	self.numCells = #self.cells
	self.numCellFaceIndexes = #self.cellFaceIndexes
	self.numCellVtxIndexes = #self.cellVtxIndexes
	self.numFaceVtxIndexes = #self.faceVtxIndexes
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

	local nx = tonumber(self.size.x) + 1
	local ny = tonumber(self.size.y) + 1
	
	local stepx = 1
	local stepy = nx
	
	local vtxsize = nx * ny
	if self.capmin.x then vtxsize = vtxsize + 1 end

	mesh.vtxs:resize(vtxsize)
	for iy=0,ny-1 do
		for ix=0,nx-1 do
			local x = vec2d(
				tonumber(ix / self.size.x * (self.maxs.x - self.mins.x) + self.mins.x),
				tonumber(iy / self.size.y * (self.maxs.y - self.mins.y) + self.mins.y))
			local u = self:coordChart(x)
			mesh.vtxs.v[ix * stepx + iy * stepy] = self.real3(u:unpack())
		end
	end
	
	local capindex = nx * ny
	if self.capmin.x then
		local sum = mesh.real3()
		for j=0,ny-1 do
			sum = sum + mesh.vtxs.v[0 + nx * j]
		end
		mesh.vtxs.v[capindex].pos = sum / ny
	end
	
	local imaxx = self.wrap.x ~= 0 and nx or nx - 1
	local imaxy = self.wrap.y ~= 0 and ny or ny - 1
	
	for iy=0,imaxy-1 do
		local niy = (iy + 1) % ny
		for ix=0,imaxx-1 do
			local nix = (ix + 1) % nx
			mesh:addCell(vector('int',{
				ix + nx * iy,
				nix + nx * iy,
				nix + nx * niy,
			}))
			mesh:addCell(vector('int',{
				nix + nx * niy,
				ix + nx * niy,
				ix + nx * iy,
			}))
		end
	end
	
	mesh:calcAux()
	return mesh
end


local Quad2DMeshFactory = class(Chart2DMeshFactory)
Quad2DMeshFactory.name = 'Quad2DMesh'
function Quad2DMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)

	local n = self.size + 1
	local step = vec2i(1, n.x)
	local vtxsize = n:volume()
	if self.capmin.x ~= 0 then vtxsize = vtxsize + 1 end
	mesh.vtxs:resize(vtxsize)

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
			mesh.vtxs.v[tonumber(i:dot(step))] = mesh.real3(u:unpack())
		end
	end
	
	local capindex = n:volume()
	if self.capmin.x ~= 0 then
		local sum = mesh.real3()
		for j=0,n.y-1 do
			sum = sum + mesh.vtxs.v[tonumber(0 + n.x * j)]
		end
		mesh.vtxs.v[tonumber(capindex)] = sum / tonumber(n.y)
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
			mesh:addCell(vector('int',{
				tonumber(i.x + n.x * i.y),
				tonumber(ni.x + n.x * i.y),
				tonumber(ni.x + n.x * ni.y),
				tonumber(i.x + n.x * ni.y),
			}))
		end
	end

	if self.capmin.x ~= 0 then
		for j=0,imax.y-1 do
			local jn = (j + 1) % n.y
			mesh.addCell(vector('int',{
				tonumber(0 + n.x * j), 
				tonumber(0 + n.x * jn), 
				tonumber(capindex),
			}))
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
	self.showVertexes = false
	self.showFaces = false
	self.showCells = true
	self.drawCellScale = 1

	MeshSolver.super.initL1(self, args)


	local meshType = assert(args.mesh.type, "expected mesh type")
	local meshFactoryClass = assert(meshFactoryClassForType[meshType], "failed to find mesh factory for type "..meshType)
	local meshFactory = meshFactoryClass(args.mesh)

--[[ script for testing mesh build times:
require 'ext.range'(1,6,.25):mapi(function(x)
	return true, math.floor(2^x)
end):keys():sort():mapi(function(i)
	local cmd = './run.lua sys=console meshsize='..i
	print('>'..cmd)
	print(os.execute(cmd))
end)
--]]
	time('creating mesh', function()
		self.mesh = meshFactory:createMesh(self)
	end)
--os.exit()

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

	-- TODO put these in solver_t
	self.numCells = self.mesh.numCells
	self.numFaces = self.mesh.numFaces
	self.numVtxs = self.mesh.numVtxs
	self.numCellFaceIndexes = self.mesh.numCellFaceIndexes
	self.numCellVtxIndexes = self.mesh.numCellVtxIndexes
	self.numFaceVtxIndexes = self.mesh.numFaceVtxIndexes
	
	--[[
	next issue ... how to render this?
	how to store it?
	where do textures come into play?
	
	UBuf will be size of # of elements
	
	how about texture, and how about rendering?
	texsize will have to be > #elems ... using size is fine as long as we need a gridsize
	--]]


	if self.app.targetSystem ~= 'console' then


	self.drawShader = GLProgram{
		vertexCode = [[
#version 460
uniform float drawCellScale;
uniform mat4 modelViewProjectionMatrix;
attribute vec3 vtx;
attribute vec3 vtxcenter;
attribute float value;
varying float valuev;
void main() {
	vec3 v = (vtx - vtxcenter) * drawCellScale + vtxcenter;
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
	valuev = value;
}
]],
		fragmentCode = [[
#version 460
uniform sampler1D gradTex;
varying float valuev;
out vec4 fragColor;
void main() {
	fragColor = texture1D(gradTex, valuev);
}
]],
		uniforms = {
			gradTex = 0,
			drawCellScale = 1,
		},
	}
	
	self.modelViewMatrix = matrix_ffi.zeros(4,4)
	self.projectionMatrix = matrix_ffi.zeros(4,4)
	self.modelViewProjectionMatrix = matrix_ffi.zeros(4,4)

	-- make a GPU version of the mesh
	-- TODO triangulate ... then include a mapping from triangle to source cell,
	-- and then just copy from display buf to a texture,
	-- then store per vertex the lookup to the texture.
	-- Also, if you are going to scale cells, then you must store a unique vertex per cell here.
	-- this means no dynamic mesh (without more coding).
	local glvtxs = vector'vec3f_t'			-- vertex position
	local glvtxcenters = vector'vec3f_t'	-- center of cell for this vertex
	local glcellindex = vector'float'		-- 0-based index of cell for this vertex
	local glvalue = vector'float'			-- value to draw the cell
	local function addTri(va,vb,vc, ci,c)
		glvtxs:push_back(vec3f(va:unpack()))
		glvtxs:push_back(vec3f(vb:unpack()))
		glvtxs:push_back(vec3f(vc:unpack()))
		for j=0,2 do
			glvtxcenters:push_back(vec3f(c.pos:unpack()))
			glcellindex:push_back(ci)
			glvalue:push_back(0)
		end
	end
	local mesh = self.mesh
	if self.dim == 2 then
		for ci,c in ipairs(mesh.cells) do
			local va = mesh.vtxs.v[mesh.cellVtxIndexes.v[0 + c.vtxOffset]]
			local vb = mesh.vtxs.v[mesh.cellVtxIndexes.v[1 + c.vtxOffset]]
			for vi=2,c.vtxCount-1 do
				local vc = mesh.vtxs.v[mesh.cellVtxIndexes.v[vi + c.vtxOffset]]
				addTri(va,vb,vc, ci, c)
				vb = vc
			end
		end
	elseif self.dim == 3 then
		for ci,c in ipairs(mesh.cells) do
			for i=0,c.faceCount-1 do
				local f = mesh.faces.v[mesh.cellFaceIndexes.v[i + c.faceOffset]]
				local va = mesh.vtxs.v[mesh.cellFaceIndexes.v[0 + f.vtxOffset]]
				local vb = mesh.vtxs.v[mesh.cellFaceIndexes.v[1 + f.vtxOffset]]
				for vi=2,f.vtxCount-1 do
					local vc = mesh.vtxs.v[mesh.cellFaceIndexes.v[vi + f.vtxOffset]]
					addTri(va,vb,vc, ci, c)
					vb = vc
				end
			end
		end
	end
	self.glvtxs = glvtxs
	self.glvtxcenters = glvtxcenters
	self.glcellindex = glcellindex
	self.glvalue = glvalue
	self.numGlVtxs = #glvtxs


	self.glVtxBuffer = GLArrayBuffer{
		data = glvtxs.v,
		size = #glvtxs * ffi.sizeof(glvtxs.type), 
	}
	self.glVtxCenterBuffer = GLArrayBuffer{
		data = glvtxcenters.v,
		size = #glvtxcenters * ffi.sizeof(glvtxcenters.type), 
	}
	self.glValueBuffer = GLArrayBuffer{
		data = glvalue.v,
		size = #glvalue * ffi.sizeof(glvalue.type),
	}


-- [[ uses VertexArray
	self.drawShaderVtxVertexArray = GLVertexArray{
		buffer = self.glVtxBuffer,
		size = 3,
		type = gl.GL_FLOAT,
	}
	self.drawShaderVtxVertexArray:bind()
	self.drawShader:setAttr('vtx', self.drawShaderVtxVertexArray)
	self.drawShaderVtxVertexArray:unbind()

	self.drawShaderVtxCenterVertexArray = GLVertexArray{
		buffer = self.glVtxCenterBuffer,
		size = 3,
		type = gl.GL_FLOAT,
	}
	self.drawShaderVtxCenterVertexArray:bind()
	self.drawShader:setAttr('vtxcenter', self.drawShaderVtxCenterVertexArray)
	self.drawShaderVtxCenterVertexArray:unbind()

	self.drawShaderValueVertexArray = GLVertexArray{
		buffer = self.glValueBuffer,
		size = 1,
		type = gl.GL_FLOAT,
	}
	self.drawShaderValueVertexArray:bind()
	self.drawShader:setAttr('value', self.drawShaderValueVertexArray)
	self.drawShaderValueVertexArray:unbind()
--]]

-- [[ doesn't use VertexArrays, but uses EnableVertexAttribArray instead
	self.glValueBuffer:bind()
	self.drawShaderValueVertexArray:setAttr(self.drawShader.attrs.value.loc)
	self.glValueBuffer:unbind()
	self.glVtxCenterBuffer:bind()
	self.drawShaderVtxCenterVertexArray:setAttr(self.drawShader.attrs.vtxcenter.loc)
	self.glVtxCenterBuffer:unbind()
	self.glVtxBuffer:bind()
	self.drawShaderVtxVertexArray:setAttr(self.drawShader.attrs.vtx.loc)
	self.glVtxBuffer:unbind()
--]]

	end	--if self.app.targetSystem ~= 'console' then

	-- no longer is dim * numCells the number of interfaces -- it is now dependent on the mesh
	-- maybe I should rename this to numFaces?
	
	local solver = self
	local Program = class(require 'cl.obj.program')
	if ffi.os == 'Windows' then
		os.execute'mkdir cache-cl 2> nul'
	else
		os.execute'mkdir cache-cl 2> /dev/null'
	end
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
	self:clalloc('facesBuf', 'face_t', self.numFaces)
	self:clalloc('cellFaceIndexesBuf', 'int', self.numCellFaceIndexes)

	-- similar to gridsolver's display buffer stuff when glSharing isn't found
	self.calcDisplayVarToTexPtr = ffi.new(self.app.real..'[?]', self.numCells * 3)
	
	-- specific to FiniteVolumeSolver
	self:clalloc('fluxBuf', self.eqn.cons_t, self.numFaces)
end

function MeshSolver:finalizeCLAllocs()
	MeshSolver.super.finalizeCLAllocs(self)

	self.vtxBufObj:fromCPU(self.mesh.vtxs.v)
	self.cellsBufObj:fromCPU(self.mesh.cells.v)
	self.facesBufObj:fromCPU(self.mesh.faces.v)
	self.cellFaceIndexesBufObj:fromCPU(self.mesh.cellFaceIndexes.v)
	self.cmds:finish()
end

function MeshSolver:checkStructSizes_getTypes()
	local typeinfos = MeshSolver.super.checkStructSizes_getTypes(self)
	typeinfos:append{
		face_t,
		cell_t,
	}
	return typeinfos
end

function MeshSolver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	-- numCells is the number of cells
	-- maybe I should rename it to numCells
	local localSize1d = math.min(maxWorkGroupSize, self.numCells)

	-- cell domain is the default
	self.domain = self.app.env:domain{
		size = {self.numCells},
		dim = 1,
	}

	-- face domain for flux computations
	-- ... and dt ?
	self.faceDomain = self.app.env:domain{
		size = {self.numFaces},
		dim = 1,
	}

	return {
		localSize1d = localSize1d, 
	}
end

function MeshSolver:resetState()
	MeshSolver.super.resetState(self)
--self:printBuf(self.UBufObj)
end

function MeshSolver:refreshSolverProgram()
	MeshSolver.super.refreshSolverProgram(self)

	self.calcFluxKernelObj = self.solverProgramObj:kernel{name='calcFlux', domain=self.faceDomain}
	self.calcFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcFluxKernelObj.obj:setArg(1, self.fluxBuf)
	self.calcFluxKernelObj.obj:setArg(2, self.UBuf)
	self.calcFluxKernelObj.obj:setArg(4, self.cellsBuf)
	self.calcFluxKernelObj.obj:setArg(5, self.facesBuf)
	self.calcFluxKernelObj.obj:setArg(6, self.cellFaceIndexesBuf)

	self.calcDerivFromFluxKernelObj = self.solverProgramObj:kernel'calcDerivFromFlux'
	self.calcDerivFromFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(3, self.cellsBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(4, self.facesBuf)
	self.calcDerivFromFluxKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
end

function MeshSolver:getSolverCode()
	return table{
		MeshSolver.super.getSolverCode(self),
	
		-- TODO flux scheme HERE
		-- this is the 'calcFluxKernel' code part of all the fvsolver subclasses
		-- so maybe make that its own module that plugs into both fvsolver and this
		template(file['hydro/solver/roe.cl'], {
			solver = self,
			eqn = self.eqn,
		}),

		-- finite volume integration
		template(file['hydro/solver/calcDerivFV.cl'], {
			solver = self,
			eqn = self.eqn,
		}),
	}:concat'\n'
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
	-- because I'm going to need more than just these...
	self.calcDTKernelObj.obj:setArg(3, self.cellsBuf)
	self.calcDTKernelObj.obj:setArg(4, self.facesBuf)
	self.calcDTKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
end

function MeshSolver:calcDT()
	if not self.useFixedDT then
		self.calcDTKernelObj.obj:setArg(3, self.cellsBuf)
		self.calcDTKernelObj.obj:setArg(4, self.facesBuf)
		self.calcDTKernelObj.obj:setArg(5, self.cellFaceIndexesBuf)
	end
	return MeshSolver.super.calcDT(self)
end

-- same as hydro/solver/fvsolver without PLM or CTU
-- TODO should MeshSolver only be a finite volume, or should we support finite difference as well?
-- or should we only support finite difference via flux calculations (like I do in fdsolver)?
function MeshSolver:calcDeriv(derivBufObj, dt)
	local dtArg = real(dt)

	self.calcFluxKernelObj.obj:setArg(3, dtArg)
	self.calcFluxKernelObj()

	self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFromFluxKernelObj()
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
	kernel:setArg(5, self.solver.facesBuf)
end


-- FPS of a 64x64 quad mesh:
-- sys=console: ~900
-- sys=imgui: ~2
-- so this function runs 500x slower
function MeshSolver:display(varName, ar)
	if self.app.targetSystem == 'console' then return end	

	local app = self.app

	local var = self.displayVarForName[varName]

-- [[ this is from the draw/*.lua
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = self:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end
--]]	

-- [[ matches GridSolver.calcDisplayVarToTex
	local component = self.displayComponentFlatList[var.component]
	local vectorField = self:isVarTypeAVectorField(component.type)
	
	local ptr = self.calcDisplayVarToTexPtr
	
	var:setToBufferArgs()
	self:calcDisplayVarToBuffer(var)

	local channels = vectorField and 3 or 1
	self.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * self.numCells * channels, ptr=ptr}
--]]

	local view = app.view
	view:projection(ar)
	view:modelview()

	
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projectionMatrix.ptr)
	self.modelViewProjectionMatrix:mul(self.projectionMatrix, self.modelViewMatrix)

	-- draw the mesh 
	local mesh = self.mesh

	-- if show cells ...
	-- TODO in my mesh solver, here I pick the color according to the display value
	-- but in my grid solvers I just overlay the whole thing with a giant texture
	-- so how to combine the two? use gradient tex.  
	-- in gridshader I have a kernel for updating tex/buffer to display values based on the cell state values
	-- i can do the same, and just update the 'displayValue' inside each cell for this.
	-- in fact, before drawing this, calcDisplayVarToBuffer will fill reduceBuf with the values (associated 1-1 with each cell?)
	if self.showCells then
		
		-- copy display values from ptr (reduceBuf) to each gl-vertex's value vertex attr
		for i=0,self.numGlVtxs-1 do
			local ci = self.glcellindex.v[i]
			local displayValue = ptr[channels * ci]
			self.glvalue.v[i] = displayValue
		end
		self.glValueBuffer:updateData()

		local gradientTex = app.gradientTex
		self.drawShader:use()
		gl.glUniform1f(self.drawShader.uniforms.drawCellScale.loc, self.drawCellScale)
		gl.glUniformMatrix4fv(self.drawShader.uniforms.modelViewProjectionMatrix.loc, 1, 0, self.modelViewProjectionMatrix.ptr)
		gradientTex:bind()
--[[ doesn't work yet
		self.drawShaderVtxVertexArray:bind()
		self.drawShaderVtxCenterVertexArray:bind()
		self.drawShaderValueVertexArray:bind()
		gl.glDrawArrays(gl.GL_TRIANGLES, 0, self.numGlVtxs * 3)
		self.drawShaderVtxVertexArray:unbind()
		self.drawShaderVtxCenterVertexArray:unbind()
		self.drawShaderValueVertexArray:unbind()
--]]
-- [[ 180 fps @ 64x64 quad grid
		gl.glEnableVertexAttribArray(self.drawShader.attrs.value.loc)
		gl.glEnableVertexAttribArray(self.drawShader.attrs.vtxcenter.loc)
		gl.glEnableVertexAttribArray(self.drawShader.attrs.vtx.loc)
		gl.glDrawArrays(gl.GL_TRIANGLES, 0, self.numGlVtxs * 3)
		gl.glDisableVertexAttribArray(self.drawShader.attrs.value.loc)
		gl.glDisableVertexAttribArray(self.drawShader.attrs.vtxcenter.loc)
		gl.glDisableVertexAttribArray(self.drawShader.attrs.vtx.loc)
--]]
--[[ 115 fps @ 64x64 quad grid
		gl.glBegin(gl.GL_TRIANGLES)
		for i=0,self.numGlVtxs-1 do
			gl.glVertexAttrib1f(self.drawShader.attrs.value.loc, self.glvalue.v[i])
			gl.glVertexAttrib3f(self.drawShader.attrs.vtxcenter.loc, self.glvtxcenters.v[i]:unpack())
			gl.glVertexAttrib3f(self.drawShader.attrs.vtx.loc, self.glvtxs.v[i]:unpack())
		end
		gl.glEnd()
--]]
		gradientTex:unbind()
		self.drawShader:useNone()
	end

	if self.showVertexes then
		gl.glPointSize(3)
		gl.glColor3f(1,1,1)
		gl.glBegin(gl.GL_POINTS)
		for _,v in ipairs(mesh.vtxs) do
			gl.glVertex3d(v:unpack())
		end
		gl.glEnd()
		gl.glPointSize(1)
	end

	if self.showFaces then
		for fi,f in ipairs(mesh.faces) do
			gl.glBegin(gl.GL_LINE_LOOP)
			for vi=0,f.vtxCount-1 do
				local v = mesh.vtxs.v[mesh.faceVtxIndexes.v[vi + f.vtxOffset]]
				gl.glVertex3d(v:unpack())
			end
			gl.glEnd()
		end
	end

-- [[ also in draw/*.lua
	local gradientValueMin = valueMin
	local gradientValueMax = valueMax
	local showName = varName
	if var.showInUnits and var.units then
		local unitScale = self:convertToSIUnitsCode(var.units).func()
		gradientValueMin = gradientValueMin * unitScale
		gradientValueMax = gradientValueMax * unitScale
		showName = showName..' ('..var.units..')'
	end
	app:drawGradientLegend(ar, showName, gradientValueMin, gradientValueMax)
--]]

	glreport'here'
end

function MeshSolver:updateGUIParams()
	MeshSolver.super.updateGUIParams(self)

	tooltip.checkboxTable('show vertexes', self, 'showVertexes')
	ig.igSameLine()
	tooltip.checkboxTable('show faces', self, 'showFaces')
	ig.igSameLine()
	tooltip.checkboxTable('show cells', self, 'showCells')
	
	tooltip.numberTable('cell scale', self, 'drawCellScale')
end

return MeshSolver 