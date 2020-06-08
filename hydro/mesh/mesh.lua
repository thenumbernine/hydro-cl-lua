--[[
parent class of meshes used by hydro/solver/meshsolver
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local vec2i = require 'vec-ffi.vec2i'
local vector = require 'hydro.util.vector'
-- one of these is bound to be real3, right?
local vec2d = require 'vec-ffi.vec2d'
local vec3f = require 'vec-ffi.vec3f'
local vec3d = require 'vec-ffi.vec3d'
local struct = require 'hydro.struct.struct'


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
local meshTypeCode = table{
	vec2i.typeCode,
	face_t:getTypeCode(),
	cell_t:getTypeCode()
}:concat'\n'

-- TODO what if real3 isn't defined yet?
ffi.cdef(meshTypeCode)



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

Mesh.meshTypeCode = meshTypeCode

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
		com = com + (a + b + c) * normal:dot((c - a):cross(c - b))
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
		for i=3,#face do
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
		local delta = self.real3()
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
			
			local fi = self:addFace(thisFaceVtxIndexes, #thisFaceVtxIndexes, ci)
			do -- if fi ~= -1 then
				local f = self.faces.v[fi]
				
				--if face area is zero then don't add it to cell's faces
				-- and don't use it later for determining cell's volume
				if f.area <= 1e-7 then
				else
					self.cellFaceIndexes:push_back(fi)

					cubeVtxs:insert(
						table.mapi(thisFaceVtxIndexes, function(i)
							return self.vtxs.v[i]
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


return Mesh
