local class = require 'ext.class'
local file = require 'ext.file'
local template = require 'template'
local gl = require 'ffi.OpenGL'
local ffi = require 'ffi'
local table = require 'ext.table'
local Draw = require 'hydro.draw.draw'


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

local function reverse(t)
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
local edgesForInside = table(single)

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


local Draw3DIso = class(Draw)

function Draw3DIso:showDisplayVar(var)
	local solver = self.solver
	local app = solver.app
	local shader = app.isobarShader

	shader:use()
	gl.glBegin(gl.GL_TRIANGLES)
	
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
	
	self:setupDisplayVarShader(shader, var, valueMin, valueMax)
	
	assert(not app.useGLSharing, "I still need to code in the GL sharing version")
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
							cornerOfs = cornerOfs + solver.stepSize.s[n]
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
											v.s[m] = v.s[m] + 1
										end
									end
									
									local frac = (isoValue - cornerValues[corner]) / (cornerValues[nextCorner] - cornerValues[corner])
									v.s[n] = v.s[n] + sign * frac
								
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
								
								x = (x - solver.numGhost) / tonumber(solver.gridSize.x - 2 * solver.numGhost) * (solver.maxs.x - solver.mins.x) + solver.mins.x
								y = (y - solver.numGhost) / tonumber(solver.gridSize.y - 2 * solver.numGhost) * (solver.maxs.y - solver.mins.y) + solver.mins.y
								z = (z - solver.numGhost) / tonumber(solver.gridSize.z - 2 * solver.numGhost) * (solver.maxs.z - solver.mins.z) + solver.mins.z
							
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
	
	gl.glEnd()
	shader:useNone()
end

function Draw3DIso:display(varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = solver.app
	app.view:setup(ar)

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
			x = x * (solver.maxs.x - solver.mins.x) + solver.mins.x
			y = y * (solver.maxs.y - solver.mins.y) + solver.mins.y
			z = z * (solver.maxs.z - solver.mins.z) + solver.mins.z
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
	

	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var)
	end
		
	gl.glDisable(gl.GL_DEPTH_TEST)
	gl.glDisable(gl.GL_CULL_FACE)
	gl.glDisable(gl.GL_BLEND)
end

-- TODO this in common with 3d_ray.lua.  subclass?
function Draw3DIso:prepareShader()
	local solver = self.solver
	if solver.volumeRayShader then return end
	
	solver.display3D_Ray_maxiter = math.max(tonumber(solver.gridSize.x), tonumber(solver.gridSize.y), tonumber(solver.gridSize.z))
	
	local volumetricCode = file['hydro/draw/volumetric.shader']
	solver.volumeRayShader = solver.GLProgram{
		name = 'volumetric',
		vertexCode = template(volumetricCode, {
			app = solver.app,
			solver = solver,
			vertexShader = true,
		}),
		fragmentCode = template(volumetricCode, {
			app = solver.app,
			solver = solver,
			fragmentShader = true,
		}),
		uniforms = {
			tex = 0,
			gradientTex = 1,
			oneOverDx = {(solver.maxs - solver.mins):unpack()},
		},
	}
end

return Draw3DIso
