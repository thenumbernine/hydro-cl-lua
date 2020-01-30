--[[
chopping up a solver
for use with multi-device
maybe also with AMR

this should have a matching equation ... or at least equation arugments and initial condition ...

matching equation object might be useful, since that means only one set of display args as well
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local vec3sz = require 'vec-ffi.vec3sz'
local vec3d = require 'vec-ffi.vec3d'

local Chopped = class()

function Chopped:init(args)
	assert(args)
	self.app = assert(args.app)

	local subsolverClass = assert(args.subsolverClass)
	self.multiSlices = assert(args.multiSlices)

	-- TODO handle like solver/gridsolver.lua
	self.dim = assert(args.dim)

	local function getmaxaxis(v)
		local maxindex = 1
		for j=2,self.dim do
			if v:ptr()[j] > v:ptr()[maxindex] then
				maxindex = j
			end
		end
		return maxindex
	end

	self.mins = args.mins or {-1, -1, -1}
	self.maxs = args.maxs or {1, 1, 1}
	
	-- TODO handle different gridSize types matching solver/gridsolver.lua
	local root = {
		size = vec3sz(table.unpack(args.gridSize)),
		mins = vec3d(table.unpack(self.mins)),
		maxs = vec3d(table.unpack(self.maxs)),
		boundary = {min={}, max={}},
	}
	root.maxaxis = getmaxaxis(root.size)
	root.maxsize = root.size:ptr()[root.maxaxis]
	
	local blocks = table{root}
	for i=2,self.multiSlices do
		local longestBlock = blocks:remove(
			(select(2, blocks:sup(function(a,b) 
				return a.maxsize > b.maxsize 
			end)))
		)
		local axis = longestBlock.maxaxis
		local newsize = vec3sz(longestBlock.size:unpack())
		local mid = .5 * (
			longestBlock.mins:ptr()[axis]
			+ longestBlock.maxs:ptr()[axis]
		)
		newsize:ptr()[axis] = bit.rshift(newsize:ptr()[axis], 1)
		if newsize:ptr()[axis] == 0 then
			blocks:insert(longestBlock)
			break	-- done
		end
		
		local blockL = {
			size = vec3sz(newsize:unpack()),
			mins = vec3d(longestBlock.mins:unpack()),
			maxs = vec3d(longestBlock.maxs:unpack()),
			boundary = {
				min = table(longestBlock.boundary.min),
				max = table(longestBlock.boundary.max),
			},
		}
		blockL.maxs:ptr()[axis] = mid 
		
		blockL.maxaxis = getmaxaxis(blockL.size)
		blockL.maxsize = blockL.size:ptr()[blockL.maxaxis]


		local blockR = {
			size = vec3sz(newsize:unpack()),
			mins = vec3d(longestBlock.mins:unpack()),
			maxs = vec3d(longestBlock.maxs:unpack()),
			boundary = {
				min = table(longestBlock.boundary.min),
				max = table(longestBlock.boundary.max),
			},
		}
		blockR.mins:ptr()[axis] = mid
		
		blockR.maxaxis = getmaxaxis(blockR.size)
		blockR.maxsize = blockR.size:ptr()[blockR.maxaxis]
	
	
		
		-- blockL and blockR inherit boundary from longestBlock
		-- ... and simulatenously take everything that points to longestBlock and change it to either point to blockL or blockR
		-- then point:
		-- blockL boundary.max[maxaxis] = blockR
		-- blockR boundary.min[maxaxis] = blockL
		blockL.boundary.max[axis+1] = blockR
		blockR.boundary.min[axis+1] = blockL
		
		
		for _,other in ipairs(blocks) do
			for _,minmax in ipairs{'min', 'max'} do
				for j=1,3 do
					if other.boundary[minmax][j] == longestBlock then
						-- redirect it to the new children
						-- make sure to match up the sections of the solvers
					end
				end
			end
		end

		blocks:insert(blockL)
		blocks:insert(blockR)
	end
	
	print'blocks:'
	for _,block in ipairs(blocks) do
		print(block.size, block.mins, block.maxs)
	end
	print()

	-- now make subsolvers for each of these solvers
	self.blocks = blocks
	self.solvers = blocks:mapi(function(block,i)
		local subsolver = subsolverClass(table(args, {
			mins = {block.mins:unpack()},
			maxs = {block.maxs:unpack()},
			gridSize = {block.size:unpack()},
			device = self.app.env.devices[(i - 1) % #self.app.env.devices + 1],
			cmds = self.app.env.cmds[(i - 1) % #self.app.env.cmds + 1],
			-- this overrides bounds in the initial condition ... it's only really used with the 'lhs' flag of Sod etc
			initCondMins = self.mins,
			initCondMaxs = self.maxs,
		}))
		subsolver.choppedupBoundaryInfo = block.boundary
		block.solver = subsolver
		return subsolver 
	end)
	for _,solver in ipairs(self.solvers) do
		for _,minmax in ipairs{'min', 'max'} do
			for j=1,3 do
				if solver.choppedupBoundaryInfo[minmax][j] then
					solver.choppedupBoundaryInfo[minmax][j] = solver.choppedupBoundaryInfo[minmax][j].solver
				end
			end
		end
	end
end

--[[
things app needs to run a solver:
	displayVars = table
	t = number
	color = {r,g,b}
	name = string
	update = function
	initDraw = function
	
these are needed for displaying things
	displayVars
	displayVarForName
--]]

function Chopped:initDraw()
	self.displayVars = table()
	self.t = 0
	self.color = vec3(math.random(), math.random(), math.random()):normalize()
	self.name = 'Chopped '..self.solvers[1].name
	for _,solver in ipairs(self.solvers) do
		solver:initDraw()
	end

	-- I guess this is safe for now
	self.displayVars = self.solvers[1].displayVars
	self.displayVarForName = self.solvers[1].displayVars
end

function Chopped:resetState()
	for _,solver in ipairs(self.solvers) do
		solver:resetState()
	end
end

function Chopped:update()
	for _,solver in ipairs(self.solvers) do
		solver:update()
	end
	self.t = self.solvers:mapi(function(solver) return solver.t end):inf()
end


return Chopped
