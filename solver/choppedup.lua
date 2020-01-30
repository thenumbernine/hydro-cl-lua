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
	
	-- TODO handle like solver/gridsolver.lua
	self.dim = assert(args.dim)
	
	local gridSize = vec3sz(table.unpack(args.gridSize))

	self.multiSlices = vec3sz(table.unpack(args.multiSlices))
	for j=self.dim+1,3 do
		self.multiSlices:ptr()[j-1] = 1
	end
	
	self.mins = args.mins or {-1, -1, -1}
	self.maxs = args.maxs or {1, 1, 1}

--[[ chop things up recursively ... 
-- but this makes copying borders a bit more painful, 
-- and if you assign 1 solver to 1 device then, if you only make as many chops as devices then some devices might get more work than otheres
--	... unless you chopped up solvers into too fine of pieces, then the differences between devices would be neglegible
	local function getmaxaxis(v)
		local maxindex = 1
		for j=2,self.dim do
			if v:ptr()[j] > v:ptr()[maxindex] then
				maxindex = j
			end
		end
		return maxindex
	end
	
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
--]]
-- [[ just explicitly state the number of chops along each axis
	local mins = vec3d(table.unpack(self.mins))
	local maxs = vec3d(table.unpack(self.maxs))

	self.blocks = table()	-- [i][j][k]
	self.flattenedBlocks = table() -- [index]
	for i1=0,tonumber(self.multiSlices.x-1) do
		local bi = table()
		self.blocks:insert(bi)
		for i2=0,tonumber(self.multiSlices.y-1) do
			local bj = table()
			bi:insert(bj)
			for i3=0,tonumber(self.multiSlices.z-1) do
				local i = vec3sz(i1, i2, i3)
				local fL = vec3d(i:unpack()) / vec3d(self.multiSlices:unpack())
				local fR = vec3d((i + 1):unpack()) / vec3d(self.multiSlices:unpack())
				local start = vec3sz((fL * vec3d(gridSize:unpack())):unpack())
				local finish = vec3sz((fR * vec3d(gridSize:unpack())):unpack())
				local block = {
					size = finish - start,
					mins = fL * (maxs - mins) + mins,
					maxs = fR * (maxs - mins) + mins,
				}
				bj:insert(block)
				self.flattenedBlocks:insert(block)
			end
		end
	end
--]]
	print'blocks:'
	for _,block in ipairs(self.flattenedBlocks) do
		--print(block.size, block.mins, block.maxs)
		print(require 'ext.tolua'(block))
	end
	print()
	-- now make subsolvers for each of these solvers
	self.solvers = self.flattenedBlocks:mapi(function(block,i)
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
--		subsolver.choppedupBoundaryInfo = block.boundary
		block.solver = subsolver
		return subsolver 
	end)

--[[	
	for _,solver in ipairs(self.solvers) do
		for _,minmax in ipairs{'min', 'max'} do
			for j=1,3 do
				if solver.choppedupBoundaryInfo[minmax][j] then
					solver.choppedupBoundaryInfo[minmax][j] = solver.choppedupBoundaryInfo[minmax][j].solver
				end
			end
		end
	end
--]]
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

	-- synchronize boundaries between devices
	-- should this happen inter-step for RK4, etc integrators?  or only at the end of update() ?
	for i=1,tonumber(self.multiSlices.x-1) do
		for j=1,math.max(1,tonumber(self.multiSlices.y-1)) do
			for k=1,math.max(1,tonumber(self.multiSlices.z-1)) do
				local vL = vec3sz(i,j,k)
				for side=1,self.dim do
					if vL:ptr()[side-1] < self.multiSlices:ptr()[side-1] then
						local vR = vec3sz(vL:unpack())
						vR:ptr()[side-1] = vR:ptr()[side-1] + 1
						--print('synchronizing '..vL..' and '..vR)
					end
				end
			end
		end
	end
end


return Chopped
