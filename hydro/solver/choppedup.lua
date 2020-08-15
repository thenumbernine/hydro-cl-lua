--[[
chopping up a solver
for use with multi-device
maybe also with AMR

this should have a matching equation ... or at least equation arugments and initial condition ...

matching equation object might be useful, since that means only one set of display args as well

TODO TODO TODO
	the init state code matches between each sub-solver
	the solver code matches between each sub solver
	the boundary code matches close enough
	so all we really need to do here is create separate grids, but reuse the same eqn, same eqn's initCond, and even the same solver-code (not the same kernels and CPU/GPU buffers though)
	that will save on a lot of initialization time and memory
	so TODO separate the code from the buffers, then chopped and AMR can both reuse code 
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local coroutine = require 'ext.coroutine'
local vec3sz = require 'vec-ffi.vec3sz'
local vec3d = require 'vec-ffi.vec3d'
local cl = require 'ffi.OpenCL'
local classert = require 'cl.assert'
local ffi = require 'ffi'
local tolua = require 'ext.tolua'

local Chopped = class()

function Chopped:init(args)
	assert(args)
	self.app = assert(args.app)

	local subsolverClass = assert(args.subsolverClass)
	
	-- TODO handle like hydro/solver/gridsolver.lua
	self.dim = assert(args.dim)
	
	local gridSize = vec3sz(table.unpack(args.gridSize))

	self.multiSlices = vec3sz(table.unpack(args.multiSlices))
	for j=self.dim+1,3 do
		self.multiSlices.s[j-1] = 1
	end
	
	self.mins = args.mins or {-1, -1, -1}
	self.maxs = args.maxs or {1, 1, 1}

	local mins = vec3d(table.unpack(self.mins))
	local maxs = vec3d(table.unpack(self.maxs))

	self.blocks = table()	-- [i][j][k]
	self.flatBlocks = table() -- [index]
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
				self.flatBlocks:insert(block)
			end
		end
	end

	if self.app.verbose then
		print'blocks:'
		for _,block in ipairs(self.flatBlocks) do
			print(tolua(block))
		end
		print()
	end
	-- now make subsolvers for each of these solvers
	self.solvers = self.flatBlocks:mapi(function(block,i)
		local deviceIndex = (i - 1) % #self.app.env.devices + 1
		local cmdsIndex = (i - 1) % #self.app.env.cmds + 1
		local device = self.app.env.devices[deviceIndex]
		local cmds = self.app.env.cmds[cmdsIndex]
		if self.app.verbose then
			print('assigning subsolver '..i..' to device #'..deviceIndex..': '..device:getName()..' and cmdqueue #'..cmdsIndex..': '..tostring(cmds))
			io.stdout:flush()
		end
		local subsolver = subsolverClass(table(args, {
			mins = {block.mins:unpack()},
			maxs = {block.maxs:unpack()},
			gridSize = {block.size:unpack()},
			device = device,
			cmds = cmds,
			-- this overrides bounds in the initial condition ... it's only really used with the 'lhs' flag of Sod etc
			initCondMins = self.mins,
			initCondMaxs = self.maxs,
		}))
		block.solver = subsolver
		return subsolver
	end)

	-- synchronize boundaries between devices
	-- TODO do this after all the solvers' buffers have been created
	-- ... but before the boundary programs have been created
	-- TODO TODO  right now the boundary code runs the left and right sides in the same kernel
	--   I would like to override the left and right separately -- and put the copy-between-devices code there
	--   but this might mean splitting the kernels
	--  alternatively, I could just change the kernel code to read from a neighboring device.
	--   this would mean adding extra buffer arguments to each boundary kernel 
	--
	-- for now I'll just override the ':boundary()' of each solver to after-the-fact copy stuff across
	-- TODO TODO TODO solver:boundary() is called every frame both at int/fe.lua:46 and at hydro/solver/gridsolver.lua:1569 ... I only need one of those
	local parent = self
	local sizeof_cons_t = ffi.sizeof(self.solvers[1].eqn.cons_t)	-- assumes all solver have the same size eqn_t
	for i=1,tonumber(self.multiSlices.x) do
		for j=1,math.max(1,tonumber(self.multiSlices.y)) do
			for k=1,math.max(1,tonumber(self.multiSlices.z)) do
				local v = vec3sz(i,j,k)
				local block = self.blocks[i][j][k]
				local solver = block.solver
				local numGhost = solver.numGhost -- == solverL.numGhost == solverR.numGhost
				
				local solverLs = table()
				local solverRs = table()
				-- src_origin is within either 'solverL' or 'solverR'
				local src_originLs = table()
				local src_originRs = table()
				-- dst_origin is always within 'solver'
				local dst_originLs = table()
				local dst_originRs = table()
				-- region is equal per-block, regardless of solver, solverL, solverR
				local regions = table()
				for side=0,self.dim-1 do
					local updateL = v.s[side] > 1
					local updateR = v.s[side] < self.multiSlices.s[side]
					
					local vL = vec3sz(v:unpack())
					vL.s[side] = vL.s[side] - 1
						
					local vR = vec3sz(v:unpack())
					vR.s[side] = vR.s[side] + 1
					
					local solverL = updateL and self.blocks[tonumber(vL.x)][tonumber(vL.y)][tonumber(vL.z)].solver or nil
					local solverR = updateR and self.blocks[tonumber(vR.x)][tonumber(vR.y)][tonumber(vR.z)].solver or nil
					
		
					local src_originL, dst_originL  
					if updateL then
						src_originL = vec3sz()
						src_originL.s[side] = solverL.gridSize.s[side] - 2 * numGhost
						dst_originL = vec3sz()
					end

					local src_originR, dst_originR
					if updateR then
						src_originR = vec3sz()
						src_originR.s[side] = numGhost
						dst_originR = vec3sz()
						dst_originR.s[side] = solver.gridSize.s[side] - numGhost
					end
					
					-- TODO just make one per side and reuse it for all blocks
					local region = vec3sz(solver.gridSize:unpack())
					region.s[side] = numGhost
					
					-- scale x-axis by structure size
					if updateL then
						src_originL.x = src_originL.x * sizeof_cons_t
						dst_originL.x = dst_originL.x * sizeof_cons_t
					end
					if updateR then
						src_originR.x = src_originR.x * sizeof_cons_t
						dst_originR.x = dst_originR.x * sizeof_cons_t
					end
					region.x = region.x * sizeof_cons_t


					solverLs[side+1] = solverL
					solverRs[side+1] = solverR
					src_originLs[side+1] = src_originL
					src_originRs[side+1] = src_originR
					dst_originLs[side+1] = dst_originL
					dst_originRs[side+1] = dst_originR
					regions[side+1] = region
				end

				function block:syncBorders()
if cmdline.dontSyncBordersOnMultiGPU then return end   -- debugging multi-gpu performance
					for side=0,solver.dim-1 do
						local solverL = solverLs[side+1]
						local solverR = solverRs[side+1]
						if solverL then
							classert(cl.clEnqueueCopyBufferRect(
								solver.cmds.id,					-- cl_command_queue command_queue,
								solverL.UBuf.id,				-- cl_mem src_buffer,
								solver.UBuf.id,					-- cl_mem dst_buffer,
								src_originLs[side+1].s,		-- const size_t src_origin[3],
								dst_originLs[side+1].s,		-- const size_t dst_origin[3],
								regions[side+1].s,			-- const size_t region[3],
								sizeof_cons_t * solverL.gridSize.x,	-- size_t src_row_pitch,
								sizeof_cons_t * solverL.gridSize.x * solverL.gridSize.y,	-- size_t src_slice_pitch,
								sizeof_cons_t * solver.gridSize.x,	-- size_t dst_row_pitch,
								sizeof_cons_t * solver.gridSize.x * solver.gridSize.y,	-- size_t dst_slice_pitch,
								0,								-- cl_uint num_events_in_wait_list,
								nil,							-- const cl_event *event_wait_list,
								nil								-- cl_event *event
							))
						end
						if solverR then
							classert(cl.clEnqueueCopyBufferRect(
								solver.cmds.id,					-- cl_command_queue command_queue,
								solverR.UBuf.id,				-- cl_mem src_buffer,
								solver.UBuf.id,					-- cl_mem dst_buffer,
								src_originRs[side+1].s,		-- const size_t src_origin[3],
								dst_originRs[side+1].s,		-- const size_t dst_origin[3],
								regions[side+1].s,			-- const size_t region[3],
								sizeof_cons_t * solverR.gridSize.x,	-- size_t src_row_pitch,
								sizeof_cons_t * solverR.gridSize.x * solverR.gridSize.y,	-- size_t src_slice_pitch,
								sizeof_cons_t * solver.gridSize.x,	-- size_t dst_row_pitch,
								sizeof_cons_t * solver.gridSize.x * solver.gridSize.y,	-- size_t dst_slice_pitch,
								0,								-- cl_uint num_events_in_wait_list,
								nil,							-- const cl_event *event_wait_list,
								nil								-- cl_event *event
							))
						end
					end
				end
			end
		end
	end

	--[[ how to sync all functions
	
	GridSolver.update:
		SolverBase.update
		solver:calcDT
		solver:step
		solver:boundary
	
	so I could either pick the update internals apart and recreate them in here
		but then I'd have to take care to recreate every single thing in here
	or I could just put update in a coroutine
		and then correctly :resume inside here
	--]]
	for i,solver in ipairs(self.solvers) do

		local oldUpdate = solver.update
		function solver:update(...)
			oldUpdate(self, ...)
			
			local result = coroutine.yield'S4) update done'	-- assert resume gets this to make sure our yields are lined up
--print(result)		
		end
	
		local oldCalcDT = solver.calcDT
		function solver:calcDT(...)
			local oldDT = oldCalcDT(self, ...)
			
			local result = coroutine.yield('S1) passing old dt', oldDT)
--print(result)		
			assert(result == 'M2) main got sub-solver dt', "expected 'M2) main got sub-solver dt', got "..tolua(result))
			
			local newDT = coroutine.yield'S2) getting new dt'
assert(newDT, "didn't get a dt back")
			
			coroutine.yield'S3) waiting for dt'
			
			return newDT
		end
	end
	
	for	_,block in ipairs(self.flatBlocks) do
		local solver = block.solver
		block.updateThread = coroutine.create(function(arg)
			assert(arg == 'M0) main init coroutine')
		
			local result = coroutine.yield'S0) start update'	-- let our calling coroutine keep the solver stopped at the beginning of the update routine
--print(result)		
			assert(result == 'M1) starting update', "expected 'M1) starting update', got "..tolua(result))

			while true do
				solver:update()
			end
		end)
		
		local err, result = coroutine.assertresume(block.updateThread, 'M0) main init coroutine')
--print(result)		
		assert(err)
		assert(result == 'S0) start update', "expected 'S0) start update', found "..tolua(result))
	
	end
	
	self.t = 0
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
	self.color = vec3d(math.random(), math.random(), math.random()):normalize()
	self.name = 'Chopped '..self.solvers[1].name
	for _,solver in ipairs(self.solvers) do
		solver:initDraw()
		solver.color = vec3d(self.color:unpack())
	end

	-- I guess this is safe for now
	self.displayVars = self.solvers[1].displayVars
	self.displayVarForName = self.solvers[1].displayVars
end

function Chopped:resetState()
	for _,solver in ipairs(self.solvers) do
		solver:resetState()
	end
	self.t = 0
end

function Chopped:update()

	local dt = math.huge
	for i,block in ipairs(self.flatBlocks) do
		-- run from 'S0) start update' until calc dt, return the old dt
		local err, result, solverdt = coroutine.assertresume(block.updateThread, 'M1) starting update')
		assert(err)
--print(result)		
		assert(result == 'S1) passing old dt', "expected 'S1) passing old dt', found "..tolua(result))
		assert(type(solverdt) == 'number', "expected to get back number as sub-solver dt, found "..tolua(result))
--print('main solver gets solver '..i..' dt = '..tostring(solverdt))
		dt = math.min(dt, solverdt)

		local err, result = coroutine.assertresume(block.updateThread, 'M2) main got sub-solver dt')
--print(result)		
		assert(err)
		--assert(result == 
	end
--print('chopped solver dt', dt)	

	for _,block in ipairs(self.flatBlocks) do
		local err, result = coroutine.assertresume(block.updateThread, dt, 'M3) main passing super-solver dt')
--print(result)		
		assert(err)
		assert(result == 'S3) waiting for dt', "expected 'S3) waiting for dt', found "..tolua(result))
	end

	for _,block in ipairs(self.flatBlocks) do
		-- run from calc dt until 'S4) update done'
		local err, result = coroutine.assertresume(block.updateThread, 'M1) starting update')
		assert(err)
--print(result)		
		assert(result  == 'S4) update done', "expected 'S4) update done', got "..tolua(result))
	end

	self.t = self.t + dt
	--self.t = self.solvers[1].t	
--print('chopped solver t', self.t)

--[[
print'before sync'
	for _,solver in ipairs(self.solvers) do
		solver:printBuf(solver.UBufObj)
	end
--]]
	
	-- synchronize boundaries between devices
	-- should this happen inter-step for RK4, etc integrators?  or only at the end of update() ?
	for i=1,tonumber(self.multiSlices.x) do
		for j=1,math.max(1,tonumber(self.multiSlices.y)) do
			for k=1,math.max(1,tonumber(self.multiSlices.z)) do
				--print('synchronizing '..vL..' and '..vR)
				local block = self.blocks[i][j][k]
				block:syncBorders()
			end
		end
	end

--[[
print'after sync'
	for _,solver in ipairs(self.solvers) do
		solver:printBuf(solver.UBufObj)
	end
--]]

end

function Chopped:updateGUI(...)
	return self.solvers[1]:updateGUI(...)
end

return Chopped
