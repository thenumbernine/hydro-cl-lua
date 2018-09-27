--[[
TODO a lot of this looks the same as the separate solvers -- combine them into a common class
--]]
local ffi = require 'ffi'
local table = require 'ext.table'
local class = require 'ext.class'
local vec3sz = require 'ffi.vec.vec3sz'
local ig = require 'ffi.imgui'
local time, getTime = table.unpack(require 'time')

local AMR = class()

AMR.name = 'AMR'

function AMR:init(implClass, args)
	self.implClass = implClass
	--self.amrMethod = 'dt vs 2dt'
	--self.amrMethod = 'gradient'

	local mainsolver = implClass(args)
	self.solvers = table{mainsolver}

	-- also in twofluid-emhd-separate-behavior
	self.displayVars = table():append(self.solvers:map(function(solver) 
		return solver.displayVars 
	end):unpack())

	self.solverForDisplayVars = table()
	for _,solver in ipairs(self.solvers) do
		for _,var in ipairs(solver.displayVars) do
			self.solverForDisplayVars[var] = solver
			var.name = solver.name:gsub('[%s]', '_')..'_'..var.name
		end
	end

	self.displayVarForName = self.displayVars:map(function(var)
		return var, var.name
	end)

	self.color = vec3(math.random(), math.random(), math.random()):normalize()

	self.numGhost = mainsolver.numGhost
	self.dim = mainsolver.dim
	self.gridSize = vec3sz(mainsolver.gridSize)
	self.sizeWithoutBorder = vec3sz(mainsolver.sizeWithoutBorder)
	self.mins = vec3(mainsolver.mins:unpack())
	self.maxs = vec3(mainsolver.maxs:unpack())

	self.coord = mainsolver.coord
	self.eqn = {
		numStates = self.solvers:map(function(solver) return solver.eqn.numStates end):sum(),
		numWaves = self.solvers:map(function(solver) return solver.eqn.numWaves end):sum(),
		getEigenTypeCode = function() end,
		getCodePrefix = function() end,
	}

	self.t = 0

	-- TODO here
	-- (this used to go in GridSolver:createBuffers
	-- depending on the AMR type (dt, gradient, etc), allocate extra buffers
	-- (

	if self.amrMethod == 'gradient' then
		--[[
		ok here's my thoughts on the size ...
		I'm gonna try to do AMR
		that means storing the nodes in the same buffer
		I think I'll pad the end of the UBuf with the leaves
		and then store a tree of index information somewhere else that says what leaf goes where
		(maybe that will go at the end)
		
		how should memory breakdown look?
		how big should the leaves be?

		how about do like reduce ...
		leaves can be 16x16 blocks
		a kernel can cycle through the root leaves and sum te amrError 
		-- in the same kernel as it is calculated
		-- then I don't need to store so much memory, only one value per leaf, not per grid ...
		-- then ... split or merge leaves, based on their error ...

		how to modify that?  in anothe kernel of its own ...
		bit array of whether each cell is divided ...
		-- then in one kernel we update all leaves
		-- in another kernel, populate leaf ghost cells
		-- in another kernel, 
		-- 		decide what should be merged and, based on that, copy from parent into this
		--		if something should be split then look for unflagged children and copy into it
		
		how many bits will we need?
		volume / leafSize bits for the root level
		
		then we have parameters of 
		- how big each node is
		- what level of refinement it is
		
		ex:
		nodes in the root level are 2^4 x 2^4 = 2^8 = 256 cells = 256 bits = 2^5 = 32 bytes
		

		leafs are 2^4 x 2^4 = 2^8 = 256 cells
		... but ghost cells are 2 border, so we need to allocate (2^n+2*2)^2 cells ... for n=4 this is 400 cells ... 
			so we lose (2^n+2*2)^2 - 2^(2n) = 2^(2n) + 2^(n+3) + 2^4 - 2^(2n) = 2^(n+3) + 2^4) cells are lost
		
		so leafs multipy at a factor of 2^2 x 2^2 = 2^4 = 16
		so the next level has 2^(8+4) = 2^12 = 4096 bits = 2^9 = 512 bytes

		--]]

		-- the size, in cells, which a node replaces
		self.amrNodeFromSize = ({
			vec3sz(8, 1, 1),
			vec3sz(8, 8, 1),
			vec3sz(8, 8, 8),
		})[self.dim]
	print('self.amrNodeFromSize', self.amrNodeFromSize)

		-- the size, in cells, of each node, excluding border, for each dimension
		self.amrNodeSizeWithoutBorder = ({
			vec3sz(16, 1, 1),
			vec3sz(16, 16, 1),
			vec3sz(16, 16, 16),
		})[self.dim]
	print('self.amrNodeSizeWithoutBorder', self.amrNodeSizeWithoutBorder)

		-- size of the root level, in terms of nodes ('from' size)
		self.amrRootSizeInFromSize = vec3sz(1,1,1)
		for i=0,self.dim-1 do
			self.amrRootSizeInFromSize:ptr()[i] = 
				roundup(self.sizeWithoutBorder:ptr()[i], self.amrNodeFromSize:ptr()[i]) 
					/ self.amrNodeFromSize:ptr()[i]
		end
	print('self.amrRootSizeInFromSize', self.amrRootSizeInFromSize)

		-- how big each node is
		self.amrNodeSize = self.amrNodeSizeWithoutBorder + 2 * self.numGhost

		-- how many nodes to allocate and use
		-- here's the next dilemma in terms of memory layout
		-- specifically in terms of rendering
		-- if I want to easily copy and render the texture information then I will need to package the leafs into the same texture as the root
		-- which means extending the texture buffer in some particular direction.
		-- since i'm already adding the leaf information to the end of the buffer
		-- and since appending to the end of a texture buffer coincides with adding extra rows to the texture
		-- why not just put our leafs in extra rows of -- both our display texture and of  
		self.amrMaxNodes = 1
		
		-- this will hold info on what leafs have yet been used
		self.amrLeafs = table()

		-- hmm, this is the info for the root node ...
		-- do I want to keep the root level data separate?
		-- or do I just want to represent everything as a collection of leaf nodes?
		-- I'll keep the root structure separate for now
		-- so I can keep the original non-amr solver untouched
		self.amrLayers = table()
		self.amrLayers[1] = table()	-- here's the root
	end

	--[[ this used to go after createBuffers UBufSize
	if self.amrMethod == 'gradient' then	
		UBufSize = UBufSize + self.amrMaxNodes * self.amrNodeSize:volume()
	end	
	--]]
	local app = self.app
	for _,solver in ipairs(self.solvers) do
		if self.amrMethod == 'dt vs 2dt' then
			-- here's my start at AMR, using the 1989 Berger, Collela two-small-steps vs one-big-step method
			solver:clalloc('lastUBuf', solver.numCells * ffi.sizeof(solver.eqn.cons_t))
			solver:clalloc('U2Buf', solver.numCells * ffi.sizeof(solver.eqn.cons_t))
		elseif self.amrMethod == 'gradient' then
			
			-- this is going to be a single value for each leaf
			-- that means the destination will be the number of nodes it takes to cover the grid (excluding the border)
			-- however, do I want this to be a larger buffer, and then perform reduce on it?
			solver:clalloc('amrErrorBuf', 
				-- self.volume 
				tonumber(self.amrRootSizeInFromSize:volume())
				* ffi.sizeof(app.real))
		end
	end

	self:getSolverCode()
end
	
function AMR:getConsLRTypeCode() return '' end


-- also found in solver/twofluid-emhd-separate-behavior.lua
function AMR:callAll(name, ...)
	local args = setmetatable({...}, table)
	args.n = select('#',...)
	return self.solvers:map(function(solver)
		return solver[name](solver, args:unpack(1, args.n))
	end):unpack()
end

function AMR:getSolverCode()
	return template(({
		['dt vs 2dt'] = [[
kernel void compareUvsU2(
	global <?=eqn.cons_t?>* U2Buf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?> *U2 = U2Buf + index;
	const global <?=eqn.cons_t?> *U = UBuf + index;
	
	//what to use to compare values ...
	//if we combine all primitives, they'll have to be appropriately weighted ...
	real sum = 0.;
	real tmp;
<? for i=0,eqn.numStates-1 do
?>	tmp = U2->ptr[<?=i?>] - U->ptr[<?=i?>]; sum += tmp * tmp;
<? end
?>	U2->ptr[0] = sum * 1e+5;
}
]],
		gradient = [==[
kernel void calcAMRError(
	global real* amrErrorBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	int4 nodei = globalInt4();
	if (nodei.x >= <?=solver.amrRootSizeInFromSize.x?> || 
		nodei.y >= <?=solver.amrRootSizeInFromSize.y?>) 
	{
		return;
	}

	int nodeIndex = nodei.x + <?=solver.amrRootSizeInFromSize.x?> * nodei.y;

	real dV_dx;	
	real sum = 0.;
	
	//hmm, it's less memory, but it's probably slower to iterate across all values as I build them here
	for (int nx = 0; nx < <?=solver.amrNodeFromSize.x?>; ++nx) {
		for (int ny = 0; ny < <?=solver.amrNodeFromSize.y?>; ++ny) {
			int4 Ui = (int4)(0,0,0,0);
			
			Ui.x = nodei.x * <?=solver.amrNodeFromSize.x?> + nx + numGhost;
			Ui.y = nodei.y * <?=solver.amrNodeFromSize.y?> + ny + numGhost;
			
			int Uindex = INDEXV(Ui);
			const global <?=eqn.cons_t?>* U = UBuf + Uindex;
				
	//TODO this wasn't the exact formula ...
	// and TODO make this modular.  some papers use velocity vector instead of density.  
	// why not total energy -- that incorporates everything?
<? for i=0,solver.dim-1 do
?>			dV_dx = (U[stepsize.s<?=i?>].rho - U[-stepsize.s<?=i?>].rho) / (2. * grid_dx<?=i?>);
			sum += dV_dx * dV_dx;
<? end
?>	
		}
	}
	amrErrorBuf[nodeIndex] = sum * 1e-2 * <?=clnumber(1/tonumber( solver.amrNodeFromSize:volume() ))?>;
}

//from is the position on the root level to read from
//to is which node to copy into
kernel void initNodeFromRoot(
	global <?=eqn.cons_t?>* UBuf,
	int4 from,
	int toNodeIndex
) {
	int4 i = (int4)(0,0,0,0);
	i.x = get_global_id(0);
	i.y = get_global_id(1);
	int dstIndex = i.x + numGhost + <?=solver.amrNodeSize.x?> * (i.y + numGhost);
	int srcIndex = from.x + (i.x>>1) + numGhost + gridSize_x * (from.y + (i.y>>1) + numGhost);

	global <?=eqn.cons_t?>* dstU = UBuf + <?=solver.numCells?> + toNodeIndex * <?=solver.amrNodeSize:volume()?>;
	
	//blitter srcU sized solver.amrNodeFromSize (in a patch of size solver.gridSize)
	// to dstU sized solver.amrNodeSize (in a patch of solver.amrNodeSize)
	
	dstU[dstIndex] = UBuf[srcIndex];
}
]==],
	})[self.amrMethod] or '', {
		solver = self,
		eqn = self.eqn,
		clnumber = clnumber,
	})
end

-- parallel in design but separate of the individual solvers's refreshSolverProgram
function AMR:refreshSolverProgram()
	local code = self:getSolverCode()
	
	time('compiling solver program', function()
		self.solverProgramObj = self.Program{code=code}
		self.solverProgramObj:compile()
	end)
	
	if self.amrMethod == 'dt vs 2dt' then
		self.compareUvsU2KernelObj = self.solverProgramObj:kernel('compareUvsU2', self.U2Buf, self.UBuf)
	elseif self.amrMethod == 'gradient' then
		self.calcAMRErrorKernelObj = self.solverProgramObj:kernel('calcAMRError', self.amrErrorBuf, self.UBuf)
		self.initNodeFromRootKernelObj = self.solverProgramObj:kernel('initNodeFromRoot', self.UBuf)
	end
end

function AMR:addDisplayVars()
	if self.amrMethod == 'dt vs 2dt' then
		self:addDisplayVarGroup{
			name = 'U2',
			type = self.eqn.cons_t,
			vars = {
				{[0] = '*value = buf[index].ptr[0];'},
			}
		}
	elseif self.amrMethod == 'gradient' then
		self:addDisplayVarGroup{
			name = 'amrError',
			type = 'real',
			vars = {
				{[0] = '*value = buf[index];'},
			}
		}
	end
end


function AMR:resetState()
	self:callAll'resetState'
	self.t = self.ion.t
end

function AMR:calcDT()
	return math.min(self:callAll'calcDT')
end

function AMR:step(dt)
	self:callAll('step', dt)
end

local SolverBase = require 'solver.solverbase'
function AMR:update()
	-- do fps / nan checks
	for _,solver in ipairs(self.solvers) do
		SolverBase.update(solver)
	end
	local dt = self:calcDT()
	
	if self.amrMethod == 'dt vs 2dt' then
		-- back up the last buffer
		self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.lastUBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
	end
	
	self:step(dt)

--[[ TODO here
for dt based, do two half steps and compare buffers
for gradient based, calc max gradient
if either signals to split then we do so..
by creating a new solver?
but I don't need to duplicate kernels ...
or maybe I do, so long as my grid sizes are hard-coded ...
but that would mean building new code for each new node in the tree
so to get around this, the AMR is going to require me to redesign the cell_x() CL function ...
but how could I possibly vary cell_x()'s bounds without regen, unless I changed some pointer passed to the kernels
and this is starting to sound like I'll have to now pass an extra global pointer to each kernel
oh well, unstructured grids needed it too, it was bound to happen
--]]
	-- now copy it to the backup buffer
	if self.amrMethod == 'dt vs 2dt' then
		-- TODO have step() provide a target, and just update directly into U2Buf?
		self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.U2Buf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
		self.app.cmds:enqueueCopyBuffer{src=self.lastUBuf, dst=self.UBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}

		local t = self.t
		self:step(.5 * dt)
		self.t = t + .5 * dt
		self:step(.5 * dt)
		self.t = t + dt

		-- now compare UBuf and U2Buf, store in U2Buf in the first real of cons_t
		self.compareUvsU2KernelObj()
	elseif self.amrMethod == 'gradient' then
		
		-- 1) compute errors from gradient, sum up errors in each root node, and output on a per-node basis
		local amrRootSizeInFromGlobalSize = vec3sz(
			roundup(self.amrRootSizeInFromSize.x, self.localSize.x),
			roundup(self.amrRootSizeInFromSize.y, self.localSize.y),
			roundup(self.amrRootSizeInFromSize.z, self.localSize.z))
		
		self.app.cmds:enqueueNDRangeKernel{
			kernel = self.calcAMRErrorKernelObj.obj, 
			dim = self.dim, 
			globalSize = amrRootSizeInFromGlobalSize:ptr(), 
			localSize = self.localSize:ptr(),
		}

		-- 2) based on what nodes' errors are past some value, split or merge...
		--[[
		1) initial tree will have nothing flagged as split
		2) then we get some split data - gradients -> errors -> thresholds -> flags 
			... which are lined up with the layout of the patches ...
			... which doesn't necessarily match the tree structure ...
		3) look through all used patches' error thresholds, min and max
			if it says to split ... 
				then look and see if we have room for any more free leafs in our state buffer
			
				the first iteration will request to split on some cells
				so go through the error buffer for each (root?) node,
				see if the error is bigger than some threshold then this node needs to be split
					then we have to add a new leaf node
				
				so i have to hold a table of what in the U extra leaf buffer is used
				which means looking
			
			if it says to merge ...
				clear the 'used' flag in the overall tree / in the layout of leafs in our state buffer
		--]]
		local vol = tonumber(self.amrRootSizeInFromSize:volume())
		local ptr = ffi.new('real[?]', vol)
		self.app.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * vol, ptr=ptr}
	-- [[
	print'amrErrors:'
	for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
		for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
			local i = nx + self.amrRootSizeInFromSize.x * ny
			io.write('\t', ('%.5f'):format(ptr[i]))
		end
		print()
	end
	--]]

	for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
		for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
			local i = nx + self.amrRootSizeInFromSize.x * ny
			local nodeErr = ptr[i]
			if nodeErr > .2 then
				print('root node '..tostring(i)..' needs to be split')
				
				-- flag for a split
				-- look for a free node to allocate in the buffer
				-- if there's one available then ...
				-- store it in a map

				local amrLayer = self.amrLayers[1]
				
				-- see if there's an entry in this layer
				-- if there's not then ...
				-- allocate a new patch and make an entry
				if not amrLayer[i+1] then
				
					-- next: find a new unused leaf
					-- for now, just this one node
					local leafIndex = 0
					
					if not self.amrLeafs[leafIndex+1] then
				
						print('splitting root node '..tostring(i)..' and putting in leaf node '..leafIndex)
			
						-- create info about the leaf
						self.amrLeafs[leafIndex+1] = {
							level = 0,	-- root
							layer = amrLayer,
							layerX = nx,		-- node x and y in the root
							layerY = ny,
							leafIndex = leafIndex,	-- which leaf we are using
						}
				
						-- tell the root layer table which node is used
						--  by pointing it back to the table of the leaf nodes 
						amrLayer[i+1] = self.amrLeafs[1]
					
						-- copy data from the root node location into the new node
						-- upsample as we go ... by nearest?
					
						-- TODO setup kernel args
						self.initNodeFromRootKernelObj.obj:setArg(1, ffi.new('int[4]', {nx, ny, 0, 0}))
						self.initNodeFromRootKernelObj.obj:setArg(2, ffi.new('int[1]', 0))
						self.app.cmds:enqueueNDRangeKernel{
							kernel = self.initNodeFromRootKernelObj.obj,
							dim = self.dim, 
							globalSize = self.amrNodeSizeWithoutBorder:ptr(),
							localSize = self.amrNodeSizeWithoutBorder:ptr(),
						}
					end
				end
			end
		end
	end
	
	self.t = self.t + dt
	self.dt = dt

--END OF THE CODE I HAVE TO FIX
	
	
	for _,solver in ipairs(self.solvers) do
		solver.t = solver.t + dt
	end
	
	self.t = self.solvers[1].t
end

function AMR:getTex(var) 
	return self.solverForDisplayVars[var].tex
end

function AMR:getHeatMap2DShader(var) 
	return self.solverForDisplayVars[var].heatMap2DShader
end

function AMR:calcDisplayVarToTex(var)
	return self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
end

function AMR:calcDisplayVarRange(var)
	return self.solverForDisplayVars[var]:calcDisplayVarRange(var)
end

function AMR:initDraw()
	self:callAll'initDraw'
end

function AMR:updateGUI()
	for i,solver in ipairs(self.solvers) do
		ig.igPushIDStr('subsolver '..i)
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			solver:updateGUI()
		end
		ig.igPopID()
	end
end

return AMR
