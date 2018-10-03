--[[
behavior that makes any GridSolver AMR-friendly

now to change solver/gridsolver so I can somehow modify the mins/maxs without reloading any kernels
...
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local vec3sz = require 'ffi.vec.vec3sz'
local roundup = require 'roundup'

local int4ptr = ffi.new'int[4]'
local function int4(a,b,c,d)
	int4ptr[0] = a
	int4ptr[1] = b
	int4ptr[2] = c
	int4ptr[3] = d
	return int4ptr
end

return function(cl)
	cl = class(cl)

	function cl:preInit(args)
		cl.super.preInit(self, args)
		
		self.amr = args.amr
		if not self.amr then
			-- in this case, we are the root:
			self.amr = {
				depth = 0,
				-- shared by all in the tree
				ctx = {
					--method = 'dt vs 2dt',
					method = 'gradient',
					maxdepth = 1,
					splitThreshold = .001,
					mergeThreshold = .0001,
				
					-- the size, in cells, which a node replaces
					nodeFromSize = ({
						vec3sz(8, 1, 1),
						vec3sz(8, 8, 1),
						vec3sz(8, 8, 8),
					})[self.dim],

					-- the size, in cells, of each node, excluding border, for each dimension
					--[[ I'm going to cheat and make this equal to the root node size
					-- this way I don't have to move the gridSize into the solver_t, and can reuse the same CL code for all nodes
					nodeSizeWithoutBorder = ({
						vec3sz(16, 1, 1),
						vec3sz(16, 16, 1),
						vec3sz(16, 16, 16),
					})[self.dim],
					--]]
					nodeSizeWithoutBorder = vec3sz(self.sizeWithoutBorder),
				},
			}
			
			-- size of the root level, in terms of nodes ('from' size)
			self.amr.ctx.rootSizeInFromSize = vec3sz(1,1,1)
			for i=0,self.dim-1 do
				self.amr.ctx.rootSizeInFromSize:ptr()[i] = 
					roundup(self.sizeWithoutBorder:ptr()[i], self.amr.ctx.nodeFromSize:ptr()[i]) 
						/ self.amr.ctx.nodeFromSize:ptr()[i]
			end
		end

		self.initArgs = table(args)
	end

	function cl:createBuffers()
		cl.super.createBuffers(self)

		if self.amr.ctx.method == 'gradient' then
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



			-- how big each node is
			self.amrNodeSize = self.amr.ctx.nodeSizeWithoutBorder + 2 * self.numGhost

			self.amr.child = table()
		end
			
		-- TODO UBufSize used to go here

		--[[ this used to go after createBuffers UBufSize
		if self.amr.ctx.method == 'gradient' then	
			UBufSize = UBufSize + self.amrMaxNodes * self.amrNodeSize:volume()
		end	
		--]]

		if self.amr.ctx.method == 'dt vs 2dt' then
			-- here's my start at AMR, using the 1989 Berger, Collela two-small-steps vs one-big-step method
			self:clalloc('lastUBuf', self.numCells * ffi.sizeof(self.eqn.cons_t))
			self:clalloc('U2Buf', self.numCells * ffi.sizeof(self.eqn.cons_t))
		elseif self.amr.ctx.method == 'gradient' then
			-- this is going to be a single value for each leaf
			-- that means the destination will be the number of nodes it takes to cover the grid (excluding the border)
			-- however, do I want this to be a larger buffer, and then perform reduce on it?
			self:clalloc('amrErrorBuf', 
				-- self.volume 
				tonumber(self.amr.ctx.rootSizeInFromSize:volume())
				* ffi.sizeof(self.app.real),
				assert(self.amr.ctx.rootSizeInFromSize))
		end
	end

	function cl:getSolverCode()
		return table{
			cl.super.getSolverCode(self),
		
			template(({
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
<? local clnumber = require 'cl.obj.number' ?>
kernel void calcAMRError(
	constant <?=solver.solver_t?>* solver,
	global real* amrErrorBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	int4 nodei = globalInt4();
	if (nodei.x >= <?=solver.amr.ctx.rootSizeInFromSize.x?> || 
		nodei.y >= <?=solver.amr.ctx.rootSizeInFromSize.y?>) 
	{
		return;
	}

	int nodeIndex = nodei.x + <?=solver.amr.ctx.rootSizeInFromSize.x?> * nodei.y;

	real dV_dx;	
	real sum = 0.;
	
//	for (int nx = 0; nx < <?=solver.amr.ctx.nodeFromSize.x?>; ++nx) {
//		for (int ny = 0; ny < <?=solver.amr.ctx.nodeFromSize.y?>; ++ny) {
	<?	-- without unrolling these for-loops, intel compiler takes 30 seconds (and produces bad asm)
for nx=0,tonumber(solver.amr.ctx.nodeFromSize.x)-1 do
	for ny=0,tonumber(solver.amr.ctx.nodeFromSize.y)-1 do
?>{	
		const int nx = <?=nx?>;
		const int ny = <?=ny?>;
			
			int4 Ui = (int4)(0,0,0,0);
			
			Ui.x = nodei.x * <?=solver.amr.ctx.nodeFromSize.x?> + nx + numGhost;
			Ui.y = nodei.y * <?=solver.amr.ctx.nodeFromSize.y?> + ny + numGhost;
			
			int Uindex = INDEXV(Ui);
			const global <?=eqn.cons_t?>* U = UBuf + Uindex;
				
	//TODO this wasn't the exact formula ...
	// and TODO make this modular.  some papers use velocity vector instead of density.  
	// why not total energy -- that incorporates everything?
<? for i=0,solver.dim-1 do
?>			dV_dx = (U[stepsize.s<?=i?>].rho - U[-stepsize.s<?=i?>].rho) / (2. * solver->grid_dx.s<?=i?>);
			sum += dV_dx * dV_dx;
<? end
?>

//		}
//	}
	}<?
	end
end
?>

	amrErrorBuf[nodeIndex] = sum * 1e-2 * <?=clnumber(1/tonumber( solver.amr.ctx.nodeFromSize:volume() ))?>;
}

kernel void initNodeFromRoot(
	global <?=eqn.cons_t?>* childUBuf,
	const global <?=eqn.cons_t?>* parentUBuf,
	int4 from	//where in the child tree
) {
	//'i' is the dest in the child node to write
	int4 i = (int4)(0,0,0,0);
	i.x = get_global_id(0);
	i.y = get_global_id(1);
	if (i.x >= <?=solver.amr.ctx.nodeSizeWithoutBorder.x?> || 
		i.y >= <?=solver.amr.ctx.nodeSizeWithoutBorder.y?>) 
	{
		return;
	}
	
	int dstIndex = i.x + numGhost + <?=solver.amrNodeSize.x?> * (i.y + numGhost);

	//'srci' is the coords within the parent node to read, relative to the child's upper-left
	int4 srci = (int4)(0,0,0,0);
	srci.x = i.x / <?= solver.amr.ctx.nodeSizeWithoutBorder.x / solver.amr.ctx.nodeFromSize.x ?>;
	srci.y = i.y / <?= solver.amr.ctx.nodeSizeWithoutBorder.y / solver.amr.ctx.nodeFromSize.y ?>;

	int srcIndex = numGhost + srci.x + from.x * <?= solver.amr.ctx.nodeSizeWithoutBorder.x / solver.amr.ctx.nodeFromSize.x ?>
		+ gridSize_x * (
			numGhost + srci.y + from.y * <?= solver.amr.ctx.nodeSizeWithoutBorder.y / solver.amr.ctx.nodeFromSize.y ?>
		);

	//blitter srcU sized solver.amr.ctx.nodeFromSize (in a patch of size solver.gridSize)
	// to dstU sized solver.amrNodeSize (in a patch of solver.amrNodeSize)
	
	childUBuf[dstIndex] = parentUBuf[srcIndex];
}
]==],
			})[self.amr.ctx.method] or '', {
				solver = self,
				eqn = self.eqn,
			}),
		}:concat'\n'
	end

	function cl:refreshSolverProgram()
		cl.super.refreshSolverProgram(self)

		if self.amr.ctx.method == 'dt vs 2dt' then
			self.compareUvsU2KernelObj = self.solverProgramObj:kernel('compareUvsU2', self.U2Buf, self.UBuf)
		elseif self.amr.ctx.method == 'gradient' then
			self.calcAMRErrorKernelObj = self.solverProgramObj:kernel('calcAMRError', self.solverBuf, self.amrErrorBuf, self.UBuf)
			self.initNodeFromRootKernelObj = self.solverProgramObj:kernel'initNodeFromRoot'
			self.initNodeFromRootKernelObj.obj:setArg(1, self.UBuf)
		end
	end

	--[[ maybe something in here is messing me up?
	-- after all, most display vars are typically all the same size
	-- the AMR ones are the only ones that are different
	-- nope...still got a crash...
	function cl:addDisplayVars()
		cl.super.addDisplayVars(self)

		if self.amr.ctx.method == 'dt vs 2dt' then
			self:addDisplayVarGroup{
				name = 'U2',
				bufferField = 'U2Buf',
				type = self.eqn.cons_t,
				vars = {
					{[0] = '*value = buf[index].ptr[0];'},
				}
			}
		elseif self.amr.ctx.method == 'gradient' then
			self:addDisplayVarGroup{
				name = 'amrError',
				bufferField = 'amrErrorBuf',
				type = 'real',
				vars = {
					{[0] = '*value = buf[index];'},
				}
			}
		end
	end
	--]]

	function cl:update()
		-- NOTICE this used to go after boundary() and before step()
		local t
		if self.amr.ctx.method == 'dt vs 2dt' then
			t = self.t
			-- back up the last buffer
			self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.lastUBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
		end
		
		cl.super.update(self)
	
		-- now copy it to the backup buffer
		if self.amr.ctx.method == 'dt vs 2dt' then
			-- TODO have step() provide a target, and just update directly into U2Buf?
			self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.U2Buf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
			self.app.cmds:enqueueCopyBuffer{src=self.lastUBuf, dst=self.UBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}

			self:step(.5 * dt)
			self.t = t + .5 * dt
			self:step(.5 * dt)
			self.t = t + dt

			-- now compare UBuf and U2Buf, store in U2Buf in the first real of cons_t
			self.compareUvsU2KernelObj()
		elseif self.amr.ctx.method == 'gradient' then
			
			-- 1) compute errors from gradient, sum up errors in each root node, and output on a per-node basis
			local amrRootSizeInFromGlobalSize = vec3sz(
				roundup(self.amr.ctx.rootSizeInFromSize.x, self.localSize.x),
				roundup(self.amr.ctx.rootSizeInFromSize.y, self.localSize.y),
				roundup(self.amr.ctx.rootSizeInFromSize.z, self.localSize.z))
--print('self.amr.ctx.rootSizeInFromSize', self.amr.ctx.rootSizeInFromSize) 
--print('amrRootSizeInFromGlobalSize', amrRootSizeInFromGlobalSize) 
--print('self.localSize', self.localSize)			
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
			local volume = tonumber(self.amr.ctx.rootSizeInFromSize:volume())
			local ptr = ffi.new('real[?]', volume)
			self.app.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * volume, ptr=ptr}
		
			-- [[
			print'amrErrors:'
			for ny=0,tonumber(self.amr.ctx.rootSizeInFromSize.y)-1 do
				for nx=0,tonumber(self.amr.ctx.rootSizeInFromSize.x)-1 do
					local i = nx + self.amr.ctx.rootSizeInFromSize.x * ny
					io.write('\t', ('%.5f'):format(ptr[i]))
				end
				print()
			end
			--]]

			-- [[
			if self.amr.depth < self.amr.ctx.maxdepth then
				for ny=0,tonumber(self.amr.ctx.rootSizeInFromSize.y)-1 do
					for nx=0,tonumber(self.amr.ctx.rootSizeInFromSize.x)-1 do
						local i = tonumber(nx + self.amr.ctx.rootSizeInFromSize.x * ny)
						local nodeErr = ptr[i]
						if nodeErr > self.amr.ctx.splitThreshold then
							if not self.amr.child[i+1] then
print("creating depth "..tonumber(self.amr.depth).." child "..tonumber(i))
								-- [==[
								-- lazy way: just recreate the whole thing
								-- except remap the mins/maxs
								local dx = 1/tonumber(self.amr.ctx.rootSizeInFromSize.x)
								local dy = 1/tonumber(self.amr.ctx.rootSizeInFromSize.y)
								local ux = nx*dx
								local uy = ny*dy
							
								local newmins = {
									(self.maxs[1] - self.mins[1]) * nx * dx + self.mins[1],
									(self.maxs[2] - self.mins[2]) * ny * dy + self.mins[2],
									1}
								local newmaxs = {
									(self.maxs[1] - self.mins[1]) * (nx+1) * dx + self.mins[1],
									(self.maxs[2] - self.mins[2]) * (ny+1) * dy + self.mins[2],
									1}	
								local _self = self
						
								-- maybe I shouldn't make a subclass
								-- maybe a patchlevel class is a better idea
								local subcl = class(cl)

								function subcl:preInit(args)
									self.maxWorkGroupSize = _self.maxWorkGroupSize
									local sizeProps = self:getSizePropsForWorkGroupSize(self.maxWorkGroupSize)
									for k,v in pairs(sizeProps) do
										self[k] = v
									end
								
									self:createEqn()
									
									--subcl.super.preInit(self, args)
									self.solver_t = _self.solver_t
								end
								
								function subcl:createEqn()
									self.eqn = _self.eqn
								end

								function subcl:createSolverBuf() 
									self.solverBuf = _self.solverBuf
								end

								function subcl:refreshEqnInitState() end
								function subcl:refreshCommonProgram() end
								function subcl:resetState() end

								function subcl:refreshSolverProgram()
									local ks = table.keys(_self):sort()
									for _,k in ipairs(ks) do
										if k:match'KernelObj$'
										or k:match'KernelObjs$'
										or k:match'ProgramObj$'
										or k:match'Shader$'
										then
											self[k] = _self[k]
										end
									end
									
									-- still needs to call ops refreshSolverProgram (or at least steal their kernels too)
									for _,op in ipairs(self.ops) do
										op:refreshSolverProgram()
									end
									-- and the display vars ...
								end

								function subcl:createBuffers()
									local ks = table.keys(_self):sort()
									for _,k in ipairs(ks) do
										-- as long as we have matching size, we can reuse all buffers (except the state buffers)
										-- TODO will need to save primBufs as well
										if (k:match'Buf$' and k ~= 'UBuf' and k ~= 'solverBuf')
										or (k:match'BufObj$' and k ~= 'UBufObj')
										then
											self[k] = _self[k]
										end
									end
									self.tex = _self.tex
									self.texCLMem = _self.texCLMem
									self.calcDisplayVarToTexPtr = _self.calcDisplayVarToTexPtr
								
									-- now for our own allocations...
									-- this is copied from GridSolver:createBuffers
print('creating subsolver ubuffer...')									
									local UBufSize = self.numCells * ffi.sizeof(self.eqn.cons_t)
									self:clalloc('UBuf', UBufSize)
								end
							
								local subsolver = subcl{
									app = self.app,
									dim = self.dim,
									gridSize = self.sizeWithoutBorder,
									coord = self.coord,
									mins = newmins,
									maxs = newmaxs,
									amr = {
										depth = self.amr.depth + 1,
										ctx = self.amr.ctx,
									},
								}
								
								-- tell the root layer table which node is used
								--  by pointing it back to the table of the leaf nodes 
								self.amr.child[i+1] = assert(subsolver)
							
								-- copy data from the root node location into the new node
								-- upsample as we go ... by nearest?
							
								-- setup kernel args
								self.initNodeFromRootKernelObj.obj:setArg(0, subsolver.UBuf)
								self.initNodeFromRootKernelObj.obj:setArg(2, int4(nx,ny,0,0))
								-- so long as node size = root size, we can use the solver.globalSize, and not --self.amr.ctx.nodeSizeWithoutBorder:ptr(),
								self.initNodeFromRootKernelObj()
								--]==]
							end
						elseif nodeErr < self.amr.ctx.mergeThreshold then
--							local subsolver= self.amr.child[i+1]
--print('node '..nx..','..ny..' needs to be merged')						
--							self.amr.child[i+1] = nil
						end
					end
				end
			end
			--]]
		end

cl.update = cl.super.update
	end

	return cl
end
