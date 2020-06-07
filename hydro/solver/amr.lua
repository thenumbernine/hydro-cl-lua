--[[
behavior that makes any GridSolver AMR-friendly

now to change hydro/solver/gridsolver so I can somehow modify the mins/maxs without reloading any kernels
...
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local vec3sz = require 'vec-ffi.vec3sz'
local roundup = require 'hydro.util.roundup'

local int4ptr = ffi.new'int[4]'
local function int4(a,b,c,d)
	int4ptr[0] = a
	int4ptr[1] = b
	int4ptr[2] = c
	int4ptr[3] = d
	return int4ptr
end

local function preInitAMR(self, args)
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
				splitThreshold = .3,
				mergeThreshold = .0001,
			
				-- the size, in cells, which a node replaces
				nodeFromSize = ({
					vec3sz(32, 1, 1),
					vec3sz(32, 32, 1),
					vec3sz(32, 32, 32),
				})[self.dim],

				nodeSizeWithoutBorder = vec3sz(self.sizeWithoutBorder),
			},
		}
		
		-- size of the root level, in terms of nodes ('from' size)
		self.amr.ctx.parentSizeInFromSize = vec3sz(1,1,1)
		for i=0,self.dim-1 do
			self.amr.ctx.parentSizeInFromSize.s[i] = 
				roundup(self.sizeWithoutBorder.s[i], self.amr.ctx.nodeFromSize.s[i]) 
					/ self.amr.ctx.nodeFromSize.s[i]
		end

-- number of cells in the node. 
-- I'm matching this with the root size.  if you want to change it then be prepared to change gridSize to be members of solver_t
print('nodeSizeWithoutBorder', self.amr.ctx.nodeSizeWithoutBorder)

-- number of cells that the node replaces
print('nodeFromSize', self.amr.ctx.nodeFromSize)

-- how many nodes cover the parent
print('parentSizeInFromSize', self.amr.ctx.parentSizeInFromSize)

		-- make sure that the child boundary is at least as big as one cell in the parent
		for i=0,2 do
			assert(self.numGhost >= self.amr.ctx.parentSizeInFromSize.s[i],
				require 'ext.tolua'{
					['self.numGhost'] =self.numGhost,
					['self.amr.ctx.parentSizeInFromSize.s[i]'] = self.amr.ctx.parentSizeInFromSize.s[i],
				})
		end

		local volume = tonumber(self.amr.ctx.parentSizeInFromSize:volume())
		self.amrErrorPtr = ffi.new('real[?]', volume)
	end

	self.initArgs = table(args)
end

local function createBuffersAMR(self)
	if self.amr.ctx.method == 'gradient' then

		-- how big each node is
		self.amr.nodeSize = self.amr.ctx.nodeSizeWithoutBorder + 2 * self.numGhost

		self.amr.child = table()
	end
		
	-- TODO UBufSize used to go here

	--[[ this used to go after createBuffers UBufSize
	if self.amr.ctx.method == 'gradient' then	
		UBufSize = UBufSize + self.amrMaxNodes * self.amr.nodeSize:volume()
	end	
	--]]

	if self.amr.ctx.method == 'dt vs 2dt' then
		-- here's my start at AMR, using the 1989 Berger, Collela two-small-steps vs one-big-step method
		self:clalloc('lastUBuf', self.eqn.cons_t, self.numCells)
		self:clalloc('U2Buf', self.eqn.cons_t, self.numCells)
	elseif self.amr.ctx.method == 'gradient' then
		-- this is going to be a single value for each leaf
		-- that means the destination will be the number of nodes it takes to cover the grid (excluding the border)
		-- however, do I want this to be a larger buffer, and then perform reduce on it?
		self:clalloc('amrErrorBuf', 'real',
			-- self.volume 
			tonumber(self.amr.ctx.parentSizeInFromSize:volume()),
			-- here is our one use of buffer 'sizevec':
			assert(self.amr.ctx.parentSizeInFromSize)
		)
	end
end

return function(cl)
	cl = class(cl)

	local subcl

	function cl:init(args)
		-- overriding the gridsize ...
		args = table(args)
		args.gridSize = {64,64,1}
		cl.super.init(self, args)
	end

	function cl:preInit(args)
		cl.super.preInit(self, args)
		preInitAMR(self, args)
	end

	function cl:createBuffers()
		cl.super.createBuffers(self)
		createBuffersAMR(self)
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
	constant <?=solver.solver_t?>* solver,	//parent node's solver_t
	global real* amrErrorBuf,
	const global <?=eqn.cons_t?>* UBuf		//parent node's UBuf
) {
	int4 nodei = globalInt4();
	if (nodei.x >= <?=solver.amr.ctx.parentSizeInFromSize.x?> || 
		nodei.y >= <?=solver.amr.ctx.parentSizeInFromSize.y?>) 
	{
		return;
	}

	int nodeIndex = nodei.x + <?=solver.amr.ctx.parentSizeInFromSize.x?> * nodei.y;

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
		
		int4 Ui = (int4)(
			numGhost + nx + <?=solver.amr.ctx.nodeFromSize.x?> * nodei.x,
			numGhost + ny + <?=solver.amr.ctx.nodeFromSize.y?> * nodei.y,
			0,0);
		
		
		int Uindex = INDEXV(Ui);
		const global <?=eqn.cons_t?>* U = UBuf + Uindex;
				
	//TODO this wasn't the exact formula ...
	// and TODO make this modular.  some papers use velocity vector instead of density.  
	// why not total energy -- that incorporates everything?
<? for i=0,solver.dim-1 do
?>		dV_dx = (U[solver->stepsize.s<?=i?>].rho - U[-solver->stepsize.s<?=i?>].rho) / (2. * solver->grid_dx.s<?=i?>);
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
	if (i.x >= <?=solver.amr.nodeSize.x?> || 
		i.y >= <?=solver.amr.nodeSize.y?>) 
	{
		return;
	}
	
	int dstIndex = i.x + <?=solver.amr.nodeSize.x?> * i.y;

	//'srci' is the coords within the parent node to read, relative to the child's upper-left
	int4 srci = (int4)(0,0,0,0);
	srci.x = i.x / <?= solver.amr.ctx.parentSizeInFromSize.x ?>;
	srci.y = i.y / <?= solver.amr.ctx.parentSizeInFromSize.y ?>;

	int srcIndex = numGhost + srci.x + from.x * <?= solver.amr.ctx.parentSizeInFromSize.x ?>
		+ gridSize_x * (
			numGhost + srci.y + from.y * <?= solver.amr.ctx.parentSizeInFromSize.y ?>
		);

	//blitter srcU sized solver.amr.ctx.nodeFromSize (in a patch of size solver.gridSize)
	// to dstU sized solver.amr.nodeSize (in a patch of solver.amr.nodeSize)
	
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
			self.calcAMRErrorKernelObj = self.solverProgramObj:kernel'calcAMRError' 
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
				bufferType = self.eqn.cons_t,
				vars = {
					{name='0', code='value.vreal = buf[index].ptr[0];'},
				}
			}
		elseif self.amr.ctx.method == 'gradient' then
			self:addDisplayVarGroup{
				name = 'amrError',
				bufferField = 'amrErrorBuf',
				bufferType = 'real',
				vars = {
					{name='0', code='value.vreal = buf[index];'},
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
			self.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.lastUBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
		end
		
		-- update children ... twice as many times, at half the timestep
		local childDT = self.fixedDT * .5
		for _,child in pairs(self.amr.child) do
			child.fixedDT = childDT
			child.useFixedDT = true
			child:update()
		end

		-- now copy it to the backup buffer
		if self.amr.ctx.method == 'dt vs 2dt' then
			-- TODO have step() provide a target, and just update directly into U2Buf?
			self.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.U2Buf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
			self.cmds:enqueueCopyBuffer{src=self.lastUBuf, dst=self.UBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}

			self:step(.5 * dt)
			self.t = t + .5 * dt
			self:step(.5 * dt)
			self.t = t + dt

			-- now compare UBuf and U2Buf, store in U2Buf in the first real of cons_t
			self.compareUvsU2KernelObj()
		elseif self.amr.ctx.method == 'gradient' then
			
			-- 1) compute errors from gradient, sum up errors in each root node, and output on a per-node basis
			local amrRootSizeInFromGlobalSize = vec3sz(
				roundup(self.amr.ctx.parentSizeInFromSize.x, self.localSize.x),
				roundup(self.amr.ctx.parentSizeInFromSize.y, self.localSize.y),
				roundup(self.amr.ctx.parentSizeInFromSize.z, self.localSize.z))
--print('self.amr.ctx.parentSizeInFromSize', self.amr.ctx.parentSizeInFromSize) 
--print('amrRootSizeInFromGlobalSize', amrRootSizeInFromGlobalSize) 
--print('self.localSize', self.localSize)			
			self.calcAMRErrorKernelObj.obj:setArgs(self.solverBuf, self.amrErrorBuf, self.UBuf)
			self.cmds:enqueueNDRangeKernel{
				kernel = self.calcAMRErrorKernelObj.obj, 
				dim = self.dim, 
				globalSize = amrRootSizeInFromGlobalSize.s, 
				localSize = self.localSize.s,
			}

			local volume = tonumber(self.amr.ctx.parentSizeInFromSize:volume())
			local ptr = assert(self.amrErrorPtr)
			self.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * volume, ptr=ptr}
		
			-- [[
print('depth '..self.amr.depth..' amrErrors:')
			for ny=0,tonumber(self.amr.ctx.parentSizeInFromSize.y)-1 do
				for nx=0,tonumber(self.amr.ctx.parentSizeInFromSize.x)-1 do
					local i = nx + self.amr.ctx.parentSizeInFromSize.x * ny
					io.write('\t', ('%.5f'):format(ptr[i]))
				end
				print()
			end
			--]]

			-- [[
			if self.amr.depth < self.amr.ctx.maxdepth then
				for ny=0,tonumber(self.amr.ctx.parentSizeInFromSize.y)-1 do
					for nx=0,tonumber(self.amr.ctx.parentSizeInFromSize.x)-1 do
						local i = tonumber(nx + self.amr.ctx.parentSizeInFromSize.x * ny)
						local nodeErr = ptr[i]
						if nodeErr > self.amr.ctx.splitThreshold then
							if not self.amr.child[i+1] then
print("creating depth "..tonumber(self.amr.depth).." child "..tonumber(i))
								-- [==[
								-- lazy way: just recreate the whole thing
								-- except remap the mins/maxs
								local dx = 1/tonumber(self.amr.ctx.parentSizeInFromSize.x)
								local dy = 1/tonumber(self.amr.ctx.parentSizeInFromSize.y)
								local ux = nx*dx
								local uy = ny*dy
							
								local newmins = {
									(self.maxs.x - self.mins.x) * nx * dx + self.mins.x,
									(self.maxs.y - self.mins.y) * ny * dy + self.mins.y,
									0}
								local newmaxs = {
									(self.maxs.x - self.mins.x) * (nx+1) * dx + self.mins.x,
									(self.maxs.y - self.mins.y) * (ny+1) * dy + self.mins.y,
									0}	
					
-- when running the following, upon the next update, I start to get nans in the amr error buffer
								local subsolver = subcl{
									app = self.app,
									parent = self,
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
								self.initNodeFromRootKernelObj.obj:setArgs(subsolver.UBuf, self.UBuf, int4(nx,ny,0,0))
								-- so long as node size = root size, we can use the solver.globalSize, and not --self.amr.ctx.nodeSizeWithoutBorder.s,
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

	end

	function cl:resetState()
		cl.super.resetState(self)

		-- TODO dealloc safely?
		self.amr.child = table()
	end



	-- maybe I shouldn't make a subclass
	-- maybe a patchlevel class is a better idea
	subcl = class(cl)

	function subcl:preInit(args)
		self.parent = assert(args.parent)
		self.maxWorkGroupSize = self.parent.maxWorkGroupSize
		local sizeProps = self:getSizePropsForWorkGroupSize(self.maxWorkGroupSize)
		for k,v in pairs(sizeProps) do
			self[k] = v
		end
	
		self:createEqn()
		
		--subcl.super.preInit(self, args)
		self.solver_t = self.parent.solver_t
		self.solverPtr = ffi.new(self.solver_t)
	
		preInitAMR(self, args)
	end
	
	function subcl:createEqn()
		self.eqn = self.parent.eqn
	end

	function subcl:refreshEqnInitState() 
		-- don't need to create init state ...
		-- refreshCodePrefix is called next
		-- and that calls refreshIntegrator ...
		self:refreshCodePrefix()
	end

	-- refreshCodePrefix calls:
	function subcl:createCodePrefix() end
	function subcl:refreshIntegrator()
		self.integrator = self.parent.integrator
	end
	function subcl:refreshInitStateProgram() end

	local template = require 'template'
	function subcl:refreshBoundaryProgram() 
--[=[
		self.boundaryProgramObj, self.boundaryKernelObjs = self:createBoundaryProgramAndKernel(
			table(
				self:getBoundaryProgramArgs(),	-- should have the buffer and type
				{
					methods = {
						xmin = 'freeflow',
						xmax = 'freeflow',
						ymin = 'freeflow',
						ymax = 'freeflow',
						zmin = 'freeflow',
						zmax = 'freeflow',
					},
				}
			)
		)

--]=]
--[=[
		local function copyBorder(args)
			return template([[
]])
		end
		self.boundaryProgramObj, self.boundaryKernelObjs = self:createBoundaryProgramAndKernel(
			table(
				self:getBoundaryProgramArgs(),	-- should have the buffer and type
				methods = {
					xmin = copyBorder,
					xmax = copyBorder,
					ymin = copyBorder,
					ymax = copyBorder,
					zmin = copyBorder,
					zmax = copyBorder,
				}
			)
		)
		-- next comes the op's boundary programs
		-- next comes the CTU's boundary programs
--]=]
	end
	
	function subcl:refreshSolverProgram()
		local ks = table.keys(self.parent):sort()
		for _,k in ipairs(ks) do
			if k:match'KernelObj$'
			or k:match'KernelObjs$'
			or k:match'ProgramObj$'
			or k:match'Shader$'
			then
				self[k] = self.parent[k]
			end
		end
		
		-- still needs to call ops refreshSolverProgram (or at least steal their kernels too)
		for _,op in ipairs(self.ops) do
			op:refreshSolverProgram()
		end
		-- and the display vars ...
	end
	
	-- called by refreshGridSize:
	function subcl:refreshCommonProgram() end
	function subcl:resetState() end

	function subcl:createBuffers()
		local ks = table.keys(self.parent):sort()
		for _,k in ipairs(ks) do
			-- as long as we have matching size, we can reuse all buffers (except the state buffers)
			-- TODO will need to save primBufs as well
			if (k:match'Buf$' and k ~= 'UBuf' and k ~= 'solverBuf')
			or (k:match'BufObj$' and k ~= 'UBufObj')
			then
				self[k] = self.parent[k]
			end
		end
		self.tex = self.parent.tex
		self.texCLMem = self.parent.texCLMem
		self.calcDisplayVarToTexPtr = self.parent.calcDisplayVarToTexPtr
		self.amrErrorPtr = self.parent.amrErrorPtr
		
		-- now for our own allocations...
		-- this is copied from GridSolver:createBuffers
print('creating subsolver ubuffer...')									
		self:clalloc('UBuf', self.eqn.cons_t, self.numCells)
		
		createBuffersAMR(self)
	end

	function subcl:boundary()
  	-- TODO give us a better boundary kernel, then this function can stay intact
	end

	return cl
end
