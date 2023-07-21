--[[
right now nodiv is a collection of a few classes, depending on what underlying solver you use
for the jacobi version:
nodiv.lua
	-> poisson_jacobi.lua
		-> poisson.clcpp (shared with poisson_krylov)
		-> poisson_jacobi.clcpp
		-> relaxation.lua

for the krylov version:
nodiv.lua
	-> poisson.clcpp (shared with poisson_jacobi)
	-> solver/cl/$solverClass.lua

what do the classes fulfill externally?
- their :step(dt) function removes divergence from the field specified in cmdline args.

what goes on within the jacobi version?
NoDiv:step()
	call NoDiv:relax()
		Relaxation:relax()
		call solveJacobiKernelObj()
		call copyWriteToPotentialNoGhsotKernelObj
		call self:potentialBoundary()
		do stop on epsilon stuff maybe
--]]
local class = require 'ext.class'
local path = require 'ext.path'
local table = require 'ext.table'
local CLBuffer = require 'cl.obj.buffer'

local NoDiv = class()

NoDiv.name = 'NoDiv'
NoDiv.scalar = 'real'

-- which cons_t field to store the solved potential value in
NoDiv.vectorField = 'B'
NoDiv.potentialField = 'psi'
NoDiv.chargeField = nil	-- nil means zero
NoDiv.chargeCode = nil
NoDiv.maxIters = cmdline.selfGravPoissonMaxIter or 20	-- 10 * numreals


function NoDiv:init(args)
	local solver = assert(args.solver)
	self.solver = solver
	
	require 'hydro.code.symbols'(self, self:getSymbolFields())

	self.writeBufObj = CLBuffer{
		env = solver.app.env,
		name = 'writeBuf',
		type = solver.app.real,
		count = solver.numCells,
	}

	self.scalar = args.scalar
	self.vectorField = args.vectorField
	self.readVectorField = args.readVectorField	
	self.writeVectorField = args.writeVectorField
	self.potentialField = args.potentialField
	self.chargeField = args.chargeField
	self.chargeCode = args.chargeCode
end

function NoDiv:getSymbolFields()
	return table{
		'initPotential',
		'copyWriteToPotential',
		'copyPotentialToReduce',
		'noDiv',
		'solveJacobi',
	}
end
	
-- TODO return BufObj
function NoDiv:getPotBuf() 
	return self.solver.UBuf 
end

-- TODO return getPotBufObj().type
function NoDiv:getPotBufType() 
	return self.solver.eqn.symbols.cons_t
end

function NoDiv:readVectorField(offset)
	return 'U['..offset..'].'..self.vectorField
end

-- returns code that performs U->$v = U->$v - $dv
function NoDiv:writeVectorField(dv)
	local scalar = self.scalar
	local UField = 'U->'..self.vectorField
	local sub = scalar..'3_sub'
	return UField..' = '..sub..'('..UField..', '..dv..');'
end


-- inserted into nodiv-unstructured.cl
function NoDiv:getPoissonDivCode()
	return self.solver.eqn:template([[

	real3 const x = cell->pos;
	for (int i = 0; i < cell->faceCount; ++i) {
		global <?=face_t?> const * const face = faces + cellFaceIndexes[i + (cell)->faceOffset];
		real const dx = face->area;		// face->cellDist?
	
		int const iL = face->cells.s0;
		int const iR = face->cells.s1;

		if (iL != -1 && iR != -1) {

			real3 dx = real3_sub(cellBuf[iR].pos, cellBuf[iL].pos);
			real3 vR = <?=op:readVectorField('iR')?>;
			real3 vL = <?=op:readVectorField('iL')?>;

<? for j=0,solver.dim-1 do ?>
			source += (vR.s<?=j?> - vL.s<?=j?>) * dx.s<?=j?>;
<? end ?>

	}

<? if op.chargeField then 
?>	source += U-><?=op.chargeField?>;
<? elseif op.chargeCode then
?><?=op.chargeCode?><? 
end 
?>
]], 
		{
			op = self,
		})
end

function NoDiv:initCodeModules()
	local solver = self.solver
	solver.modules:addFromMarkup{
		code = solver.eqn:template(
			path'hydro/op/nodiv-unstructured.cl':read(),
			table(self.symbols, {op = self})
		),
	}
	solver.solverModulesEnabled[self.symbols.noDiv] = true
end

function NoDiv:refreshSolverProgram()
	NoDiv.super.refreshSolverProgram(self)
	local solver = self.solver
	

	self.copyWriteToPotentialKernelObj = solver.solverProgramObj:kernel(self.symbols.copyWriteToPotential, solver.solverBuf, solver.UBuf, self.writeBufObj)

	self.noDivKernelObj = solver.solverProgramObj:kernel(self.symbols.noDiv, solver.solverBuf, solver.UBuf)
	self.copyPotentialToReduceKernelObj = solver.solverProgramObj:kernel(self.symbols.copyPotentialToReduce, solver.solverBuf, solver.reduceBuf, solver.UBuf)
	self.solveJacobiKernelObj = solver.solverProgramObj:kernel(self.symbols.solveJacobi, solver.solverBuf, solver.reduceBuf, solver.UBuf)
end

function NoDiv:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initPotentialKernelObj()
	self:relax()
end

function NoDiv:relax()
	local solver = self.solver
	for i=1,self.maxIters do
		self.lastIter = i
	
		-- writes new potentialField to writeBuf
		-- stores deltas in reduceBuf
		self.solveJacobiKernelObj()

		-- copy new values back from writeBuf to UBuf potentialField
		-- copies square of increment into reduceBuf
		self.copyWriteToPotentialKernelObj()

		-- apply boundary to UBuf.potentialField
		self:potentialBoundary()

		if self.stopOnEpsilon then
			local residual = math.sqrt(fromreal(solver.reduceSum()) / tonumber(solver.volumeWithoutBorder))
			local lastResidual = self.lastResidual	
			self.lastResidual = residual
			if self.verbose then
				self.setReduceToPotentialSquaredKernelObj()
				local xNorm = math.sqrt(fromreal(solver.reduceSum()) / tonumber(solver.volumeWithoutBorder))
				self.copyPotentialToReduceKernelObj()
				local xmin = fromreal(solver.reduceMin())
				self.copyPotentialToReduceKernelObj()
				local xmax = fromreal(solver.reduceMax())
				io.stderr:write(table{i-1, residual, xNorm, xmin, xmax}:mapi(tostring):concat'\t','\n')
			end
			
			-- TODO compare residual
			if math.abs(residual) <= self.stopEpsilon then break end
		end
	end
end


function NoDiv:step(dt)
	local solver = self.solver
	self:relax()
	self.noDivKernelObj()
end

return NoDiv
