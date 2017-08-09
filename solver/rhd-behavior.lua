--[[
behavior for allocating a primbuf
and updating it using newton descent every iteration
and passing it as an extra arg to all those functions that need it

this is different from Euler in that Euler doesn't hold a prim buf
but maybe it should -- because this is running faster than Euler
--]]

--[[
instead of keeping two separate buffers
 and always passing two functions together
 and overriding *every single thing* that uses the UBuf
I'm just going to combine the structures
and only integrate over the first 5 terms

because for RK4 I need to save the primBuf as well as the consBuf ...
but this forces the integrators to allocate numStates instead of just numIntStates...

maybe it is the integrators that need to be overridden...
or maybe numIntStates should include the primBuf as well?

so the RK4 on SRHD has to push both UBuf and primBuf upon pushing the state
and after each update it has to re-converge the primBuf

ok so there's a few ways we can go about this ...
0) like we're doing, which is wrong
1) combine UBuf and primBuf, so when the RK pushes/pops state the UBuf and primBuf are together
	the downside is the # of allocations 
	and the necessity of rk4 to push non-integratable data 
	which it has to do anyways, right, so it's not that big of a risk ... right?
	do buffers with non-integratable variables need those variables in calcDeriv?  yes...?
	so this method would work.  and require less setArgs() to be done everywhere.
	(more reworking, but better in the end? and more memory used by RK4)
2) converge primBuf at the start of calcDeriv 
	so it hopefully matches up with U's values
	but this means it will be converging *from* different U's depending on what step in RK4 has been run,
	which means it will convert *to* possibly different prim values even for same U's ... which could be disastrous
	(the easy way)
--]]

local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'

return function(parent)
	local template = class(parent)

	function template:createBuffers()
		template.super.createBuffers(self)
		self:clalloc('primBuf', self.volume * self.dim * ffi.sizeof(self.eqn.prim_t))
	end

	template.ConvertToTex_U = class(template.ConvertToTex_U)
	function template.ConvertToTex_U:setArgs(kernel, var)
		kernel:setArg(1, self.solver.UBuf)
		kernel:setArg(2, self.solver.primBuf)
	end

	-- replace the U convertToTex with some custom code 
	function template:getAddConvertToTexUBufArgs()
		return table(template.super.getAddConvertToTexUBufArgs(self),
			{extraArgs = {'const global '..self.eqn.prim_t..'* primBuf'}})
	end

	function template:addConvertToTexs()
		template.super.addConvertToTexs(self)

		self:addConvertToTex{
			name = 'prim', 
			type = self.eqn.prim_t,
			varCodePrefix = self.eqn:getPrimDisplayVarCodePrefix(),
			vars = self.eqn.primDisplayVars,
		}
	end

	function template:refreshInitStateProgram()
		template.super.refreshInitStateProgram(self)
		self.initStateKernel:setArg(1, self.primBuf)
	end

	--[[
	calcDT, calcEigenBasis use primBuf
	so for the Roe implicit linearized solver,
	(TODO?) primBuf must be push/pop'd as well as UBuf
	(or am I safe just using the last iteration's prim values, and doing the newton descent to update them?)
	--]]
	function template:refreshSolverProgram()
		-- createKernels in particular ...
		template.super.refreshSolverProgram(self)

		self.calcDTKernel:setArg(1, self.primBuf)
		self.calcEigenBasisKernel:setArg(2, self.primBuf)
		
		-- grhd has one of these, srhd doesn't
		if self.addSourceKernel then
			self.addSourceKernel:setArg(2, self.primBuf)
		end

		self.updatePrimsKernel = self.solverProgram:kernel('updatePrims', self.primBuf, self.UBuf)
	end

	--[[ method 1: update prims after step() overall is called
	-- this leaves bad prims throughout RK4
	function template:step(dt)
		template.super.step(self, dt)

		self.app.cmds:enqueueNDRangeKernel{kernel=self.updatePrimsKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	end
	--]]
	-- [[ method 2: update the before calcDeriv
	-- this might take some extra calculations
	-- and prims could converge *from* different prevoius values even when converging *to* the same destination UBuf 
	function template:calcDeriv(derivBuf, dt)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.updatePrimsKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}

		template.super.calcDeriv(self, derivBuf, dt)
	end
	--]]
	
	function template:boundary()
		-- U boundary
		template.super.boundary(self)
		-- prim boundary
		self.boundaryKernel:setArg(0, self.primBuf)
		self:applyBoundaryToBuffer(self.boundaryKernel)
		self.boundaryKernel:setArg(0, self.UBuf)
	end

	return template
end
