local class = require 'ext.class'
local table = require 'ext.table'
local GridSolver = require 'hydro.solver.gridsolver'

local EinsteinFiniteDifferenceSolver = class(GridSolver)

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
EinsteinFiniteDifferenceSolver.numGhost = 3

EinsteinFiniteDifferenceSolver.name = 'EinsteinFiniteDifference'

function EinsteinFiniteDifferenceSolver:init(...)
	EinsteinFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function EinsteinFiniteDifferenceSolver:refreshSolverProgram()
	EinsteinFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(0, assert(self.solverBuf))
	self.calcDerivKernelObj.obj:setArg(2, assert(self.UBuf))
	self.calcDerivKernelObj.obj.setArg(3, assert(self.cellBuf))
end

function EinsteinFiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
-- on my Intel HD 520 I have to set these arguments twice.  only in einstein-fd, not in fvsolver.  otherwise I get CL_INVALID_KERNEL_ARGS
-- on NVIDIA GTX 1080 Ti I'm getting a CL_INVALID_KERNEL_ARGS error here even after setting it again.
-- https://stackoverflow.com/questions/20562637/opencl-code-works-on-a-machine-but-i-am-getting-cl-invalid-kernel-args-on-anothe/20566270#20566270
-- ... says it could be an issue with constant-allocated data being limited by CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE 
-- on linux for my Neo drivers CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE is 1601060864 bytes, ~1.5 GB, which is equal to the max mem alloc size, and half of the global memory
self.calcDerivKernelObj.obj:setArg(0, assert(self.solverBuf))
self.calcDerivKernelObj.obj:setArg(2, assert(self.UBuf))
self.calcDerivKernelObj.obj.setArg(3, assert(self.cellBuf))
	
	self.calcDerivKernelObj.obj:setArg(1, assert(derivBufObj.obj))
	self.calcDerivKernelObj()
end

-- only set these for certain types ... 
function EinsteinFiniteDifferenceSolver:createDisplayComponents()
	EinsteinFiniteDifferenceSolver.super.createDisplayComponents(self)
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted',
		code = [[
	const global <?=eqn.cons_t?>* U = buf + index;
	value->vreal = real3_weightedLen(value->vreal3, calc_gamma_ll(U, x));
]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted',
		code = [[
	const global <?=eqn.cons_t?>* U = buf + index;
	value->vreal = sym3_dot(value->vsym3, calc_gamma_uu(U, x));]],
	})
end

return EinsteinFiniteDifferenceSolver
