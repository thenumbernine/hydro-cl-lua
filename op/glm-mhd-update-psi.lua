local class = require 'ext.class'
local template = require 'template'
local ffi = require 'ffi'

local GLM_MHD_UpdatePsi = class()

function GLM_MHD_UpdatePsi:init(args)
	self.solver = assert(args.solver)
end

function GLM_MHD_UpdatePsi:getSolverCode()
	local solver = self.solver
	return template([[
kernel void updatePsi(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	real dt
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	//TODO don't need the whole eigen here, just the Ch	
<? if not eqn.useFixedCh then ?>
	real Ch = 0;
	<? for side=0,solver.dim-1 do ?>{
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(*U, x);
		Ch = max(Ch, eig.Ch);
	}<? end ?>
<? end ?>
	
	U->psi *= exp(-dt * Ch * Ch / (Cp * Cp) * U->psi);
}
]], {
		solver = solver,
		eqn = solver.eqn,
	})
end

function GLM_MHD_UpdatePsi:refreshSolverProgram()
	local solver = self.solver
	self.updatePsiKernelObj = solver.solverProgramObj:kernel('updatePsi', solver.solverBuf, solver.UBuf)
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

function GLM_MHD_UpdatePsi:step(dt)
	local dtArg = real(dt)
	self.updatePsiKernelObj.obj:setArg(2, dtArg)
end

return GLM_MHD_UpdatePsi
