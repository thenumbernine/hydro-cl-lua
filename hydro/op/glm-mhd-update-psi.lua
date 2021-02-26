local class = require 'ext.class'
local template = require 'template'
local ffi = require 'ffi'
local real = require 'hydro.real'

local GLM_MHD_UpdatePsi = class()

function GLM_MHD_UpdatePsi:init(args)
	self.solver = assert(args.solver)
end

function GLM_MHD_UpdatePsi:initCodeModules()
	local solver = self.solver
	solver.modules:add{
		name = 'op.GLM_MHD_UpdatePsi',
		depends = {
			self.solver.eqn.symbols.eqn_guiVars_compileTime,	-- Cp
		},
		code = solver.eqn:template([[
kernel void updatePsi(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	real const dt,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(0,0);

	global <?=cons_t?> * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	
	//TODO don't need the whole eigen here, just the Ch	
<? if not eqn.useFixedCh then ?>
	real Ch = 0;
	<? for side=0,solver.dim-1 do ?>{
		<?=eigen_t?> eig;
		<?=eigen_forCell?>(&eig, solver, U, cell, normal_fromSide<?=side?>(cell->pos));
		Ch = max(Ch, eig.Ch);
	}<? end ?>
<? else ?>
	real Ch = solver->Ch;
<? end ?>
	
	U->psi *= exp(-dt * Ch * Ch / (Cp * Cp) * U->psi);
}
]], 	{
			op = self,
		}),
	}
	solver.solverModulesEnabled['op.GLM_MHD_UpdatePsi'] = true
end

function GLM_MHD_UpdatePsi:refreshSolverProgram()
	local solver = self.solver
	self.updatePsiKernelObj = solver.solverProgramObj:kernel('updatePsi', solver.solverBuf, solver.UBuf, solver.cellBuf)
end

function GLM_MHD_UpdatePsi:step(dt)
	local dtArg = real(dt)
	self.updatePsiKernelObj.obj:setArg(2, dtArg)
end

return GLM_MHD_UpdatePsi
