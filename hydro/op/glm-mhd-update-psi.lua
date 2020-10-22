local class = require 'ext.class'
local template = require 'template'
local ffi = require 'ffi'
local real = require 'hydro.real'

local GLM_MHD_UpdatePsi = class()

function GLM_MHD_UpdatePsi:init(args)
	self.solver = assert(args.solver)
end

function GLM_MHD_UpdatePsi:initCodeModules(solver)
	solver.modules:add{
		name = 'op.GLM_MHD_UpdatePsi',
		depends = {
			'cell_x',
			'eqn.guiVars.compileTime',	-- Cp
		},
		code = template([[
kernel void updatePsi(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	real dt,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);

	global <?=eqn.cons_t?>* U = UBuf + index;
	
	//TODO don't need the whole eigen here, just the Ch	
<? if not eqn.useFixedCh then ?>
	real Ch = 0;
	<? for side=0,solver.dim-1 do ?>{
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(solver, *U, x);
		Ch = max(Ch, eig.Ch);
	}<? end ?>
<? else ?>
	real Ch = solver->Ch;
<? end ?>
	
	U->psi *= exp(-dt * Ch * Ch / (Cp * Cp) * U->psi);
}
]], 	{
			op = self,
			solver = self.solver,
			eqn = self.solver.eqn,
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
