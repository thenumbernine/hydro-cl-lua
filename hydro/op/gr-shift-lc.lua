--[[
Why not just look back along the shift vector and advect?
I'm sure there's something wrong with this
--]]

local ffi = require 'ffi'
local class = require 'ext.class'
local template = require 'template'
local real = require 'hydro.real'

local LagrangianCoordinateShift = class()

function LagrangianCoordinateShift:init(args)
	self.solver = assert(args.solver)
end

function LagrangianCoordinateShift:initCodeModules()
	local solver = self.solver
	solver.modules:add{
		name = 'op.LagrangianCoordinateShift',
		code = solver.eqn:template([[
kernel void lagrangianCoordinateAdvect(
	<?=cons_t?> const * const dstUBuf,
	global <?=cons_t?> const * const UBuf,
	real dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	// U[i] = U[i - beta * dt]
	global <?=cons_t?> const * const U = UBuf + index;

	real4 src = (real4)(
<? 
for j=0,solver.dim-1 do
?>		(real)i.s<?=j?> - dt * U->beta_u.s<?=j?> / grid_dx<?=j?>,
<? 
end
for j=solver.dim,2 do
?>		0.,
<? 
end
?>		0.);

	real4 srcf = floor(src);
	real4 f = src - srcf;
	real4 nf = 1. - f;

	//TODO boundary conditions... for now just freeflow
	isrc = clamp((int4)srcf, (int4)(0,0,0,0), gridSize - 1);

	//here's the next cell
	int4 isrc2 = clamp(isrc + 1, gridSize - 1);

	global <?=cons_t?> *dstU = dstUBuf + index;
<? if solver.dim == 1 then
?>	global <?=cons_t?> const * const srcUL = UBuf[INDEX(solver, isrc.x, isrc.y, isrc.z)];
	global <?=cons_t?> const * const srcUR = UBuf[INDEX(solver, isrc2.x, isrc.y, isrc.z)];
	for (int j = 0; j < numStates; ++j) {
		dstU.ptr[j] = (1. - f.x) * srcUL->ptr[j] + f.x * srcUR->ptr[j];
	}
<? elseif solver.dim == 2 then
?>	global <?=cons_t?> const * const srcULL = UBuf[INDEX(solver, isrc.x, isrc.y, isrc.z)];
	global <?=cons_t?> const * const srcULR = UBuf[INDEX(solver, isrc.x, isrc2.y, isrc.z)];
	global <?=cons_t?> const * const srcURL = UBuf[INDEX(solver, isrc2.x, isrc.y, isrc.z)];
	global <?=cons_t?> const * const srcURR = UBuf[INDEX(solver, isrc2.x, isrc2.y, isrc.z)];
	for (int j = 0; j < numStates; ++j) {
		real UL = nf.x * srcULL->ptr[j] + f.x * srcURL->ptr[j];
		real UR = nf.x * srcULR->ptr[j] + f.x * srcURR->ptr[j];
		dstU.ptr[j] = (1. - f.y) * UL + f.y * UR;
	}
<? elseif solver.dim == 3 then
?>	global <?=cons_t?> const * const srcULLL = UBuf[INDEX(solver, isrc.x, isrc.y, isrc.z)];
	global <?=cons_t?> const * const srcULRL = UBuf[INDEX(solver, isrc.x, isrc2.y, isrc.z)];
	global <?=cons_t?> const * const srcURLL = UBuf[INDEX(solver, isrc2.x, isrc.y, isrc.z)];
	global <?=cons_t?> const * const srcURRL = UBuf[INDEX(solver, isrc2.x, isrc2.y, isrc.z)];
	global <?=cons_t?> const * const srcULLR = UBuf[INDEX(solver, isrc.x, isrc.y, isrc2.z)];
	global <?=cons_t?> const * const srcULRR = UBuf[INDEX(solver, isrc.x, isrc2.y, isrc2.z)];
	global <?=cons_t?> const * const srcURLR = UBuf[INDEX(solver, isrc2.x, isrc.y, isrc2.z)];
	global <?=cons_t?> const * const srcURRR = UBuf[INDEX(solver, isrc2.x, isrc2.y, isrc2.z)];
	for (int j = 0; j < numStates; ++j) {
		real ULL = nf.x * srcULLL->ptr[j] + f.x * srcURLL->ptr[j];
		real ULR = nf.x * srcULLR->ptr[j] + f.x * srcURLR->ptr[j];
		real URL = nf.x * srcULRL->ptr[j] + f.x * srcURRL->ptr[j];
		real URR = nf.x * srcULRR->ptr[j] + f.x * srcURRR->ptr[j];
		real UL = nf.y * ULL + f.y * URL;
		real UR = nf.y * ULR + f.y * URR;
		dstU.ptr[j] = nf.z * UL + f.z * UR;
	}
<? end
?>
}
]], 	{
			op = self,
		}),
	}
	solver.solverModulesEnabled['op.LagrangianCoordinateShift'] = true
end

function LagrangianCoordinateShift:refreshSolverProgram()
	local solver = self.solver
	self.lagrangianCoordinateAdvectKernelObj = solver.solverProgramObj:kernel(
		'lagrangianCoordinateAdvect',
		solver.integrator.derivBuf,	-- used as a temp
		solver.UBuf)
end

function LagrangianCoordinateShift:step(dt)
	local solver = self.solver
	self.lagrangianCoordinateAdvectKernelObj.obj:setArg(2, real(dt))
	self.lagrangianCoordinateAdvectKernelObj()
	self.cmds:enqueueCopyBuffer{
		src = solver.integrator.derivBuf,
		dst = solver.UBuf,
		size = solver.UBufObj.size * ffi.sizeof(solver.cons_t),
	}
end

return LagrangianCoordinateShift 
