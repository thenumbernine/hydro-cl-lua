local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local math = require 'ext.math'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'hydro.tooltip'
local template = require 'template'
local CLBuffer = require 'cl.obj.buffer'

local half = require 'hydro.half'
local toreal, fromreal = half.toreal, half.fromreal


local PoissonKrylov = class()

PoissonKrylov.name = 'PoissonKrylov'
PoissonKrylov.scalar = 'real'
PoissonKrylov.potentialField = 'ePot'

--[[
TODO how to make this universal across all linear solver (including Jacobi).
Right now I'm comparing it to the norm of x.  I should be using the residual.
--]]
PoissonKrylov.stopEpsilon = 1e-10

PoissonKrylov.verbose = false
PoissonKrylov.linearSolverType = 'conjres'

function PoissonKrylov:getPotBufType()
	return self.solver.UBufObj.type
end

function PoissonKrylov:getPotBuf()
	return self.solver.UBuf
end

function PoissonKrylov:init(args)
	local solver = assert(args.solver)
	self.solver = solver
	
	self.potentialField = args.potentialField
	self.verbose = args.verbose
	self.linearSolverType = args.linearSolver

	-- matches hydro/op/relaxation.lua
	self.name = solver.app:uniqueName(self.name)
end

function PoissonKrylov:initSolver()
	local solver = self.solver


	local CLKrylov = require('solver.cl.'..self.linearSolverType)

	-- can I combine this with the CLKrylov in int/be.lua somehow?
	local ThisKrylov = class(CLKrylov)

	function ThisKrylov:newBuffer(name)
		if not self.cache then self.cache = {} end
		local cached = self.cache[name]
		if cached then return cached end
		cached = ThisKrylov.super.newBuffer(self, name)
		cached:fill()
		self.cache[name] = cached
		return cached
	end

	
	
	for _,name in ipairs{
		'krylov_b',
		'krylov_x',
	} do
		self[name..'Obj'] = CLBuffer{
			env = solver.app.env,
			name = name,
			type = 'real',
			count = solver.numCells,	-- without border?
		}
	end

	local codePrefix = table{
		solver.modules:getHeader(solver.solverModuleNames:unpack()),
		solver.modules:getCode(solver.solverModuleNames:unpack()),
	}:concat'\n'

	local mulWithoutBorderKernelObj = solver.domain:kernel{
		name = 'Poisson_mulWithoutBorder'..self.name,
		header = codePrefix,
		argsOut = {
			{name='y', type=solver.app.real, obj=true},
		},
		argsIn = {
			solver.solverBuf,
			{name='a', type=solver.app.real, obj=true},
			{name='b', type=solver.app.real, obj=true},
		},
		body = [[	
	if (OOB(numGhost, numGhost)) {
		y[index] = 0;
		return;
	}
	y[index] = a[index] * b[index];
]],
	}

	local squareKernelObj = solver.domain:kernel{
		name = 'Poisson_square'..self.name,
		header = codePrefix,
		argsOut = {
			{name = 'y', type=solver.app.real, obj=true},
		},
		argsIn = {
			solver.solverBuf,
			{name = 'x', type=solver.app.real, obj=true},
		},
		body = [[
	if (OOB(0,0)) return;
	y[index] = x[index] * x[index];
]],
	}

	local numreals = solver.numCells
	local volumeWithoutBorder = solver.volumeWithoutBorder
	local numRealsWithoutBorder = volumeWithoutBorder

	local sum = solver.app.env:reduce{
		count = numreals,
		op = function(x,y) return x..' + '..y end,
	}
	local dotWithoutBorder = function(a,b)
		mulWithoutBorderKernelObj(sum.buffer, self.solverBuf, a, b)
		return sum()
	end
	
	local restart = args and args.restart or 10
	self.lastResidual = 0
	self.lastIter = 0
	local linearSolverArgs = {
		env = solver.app.env,
		x = self.krylov_xObj,
		count = numreals,
		epsilon = self.stopEpsilon,
		restart = restart,
		maxiter = cmdline.selfGravPoissonMaxIter or 20,	--restart * numreals,
		-- logging:
		errorCallback = function(residual, iter, x, rLenSq)
			local lastResidual, lastIter = self.lastResidual, self.lastIter
			self.lastResidual, self.lastIter = residual, iter
			if self.verbose then
				-- square from x to reduceBuf
				squareKernelObj(x, solver.solverBuf, x)
				local xNorm = math.sqrt(solver.reduceSum() / numRealsWithoutBorder)
				
				solver.cmds:enqueueCopyBuffer{
					src=assert(x.obj),
					dst=assert(solver.reduceBuf),
					size = ffi.sizeof(solver.app.real) * numreals,
				}
				local xmin = fromreal(solver.reduceMin())
				solver.cmds:enqueueCopyBuffer{
					src=assert(x.obj),
					dst=assert(solver.reduceBuf),
					size = ffi.sizeof(solver.app.real) * numreals,
				}
				local xmax = fromreal(solver.reduceMax())
				io.stderr:write(table{iter, residual, xNorm, xmin, xmax}:map(tostring):concat'\t','\n')
			end
			if math.abs(residual) < self.stopEpsilon then
				return true
			end
		end,
		dot = function(a,b)
			return dotWithoutBorder(a,b) / numRealsWithoutBorder
		end,
	}

	linearSolverArgs.b = self.krylov_bObj
	
	linearSolverArgs.A = function(UNext, U)
		-- A(x) = div x
		-- but don't use the routine in hydro/op/poisson.cl, that's for Jacobi
		self.poissonKrylovLinearFuncKernelObj(solver.solverBuf, UNext, U)
	end	
	self.linearSolver = ThisKrylov(linearSolverArgs)
end

local poissonKrylovCode = [[
<?
local solver = op.solver
local eqn = op.solver.eqn
?>

kernel void poissonKrylovLinearFunc<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global real* Y,
	global real* X
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost, numGhost)) {
		Y[index] = 0.;
		return;
	}
	real3 x = cell_x(i);

<? for j=0,solver.dim-1 do ?>
	real dx<?=j?> = cell_dx<?=j?>(x);
<? end ?>

	real3 xInt = x;
	real3 volL, volR;
<? for j=0,solver.dim-1 do 
?>	xInt.s<?=j?> = x.s<?=j?> - .5 * solver->grid_dx.s<?=j?>;
	volL.s<?=j?> = cell_sqrt_det_g(solver, xInt);
	xInt.s<?=j?> = x.s<?=j?> + .5 * solver->grid_dx.s<?=j?>;
	volR.s<?=j?> = cell_sqrt_det_g(solver, xInt);
	xInt.s<?=j?> = x.s<?=j?>;
<? end 
?>	real volAtX = cell_sqrt_det_g(solver, x);

	real sum = (0.
<? for j=0,solver.dim-1 do ?>
		+ volR.s<?=j?> * X[index + solver->stepsize.s<?=j?>] / (dx<?=j?> * dx<?=j?>)
		+ volL.s<?=j?> * X[index - solver->stepsize.s<?=j?>] / (dx<?=j?> * dx<?=j?>)
<? end 
?>	) / volAtX;

	sum += X[index] * (0.
<? for j=0,solver.dim-1 do ?>
		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end ?>
	) / volAtX;

	Y[index] = sum;
}

kernel void copyPotentialFieldToVecAndInitB<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global real* x,
	global real* b,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0, 0);
	
	global const <?=eqn.cons_t?>* U = UBuf + index;
	x[index] = U-><?=op.potentialField?>;
	
	real source = 0.;
<?=op:getPoissonDivCode() or ''?>
	b[index] = source;
}

kernel void copyVecToPotentialField<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	global const real* x
) {
	SETBOUNDS(0, 0);
	UBuf[index].<?=op.potentialField?> = x[index];
}
]]

function PoissonKrylov:initCodeModules(solver)
	solver.modules:add{
		name ='op.PoissonKrylov',
		code = template(
			file['hydro/op/poisson.cl']..'\n'
			..poissonKrylovCode,
			{
				op = self,
			}
		)..'\n'
		..(self:getPoissonCode() or ''),
	}
	solver.solverModuleNames:insert'op.PoissonKrylov'
end

function PoissonKrylov:refreshSolverProgram()
	local solver = self.solver
	self:initSolver()
	self.initPotentialKernelObj = solver.solverProgramObj:kernel('initPotential'..self.name, solver.solverBuf, self:getPotBuf())
	self.copyPotentialFieldToVecAndInitBKernelObj = solver.solverProgramObj:kernel('copyPotentialFieldToVecAndInitB'..self.name, solver.solverBuf, assert(self.krylov_xObj.obj), self.krylov_bObj.obj, self:getPotBuf())
	self.copyVecToPotentialFieldKernelObj = solver.solverProgramObj:kernel('copyVecToPotentialField'..self.name, solver.solverBuf, self:getPotBuf(), self.krylov_xObj.obj)
	self.poissonKrylovLinearFuncKernelObj = solver.solverProgramObj:kernel('poissonKrylovLinearFunc'..self.name)
end

-- matches hydro/op/relaxation.lua
function PoissonKrylov:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to PoissonKrylov:potentialField
	self.potentialBoundaryProgramObj, self.potentialBoundaryKernelObjs =
		solver:createBoundaryProgramAndKernel{
			type = self:getPotBufType(),
			methods = solver.boundaryMethods,
			fields = {self.potentialField},
			programNameSuffix = '-'..self.name,	-- self.name is used to uniquely name the kernel as well as this boundary program
		}
	for _,obj in ipairs(self.potentialBoundaryKernelObjs) do
		obj.obj:setArg(1, self:getPotBuf())
	end
end

-- matches hydro/op/relaxation.lua
function PoissonKrylov:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initPotentialKernelObj()
	self:potentialBoundary()
	self:relax()
end

function PoissonKrylov:relax()
	local solver = self.solver
	-- apply boundary conditions to UBuf.potentialField before our iterative solution	
	-- (TODO should we be applying it each iteration as well?)
	self:potentialBoundary()
	-- copy potential field into krylov_x
	-- calculate krylov_b by PoissonKrylov:getCodeParams():getPoissonDivCode()
	self.copyPotentialFieldToVecAndInitBKernelObj()

	-- solve
	self.linearSolver()
	
	-- copy krylov_x back to potential field
	self.copyVecToPotentialFieldKernelObj()
end

-- matches to hydro/op/relaxation.lua
function PoissonKrylov:potentialBoundary()
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernelObjs)
end

-- similar to hydro/op/relaxation.lua
function PoissonKrylov:updateGUI()
	ig.igPushIDStr(self.name..' solver')
	-- TODO name from 'field' / 'enableField', though those aren't properties of PoissonKrylov
	if ig.igCollapsingHeader'PoissonKrylov solver' then
		tooltip.numberTable('Krylov epsilon', self.linearSolver.args, 'epsilon')
		tooltip.intTable('GMRES restart', self.linearSolver.args, 'restart')
		tooltip.intTable('maxiter', self.linearSolver.args, 'maxiter')	-- typically restart * number of reals = restart * numCells * number of states
		ig.igText('residual = '..self.lastResidual)
		ig.igText('iter = '..self.lastIter)
	end
	ig.igPopID()
end

return PoissonKrylov
