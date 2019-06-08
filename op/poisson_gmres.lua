local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local math = require 'ext.math'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local template = require 'template'

local CLBuffer = require 'cl.obj.buffer'

--local CLKrylov = require 'solver.cl.conjgrad'
local CLKrylov = require 'solver.cl.conjres'
--local CLKrylov = require 'solver.cl.gmres'

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

local Poisson = class()

Poisson.name = 'Poisson'

Poisson.potentialField = 'ePot'

Poisson.scalar = 'real'

function Poisson:getPotBufType()
	return self.solver.UBufObj.type
end

function Poisson:getPotBuf()
	return self.solver.UBuf
end

function Poisson:init(args)
	local solver = assert(args.solver)
	self.solver = solver
	self.potentialField = args.potentialField

	-- matches op/relaxation.lua
	self.name = solver.app:uniqueName(self.name)
end

function Poisson:initSolver()
	local solver = self.solver
	
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

	local mulWithoutBorder = solver.domain:kernel{
		name = 'Poisson_mulWithoutBorder',
		header = solver.codePrefix,
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
	local volumeWithoutBorder = tonumber(solver.sizeWithoutBorder:volume())

	local sum = solver.app.env:reduce{
		count = volumeWithoutBorder,
		op = function(x,y) return x..' + '..y end,
	}
	local dotWithoutBorder = function(a,b)
		mulWithoutBorder(sum.buffer, self.solverBuf, a, b)
		return sum()
	end
	
	local restart = args and args.restart or 10
	local epsilon = 1e-10
	self.lastResidual = 0
	self.lastIter = 0
	local linearSolverArgs = {
		env = solver.app.env,
		x = self.krylov_xObj,
		count = volumeWithoutBorder,
		epsilon = epsilon,
		--maxiter = 1000,
		restart = restart,
		maxiter = restart * volumeWithoutBorder,
		-- logging:
		errorCallback = function(residual, iter, x, rLenSq)
			local lastResidual, lastIter = self.lastResidual, self.lastIter
			self.lastResidual, self.lastIter = residual, iter
print('krylov iter', iter, 'residual', residual)			
			if not math.isfinite(residual) then
				print("got non-finite residual: "..residual)	-- error?
				return true	-- fail
			end
			if math.abs(residual - lastResidual) < epsilon then
				return true
			end
		end,
		dot = function(a,b)
			return dotWithoutBorder(a,b) / volumeWithoutBorder 
		end,
	}

	linearSolverArgs.b = self.krylov_bObj
	
	linearSolverArgs.A = function(UNext, U)
		-- A(x) = div x
		-- but don't use the poisson.cl one, that's for in-place Gauss-Seidel
		self.poissonKrylovLinearFuncKernelObj(solver.solverBuf, UNext, U)
	end	
	self.linearSolver = ThisKrylov(linearSolverArgs)
end

local poissonKrylovCode = [[

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
	volL.s<?=j?> = cell_volume(solver, xInt);
	xInt.s<?=j?> = x.s<?=j?> + .5 * solver->grid_dx.s<?=j?>;
	volR.s<?=j?> = cell_volume(solver, xInt);
	xInt.s<?=j?> = x.s<?=j?>;
<? end 
?>	real volAtX = cell_volume(solver, x);

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
	b[index] = -source;
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

-- TODO rename to 'getCode'
function Poisson:getSolverCode()
	return table{
		template(
			table{
				file['op/poisson.cl'],
				poissonKrylovCode,
				self:getPoissonCode() or '',
			}:concat'\n', {
				op = self,
				solver = self.solver,
				eqn = self.solver.eqn,
			}
		),
	}:concat'\n'
end

function Poisson:refreshSolverProgram()
	local solver = self.solver
	self:initSolver()
	self.initPotentialKernelObj = solver.solverProgramObj:kernel('initPotential'..self.name, solver.solverBuf, self:getPotBuf())
	self.copyPotentialFieldToVecAndInitBKernelObj = solver.solverProgramObj:kernel('copyPotentialFieldToVecAndInitB'..self.name, solver.solverBuf, assert(self.krylov_xObj.obj), self.krylov_bObj.obj, self:getPotBuf())
	self.copyVecToPotentialFieldKernelObj = solver.solverProgramObj:kernel('copyVecToPotentialField'..self.name, solver.solverBuf, self:getPotBuf(), self.krylov_xObj.obj)
	self.poissonKrylovLinearFuncKernelObj = solver.solverProgramObj:kernel('poissonKrylovLinearFunc'..self.name)
end

-- matches op/relaxation.lua
function Poisson:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to Poisson:potentialField
	self.potentialBoundaryProgramObj, self.potentialBoundaryKernelObjs =
		solver:createBoundaryProgramAndKernel{
			type = self:getPotBufType(),
			methods = table.map(solver.boundaryMethods, function(v)
				return (select(2, next(solver.boundaryOptions[v])))
			end),
			assign = function(a,b)
				return a..'.'..self.potentialField..' = '..b..'.'..self.potentialField
			end,
		}
	for _,obj in ipairs(self.potentialBoundaryKernelObjs) do
		obj.obj:setArg(1, self:getPotBuf())
	end
end

-- matches op/relaxation.lua
function Poisson:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initPotentialKernelObj()
	self:potentialBoundary()
	self:relax()
end

function Poisson:relax()
	local solver = self.solver
	
	-- apply boundary conditions to UBuf.potentialField before our iterative solution	
	-- (TODO should we be applying it each iteration as well?)
	self:potentialBoundary()
	-- copy potential field into krylov_x
	-- calculate krylov_b by Poisson:getCodeParams():getPoissonDivCode()
	self.copyPotentialFieldToVecAndInitBKernelObj()
	-- solve
	self.linearSolver()
	-- copy krylov_x back to potential field
	self.copyVecToPotentialFieldKernelObj()
end

-- matches to op/relaxation.lua
function Poisson:potentialBoundary()
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernelObjs)
end

-- similar to op/relaxation.lua
function Poisson:updateGUI()
	ig.igPushIDStr(self.name..' solver')
	-- TODO name from 'field' / 'enableField', though those aren't properties of Poisson
	if ig.igCollapsingHeader'Poisson solver' then
		tooltip.numberTable('Krylov epsilon', self.linearSolver.args, 'epsilon')
		tooltip.intTable('GMRES restart', self.linearSolver.args, 'restart')
		tooltip.intTable('maxiter', self.linearSolver.args, 'maxiter')	-- typically restart * number of reals = restart * numCells * number of states
		ig.igText('residual = '..self.lastResidual)
		ig.igText('iter = '..self.lastIter)
	end
	ig.igPopID()
end

return Poisson
