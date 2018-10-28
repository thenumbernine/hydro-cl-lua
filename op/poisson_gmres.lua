local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local math = require 'ext.math'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local template = require 'template'

local CLBuffer = require 'cl.obj.buffer'
local CLGMRES = require 'solver.cl.gmres'

-- can I combine this with the CLGMRES in int/be.lua somehow?
local ThisGMRES = class(CLGMRES)

function ThisGMRES:newBuffer(name)
	if not self.cache then self.cache = {} end
	local cached = self.cache[name]
	if cached then return cached end
	cached = ThisGMRES.super.newBuffer(self, name)
	cached:fill()
	self.cache[name] = cached
	return cached
end

local Poisson = class()

Poisson.potentialField = 'ePot'

function Poisson:init(args)
	self.solver = assert(args.solver)
	self.potentialField = args.potentialField
end

function Poisson:getPotBufType()
	return self.solver.eqn.cons_t
end

function Poisson:getPotBuf()
	return self.solver.UBuf
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
			size = solver.numCells,	-- without border?
		}
	end

	local mulWithoutBorder = solver.domain:kernel{
		name = 'Poisson_mulWithoutBorder',
		header = solver.codePrefix,
		argsOut = {
			{name='y', type=solver.app.real, obj=true},
		},
		argsIn = {
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
		size = volumeWithoutBorder,
		op = function(x,y) return x..' + '..y end,
	}
	local dotWithoutBorder = function(a,b)
		mulWithoutBorder(sum.buffer, a, b)
		return sum()
	end

	self.last_err = 0
	self.last_iter = 0
	local linearSolverArgs = {
		env = solver.app.env,
		x = self.krylov_xObj,
		size = volumeWithoutBorder,
		epsilon = 1e-10,
		--maxiter = 1000,
		restart = 10,
		maxiter = 10 * volumeWithoutBorder,
		-- logging:
		errorCallback = function(err, iter, x, rLenSq)
			self.last_err = err
			self.last_iter = iter
			
			if not math.isfinite(err) then
				print("got non-finite err: "..err)	-- error?
				return true	-- fail
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
		self.poissonGMRESLinearFuncKernelObj(solver.solverBuf, UNext, U)
	end	
	self.linearSolver = ThisGMRES(linearSolverArgs)
end

local poissonGMRESCode = [[

kernel void poissonGMRESLinearFunc(
	constant <?=solver.solver_t?>* solver,
	global real* y,
	global real* x
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost, numGhost)) {
		y[index] = 0.;
		return;
	}

<? for j=0,solver.dim-1 do ?>
	real dx<?=j?> = dx<?=j?>_at(i);
<? end ?>

	real3 intIndex = _real3(i.x, i.y, i.z);
	real3 volL, volR;
<? for j=0,solver.dim-1 do ?>
	intIndex.s<?=j?> = i.s<?=j?> - .5;
	volL.s<?=j?> = volume_at(solver, cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?> + .5;
	volR.s<?=j?> = volume_at(solver, cell_x(intIndex));
	intIndex.s<?=j?> = i.s<?=j?>;
<? end ?>
	real volAtX = volume_at(solver, cell_x(i));

	real sum = (0.
<? for j=0,solver.dim-1 do ?>
		+ volR.s<?=j?> * x[index + solver->stepsize.s<?=j?>] / (dx<?=j?> * dx<?=j?>)
		+ volL.s<?=j?> * x[index - solver->stepsize.s<?=j?>] / (dx<?=j?> * dx<?=j?>)
<? end 
?>	) / volAtX;

	sum += x[index] * (0.
<? for j=0,solver.dim-1 do ?>
		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end ?>
	) / volAtX;

	y[index] = sum;
}

kernel void copyPotentialFieldToVecAndInitB(
	global real* x,
	global real* b,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0, 0);

	global const <?=eqn.cons_t?>* U = UBuf + index;

	x[index] = U-><?=poisson.potentialField?>;
	
	real rho = 0.;
	<?=calcRho?>
	b[index] = rho;
}

kernel void copyVecToPotentialField(
	global <?=eqn.cons_t?>* UBuf,
	global const real* x
) {
	SETBOUNDS(0, 0);
	
	UBuf[index].<?=poisson.potentialField?> = x[index];
}
]]

-- TODO rename to 'getCode'
function Poisson:getSolverCode()
	return table{
		template(
			table{
				file['op/poisson.cl'],
				poissonGMRESCode,
				self:getPoissonCode() or '',
			}:concat'\n',
			table(self:getCodeParams(), {
				op = self,
				solver = self.solver,
				eqn = self.solver.eqn,
			})),
	}:concat'\n'
end

function Poisson:refreshSolverProgram()
	local solver = self.solver
	self:initSolver()
	self.initRelaxationPotentialKernelObj = solver.solverProgramObj:kernel('initRelaxationPotential', self:getPotBuf())
	self.copyPotentialFieldToVecAndInitBKernelObj = solver.solverProgramObj:kernel('copyPotentialFieldToVecAndInitB', assert(self.krylov_xObj.obj), self.krylov_bObj.obj, self:getPotBuf())
	self.copyVecToPotentialFieldKernelObj = solver.solverProgramObj:kernel('copyVecToPotentialField', self:getPotBuf(), self.krylov_xObj.obj)
	self.poissonGMRESLinearFuncKernelObj = solver.solverProgramObj:kernel'poissonGMRESLinearFunc'
end

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

-- TODO
-- for Euler, add potential energy into total energy
-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
function Poisson:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initRelaxationPotentialKernelObj()
	solver:potentialBoundary()
	self:relax()
end

function Poisson:relax()
	local solver = self.solver
	-- copy potential field into krylov_x
	-- calculate krylov_b by Poisson:getCodeParams().calcRho
	self.copyPotentialFieldToVecAndInitBKernelObj()
	-- solve
	self.linearSolver()
	-- copy krylov_x back to potential field
	self.copyVecToPotentialFieldKernelObj()
end

function Poisson:potentialBoundary()
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernelObjs)
end

function Poisson:updateGUI()
	-- TODO unique name for other Poisson solvers?
	ig.igPushIDStr'Poisson GMRES solver'
	-- TODO name from 'field' / 'enableField', though those aren't properties of Poisson
	if ig.igCollapsingHeader'Poisson solver' then
		tooltip.numberTable('Krylov epsilon', self.linearSolver.args, 'epsilon')
		tooltip.intTable('GMRES restart', self.linearSolver.args, 'restart')
		tooltip.intTable('Krylov maxiter', self.linearSolver.args, 'maxiter')	-- typically restart * number of reals = restart * numCells * number of states
		-- read-only:
		ig.igText('err = '..self.last_err)
		ig.igText('iter = '..self.last_iter)
	end
	ig.igPopID()
end

return Poisson
