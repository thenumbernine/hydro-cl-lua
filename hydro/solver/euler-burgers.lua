local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local real = require 'hydro.real'
local FiniteVolumeSolver = require 'hydro.solver.fvsolver'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


-- TODO make this work with ops, specifically Euler's SelfGrav
local EulerBurgers = class(FiniteVolumeSolver)
EulerBurgers.name = 'EulerBurgers'

function EulerBurgers:init(...)
	EulerBurgers.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

-- override this -- don't expect a flux argument
function EulerBurgers:createFlux(fluxName, fluxArgs)
	self.flux = {
		initCodeModules = function()
			self.modules:add{name = 'calcFlux'}
		end,
	}
end

function EulerBurgers:createEqn()
	-- TODO put this in its own eqn file? eqn/euler-burgers.lua?
	local EulerEqn = require 'hydro.eqn.euler'
	
	local EulerBurgersEqn = class(EulerEqn)
	
	function EulerBurgersEqn:initCodeModule_calcDT()
		local solver = self.solver
		solver.modules:add{
			name = 'calcDT',
			depends = {
				'solver_t',
				'eqn.prim-cons',
				'eqn.guiVars.compileTime',
			},
			code = self:template[[
<? if require 'hydro.solver.gridsolver'.is(solver) then ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf			//[numCells]
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cellBuf[index].pos;

	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real Cs = calc_Cs(solver, &W);

	real dt = INFINITY;
<? for side=0,solver.dim-1 do 
?>	dt = min(dt, (real)solver->grid_dx.s<?=side?> / (Cs + fabs(W.v.s<?=side?>)));
<? end
?>	dtBuf[index] = dt;
}

<? else -- mesh solver ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,					//[numCells]
	const global <?=eqn.cons_t?>* UBuf,	//[numCells]
	const global <?=solver.coord.cell_t?>* cellBuf,			//[numCells]
	const global <?=solver.coord.face_t?>* faceBuf,			//[numFaces]
	const global int* cellFaceIndexes	//[numCellFaceIndexes]
) {
	int cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	const global cell_t* cell = cellBuf + cellIndex;
	real3 x = cell->pos;
	
	const global <?=eqn.cons_t?>* U = UBuf + cellIndex;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real Cs = calc_Cs(solver, &W);

	real dt = INFINITY;
	for (int i = 0; i < cell->faceCount; ++i) {
		const global face_t* face = faceBuf + cellFaceIndexes[i + cell->faceOffset];
		normal_t n = normal_forFace(face);
		real v_n = normal_vecDotN1(n, W.v);
		real dx = face->area;
		dt = (real)min(dt, dx / (Cs + fabs(v_n)));
	}
	dtBuf[cellIndex] = dt;
}

<? end -- mesh solver ?>
]],
		}
	end
	
	self.eqn = EulerBurgersEqn(table(self.eqnArgs, {solver=self}))
end

function EulerBurgers:createBuffers()
	EulerBurgers.super.createBuffers(self)
	
	self:clalloc('intVelBuf', self.app.real, self.numCells * self.dim)
	self:clalloc('PBuf', self.app.real, self.numCells)
end

function EulerBurgers:initCodeModules()
	EulerBurgers.super.initCodeModules(self)

	self.modules:add{
		name = 'EulerBurgers.solver',
		depends = {
			'eqn.prim-cons',
			'fluxLimiter',
			'eigen_forInterface',
		},
		code = template(file['hydro/solver/euler-burgers.cl'], {solver=self, eqn=self.eqn}),
	}
	self.solverModulesEnabled['EulerBurgers.solver'] = true
end

function EulerBurgers:refreshSolverProgram()
	EulerBurgers.super.refreshSolverProgram(self)

	-- no mention of ULR just yet ...

	self.calcIntVelKernelObj = self.solverProgramObj:kernel'calcIntVel'
	self.calcFluxKernelObj = self.solverProgramObj:kernel'calcFlux'

	self.computePressureKernelObj = self.solverProgramObj:kernel('computePressure', self.solverBuf, self.PBuf, self.UBuf, self.cellBuf)
	
	self.diffuseMomentumKernelObj = self.solverProgramObj:kernel{name='diffuseMomentum', domain=self.domainWithoutBorder}
	self.diffuseMomentumKernelObj.obj:setArg(0, self.solverBuf)
	self.diffuseMomentumKernelObj.obj:setArg(2, self.PBuf)
	
	self.diffuseWorkKernelObj = self.solverProgramObj:kernel'diffuseWork'
	self.diffuseWorkKernelObj.obj:setArg(0, self.solverBuf)
	self.diffuseWorkKernelObj.obj:setArg(2, self.UBuf)
	self.diffuseWorkKernelObj.obj:setArg(3, self.PBuf)
end

function EulerBurgers:addDisplayVars()
	EulerBurgers.super.addDisplayVars(self)

	self:addDisplayVarGroup{
		name = 'P', 
		bufferType = 'real',
		bufferField = 'PBuf',
		vars = {
			{name='P', code='value.vreal = buf[index];'},
		},
	}

	for j,xj in ipairs(xNames) do
		self:addDisplayVarGroup{
			name = 'intVel', 
			bufferType = 'real',
			bufferField = 'intVelBuf',
			codePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
]],
			vars = range(0,self.dim-1):map(function(i)
				return {name=xj..'_'..i, code='value.vreal = buf['..i..' + indexInt];'}
			end),
		}
	end
end

function EulerBurgers:step(dt)
	-- calc deriv here
	self.integrator:integrate(dt, function(derivBufObj)
		self.calcIntVelKernelObj(self.solverBuf, self.intVelBuf, self:getULRBuf(), self.cellBuf)

		self.calcFluxKernelObj(self.solverBuf, self.fluxBuf, self:getULRBuf(), self.intVelBuf, self.cellBuf, real(dt))
	
		self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
		self.calcDerivFromFluxKernelObj()
	
		if self.eqn.useSourceTerm then
			self.addSourceKernelObj(self.solverBuf, derivBufObj.obj, self:getULRBuf(), self.cellBuf)
		end
	end)

	self:boundary()

	-- TODO potential update here
	
	self.integrator:integrate(dt, function(derivBufObj)
		self.computePressureKernelObj()
	
		self.diffuseMomentumKernelObj.obj:setArg(1, derivBufObj.obj)
		self.diffuseMomentumKernelObj()
	end)
	
	self:boundary()
	
	self.integrator:integrate(dt, function(derivBufObj)
		self.diffuseWorkKernelObj.obj:setArg(1, derivBufObj.obj)
		self.diffuseWorkKernelObj()
	end)
end

return EulerBurgers
