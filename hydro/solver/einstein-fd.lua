local class = require 'ext.class'
local table = require 'ext.table'
local GridSolver = require 'hydro.solver.gridsolver'

-- only needed for cmdline.renderBssnVideo:
local vec2i = cmdline.renderBssnVideo and require 'vec-ffi.vec2i'
local vec3i = cmdline.renderBssnVideo and require 'vec-ffi.vec3i'
local vec3d = cmdline.renderBssnVideo and require 'vec-ffi.vec3d'
local Image = cmdline.renderBssnVideo and require 'image'


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


	-- HACK FOR RENDERING SPHERE-LOG-RADIAL COORD SIMS TO CARTESIAN FOR MOVIE DUMPING
	-- right now the rest of my movie-making is done by capturing OpenGL-rendered content
	-- so maybe this will morph into a more general movie dumping for all simulations
	if cmdline.renderBssnVideo then
		self.xcmin = vec3d(-1.5, 0, -1.5)
		self.xcmax = vec3d(1.5, 0, 1.5)
		self.csize = vec2i(80, 80)

		self.writeOutImage = Image(csize.x, csize.y, 1, 'real')

		local coord = self.coord
		
		assert(require 'hydro.coord.sphere-log-radial'.is(coord))
		local symmath = require 'symmath'
		local rho, theta, phi = coord.baseCoords:unpack()
		
		local env = self.app.env

		self.writeResultBuf = env:buffer{name='writeResult', type='real', count=csize:volume()}

		self.writeResultProgram = env:program{
			name = 'writeResult',
			code = self.modules:getCodeAndHeader(
				self.sharedModulesEnabled:keys():unpack()
			)
			..template([[
// right now I could use coordMapInv but it is hardcoded to GLSL
// but I could ...
// 1) make it a module
// 2) move out sinh, cosh, asinh to dependency modules
// 3) use the 'toGLSL' conversion to provide the old getCoordMapInvGLSLCode()

real3 getCoordMapInv(real3 pt) {
	real r = real3_len(pt);
	real theta = fmod(acos(pt.z / r), M_PI
	real phi = fmod(math.atan2(pt.y, pt.x), 2. * M_PI
	real rho = calcRho(r, theta, phi, coord.amplitude, coord.sinh_w, math)
	return rho, theta, phi
}

kernel void writeResult(
	constant <?=solver.solver_t?>* solver,
	global real* writeResultBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	int dstx = get_global_id(0);
	int dsty = get_global_id(1);
	if (dstx >= <?=csize.x?> || dsty >= <?=csize.y?>) return;
	// write out a blitted frame of [-1.5, 1.5]^2 XZ
	real3 xc = _real3(
		((real)dstx+.5)/<?=csize.x?> * <?=xcmax.x - xcmin.x?> + <?=xcmin.x?>,
		<?= .5 * (xcmin.y + xcmax.y) ?>,
		((real)dsty+.5)/<?=csize.y?> * <?=xcmax.z - xcmin.z?> + <?=xcmin.z?>);
	real3 x = getCoordMapInv(xc);
	real4 xf = (real4)(
		(x.x - solver->mins.x) / (solver->maxs.x - solver->mins.x) * (real)(solver->gridSize.x - 2 * numGhost) + numGhost,
		(x.y - solver->mins.y) / (solver->maxs.y - solver->mins.y) * (real)(solver->gridSize.y - 2 * numGhost) + numGhost,
		(x.z - solver->mins.z) / (solver->maxs.z - solver->mins.z) * (real)(solver->gridSize.z - 2 * numGhost) + numGhost,
		0);
	int4 i = (int4)((int)xf.x, (int)xf.y, (int)xf.z, (int)xf.w);
	real4 s = xf - (real4)((real)i.x, (real)i.y, (real)i.z, (real)i.w);
	real4 t = (real4)(1,1,1,1) - s;	

	real value = 0;	// or some other 'non-mapped cell' value
	if (!OOB(numGhost,numGhost)) {
		int index = INDEXV(i);
		value = t.x * t.y * UBuf[index].W
			  + s.x * t.y * UBuf[index + solver->stepsize.x].W
			  + t.x * s.y * UBuf[index + solver->stepsize.y].W
			  + s.x * s.y * UBuf[index + solver->stepsize.x + solver->stepsize.y].W;
	}
	writeResultBuf[dstx + <?=csize.x?> * dsty] = value;
}
]], 		{
				solver = self,
				eqn = self.eqn,
				coord = self.coord,
				csize = csize,
				xcmin = xcmin,
				xcmax = xcmax,
			}),
		}
		-- crashing without warning when i compile without 'kernel'
		self.writeResultProgram:compile()

		self.writeResultKernel = self.writeResultProgram:kernel{
			name = 'writeResult',
			domain = env:domain{size={csize.x, csize.y}},
		}
		self.writeResultKernel.obj:setArgs(
			self.solverBuf,
			self.writeResultBuf,
			self.UBufObj
		)

	end
end

function EinsteinFiniteDifferenceSolver:refreshSolverProgram()
	EinsteinFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivKernelObj.obj:setArg(2, self.UBuf)
	self.calcDerivKernelObj.obj:setArg(3, self.cellBuf)
end

function EinsteinFiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivKernelObj()
end

-- [=[ commenting out to try to reduce build time
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
--]=]

if cmdline.renderBssnVideo then
	function EinsteinFiniteDifferenceSolver:update()
		EinsteinFiniteDifferenceSolver.super.update(self)

		self.epochCounter = (self.epochCounter or -1) + 1
		if self.epochCounter % 100 == 0 then
			self.writeResultKernel()
			self.writeResultBuf:toCPU(self.writeOutImage.buffer)
			
			self.imageCounter = (self.imageCounter or -1) + 1
			local fn = rundir..'/frame-'..('%05d'):format(self.imageCounter)..'.png'
			print('writing '..fn)
			self.writeOutImage:normalize():rgb():save(fn)
		end
	end
end

return EinsteinFiniteDifferenceSolver
