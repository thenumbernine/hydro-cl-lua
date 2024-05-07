local table = require 'ext.table'
local InitCond = require 'hydro.init.init'

local LatticeBoltzmannInitCond = InitCond:subclass()

local initConds = table{
	{
		name = 'cylinder',
		mins = {-1, -.25, -1},
		maxs = {1, .25, 1},
		guiVars = {
			{name='r', value=.1},
			{name='cx', value=-.8},
			{name='cy', value=0},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			-- TODO custom boundary.  rhs is set to zero.  lhs is U[-2] = U[2], U[-1] = U[1], and U[0] is not modified
			solver:setBoundaryMethods'freeflow'
			-- solver_macros has M_PI
			return solver.eqn:template[[
//// MODULE_DEPENDS: <?=solver_macros?>
// TODO put this in lua ext.math?
#define DBL_EPSILON 2.220446049250313080847e-16
#define DBL_EPS_COMP (1. - DBL_EPSILON)

solid = (real)(real3_len(real3_sub(xc, _real3(initCond->cx, initCond->cy, 0.))) < initCond->r);

<?
for i,ofs in ipairs(solver.offsets) do 
	local c = ofs.c
	local ofsindex = i-1
?>{
	real randunit = U->F<?=ofsindex?>;
	real const s = randunit > .5 ? -1 : 1;
	randunit = fmod(randunit * 2., 1.);
	real const mu = 0.;
	real const sigma = 1.;
	real const randn = mu + s * sigma * sqrt(-log(DBL_EPS_COMP * randunit));

	F[<?=ofsindex?>] = 1. + .01 * randn;
	if (<?=c.x?> == 1 && <?=c.y?> == 0 && <?=c.z?> == 0) {
		F[<?=ofsindex?>] += 2 * (1 + .2 * cos(2. * M_PI * xc.x));
	}

	rho += F[<?=ofsindex?>];
}
<?
end
?>

real const rho0 = 100.;

<?
for i,ofs in ipairs(solver.offsets) do 
	local c = ofs.c
	local ofsindex = i-1
?>	F[<?=ofsindex?>] *= rho0 / rho;
<?
end
?>
rho = rho0;
]]
		end,
	},
}:mapi(function(cl)
	return LatticeBoltzmannInitCond:subclass(cl)
end)

function LatticeBoltzmannInitCond:getList()
	return initConds
end

return LatticeBoltzmannInitCond 
