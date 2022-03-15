local class = require 'ext.class'
local table = require 'ext.table'
local InitCond = require 'hydro.init.init'

local TwoFluidInitCond = class(InitCond)

local initConds = table{
	{
		name = 'two-fluid Brio-Wu',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 2
			end
			return [[
	ion_rho = lhs ? 1 : .125;
	ion_P = lhs ? 1 : .1;
	elec_rho = lhs ? 1 : .125;
	elec_P = lhs ? 1 : .1;
	B.x = .75;
	B.y = lhs ? 1 : -1;
	B.z = 0;
]]
		end,
	},

	-- 2008 Johnson, "The Geospace Environmental Modeling (GEM) Magnetic Reconnection Challenge Problem"
	-- section 2.4
	-- electron pressure diverges in just a single point, then slowly moves to the south boundary, and creates a singularity
	{
		name = 'GEM challenge',
		mins = {-4*math.pi, -2*math.pi, -1},
		maxs = {4*math.pi, 2*math.pi, 1},
		depends = {'sech'},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			solver:setBoundaryMethods{
				xmin = 'periodic',
				xmax = 'periodic',
				ymin = 'mirror',	-- TODO boundary for zero derivative
				ymax = 'mirror',
				zmin = 'freeflow',
				zmax = 'freeflow',
			}
			-- mi/me = 25
			-- mu_0 = 1
			-- c = 10 B_0 (B_0 = Alfven speed)
			-- gyroradius of ion, Debye length, light speed, all are 1
			return [[
	real const L_x = 8. * M_PI;	//should correspond with the boundary
	real const L_y = 4. * M_PI;	//should correspond with the boundary
	real const lambda = .5;
	real const B_0 = 1.;
	real const n_0 = 1.;
	real const phi_0 = B_0 / 10.;
	B.x = B_0 * tanh(x.y / lambda);
	real phi = phi_0 * cos(2. * M_PI * x.x / L_x) * cos(2. * M_PI * x.y / L_y);
	real dphi_dx = -2. * M_PI * phi_0 / L_x * sin(2. * M_PI / L_x * x.x) * cos(2. * M_PI / L_y * x.y);
	real dphi_dy = -2. * M_PI * phi_0 / L_y * cos(2. * M_PI / L_x * x.x) * sin(2. * M_PI / L_y * x.y);
	//delta B = -e_z cross nabla phi = [dphi/dy, -dphi/dx, 0]
	B.x += dphi_dy;
	B.y -= dphi_dx;
	real n = n_0 * (.2 + sech(x.y/lambda) * sech(x.y/lambda));
	real ion_n = n_0;
	real elec_n = n_0;
	//n = rho * dV ?
	ion_rho = 1. / ion_n;
	elec_rho = 1. / elec_n;
	real P_ = B_0 * B_0 * n / (2 * n_0);
	elec_P = P_ / 6.;
	ion_P = P_ * 5. / 6.;
]]
		end,
	},

}:mapi(function(cl)
	return class(TwoFluidInitCond, cl)
end)

function TwoFluidInitCond:getList()
	return initConds
end

return TwoFluidInitCond 
