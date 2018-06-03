local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local Geometry = require 'geom.geom'

local Tensor = symmath.Tensor

local _1D_Radial = class(Geometry)
_1D_Radial.name = '1d_radial'
_1D_Radial.coords = {'r', 'θ', 'φ'}

--[[
fakedim = what dim to simulate the radially symmetric 
--]]
function _1D_Radial:init(args)
	local fakedim = (args.fakedim or 3) - 1

	assert(args.solver.dim == 1)
	
	args.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, theta, phi = symmath.vars('r', 'θ', 'φ')
	
	-- n-D sphere radial only
	-- holonomic or anholonomic
	args.coords = table{r, theta, phi}

	--[[
	nD radial metric has the following properties ...

	g_rr = 1
	g_{theta_i,theta_j} = delta_ij * r^2 * sin theta_1 * ... * sin theta_i
	
	...fixed along theta_i = 0 ...
	g_rr = 1
	g_{theta_i,theta_j = delta_ij r^2
	--]]
	args.chart = function() return Tensor('^I', r^fakedim, 0, 0) end

	_1D_Radial.super.init(self, args)
end

return _1D_Radial 
