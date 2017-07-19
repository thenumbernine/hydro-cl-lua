--[[
this is the GR Euler Fluid Equation solver 
with the BSSN or ADM3D (or some other NR solver) plugged into it
--]]

local class = require 'ext.class'

--[[
the parent class is your scheme... just Roe right now
--]]
local function GRPlusHDBehavior(parent)
	local templateClass = class()

	templateClass.name = 'GR+HD '..parent.name

	function templateClass:init(args)
		self.app = assert(args.app)
	
		-- TODO copy ion + electron two-fluid solver
		-- TODO better yet, make a behavior for composite solvers
	end
end

return GRPlusHDBehavior
