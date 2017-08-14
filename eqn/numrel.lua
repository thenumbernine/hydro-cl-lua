--[[
common functions for all num rel equations
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local Equation = require 'eqn.eqn'

local NumRelEqn = class(Equation)

local xNames = table{'x', 'y', 'z'}

NumRelEqn.initStates = require 'init.numrel'

function NumRelEqn:createInitState()
	NumRelEqn.super.createInitState(self)
	self:addGuiVar{
		type = 'combo',
		name = 'f',
		options = {
			'2/alpha',	-- 1+log slicing
			'1 + 1/alpha^2', 	-- Alcubierre 10.2.24: "shock avoiding condition" for Toy 1+1 spacetimes 
			'1', 		-- Alcubierre 4.2.50 - harmonic slicing
			'0', '.49', '.5', '1.5', '1.69',
		}
	}
end

-- add an option for fixed Minkowsky boundary spacetime
function NumRelEqn:createBoundaryOptions()
	self.solver.boundaryOptions:insert{
		fixed = function(args)
			local lines = table()
			local gridSizeSide = 'gridSize_'..xNames[args.side]
			for _,U in ipairs{
				'buf['..args.index'j'..']',
				'buf['..args.index(gridSizeSide..'-numGhost+j')..']',
			} do
				lines:insert(template([[
	setFlatSpace(&<?=U?>);
]], {eqn=eqn, U=U}))
			end
			return lines:concat'\n'
		end,
	}
end

return NumRelEqn
