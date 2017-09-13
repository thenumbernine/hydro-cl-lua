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

-- and now for fillRandom ...
local ffi = require 'ffi'
local function crand() return 2 * math.random() - 1 end
function NumRelEqn:fillRandom(epsilon)
	local solver = self.solver
	local ptr = ffi.new(self.cons_t..'[?]', solver.volume)
	ffi.fill(ptr, 0, ffi.sizeof(ptr))
	for i=0,solver.volume-1 do
		for _,var in ipairs(intVars) do
			local name, ctype = next(var)
			if ctype == 'real' then
				ptr[i][name] = epsilon * crand()
			elseif ctype == 'real3' then
				for j=0,2 do
					ptr[i][name].s[j] = epsilon * crand()
				end
			elseif ctype == 'sym3' then
				for jk=0,5 do
					ptr[i][name].s[jk] = epsilon * crand()
				end
			else
				error("don't know how to handle ctype "..ctype.." for field "..name)
			end
		end
		
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gammaBar_ll.xx = ptr[i].gammaBar_ll.xx + 1
		ptr[i].gammaBar_ll.yy = ptr[i].gammaBar_ll.yy + 1
		ptr[i].gammaBar_ll.zz = ptr[i].gammaBar_ll.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end


return NumRelEqn
