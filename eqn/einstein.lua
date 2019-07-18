--[[
common functions for all Einstein field equation solvers
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local Equation = require 'eqn.eqn'

local common = require 'common'
local xNames = common.xNames
local symNames = common.symNames
local from6to3x3 = common.from6to3x3

local EinsteinEquation = class(Equation)

EinsteinEquation.initStates = require 'init.einstein'

function EinsteinEquation:createInitState()
	EinsteinEquation.super.createInitState(self)
	self:addGuiVars{
		{
			type = 'combo',
			name = 'f_eqn',
			options = {
				'2/alpha',							-- 1+log slicing
				'1 + 1/alpha^2', 					-- Alcubierre 10.2.24: "shock avoiding condition" for Toy 1+1 spacetimes 
				'1', 								-- Alcubierre 4.2.50 - harmonic slicing
				'0', '.49', '.5', '1.5', '1.69',
			}
		},
	}
end

-- add an option for fixed Minkowsky boundary spacetime
-- TODO now there is already a BoundaryFixed in solver/gridsolver, but no easy way to parameterize how to set what fixed values it is
function EinsteinEquation:createBoundaryOptions()
	local Boundary = self.solver.Boudary
	local BoundaryFixed = class(Boundary)
	BoundaryFixed.name = 'fixed'
	function BoundaryFixed:getCode(args)
		local lines = table()
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		for _,j in ipairs{'j', gridSizeSide..'-numGhost+j'} do
			local index = args.indexv(j)
			local U = 'buf[INDEX('..index..')]'
			lines:insert(template([[
	setFlatSpace(solver, &<?=U?>, cell_x((int4)(<?=index?>, 0)));
]], 		{
				eqn = eqn,
				U = U,
				index = index,
			}))
		end
		return lines:concat'\n'
	end
	
	self.solver:addBoundaryOption(BoundaryFixed)
end

-- used by bssnok-fd-num and bssnok-fd-sym
-- useful with spherical grids
-- (which no other eqn has attempted to implement yet)
-- signs = array of 1 or -1 based on the index parity wrt reflection about the boundary condition
-- see table III in 2017 Ruchlin
function EinsteinEquation:getParityVars(...)
	local sign = {...}
	local vars = table()
	for _,var in ipairs(self.consVars) do
		if var.type == 'real' then
		elseif var.type == 'real3' then
			for i,xi in ipairs(xNames) do
				if sign[i] == -1 then
					vars:insert(var.name..'.'..xi)
				end
			end
		elseif var.type == 'sym3' then
			for ij,xij in ipairs(symNames) do
				local i,j = from6to3x3(ij)
				if sign[i] * sign[j] == -1 then
					vars:insert(var.name..'.'..xij)
				end
			end
		else
			error"you are here"
		end
	end
	return vars
end

-- and now for fillRandom ...
local ffi = require 'ffi'
local function crand() return 2 * math.random() - 1 end
function EinsteinEquation:fillRandom(epsilon)
	local solver = self.solver
	local ptr = ffi.new(self.cons_t..'[?]', solver.numCells)
	ffi.fill(ptr, 0, ffi.sizeof(ptr))
	for i=0,solver.numCells-1 do
		for j=0,self.numStates-1 do
			ptr[i].ptr[j] = epsilon * crand()
		end
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return EinsteinEquation
