local class = require 'ext.class'
local table = require 'ext.table'
local FiniteVolumeSolver = require 'hydro.solver.fvsolver'


local Roe = class(FiniteVolumeSolver)
Roe.name = 'Roe'

function Roe:initL1(args)
	Roe.super.initL1(self, args)

	local RoeFlux = require 'hydro.flux.roe'
	self.flux = RoeFlux(self)
end

function Roe:getSolverCode()
	return table{
		Roe.super.getSolverCode(self),
		self.flux:getSolverCode(),
	}:concat'\n'
end

return Roe
