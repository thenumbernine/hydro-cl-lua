local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local template = require 'template'
local Relaxation = require 'op.relaxation'

local Poisson = class(Relaxation)

Poisson.name = 'Poisson'

Poisson.solverCodeFile = 'op/poisson.cl'

function Poisson:getSolverCode()
	return table{
		Poisson.super.getSolverCode(self),
		self:getPoissonCode() or '',
	}:concat'\n'
end

return Poisson
