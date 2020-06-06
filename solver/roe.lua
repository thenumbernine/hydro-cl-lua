local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local math = require 'ext.math'
local file = require 'ext.file'
local template = require 'template'
local FiniteVolumeSolver = require 'solver.fvsolver'

local common = require 'common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local Roe = class(FiniteVolumeSolver)
Roe.name = 'Roe'

function Roe:createBuffers()
	Roe.super.createBuffers(self)

	-- to get sizeof
	ffi.cdef(self.eqn:getEigenTypeCode())
end

function Roe:getSolverCode()
	return table{
		Roe.super.getSolverCode(self),
	
		-- before this went above solver/plm.cl, now it's going after it ...
		template(file['solver/roe.cl'], {
			solver = self,
			eqn = self.eqn,
			clnumber = require 'cl.obj.number',
		}),
	}:concat'\n'
end

return Roe
