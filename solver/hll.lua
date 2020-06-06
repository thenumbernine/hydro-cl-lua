local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local template = require 'template'
local FiniteVolumeSolver = require 'solver.fvsolver'


local HLL = class(FiniteVolumeSolver)
HLL.name = 'HLL'

HLL.solverCodeFile = 'solver/hll.cl'


--HLL.calcWaveMethod = 'Davis direct'
HLL.calcWaveMethod = 'Davis direct bounded'

function HLL:init(args)
	self.calcWaveMethod = args.calcWaveMethod
	HLL.super.init(self, args)
end

function HLL:getSolverCode()
	return table{
		HLL.super.getSolverCode(self),
	
		-- before this went above solver/plm.cl, now it's going after it ...
		template(file[self.solverCodeFile], {
			solver = self,
			eqn = self.eqn,
			clnumber = require 'cl.obj.number',
		}),
	}:concat'\n'
end

return HLL
