local class = require 'ext.class'
local template = require 'template'
local Relaxation = require 'hydro.op.relaxation'

local MinimalDistortionEllipticShift = class(Relaxation)

MinimalDistortionEllipticShift.name = 'MDEShift'

-- solver-mde or shift-mde?
MinimalDistortionEllipticShift.solverCodeFile = 'hydro/op/gr-solver-mde.cl'

MinimalDistortionEllipticShift.potentialField = 'betaLap_u'

function MinimalDistortionEllipticShift:getPotBufType() return 'real3' end
function MinimalDistortionEllipticShift:getPotBuf() return self.solver.UBuf end

function MinimalDistortionEllipticShift:getSymbolFields()
	return MinimalDistortionEllipticShift.super.getSymbolFields(self):append{
		'solveMinimalDistortionEllipticShift',
	}
end

return MinimalDistortionEllipticShift 
