local class = require 'ext.class'
local file = require 'ext.file'
local template = require 'template'
local Relaxation = require 'hydro.op.relaxation'

local MinimalDistortionEllipticShift = class(Relaxation)

MinimalDistortionEllipticShift.name = 'MDEShift'

MinimalDistortionEllipticShift.solverCodeFile = 'op/gr-solver-mde.cl'

MinimalDistortionEllipticShift.potentialField = 'betaLap_u'
function MinimalDistortionEllipticShift:getPotBufType() return 'real3' end
function MinimalDistortionEllipticShift:getPotBuf() return self.solver.UBuf end

return MinimalDistortionEllipticShift 
