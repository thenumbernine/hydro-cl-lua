--[[
this will use an inverse Laplace equation solver to converge the shift 

using eqns Alubierre's book 4.3.14, 4.3.15 ...

D^2 beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j = 2 D_j (alpha A^ij)
A^ij = K^ij - 1/3 K gamma^ij
solve for beta^i

here's my thoughts ...
1) calculate Gamma^i_jk
2) calculate R^i_j
3) solve inverse linear system

this means the 

looks like I'll need a grid of gamma^ij, Gamma^i_jk, and R_ij
--]]
local class = require 'ext.class'

local MinimalDistortionEllipticShift = class()

function MinimalDistortionEllipticShift:init(args)
	self.solver = args.solver

	--
end

function MinimalDistortionEllipticShift:getSolverCode()

end

function MinimalDistortionEllipticShift:refreshSolverProgram()

end

function MinimalDistortionEllipticShift:refreshBoundaryProgram()

end

function MinimalDistortionEllipticShift:resetState()

end

return MinimalDistortionEllipticShift 
