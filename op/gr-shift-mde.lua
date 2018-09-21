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

this means the equation will look like ...
gamma^jk D_j D_k beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j = 2 D_j (alpha A^ij)
gamma^jk D_j (beta^i_,k + Gamma^i_lk beta^l) + 1/3 gamma^ik D_k (beta^j_,j + Gamma^j_lj beta^l) + R^i_j beta^j = 2 (D_j alpha A^ij + alpha D_j A^ij)
gamma^jk (D_j beta^i_,k + D_j Gamma^i_lk beta^l + Gamma^i_lk D_j beta^l) + 1/3 gamma^ik (D_k beta^j_,j + D_k Gamma^j_lj beta^l + Gamma^j_lj D_k beta^l) + R^i_j beta^j = 2 (alpha_,j A^ij + alpha (A^ij_,j + Gamma^i_kj A^kj + Gamma^j_kj A^ik))

(
	+ gamma^jk beta^i_,jk 
	+ 1/3 gamma^ik beta^j_,jk 
	+ gamma^jk Gamma^i_lk beta^l_,j 
	- gamma^jk Gamma^l_kj beta^i_,l 
	+ 1/3 gamma^ik Gamma^j_lj beta^l_,k 
	+ 1/3 gamma^ik beta^l Gamma^j_lj,k 
	+ gamma^jk beta^l Gamma^i_lk,j 
	+ gamma^jk Gamma^i_lj beta^l_k 
	+ gamma^jk beta^l Gamma^i_mj Gamma^m_lk 
	- gamma^jk beta^l Gamma^m_lj Gamma^i_mk 
	- gamma^jk beta^l Gamma^m_kj Gamma^i_lm 
	+ gamma^jk Gamma^i_lk Gamma^l_mj beta^m
	+ 1/3 gamma^ik beta^l Gamma^j_mk Gamma^m_lj 
	- 1/3 gamma^ik beta^l Gamma^m_lk Gamma^j_mj 
	- 1/3 gamma^ik beta^l Gamma^m_jk Gamma^j_lm
	+ 1/3 gamma^ik Gamma^j_lj Gamma^l_mk beta^m
	+ R^i_j beta^j 
) = 2 alpha_,j A^ij + 2 alpha A^ij_,j + 2 alpha Gamma^i_kj A^kj + 2 alpha Gamma^j_kj A^ik


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
