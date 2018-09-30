--[[
this looks a lot like GRHDBehavior
it does the same thing -- adds 'getADMArgs' and 'getADMVarCode' for the eqn to use
make it a behavior for adding ADM args to eqns
--]]
local class = require 'ext.class'
local template = require 'template'

return function(parent)
	local templateClass = class(parent)
	
	-- eqn/grhd.cl needs this implemented for deducing alpha, beta, gamma
	function templateClass:getADMArgs()
		return ''	
	end

	--[[
	args:
		alpha = (optional) alpha var name
		beta = (optional) beta var name
		gamma = (optional) gamma var name
		suffix = (optional) suffix for default variable names
		index = (optional) index in grid associated with metric. default 'index'
	I made args volatile because I was lazy
	--]]
	function templateClass:getADMVarCode(args)
		args = args or {}
		args.suffix = args.suffix or ''
		args.alpha = args.alpha or ('alpha'..args.suffix)
		args.beta = args.beta or ('beta'..args.suffix)
		args.gamma = args.gamma or ('gamma'..args.suffix)
		return template([[
	real <?=args.alpha?> = 1;
	real3 <?=args.beta?> = real3_zero;
	sym3 <?=args.gamma?> = sym3_ident;
]], {args=args})
	end

	return templateClass
end
