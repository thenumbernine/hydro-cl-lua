local class = require 'ext.class'

local InitState = class()

function InitState:init(args)
	for k,v in pairs(args) do self[k] = v end
end

return InitState
