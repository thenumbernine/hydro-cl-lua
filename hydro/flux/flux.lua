local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local string = require 'ext.string'
local template = require 'template'
local clnumber = require 'cl.obj.number'

local Flux = class()

function Flux:init(args)
	self.solver = assert(args.solver)
end

function Flux:initCodeModules()
	local code = self.solver.eqn:template(file(self.solverCodeFile):read(), {
		clnumber = clnumber,
		flux = self,
	})

	-- in case any MODULE_* cmds were in there, align them back with the lhs
	--[[
	and in case anyone put a trailing \ ...
	TODO think about this ...
	do I want the MODULE_* cmds to respect \ and search for the next line?
	that'll make it tough for the auto-inline to #define code (like in fvsolver.cl)
	that just replaces all \n's with \\\n's ... those will now have to filter-out MODULE_ cmds.
	maybe I do want that later?
	for now I will just remove them here.
	--]]
	code = string.split(code, '\n'):mapi(function(l)
		return (l:gsub('^%s*(//// MODULE_)', '%1'))
	end):mapi(function(l)
		return (l:gsub('^(//// MODULE_.*)\\$', '%1'))
	end):concat'\n'

	self.solver.modules:addFromMarkup(code)
end

return Flux
