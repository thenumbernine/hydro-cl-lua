#!/usr/bin/env luajit
-- plot the var-ranges.txt created by running the app.lua with 'output variable ranges' set

require 'ext'
local matrix = require 'matrix'
local _ = range
local ls = file['var-ranges.txt']:trim():split'\n'

local cols = ls:remove(1):split'\t'

local m = matrix(ls:map(function(l)
 	return l:split'\t':map(function(w)
		return w == 'inf' and math.huge 
		or (w == '-inf' and -math.huge 
		or assert(tonumber(w), "failed to convert to number "..tostring(w))) 
	end) 
end)):transpose()

--[[ using gnuplot

cols = range((#cols-1)/2):map(function(i)
	local name = cols[2*i-1]
	if i > 1 then
		name = assert(name:match('(.*)_max')):gsub('_', ' ')
	end
	return name
end)

require 'gnuplot'(table({
	output = 'gnuplot.png', 
	--style = 'data lines', 
	--style = 'fill transparent solid 0.2 noborder', 
	style = 'fill transparent pattern 4 bo', 
	data = m(_,_(1,3000)), 
		
}, cols:sub(2):map(function(title, i)
	return {
		using = table{1,2*i,2*i+1}:concat':',
		title = title,
		with = 'filledcurves',
	}
end)))

--]]
-- [[ using interactive plot

m = m(_, _(1, 3199))

local t = m[1]
require 'plot2d'(cols:sub(2):map(function(title,i)
	local s = matrix{t, m[i+1]}
	s.enabled = true
	return s, title
end))

--]]
