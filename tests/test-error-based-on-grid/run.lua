#!/usr/bin/env lua
-- run some test for {2^n, 2^n, 1}
-- plot error

require 'ext'
local exec = require 'exec'

local DIR = file:cwd()
assert(file'../..':cd())

local writeOut = false
local writeFits = true

for _,i in ipairs(range(1,8):mapi(function(i) return math.floor(2^i) end)) do
	-- hmm, save on exit command?  saving fits files?
	local cmd = './run.lua sys=console "config='..DIR..'/config.lua" "gridSize={'..i..','..i..',1}" exitTime=1'

	if writeFits then
		cmd = cmd .. ' "saveOnExit=results"' -- should produce results.fits
	end

	if writeOut then
		cmd = cmd .. ' "trackvars=U rho" > out.txt'
	end
	
	exec(cmd)
	
	if writeOut then
		file'out.txt':move(DIR..'/out-'..i..'.txt')
	end
	if writeFits then
		file'results_UBuf.fits':move(DIR..'/U-'..i..'.fits')
	end
end

if writeFits then
	for f in file:dir() do
		if f:match'%.fits$' then
			exec('rm "'..f..'"')
		end
	end
end
