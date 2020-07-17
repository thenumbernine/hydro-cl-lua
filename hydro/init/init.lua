local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local time = table.unpack(require 'hydro.util.time')

--[[
name = name of the initial condition
guiVars = any gui variables that the initial conditions wants to expose
overrideGuiVars = any of the equation's gui vars that the initial conditions wants to override
initState = function(self, solver) returns the OpenCL code for the initial conditions
--]]

local InitCond = class()

function InitCond:refreshInitStateProgram(solver)
	solver:buildMathCLUnlinked()

	local initStateCode 
	time('generating init state code', function()
		initStateCode = table{
			solver.codePrefix,
			self.header and self:header(solver) or '',
			solver.eqn:getInitStateCode(),
		}:concat'\n'
	end)

if solver.useCLLinkLibraries then 
	--local code = initStateCode
	local code = initStateCode..'\n'..template(require 'ext.io'.readfile'hydro/math.cl')
	time('compiling init state program', function()
		solver.initStateUnlinkedObj = solver.Program{name='initState', code=code}
		solver.initStateUnlinkedObj:compile{dontLink=true}
	end)
	time('linking init state program', function()
		solver.initStateProgramObj = solver.Program{
			programs = {
				solver.mathUnlinkedObj, 
				solver.initStateUnlinkedObj,
			},
		}
	end)
else	-- not useCLLinkLibraries
	local file = require 'ext.file'
	time('building init state program', function()
		solver.initStateProgramObj = solver.Program{
			name = 'initState',
			code = initStateCode,
		}
		solver.initStateProgramObj:compile()
	end)
end

	-- how to give initCond access to rand()?
	-- fill UBuf with random numbers before calling it
	time('randomizing UBuf...', function()
		local ptr = solver.UBufObj:toCPU()
		for i=0,solver.numCells-1 do
			for j=0,solver.eqn.numStates-1 do
				ptr[i].ptr[j] = math.random()
			end
		end
		solver.UBufObj:fromCPU(ptr)
	end)

	solver.initStateKernelObj = solver.initStateProgramObj:kernel('initState', solver.solverBuf, solver.UBuf, solver.cellBuf)

	-- here's an ugly hack ...
	-- I need a custom init state kernel for the GLM_MHD only
	-- and it shares init conditions with a lot of other solvers
	-- so ...
	-- (don't the Einstein solvers also use initDerivs?)
	if require 'hydro.eqn.glm-mhd'.is(solver.eqn) then
		solver.initDerivsKernelObj = solver.initStateProgramObj:kernel('initDerivs', solver.solverBuf, solver.UBuf, solver.cellBuf)
	end
end

function InitCond:initState(solver)
	return '//no code from InitCond:initState() was provided'
end

-- called when the solver resets
function InitCond:resetState(solver)
	solver.initStateKernelObj()

	if cmdline.printBufs then
		print('init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
	
	if require 'hydro.eqn.glm-mhd'.is(solver.eqn) then
		solver:boundary()
		solver.initDerivsKernelObj()
	end
	solver:boundary()

	if cmdline.printBufs then
		print('post-boundary init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
end

return InitCond
