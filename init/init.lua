local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local time = table.unpack(require 'util.time')

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
	time('compiling init state program', function()
		solver.initStateUnlinkedObj = solver.Program{
			name = 'initState',
			--code = initStateCode,
			code = initStateCode..'\n'..template(require'ext.io'.readfile'math.cl'),
		}
		solver.initStateUnlinkedObj:compile{
			dontLink = true,
			buildOptions = table{
				'-Werror',
				'-cl-std=CL2.0',
				'-cl-kernel-arg-info',
				'-g',
				'-create-library',
				'-enable-link-options',
			}:concat' ',
		}
	end)
	
	time('linking init state program', function()
		solver.initStateProgramObj = solver.Program{
			programs = {
				--solver.mathUnlinkedObj, 
				solver.initStateUnlinkedObj,
			},
			buildOptions =  table{
				'-Werror',
				'-cl-std=CL2.0',
				'-cl-kernel-arg-info',
				'-g',
				'-create-library',
				'-enable-link-options',
			}:concat' ',
		}
	end)

local cl = require 'ffi.OpenCL'
--not supported?
--print('CL_PROGRAM_NUM_KERNELS', cl.CL_PROGRAM_NUM_KERNELS)
--print('CL_PROGRAM_NUM_KERNELS', solver.initStateProgramObj.obj:getInfo'CL_PROGRAM_NUM_KERNELS')
print('CL_PROGRAM_KERNEL_NAMES', cl.CL_PROGRAM_KERNEL_NAMES)
local s = solver.initStateProgramObj.obj:getInfo'CL_PROGRAM_KERNEL_NAMES'
print('CL_PROGRAM_KERNEL_NAMES', #s, s)
os.exit()
else	-- not useCLLinkLibraries
	local file = require 'ext.file'
	time('building init state program', function()
		solver.initStateProgramObj = solver.Program{
			name = 'initState',
			code = initStateCode,
		}
		solver.initStateProgramObj:compile()
	end)
--[[ works
local cl = require 'ffi.OpenCL'
print('CL_PROGRAM_KERNEL_NAMES', cl.CL_PROGRAM_KERNEL_NAMES)
local s = solver.initStateProgramObj.obj:getInfo'CL_PROGRAM_KERNEL_NAMES'
print('CL_PROGRAM_KERNEL_NAMES', #s, s)
--]]
end
	solver.initStateKernelObj = solver.initStateProgramObj:kernel('initState', solver.solverBuf, solver.UBuf)

	-- here's an ugly hack ...
	-- I need a custom init state kernel for the GLM_MHD only
	-- and it shares init conditions with a lot of other solvers
	-- so ...
	-- (don't the Einstein solvers also use initDerivs?)
	if require 'eqn.glm-mhd'.is(solver.eqn) then
		solver.initDerivsKernelObj = solver.initStateProgramObj:kernel('initDerivs', solver.solverBuf, solver.UBuf)
	end
end

function InitCond:initState(solver)
	return '//no code from InitCond:initState() was provided'
end

-- called when the solver resets
function InitCond:resetState(solver)
	solver.initStateKernelObj()
	if require 'eqn.glm-mhd'.is(solver.eqn) then
		solver:boundary()
		solver.initDerivsKernelObj()
	end
	solver:boundary()
end

return InitCond
