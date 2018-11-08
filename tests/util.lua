local table = require 'ext.table'
local range = require 'ext.range'

--[[
outer({{a=1}, {a=2}}, {{b=1}, {b=2}})
returns {{a=1,b=1}, {a=1,b=2}, {a=2,b=1}, {a=2,b=2}}
--]]
local function outer(...)
	local args = {...}
	local n = #args
	local is = table()
	for i=1,n do is[i] = 1 end
	local ts = table()
	local done
	repeat
		ts:insert(table( is:map(function(j,i) return args[i][j] end):unpack() ))
		for j=1,n do
			is[j] = is[j] + 1
			if is[j] > #args[j] then
				is[j] = 1
				if j == n then
					done = true
					break
				end
			else
				break
			end
		end
	until done
	return ts
end

local function run(...)
	print(...)
	return os.execute(...)
end

local function nameForConfig(cfg, args)
	return table(
		args.eqn and {'eqn='..args.eqn} or nil
	):append{
		'solver='..cfg.solver,
		'integrator='..args.integrator,
	}:append(args.usePLM 
		-- plm:
		and table{
			'plm='..args.usePLM,
		}:append(
			args.slopeLimiter and {'slopeLimiter='..args.slopeLimiter} or nil
		)
		-- non-plm: use flux limiter
		or (
			args.fluxLimiter and {'fluxLimiter='..args.fluxLimiter} or nil
		)
	):append{	
		'init='..args.initState,
		'gridSize='..range(args.dim):map(function(i) return args.gridSize[i] end):concat'x',
	}:append{
		cfg.movieStartTime and ('t0='..cfg.movieStartTime) or nil,
	}:append{
		cfg.movieEndTime and ('t1='..cfg.movieEndTime) or nil,
	}:append{
		cfg.movieFrameDT and ('dt='..cfg.movieFrameDT) or nil,
	}:concat', '
end

return {
	outer = outer,
	run = run,
	nameForConfig = nameForConfig,
}
