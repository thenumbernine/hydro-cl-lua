local table = require 'ext.table'
--[[
usage: processCL(code, {a=1, b=2, ...})
--]]
local function processCL(code, args)
	local argKeys, argValues = table(), table()
	for k,v in pairs(args) do
		argKeys:insert(k)
		argValues:insert(v)
	end
	local outputFunc = '__output'
	local newcode = table{
		'local '..table():append({outputFunc},argKeys):concat', '..' = ...\n',
	}
	local function addprint(from,to)
		local block = code:sub(from,to)
		local eq = ('='):rep(5)	-- TODO make sure no such [=..=[ appears in the code block
		local nl = block:find'\n' and '\n' or ''
		newcode:insert(outputFunc..' ['..eq..'['..nl..block..']'..eq..']\n')
	end
	local pos = 1
	while true do
		local start1, start2 = code:find('<%?', pos)
		if not start1 then
			addprint(pos, #code)
			break
		else
			local ret
			if code:sub(start2+1,start2+1) == '=' then 
				ret = true
				start2 = start2 + 1
			end
			local end1, end2 = code:find('%?>', start2+1)
			end1 = end1 or #code+1
			end2 = end2 or #code-1
			addprint(pos, start1-1)
			local block = code:sub(start2+1, end1-1)
			if ret then
				newcode:insert(outputFunc..'(tostring('..block..'))\n')
			else
				newcode:insert(block..'\n')
			end
			pos = end2+1
			if pos >= #code then break end
		end
	end

	newcode = newcode:concat()
	local f, msg = loadstring(newcode)
	if not f then
		print(require 'showcode'(newcode))
		error(msg)
	end
	
	local outputStrs = table()
	f(function(str)
		outputStrs:insert(str)
	end, argValues:unpack())
	return outputStrs:concat()
end

return processCL
