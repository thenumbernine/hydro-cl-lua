local os = require 'ext.os'
local function exec(s)
	print('>'..s)
	return os.execute(s)
end
return exec
