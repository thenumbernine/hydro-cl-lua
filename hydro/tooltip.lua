-- wrapper for imgui stuff to put a tooltip over it (and give it no title)
-- (and push/pop id strs so the no-title doesn't cause a problem)

local ffi = require 'ffi'
local ig = require 'imgui'
local table = require 'ext.table'
require 'ffi.c.string'	-- strlen

-- wrappers for table[key] access rather than ffi-allocated primitive access


-- unlike others, this is a bit more than simple ffi primitive read/write
-- because imgui's float formatting isn't very flexible
-- also not 'tableFloat' because of the higher accuracy of double/string than float/string serializing
local buf = ffi.new'char[256]'
local function numberTable(title, t, k, ...)
	local src = tostring(t[k])
	local len = math.min(ffi.sizeof(buf)-1, #src)
	ffi.copy(buf, src, len)
	buf[len] = 0
	-- TODO maybe ig.ImGuiInputTextFlags_EnterReturnsTrue
	if ig.tooltipInputText(title, buf, ffi.sizeof(buf), ...) then
		-- ffi.string doesn't stop at the null term when you pass a fixed size?
		local s = ffi.string(buf, math.min(ffi.sizeof(buf), tonumber(ffi.C.strlen(buf))))
		local v = tonumber(s)
		if v then
			t[k] = v
			return true
		end
	end
end

-- here's another exception: combo boxes
-- because t[k] is Lua-based, lets make our values 1-based instead of 0-based

-- tooltip wrappers

local tooltip = {
	-- non table suffix = tooltip prefix
	button = ig.tooltipButton,
	checkbox = ig.tooltipCheckbox,
	combo = ig.tooltipCombo,
	float = ig.tooltipInputFloat,
	int = ig.tooltipInputInt,
	slider = ig.tooltipSliderFloat,
	-- table suffix = luatableTooltip prefix
	sliderTable = ig.luatableTooltipSliderFloat,
	checkboxTable = ig.luatableTooltipCheckbox,
	intTable = ig.luatableTooltipInputInt,
	numberTable = numberTable,
	comboTable = ig.luatableTooltipCombo,
	textTable = ig.luatableTooltipInputText,
}

return tooltip
