-- wrapper for imgui stuff to put a tooltip over it (and give it no title)
-- (and push/pop id strs so the no-title doesn't cause a problem)

local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local table = require 'ext.table'

local function hoverTooltip(name)
	if ig.igIsItemHovered() then
		ig.igBeginTooltip()
		ig.igText(name)
		ig.igEndTooltip()
	end
end

local function makeWrapTooltip(f)
	return function(name, ...)
		ig.igPushIdStr(name)
		local result = f('', ...)
		hoverTooltip(name)
		ig.igPopId()
		return result
	end
end

local function tooltipLabel(label, str)
	ig.igPushIdStr(label)
	ig.igText(str)
	hoverTooltip(label)
	ig.igPopId()
end

-- naive wrappers of makeWrapTooltip

local wrap = table.map({
	slider = ig.igSliderFloat,
	combo = ig.igCombo,
	button = ig.igButton,
	float = ig.igInputFloat,
	int = ig.igInputInt,
	checkbox = ig.igCheckbox,
	text = ig.igInputText,
}, function(funcName, wrapName)
	return makeWrapTooltip(funcName), wrapName
end)

-- wrappers for table[key] access rather than ffi-allocated primitive access

local function makeTableAccess(prim, orig)
	local ptr = ffi.new(prim..'[1]')
	return function(title, t, k, ...)
		ptr[0] = t[k]
		if orig(title, ptr, ...) then
			t[k] = ptr[0]
			return true
		end
	end
end

local checkboxTable = makeTableAccess('bool', wrap.checkbox)
local intTable = makeTableAccess('int', wrap.int)
local sliderTable = makeTableAccess('float', wrap.slider)

-- unlike others, this is a bit more than simple ffi primitive read/write
-- because imgui's float formatting isn't very flexible
-- also not 'tableFloat' because of the higher accuracy of double/string than float/string serializing
local buf = ffi.new'char[256]'
local function numberTable(title, t, k, ...)
	local s = tostring(t[k])
	ffi.copy(buf, s, ffi.sizeof(buf))
	-- TODO maybe ig.ImGuiInputTextFlags_EnterReturnsTrue
	if wrap.text(title, buf, ffi.sizeof(buf), ...) then
		local s = ffi.string(buf, ffi.sizeof(buf))
		local v = tonumber(s)
		if v then
			t[k] = v
			return true
		end
	end
end

-- here's another exception: combo boxes
-- because t[k] is Lua-based, lets make our values 1-based instead of 0-based

local int = ffi.new'int[1]'
local function comboTable(title, t, k, ...)
	int[0] = t[k]-1
	if wrap.combo(title, int, ...) then
		t[k] = int[0]+1
		return true
	end
end

-- tooltip wrappers

local tooltip = {
	button = wrap.button,
	checkbox = wrap.checkbox,
	combo = wrap.combo,
	float = wrap.float,
	int = wrap.int,
	slider = wrap.slider,
	sliderTable = sliderTable,
	checkboxTable = checkboxTable,
	intTable = intTable,
	numberTable = numberTable,
	comboTable = comboTable,
	text = wrap.text,
	label = tooltipLabel,
}

return tooltip
