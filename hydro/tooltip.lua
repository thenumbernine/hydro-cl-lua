-- wrapper for imgui stuff to put a tooltip over it (and give it no title)
-- (and push/pop id strs so the no-title doesn't cause a problem)

local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local table = require 'ext.table'

local function hoverTooltip(name)
	if ig.igIsItemHovered(ig.ImGuiHoveredFlags_None) then
		ig.igBeginTooltip()
		ig.igText(name)
		ig.igEndTooltip()
	end
end

local function makeWrapTooltip(f)
	return function(name, ...)
		ig.igPushID_Str(name)
		local result = f('', ...)
		hoverTooltip(name)
		ig.igPopID()
		return result
	end
end

local function tooltipLabel(label, str)
	ig.igPushID_Str(label)
	ig.igText(str)
	hoverTooltip(label)
	ig.igPopID()
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
}, function(f, wrapName)
	return makeWrapTooltip(f), wrapName
end)

-- wrappers for table[key] access rather than ffi-allocated primitive access

local function makeTableAccess(prim, orig)
	local ptr = ffi.new(prim..'[1]')
	return function(title, t, k, ...)
		if t[k] == nil then
			error("failed to find value "..k.." in table "..tostring(t))
		end
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

local buf = ffi.new'char[256]'

-- unlike others, this is a bit more than simple ffi primitive read/write
-- because imgui's float formatting isn't very flexible
-- also not 'tableFloat' because of the higher accuracy of double/string than float/string serializing
local function numberTable(title, t, k, ...)
	local src = tostring(t[k])
	local len = math.min(ffi.sizeof(buf)-1, #src)
	ffi.copy(buf, src, len)
	buf[len] = 0
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
assert(t[k])
assert(type(t[k]) == 'number')
	int[0] = t[k]-1
	if wrap.combo(title, int, ...) then
		t[k] = int[0]+1
		return true
	end
end

-- TODO dynamic sized buffer?
local function textTable(title, t, k, ...)
	local src = tostring(t[k])
	local len = math.min(ffi.sizeof(buf)-1, #src)
	ffi.copy(buf, src, len)
	buf[len] = 0
	if wrap.text(title, buf, ffi.sizeof(buf), ...) then
		t[k] = ffi.string(buf)
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
	textTable = textTable,
	label = tooltipLabel,
}

return tooltip
