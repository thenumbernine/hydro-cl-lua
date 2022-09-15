require 'ext'
local gnuplot = require 'gnuplot'
local txts = table()
for fn in file:dir() do
	local name,ext = file(fn):getext()
	if ext == 'txt' then
		txts:insert(fn)
	end
end
for _,fn in ipairs(txts) do
	print(fn)
	local name,ext = file(fn):getext()
	assert(name)
	gnuplot{
		output = name..'.png',
		style = 'data lines',
		xlabel = 't',
		ylabel = 'H',
		title = fn,
		{datafile=fn, using='1:2', title='min'},
		{datafile=fn, using='1:3', title='avg'},
		{datafile=fn, using='1:4', title='max'},
	}
	gnuplot{
		output = name..'-log.png',
		style = 'data lines',
		xlabel = 't',
		ylabel = 'H',
		title = fn,
		log = 'y',
		{datafile=fn, using='1:4', title='max'},
	}
end
