-- https://www.cfd-online.com/Forums/blogs/praveen/315-flux-computation-unstructured-grids.html
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local Rusanov = Flux:subclass()

Rusanov.name = 'rusanov'
Rusanov.solverCodeFile = 'hydro/flux/rusanov.clcpp'

return Rusanov
