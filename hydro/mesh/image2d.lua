local table = require 'ext.table'
local Image = require 'image'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Image2DMeshFactory = Quad2DMeshFactory:subclass()

Image2DMeshFactory.name = 'image2d'

function Image2DMeshFactory:init(args)
	self.image = Image(assert(args.image, "expected image"))

	args = table(args)
	args.size = {self.image.width, self.image.height, args.extrude or 1}

	Image2DMeshFactory.super.init(self, args)
end

function Image2DMeshFactory:testMakeCell(ix, iy)
	return self.image(ix, iy) > .5
end

return Image2DMeshFactory 
