comparing GEM simulations of Earth using dif charts/components

Cartesian = lookin good
chart=Cylindrical component=cartesian = not too far from baseline
chart=Spherical component=cartesian = error starting to wander away ... boundary r=0 might be faulty for component=Cartesian - it's designed for (an)holonomic components.
chart=Cylindrical component=anholonomic = exploding.  I wonder if my GEM / MHD are correctly implemented for curvilinear components...
chart=Spherical component=anholonomic = exploding really badly

maybe change grids / cfls to adjust for the # of rows produced?
cfl=.1	gridSize=32,32,32	coord=cartesian	vectorComponent=cartesian	produces 2993 rows
cfl=.1	gridSize=32,32,32	coord=cylinder	vectorComponent=cartesian	produces 63307 rows
cfl=.1	gridSize=16,16,16	coord=sphere	vectorComponent=cartesian	produces 97800 rows
cfl=.1	gridSize=32,32,32	coord=cylinder	vectorComponent=anholonomic	produces 60984 rows
cfl=.1	gridSize=16,16,16	coord=sphere	vectorComponent=anholonomic	produces 98803 rows
