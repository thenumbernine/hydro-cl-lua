test the error in the scheme based on curvilinear coordinates

tests constant flow problem with wave equation, so if the outer boundary is freeflow then it should be in a steady state

things to compare

vectorComponent=holonomic
vectorComponent=anholonomic
vectorComponent=cartesian

coord == cartesian

coord == cylinder or coord == sphere or coord == sphere_sinh_radial
	- rmin != 0
	- rmin = 0
		- bounds.rmin = none
		- bounds.rmin = cylinderRMin (assumes phi is full range, steps across r=0 and flips components accordingly)
	
for coord == sphere or coord == sphere_sinh_radial:
	- θmin != 0
	- θmin = 0
		- bounds.θmin & θmax = 'none'
		- bounds.θmin & θmax = 'sphereTheta'
	- also same rmin options, either 'none' or 'sphereRMin'
