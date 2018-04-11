HydroJS was the first
then Hydro, a C++/multithread version of HydroJS
then HydroGPU, a OpenCL version of Hydro
then with Lua script config, 
then the Lua got out of hand until the C++ was doing nothing but managing strings
now this project, lua-hydro-cl, pushes the middleman (C++) out completely


Features:
- script-generated OpenCL GPGPU code regenerated on the fly as soon as you change GUI options
- automatic tensor index representation of equations / symbolic differentiation (via symmath project)
- 1D, 2D, 3D simulations and visualizations
- solvers are usually a Roe scheme, though some implementations vary
- various flux limiters
- PLM for certain solvers
- various boundary conditions
- integrators: Forward Euler, several Runge-Kutta and Runge-Kutta TVD, and implicit linearized GMRES on the GPU
- GUI-driven everything.  no more restarting the program to switch solvers or initial conditions.
- Euler equations from Toro's book (with some modifications for curved geometry) 
- Maxwell equations from Trangenstein's book with poisson solver constraints
- Maxwell equations with GLM from 2000 Munz
- ideal MHD from Stone et al 2008
- two-fluid electron/ion plasma model from 2014 Abgrall, Kumar
- SRHD from Marti & Muller 2008
- GRHD from Font 2008
- numerical relativity via Bona-Masso formalism described in Alcubierre 1997 and Alcubierre's 2008 book
- numerical relativity via finite difference BSSNOK (Baumgarte & Shapiro 2010)
- self-gravitation for some schemes (Euler equations)
- Z4c from Cao, Hilditch 2011

Example Videos:

[![Rotating Black Hole / Ergosphere Formation](http://img.youtube.com/vi/i2rCutSayXE/0.jpg)](https://www.youtube.com/watch?v=i2rCutSayXE)

[![3D Rotating Black Hole](http://img.youtube.com/vi/LrL-HKcNZLM/0.jpg)](https://www.youtube.com/watch?v=LrL-HKcNZLM)

[![3D Alcubierre warp bubble](http://img.youtube.com/vi/QFwMSj1485M/0.jpg)](https://www.youtube.com/watch?v=QFwMSj1485M)

[![3D Alcubierre warp bubble](http://img.youtube.com/vi/tfMLMxdRid8/0.jpg)](https://www.youtube.com/watch?v=tfMLMxdRid8)


TODO:
- ADM3D with shift as a hyperbolic conservation law system
- ADM3D (and BSSNOK, and any other GR solver) for minimal-distortion elliptical shift solved as a Poisson equation -- which doesn't require extra time-iterating variables. 
- GR horizon tracking / moving puncture
- FOBSSN would be nice.  Something with the analytic stability of BSSN and the algorithmic stability of finite-volume.
- Z4 ... I need to finish typing in the source terms.  I also need a shift condition.
- implement eigen-stuff code in SRHD so that PLM can work 
- PLM for BSSNOK-FD and Euler-Burgers
- PPM
- higher-order polynomial stuff - WENO or whatever
- better divergence removal (multigrid)
- finish GLM-(ideal)MHD ... especially implement the source term as a second step with an exp(dt)  (which I'm not doing at the moment)
- rename mhd to ideal-mhd
- how about a GLM method for Maxwell equations, so I can remove the divergence-free constraint
- get curved grid coordinates to work (cylindrical, sphere 1D radial, sphere 2D surface, sphere r+phi etc)
- get two-fluid-separate EMHD working, so I can update the glm-maxwell with an implicit and update the ion and electron with an explicit solver
- currently seeing errors when two solvers run simultaneously ... which makes EM+HD difficult
- add HLLC/D solvers
- implement Navier-Stokes, compressible & incompressible
- initialize NR stuff to the geometry metric ... or decide what kind of geometry metric to use (holonomic vs anholonomic) ... or just use euclidian components and calculate the normals and volumes and surfaces using geom.
- test out the GR+HD solvers
- add source terms to GRHD -- or at least plugins for 'gr-hd-separate' to fill in from the NR solver
- finish the GR+EM solver
- add EM+GR+HD by winning
- Figure out what to do with self-gravitational potential energy in the Euler simulation.  offsetting it positive makes gravitational instability stable.  offsetting it too positive explodes -- even a forward euler integrator (why).  offsetting it negative causes more instability.
- change vector field from immediate mode to buffered geometry, and gometry shaders if they're available
- coroutines to iterative solvers?  so they don't stall the app execution?
- RHD W error in >1 dimension
- GR flat space simulations make an initial wave.  but shouldn't flat space be stable?

### Sources:

* Hydrodynamics:
	* Duellemond, 2009. Lecture on Hydrodynamics II http://www.mpia-hd.mpg.de/homes/dullemon/lectures/hydrodynamicsII/ 
	* Masatsuka, I Do Like CFD.  http://www.cfdbooks.com/cfdcodes.html 
	* Toro, Eleuterio F. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction. Springer, Germany, 1999. 2nd Edition.
	* http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf
* Flux Limiters:
	* https://en.wikipedia.org/wiki/Flux_limiter
* Electromagnetics:
	* Trangenstein "Numerical Simulation of Hyperbolic Conservation Laws"
* Einstein Field Equations- ADM, BSSN, etc:
	* Alcubierre, Miguel. Introduction to 3+1 Numerical Relativity. Oxford Science Publications, Oxford, 2008.
	* Baumgarte, Shapiro. Numerical Relativity: Solving Einstein's Equations on the Computer, 2010.
	* Alcubierre (1997) "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
	* Bona, Palenzuela-Luque, Bona-Casas.  Elements of Numerical Relativity and Relativistic Hydrodynamics, 2009.
* Stellar Schwarzschild initial conditions:
	* Misner, Thorne, Wheeler. Gravitation, 1973
* SRHD:
	* Marti, J. M. and Muller, E. Numerical Hydrodynamics in Special Relativity Living Reviews in Relativity 6 (2003), 7 http://relativity.livingreviews.org/Articles/lrr-2003-7
	* Sheck, Aloy, Marti, Gomez, Muller Does the plasma composition affect the long-term evolution of relativistic jets? Monthly Notices of Royal Astronomical Society 331, 615-634 2002.
	* Anton, Luis; Zanotti, Olindo; Miralles, Juan; Marti, Jose; Ibanez, Jose; Font, Jose; Pons, Jose. Numerical 3+1 General Relativistic Magnetohydrodynamics: A Local Characteristic Approach. February 2, 2008 https://arxiv.org/abs/astro-ph/0506063
	* Font (2008) "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity"
* HLLC:
	* http://math.lanl.gov/~shenli/publications/hllc_mhd.pdf
	* http://marian.fsik.cvut.cz/~bodnar/PragueSum_2012/Toro_2-HLLC-RiemannSolver.pdf
* ideal MHD Roe:
	* Athena: A New Code for Astrophysical MHD (2008) https://arxiv.org/pdf/0804.0402v1.pdf
* two-fluid plasma model:
	* Abgrall, Kumar 2014
* MUSCL:
	* https://en.wikipedia.org/wiki/MUSCL_scheme
* Euler Initial Conditoins:
	* http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
	* http://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html
	* http://www.astro.virginia.edu/VITA/ATHENA/ot.html
	* http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
	* http://www.cfd-online.com/Wiki/Explosion_test_in_2-D
	* http://www.astro.virginia.edu/VITA/ATHENA/dmr.html
	* http://www.cfd-online.com/Wiki/2-D_laminar/turbulent_driven_square_cavity_flow
* MHD initial conditions:
	* Brio, M. & C.C. Wu, "An Upwind Differencing Scheme for the Equations of Ideal Magnetohydrodynamics", Journal of Computational Physics, 75, 400-422 (1988). The test is described in Section V.
* Runge Kutta & TVD RK:
	* https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Classic_fourth-order_method	
	* http://www.ams.org/journals/mcom/1998-67-221/S0025-5718-98-00913-2/S0025-5718-98-00913-2.pdf
	* http://lsec.cc.ac.cn/lcfd/DEWENO/paper/WENO_1996.pdf
* OpenCL reduce:
	* http://developer.amd.com/resources/documentation-articles/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
