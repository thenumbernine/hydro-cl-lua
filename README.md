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
- Maxwell equations from Trangenstein's book
- ideal MHD from Stone et al 2008
- SRHD from Marti & Muller 2008
- GRHD from Font 2008
- numerical relativity via Bona-Masso formalism described in Alcubierre 1997 and Alcubierre's 2008 book
- numerical relativity via finite difference BSSNOK (Baumgarte & Shapiro 2010)
- self-gravitation for some schemes (Euler equations)

TODO:
- ADM3D needs shift conditions, and moving puncture, etc ...
- FOBSSN would be nice.  something with the equation stability of BSSN and the algorithmic stability of finite-volume.
- for 2D SRHD, it looks like my W = U.D / prim.rho is drifting from W = 1 / sqrt(1 - prim.v^2)
- implement eigen-stuff code in SRHD so that PLM can work 
- better divergence removal
- get curved coordinates to work
- get two-fluid EMHD working (currently has nans)
- add HLL / HLLC solvers
- implement Navier-Stokes, compressible & incompressible
- rename 'adm\_' prefixes to 'nr\_' or 'gr\_' or something else ...
- PLM support that works on a wide range of equations (currently have a few vying options)
- PLM for BSSNOK-FD and Euler-Burgers
- higher-order polynomial stuff - WENO or whatever
- initialize NR stuff to the geometry metric ... or decide what kind of geometry metric to use (holonomic vs anholonomic) ... or just use euclidian components and calculate the normals and volumes and surfaces using geom.
- currently seeing errors when two solvers run simultaneously ... which makes EM+HD difficult
- finish GR+HD by taking the SRHD and giving it the metric from GR
- add EM+SRHD by mixing and matching SRHD and EMHD
- add EM+GR+HD by winning
- Figure out what to do with self-gravitational potential energy in the Euler simulation.  offsetting it positive makes gravitational instability stable.  offsetting it too positive explodes -- even a forward euler integrator (why).  offsetting it negative causes more instability.
- change vector field from immediate mode to buffered geometry, and gometry shaders if they're available
- add different vector field scales for different solvers (use valuemin and valuemax?)
- looks like RHD has multistep errors in RK, probably because it needs more (either boundary() or updatePrimsFromCons()) to be called between substeps
- add source terms to GRHD -- or at least plugins for 'gr-hd-separate' to fill in from the NR solver
- coroutines to iterative solvers?  so they don't stall the app execution?

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
* Numerical Relativity- ADM, BSSN, etc:
	* Alcubierre, Miguel. Introduction to 3+1 Numerical Relativity. Oxford Science Publications, Oxford, 2008.
	* Baumgarte, Shapiro. Numerical Relativity: Solving Einstein's Equations on the Computer, 2010.
	* Alcubierre (1997) "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
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
* MHD Roe:
	* Athena: A New Code for Astrophysical MHD (2008) https://arxiv.org/pdf/0804.0402v1.pdf
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
