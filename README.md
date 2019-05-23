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
- Euler equations from Toro's book (with some modifications for curvilinear coordinate systems) 
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
- nonlinear Schrodinger equation from Colliander, Simpso, Sulem "Numerical Simulations of the Energy-Supercritical nonlinear Schrodinger equation" 2010 
- Modular system.  Integrators go in the 'int' folder.  Schemes go in the 'solver' folder.  Equations go in the 'eqn' folder.  The idea is that you can mix and match as you choose, provided the functionality matches up.

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
- Z4 2008 Alcubierre implementation (I have the eigenvectors in 'numerical relativity codegen/run.lua') vs 2008 Yano (in 'numerical relativity codegen/verify 2008 yano')
- Implement eigen-stuff code in SRHD so that PLM can work 
- PLM for Euler-Burgers
- Any kind of slope extrapolation for finite-difference solvers? (BSSN and Z4c)
- PPM
- Better divergence removal.  multigrid GPU?
- Finish GLM-(ideal)MHD ... especially implement the source term as a second step with an exp(dt)  (which I'm not doing at the moment)
- Rename mhd to ideal-mhd
- How about a GLM method for Maxwell equations, so I can remove the divergence-free constraint
- Calculate and implement source terms for curvilinear coordinate systems (working on a tool to do this)
- Get two-fluid-separate EMHD working, so I can update the glm-maxwell with an implicit and update the ion and electron with an explicit solver
- Two-fluid plasma EMHD combined has numerical issues: the ion wavespeeds surpass the speed of light
- Currently seeing errors when two solvers run simultaneously ... which makes EM+HD difficult
- Add HLLD solver
- Finish implementing Navier-Stokes, compressible & incompressible
- BSSN connections based on difference with grid coordinate system
- Test out the GR+HD solvers
- Add source terms to GRHD -- or at least plugins for 'gr-hd-separate' to fill in from the NR solver
- Finish the GR+EM solver
- Add EM+GR+HD by winning
- Figure out what to do with self-gravitational potential energy in the Euler simulation.  offsetting it positive makes gravitational instability stable.  offsetting it too positive explodes -- even a forward euler integrator (why).  offsetting it negative causes more instability.
- Change vector field from immediate mode to buffered geometry, and gometry shaders if they're available
- Coroutines to iterative solvers?  so they don't stall the app execution?
- RHD W error in >1 dimension
- GR flat space simulations make an initial wave.  but shouldn't flat space be stable?
- BSSN, ADM, Z4C, etc still need momentum constraints?
- Not all of the Maxwell initial conditions are working with the non-GLM maxwell, especially with scalar vs complex
- Right now I'm implementing weno similar to the 1996 Jiang Shu paper: 1) calc interface flux, 2) weno reconstruct interface flux 3) finite volume integrate.  There are alternative uses of WENO (using PLM or flux-limiters to find the initial interface flux values, using WENO on left and right states and then applying a flux (HLL, Roe, etc), etc).  Try these out?
- Also, right now I am performing the PLM slope extrapolation in a separate kernel.  Not for WENO.  Combining kernels seems to run faster.  Maybe I should just inline the PLM stuff into the same kernel?
- Cram as much into a single kernel as possible.  More inlining, more redundant calculations.  This seems to run faster than separate solvers and separate buffers.  Best would be a way to automate the inlining.

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
	* 1996 Jiang, Shu, "Efficient Implementation of Weighted ENO Schemes"
* OpenCL reduce:
	* http://developer.amd.com/resources/documentation-articles/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
* WENO:
	* 1996 Jiang, Shu, "Efficient Implementation of Weighted ENO Schemes"
	* 1998 Shu "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws"
	* "A hybrid approach for the regularized long wave-Burgers equation"
	* 2016 Rathan, Raju "An improved Non-linear Weights for Seventh-Order WENO Scheme"
	* https://github.com/jzrake/Mara for weno5 examples
	* https://github.com/wme7/WENO7-Z/blob/master/WENO7ZresAdv1d.m for weno7 examples
	* https://github.com/python-hydro/hydro_examples/blob/master/compressible/weno_coefficients.py for weno7-13 coefficients

