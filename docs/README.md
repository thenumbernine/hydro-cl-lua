I've got a lot of equations in the [`hydro/eqn/`](../hydro/eqn/) folder.<br>
I'll try to sort out where the deriving work comes from.<br>
In a perfect world, I'd move all their calculation symmath notebooks into here and just automate the code generation of every equation in one fell swoop.<br>
<br>

Wave Equation:<br>
- [`wave.lua`](../hydro/eqn/wave.lua) - Wave Equation, Finite-Volume.<br>
- [`wave_metric.lua`](../hydro/eqn/wave_metric.lua) - Wave Equation with background metric tensor support (distinct of the grid coordinate chart?).<br>
- - work and notes:<br>
- - [`lua/symmath/tests/output/wave equation in spacetime.html`](https://thenumbernine.github.io/lua/symmath/tests/output/wave%20equation%20in%20spacetime.html) - symmath script to rearrange curved-space wave equation terms into a PDE.<br>
- - [`lua/symmath/tests/output/wave equation in spacetime - flux eigenvectors.html`](https://thenumbernine.github.io/lua/symmath/tests/output/wave%20equation%20in%20spacetime%20%2d%20flux%20eigenvectors.html) - has flux derivative, eigensystem wrt a metric, but not in terms of an arbitrary-normal basis.<br>
- - [`math/wave equation hyperbolic form.html`](https://thenumbernine.github.io/math/wave%20equation%20hyperbolic%20form.html) - by hand, succeeds "wave equation in spacetime".  has flux eigensystem in terms of an arbitrary coordinate system and arbitrary normal basis.<br>
<br>

Wave Equation - Finite-Difference:<br>
- [`wave-fd.lua`](../hydro/eqn/wave-fd.lua) - Wave Equation, Finite-Difference.<br>
- - work and notes:<br>
- - [`math/wave equation in curved spacetime - finite difference.html`](https://thenumbernine.github.io/math/wave%20equation%20in%20curved%20spacetime%20%2d%20finite%20difference.html) - the filenames relate, but maybe it's completely different.<br>
<br>

Quantum Physics:<br>
- [`nls.lua`](../hydro/eqn/nls.lua) - Non-Linear Schrodinger Equation, Finite-Difference (I think?)<br>
- [`schrodinger-fd.lua`](../hydro/eqn/schrodinger-fd.lua) - Schrodinger Equation, Finite-Difference.<br>
<br>

Fluid Dynamics - Shallow Water Equations:<br>
- [`shallow-water.lua`](../hydro/eqn/shallow-water.lua)<br>
- - work and notes:<br>
- - [`math/CFD/Shallow Water.html`](https://thenumbernine.github.io/math/CFD/Shallow%20Water.html) - by hand, has flux derivative, acoustic matrix, arbitrary normal basis.<br>
- - [`lua/symmath/tests/output/Shallow Water equations - flux eigenvectors.html`](https://thenumbernine.github.io/lua/symmath/tests/output/Shallow%20Water%20equations%20%2d%20flux%20eigenvectors.html) - symmath script, has flux derivative, acoustic matrix, arbitrary normal eigensystem, and code-generation.<br>
<br>

Fluid Dynamics - Euler Fluid Equations:<br>
- [`euler.lua`](../hydro/eqn/euler.lua) - Euler fluid equation, Finite-Volume.<br>
- [`euler_prim.lua`](../hydro/eqn/euler_prim.lua) - An attempt to recast the Euler equations flux-vector-split using their primitive vars, at the suggestion of a professor, which is currently very unstable and going nowhere.<br>
- - work and notes:<br>
- - [`math/MHD/0 - Euler Fluid Equations.html`](https://thenumbernine.github.io/math/MHD/0%20%2d%20Euler%20Fluid%20Equations.html) - derivation, and proly should go in the CFD folder...<br>
- - [`lua/symmath/tests/output/Euler fluid equations - flux eigenvectors.html`](https://thenumbernine.github.io/lua/symmath/tests/output/Euler%20fluid%20equations%20%2d%20flux%20eigenvectors.html) - has code-generation in it.<br>
- - [`Euler fluid equations - flux eigenvectors.symmath`](https://thenumbernine.github.io/symmath/index.html?open=/symmath/tests/Euler%20fluid%20equations%20-%20flux%20eigenvectors.symmath) - matches above, but without code-generation, and as a worksheet.<br>
- - [`math/CFD/Euler Fluid Equations - Curved Geometry - Contravariant.html`](https://thenumbernine.github.io/math/CFD/Euler%20Fluid%20Equations%20%2d%20Curved%20Geometry%20%2d%20Contravariant.html) - matches above, but by hand and not fully ported to arbitrary-normal-basis.</a><br>
- - [`lua/symmath/tests/output/Euler fluid equations - primitive form.html`](https://thenumbernine.github.io/lua/symmath/tests/output/Euler%20fluid%20equations%20%2d%20primitive%20form.html)<br>
- - [`math/CFD/Euler Fluid Equations - Entropy Function.html`](https://thenumbernine.github.io/math/CFD/Euler%20Fluid%20Equations%20%2d%20Entropy%20Function.html)<br>
<br>

Fluid Dynamics - Navier-Stokes Equations:<br>
- [`navstokes-incomp.lua`](../hydro/eqn/navstokes-incomp.lua) - Navier-Stokes, incompressible (via divergence-free projection?), maybe with viscous term as a separate solver?<br>
- [`navstokes-wilcox.lua`](../hydro/eqn/navstokes-wilcox.lua) - an attempt at Navier-Stokes-Wilcox which I haven't got working yet.<br>
<br>

Special-Relativistic Euler Fluid Equations:<br>
- [`srhd.lua`](../hydro/eqn/srhd.lua)<br>
- - work and notes:<br>
- - [`math/CFD/SRHD.html`](https://thenumbernine.github.io/math/CFD/SRHD.html) - by hand notes, and the most complete, gets as far as acoustic matrix calcuations.<br>
- - [`lua/symmath/tests/output/SRHD_1D.html`](https://thenumbernine.github.io/lua/symmath/tests/output/SRHD_1D.html) - symmath script, has half of the previous, gets as far as the flux partial wrt the primitives.<br>
- - [`lua/symmath/tests/output/SRHD.html`](https://thenumbernine.github.io/lua/symmath/tests/output/SRHD.html) - symmath script, the very start of me doing a 3D version of the above.<br>
<br>

Maxwell Equations:<br>
- [`maxwell.lua`](../hydro/eqn/maxwell.lua) - Maxwell using divergence-free projection.<br>
- - work and notes:<br>
- - [`lua/symmath/tests/output/Maxwell equations - flux eigenvectors.html`](https://thenumbernine.github.io/lua/symmath/tests/output/Maxwell%20equations%20%2d%20flux%20eigenvectors.html)<br>
- - [`math/Electromagnetism/Maxwell equations in hyperbolic form.html`](https://thenumbernine.github.io/math/Electromagnetism/Maxwell%20equations%20in%20hyperbolic%20form.html) - has flux eigensystem in arbitrary normal basis, and includes GLM calculations as well.<br>
- - [`math/Electromagnetism/Maxwell Equations - Eigenmode Analysis.html`](https://thenumbernine.github.io/math/Electromagnetism/Maxwell%20Equations%20%2d%20Eigenmode%20Analysis.html)<br>
- <br>
- [`glm-maxwell.lua`](../hydro/eqn/glm-maxwell.lua) - Maxwell using generalized-lagrangian-multiplier (G.L.M.).<br>
- - work and notes:<br>
- - [`lua/symmath/tests/output/GLM-Maxwell equations - flux eigenvectors.html`](https://thenumbernine.github.io/lua/symmath/tests/output/GLM%2dMaxwell%20equations%20%2d%20flux%20eigenvectors.html) - has flux eigensystem in arbitrary normal basis, and has code generation.<br>
- <br>
- [`maxwell_A_wave.lua`](../hydro/eqn/maxwell_A_wave.lua) - Maxwell in its wave equation form, using divergence-free projection, probably not finished.<br>
- [`glm-maxwell-A-wave.lua`](../hydro/eqn/glm-maxwell-A-wave.lua) Maxwell using G.L.M., probably not finished.<br>
<br>

Magnetohydrodynamics:<br>
- [`mhd.lua`](../hydro/eqn/mhd.lua) - Ideal MHD with divergence-free projection<br>
- [`glm-mhd.lua`](../hydro/eqn/glm-mhd.lua) - Ideal MHD with G.L.M.<br>
- - work and notes:<br>
- - [`math/MHD/1 - Ideal Microscopic Magnetohydrodynamics.html`](https://thenumbernine.github.io/math/MHD/1%20%2d%20Ideal%20Microscopic%20Magnetohydrodynamics.html)<br>
- - [`math/MHD/2 - Ideal MHD - Curved Geometry - Contravariant.html`](https://thenumbernine.github.io/math/MHD/2%20%2d%20Ideal%20MHD%20%2d%20Curved%20Geometry%20%2d%20Contravariant.html)<br>
- - [`lua/symmath/tests/output/MHD - flux eigenvectors.html`](https://thenumbernine.github.io/lua/symmath/tests/output/MHD%20%2d%20flux%20eigenvectors.html)<br>
- - [`lua/symmath/tests/output/MHD inverse.html`](https://thenumbernine.github.io/lua/symmath/tests/output/MHD%20inverse.html)<br>
- - [`lua/symmath/tests/output/MHD symmetrization.html`](https://thenumbernine.github.io/lua/symmath/tests/output/MHD%20symmetrization.html)<br>
- - [`math/MHD/MHD A-wave.html`](https://thenumbernine.github.io/math/MHD/MHD%20A%2dwave.html)<br>
<br>

Maxwell (E+B fields) + Hydrodynamics:<br>
- [`twofluid-emhd.lua`](../hydro/eqn/twofluid-emhd.lua)<br>
- [`twofluid-emhd-composite.lua`](../hydro/eqn/twofluid-emhd-composite.lua)<br>
<br>

Gravito-Electro-Magnetism, aka de-Donder Gauge Linearized General Relativity:<br>
- [`lingr.lua`](../hydro/eqn/lingr.lua)<br>
<br>

Euler Fluid Equations + Linearized GR:<br>
- [`euler-lingr.lua`](../hydro/eqn/euler-lingr.lua)<br>
<br>

Maxwell + Euler Fluid Equations + Linearized GR:<br>
- [`twofluid-emhd-lingr.lua`](../hydro/eqn/twofluid-emhd-lingr.lua)<br>
<br>

General Relativity:<Br>
- [`adm1d_v1.lua`](../hydro/eqn/adm1d_v1.lua) - ADM, 1D only<br>
- [`adm1d_v2b.lua`](../hydro/eqn/adm1d_v2b.lua) - ADM, 1D only<br>
- [`adm1d_v2.lua`](../hydro/eqn/adm1d_v2.lua) - ADM, 1D only<br>
- [`adm3d.lua`](../hydro/eqn/adm3d.lua) - ADM, 3D, Finite-Volume, violates Hamiltonian constraint for the sake of stability.<br>
- [`bssnok-fd-num.lua`](../hydro/eqn/bssnok-fd-num.lua) - BSSNOK, Finite-Difference, code generation done ahead of time.<br>
- - [`lua/symmath/tests/output/BSSN - generate - cartesian.html`](https://thenumbernine.github.io/lua/symmath/tests/output/BSSN%20%2d%20generate%20%2d%20cartesian.html)
- - [`lua/symmath/tests/output/BSSN - generate - spherical.html`](https://thenumbernine.github.io/lua/symmath/tests/output/BSSN%20%2d%20generate%20%2d%20spherical.html)
- [`bssnok-fd-sym.lua`](../hydro/eqn/bssnok-fd-sym.lua) - BSSNOK, Finite-Difference, attempt at runtime code generation + simplification (to remove 1/r's), not working last I checked<br>
- [`z4c-fd.lua`](../hydro/eqn/z4c-fd.lua) - Z4C, Finite-Difference, probably neglected and broken.<br>
- [`z4.lua`](../hydro/eqn/z4.lua) - Z4 based on eigensystem I worked out by hand.  Works.  But no shift block is included in the matrix.<br>
- - [`lua/symmath/tests/output/Z4.html`](https://thenumbernine.github.io/lua/symmath/tests/output/Z4.html)
- coordinate system code generation:<br>
- - [`lua/symmath/tests/output/Schwarzschild - isotropic.html`](https://thenumbernine.github.io/lua/symmath/tests/output/Schwarzschild%20%2d%20isotropic.html)
- <br>

- [`z4_2008yano.lua`](../hydro/eqn/z4_2008yano.lua) - Z4 based on 2008 Yano eigensystem, which I found math errors in the paper's eigensystem, and coincidentally I didn't get it to run stable.<br>
- - work and notes:<br>
- - [`numrel-codegen/flux_matrix_output/verify_2008_yano.html`](https://thenumbernine.github.io/numrel%2dcodegen/flux_matrix_output/verify_2008_yano.html)<br>
<br>

General Relativity + Maxwell:<br>
- [`gr-maxwell.lua`](../hydro/eqn/gr-maxwell.lua) - not finished<br>
<br>

General Relativity + Euler Fluid Equations:<br>
- [`grhd.lua`](../hydro/eqn/grhd.lua) - not finished<br>
<br>

General Relativity + Ideal MHD:<br>
- [`grmhd.lua`](../hydro/eqn/grmhd.lua) - not finished<br>

