//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: <?=solver_macros?> math <?=eigen_forInterface?> <?=eqn_waveCode_depends?> <?=fluxFromCons?>

//HLLC based on
//http://math.lanl.gov/~shenli/publications/hllc_mhd.pdf
//2012 Toro, "The HLLC Riemann Solver" presentation: http://marian.fsik.cvut.cz/~bodnar/PragueSum_2012/Toro_2-HLLC-RiemannSolver.pdf
<?
if not require "hydro.eqn.euler".is(solver.eqn) then 
	error("euler-hllc only works with euler eqn") 
end
?>


#define <?=calcFluxForInterface?>(\
	/*global <?=cons_t?> * const */flux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */xInt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> WL;\
	<?=primFromCons?>(&WL, solver, UL, (cellL)->pos);\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, (cellR)->pos);\
\
<? if true then ?> /* use interface waves? */\
	/*\
	get min/max lambdas of UL, UR, and interface U (based on Roe averaging)\
	TODO this in a more computationally efficient way\
	\
	This still needs the 'xInt' arg because some things (eig parameters derived from the averaged parameters, glm-mhd's calcRoeValues, navstokes-wilcox, twofluid-emhd & -lingr, etc) use it, and why recalculate it?\
	But mind you not everyone uses the xInt, so it shouldn't always be calculated right?\
	This is where fully-analytical implementations are handy (like BSSN is becoming).\
	*/\
	<?=eigen_t?> eigInt;\
	<?=eigen_forInterface?>(&eigInt, solver, UL, UR, cellL, cellR, xInt, n);\
	\
	real lambdaIntMin, lambdaIntMax;\
	{\
		<?=eqn:eigenWaveCodePrefix("n", "&eigInt", "xInt"):gsub("\n", "\\\n\t\t")?>\
		lambdaIntMin = <?=eqn:eigenMinWaveCode("n", "&eigInt", "xInt")?>;\
		lambdaIntMax = <?=eqn:eigenMaxWaveCode("n", "&eigInt", "xInt")?>;\
	}\
	\
<? 	if solver.flux.hllCalcWaveMethod == "Davis direct" then ?>\
	real sL = lambdaIntMin;\
	real sR = lambdaIntMax;\
<? 	end ?>\
<? 	if solver.flux.hllCalcWaveMethod == "Davis direct bounded" then ?>\
	real lambdaLMin;\
	{\
		<?=eqn:consWaveCodePrefix("n", "UL", "xInt"):gsub("\n", "\n\t\t"):gsub("\n", "\\\n\t\t")?>\
		lambdaLMin = <?=eqn:consMinWaveCode("n", "UL", "xInt")?>;\
	}\
\
	real lambdaRMax;\
	{\
		<?=eqn:consWaveCodePrefix("n", "UR", "xInt"):gsub("\n", "\\\n\t\t")?>\
		lambdaRMax = <?=eqn:consMaxWaveCode("n", "UR", "xInt")?>;\
	}\
	\
	real const sL = min(lambdaLMin, lambdaIntMin);\
	real const sR = max(lambdaRMax, lambdaIntMax);\
<? 	end ?>\
<? else ?> /* don't use interface waves.  looks no different than above, soo ... why even use interface state waves? */\
	real sLMin, sLMax, sRMin, sRMax;\
	{\
		<?=eqn:consWaveCodePrefix(0, "UL", "xL"):gsub("\n", "\\\n\t\t")?>\
		sLMin = <?=eqn:consMinWaveCode(0, "UL", "xL")?>;\
		sLMax = <?=eqn:consMaxWaveCode(0, "UL", "xL")?>;\
	}\
	{\
		<?=eqn:consWaveCodePrefix(0, "UR", "xR"):gsub("\n", "\\\n\t\t")?>\
		sRMin = <?=eqn:consMinWaveCode(0, "UR", "xR")?>;\
		sRMax = <?=eqn:consMaxWaveCode(0, "UR", "xR")?>;\
	}\
	real const sL = min(sLMin, sRMin);\
	real const sR = max(sLMax, sRMax);\
<? end ?>\
\
	/* notice that, for Euler, consWaveCodePrefix calculates v dot n1 ... while here we're using v dot nj */\
	real3 const vnL = normal_vecDotNs(n, WL.v);\
	real3 const vnR = normal_vecDotNs(n, WR.v);\
\
<? assert(solver.flux.hllcMethod) ?>	\
	real const sStar = (WR.rho * vnR.x * (sR - vnR.x) - WL.rho * vnL.x * (sL - vnL.x) + WL.P - WR.P) \
		/ (WR.rho * (sR - vnR.x) - WL.rho * (sL - vnL.x));\
\
	if (0 <= sL) {\
		<?=fluxFromCons?>(flux, solver, UL, cellL, n);\
<? if solver.flux.hllcMethod == 0 then ?>\
	\
	} else if (sL <= 0. && 0. <= sStar) {\
		<?=cons_t?> FL;\
		<?=fluxFromCons?>(&FL, solver, UL, cellL, n);\
		<?=cons_t?> ULStar;\
		ULStar.rho = (UL)->rho * (sL - vnL.x) / (sL - sStar);\
		\
		real3 vStar = normal_vecFromNs(n, _real3(sStar, vnL.y, vnL.z));\
		ULStar.m.x = ULStar.rho * vStar.x;\
		ULStar.m.y = ULStar.rho * vStar.y;\
		ULStar.m.z = ULStar.rho * vStar.z;\
		\
		ULStar.ETotal = ULStar.rho * (\
			(UL)->ETotal / (UL)->rho\
			+ (sStar - vnL.x) \
				* (sStar + WL.P / ((UL)->rho * (sL - vnL.x)))\
		);\
		for (int i = 0; i < numStates; ++i) {\
			flux->ptr[i] = FL.ptr[i] + sL * (ULStar.ptr[i] - (UL)->ptr[i]);\
		}\
	} else if (sStar <= 0. && 0. <= sR) {\
		<?=cons_t?> FR;\
		<?=fluxFromCons?>(&FR, solver, UR, cellR, n);\
		<?=cons_t?> URStar;\
		URStar.rho = (UR)->rho * (sR - vnR.x) / (sR - sStar);\
		\
		real3 vStar = normal_vecFromNs(n, _real3(sStar, vnR.y, vnR.z));\
		URStar.m.x = URStar.rho * vStar.x;\
		URStar.m.y = URStar.rho * vStar.y;\
		URStar.m.z = URStar.rho * vStar.z;\
		\
		URStar.ETotal = URStar.rho * (\
			(UR)->ETotal / (UR)->rho\
			+ (sStar - vnR.x) \
				* (sStar + WR.P / ((UR)->rho * (sR - vnR.x))));\
		for (int i = 0; i < numStates; ++i) {\
			flux->ptr[i] = FR.ptr[i] + sR * (URStar.ptr[i] - (UR)->ptr[i]);\
		}\
\
<? elseif solver.flux.hllcMethod == 1 then ?>\
\
	} else if (sL <= 0. && 0. <= sStar) {\
		<?=cons_t?> FL;\
		<?=fluxFromCons?>(&FL, solver, UL, cellL, n);\
		flux->rho = (sStar * (sL * (UL)->rho - FL.rho)) / (sL - sStar);\
	\
		real3 ULmn = normal_vecDotNs(n, (UL)->m);\
		real3 FLmn = normal_vecDotNs(n, FL.m);\
		flux->m = normal_vecFromNs(n, _real3(\
			(sStar * (sL * ULmn.x - FLmn.x) + sL * (WL.P + WL.rho * (sL - vnL.x) * (sStar - vnL.x))) / (sL - sStar),\
			(sStar * (sL * ULmn.y - FLmn.y)) / (sL - sStar),\
			(sStar * (sL * ULmn.z - FLmn.z)) / (sL - sStar)\
		));\
		\
		flux->ETotal = (sStar * (sL * (UL)->ETotal - FL.ETotal) + sL * (WL.P + WL.rho * (sL - vnL.x) * (sStar - vnL.x)) * sStar) / (sL - sStar);\
	} else if (sStar <= 0. && 0. <= sR) {\
		<?=cons_t?> FR;\
		<?=fluxFromCons?>(&FR, solver, UR, cellR, n);\
		flux->rho = (sStar * (sR * (UR)->rho - FR.rho)) / (sR - sStar);\
		\
		real3 URmn = normal_vecDotNs(n, (UR)->m);\
		real3 FRmn = normal_vecDotNs(n, FR.m);\
		flux->m = normal_vecFromNs(n, _real3(\
			(sStar * (sR * URmn.x - FRmn.x) + sR * (WR.P + WR.rho * (sR - vnR.x) * (sStar - vnR.x))) / (sR - sStar),\
			(sStar * (sR * URmn.y - FRmn.y)) / (sR - sStar),\
			(sStar * (sR * URmn.z - FRmn.z)) / (sR - sStar)\
		));\
\
		flux->ETotal = (sStar * (sR * (UR)->ETotal - FR.ETotal) + sR * (WR.P + WR.rho * (sR - vnR.x) * (sStar - vnR.x)) * sStar) / (sR - sStar);\
\
<? elseif solver.flux.hllcMethod == 2 then ?>\
\
	} else if (sL <= 0. && 0. <= sStar) {\
		<?=cons_t?> FL;\
		<?=fluxFromCons?>(&FL, solver, UL, cellL, n);\
		real PLR = .5 * (\
			WL.P\
			+ WR.P\
			+ WL.rho * (sL - vnL.x) * (sStar - vnL.x)\
			+ WR.rho * (sR - vnR.x) * (sStar - vnR.x)\
		);\
		flux->rho = (sL * (UL)->rho - FL.rho) * sStar / (sL - sStar);\
		\
		real3 ULmn = normal_vecDotNs(n, (UL)->m);\
		real3 FLmn = normal_vecDotNs(n, FL.m);\
		flux->m = normal_vecFromNs(n, _real3(\
			((sL * ULmn.x - FLmn.x) * sStar + sL * PLR) / (sL - sStar),\
			sStar * (sL * ULmn.y - FLmn.y) / (sL - sStar),\
			sStar * (sL * ULmn.z - FLmn.z) / (sL - sStar)\
		));\
\
		flux->ETotal = (sStar * (sL * (UL)->ETotal - FL.ETotal) + sL * PLR * sStar) / (sL - sStar);\
	} else if (sStar <= 0. && 0. <= sR) {\
		<?=cons_t?> FR;\
		<?=fluxFromCons?>(&FR, solver, UR, cellR, n);\
		real PLR = .5 * (WL.P + WR.P + WL.rho * (sL - vnL.x) * (sStar - vnL.x) + WR.rho * (sR - vnR.x) * (sStar - vnR.x));\
		flux->rho = sStar * (sR * (UR)->rho - FR.rho) / (sR - sStar);\
		\
		real3 URmn = normal_vecDotNs(n, (UR)->m);\
		real3 FRmn = normal_vecDotNs(n, FR.m);\
		flux->m = normal_vecFromNs(n, _real3(\
			/*TODO why were these assigns here:?*/\
			/*flux->m.x = */(sStar * (sR * URmn.x - FRmn.x) + sR * PLR) / (sR - sStar),\
			/*flux->m.y = */sStar * (sR * URmn.y - FRmn.y) / (sR - sStar),\
			/*flux->m.z = */sStar * (sR * URmn.z - FRmn.z) / (sR - sStar)\
		));\
		\
		flux->ETotal = (sStar * (sR * (UR)->ETotal - FR.ETotal) + sR * PLR * sStar) / (sR - sStar);\
\
<? end	--solver.flux.hllcMethod ?>\
	\
	} else if (sR <= 0) {\
		<?=fluxFromCons?>(flux, solver, UR, cellR, n);\
<? if true then ?>	/*why is this here? for when sStar is not between sL and sR*/\
	} else if (sL <= 0 && 0 <= sR) {\
		<?=cons_t?> FL;\
		<?=fluxFromCons?>(&FL, solver, UL, cellL, n);\
		<?=cons_t?> FR;\
		<?=fluxFromCons?>(&FR, solver, UR, cellR, n);\
		for (int j = 0; j < numIntStates; ++j) {\
			flux->ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * ((UR)->ptr[j] - (UL)->ptr[j])) / (sR - sL);\
		}\
<? end ?>\
	}\
	/*flux->m = real3_rotTo<?=side?>(flux->m);*/\
}
