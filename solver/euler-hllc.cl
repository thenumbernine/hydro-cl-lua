//HLLC based on
//http://math.lanl.gov/~shenli/publications/hllc_mhd.pdf
//2012 Toro, "The HLLC Riemann Solver" presentation: http://marian.fsik.cvut.cz/~bodnar/PragueSum_2012/Toro_2-HLLC-RiemannSolver.pdf

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	<?=solver.getULRArg?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - solver->stepsize.s<?=side?>;
		real3 xL = xR;
		xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		int indexInt = side + dim * index;

		<?=solver:getULRCode{suffix='_', indexL='indexL', indexR='indexR'}?>
		<?=eqn.cons_t?> UL = *UL_;
		<?=eqn.cons_t?> UR = *UR_;

		//align to x-axis
		UL.m = real3_rotFrom<?=side?>(UL.m);
		UR.m = real3_rotFrom<?=side?>(UR.m);
		
		<?=eqn.prim_t?> WL = primFromCons(solver, UL, xInt);
		<?=eqn.prim_t?> WR = primFromCons(solver, UR, xInt);

#if 1	//use interface waves?
		// get min/max lambdas of UL, UR, and interface U (based on Roe averaging)
		// TODO this in a more computationally efficient way
		<?=eqn.eigen_t?> eigInt = eigen_forInterface(solver, UL, UR, xInt, normalForSide0());
		
		real lambdaIntMin, lambdaIntMax;
		{
			<?=eqn:eigenWaveCodePrefix(0, 'eigInt', 'xInt')?>
			lambdaIntMin = <?=eqn:eigenMinWaveCode(0, 'eigInt', 'xInt')?>;
			lambdaIntMax = <?=eqn:eigenMaxWaveCode(0, 'eigInt', 'xInt')?>;
		}
		
<? if solver.calcWaveMethod == 'Davis direct' then ?>
		real sL = lambdaIntMin;
		real sR = lambdaIntMax;
<? end ?>
<? if solver.calcWaveMethod == 'Davis direct bounded' then ?>
		real lambdaLMin;
		{
			<?=eqn:consWaveCodePrefix(0, 'UL', 'xL')?>
			lambdaLMin = <?=eqn:consMinWaveCode(0, 'UL', 'xL')?>;
		}

		real lambdaRMax;
		{
			<?=eqn:consWaveCodePrefix(0, 'UR', 'xR')?>
			lambdaRMax = <?=eqn:consMaxWaveCode(0, 'UR', 'xR')?>;
		}
		
		real sL = min(lambdaLMin, lambdaIntMin);
		real sR = max(lambdaRMax, lambdaIntMax);
<? end ?>
#else	//don't use interface waves.  looks no different than above, soo ... why even use interface state waves?
		real sLMin, sLMax, sRMin, sRMax;
		{
			<?=eqn:consWaveCodePrefix(0, 'UL', 'xL')?>
			sLMin = <?=eqn:consMinWaveCode(0, 'UL', 'xL')?>;
			sLMax = <?=eqn:consMaxWaveCode(0, 'UL', 'xL')?>;
		}
		{
			<?=eqn:consWaveCodePrefix(0, 'UR', 'xR')?>
			sRMin = <?=eqn:consMinWaveCode(0, 'UR', 'xR')?>;
			sRMax = <?=eqn:consMaxWaveCode(0, 'UR', 'xR')?>;
		}
		real sL = min(sLMin,sRMin);
		real sR = max(sLMax,sRMax);
#endif

<? if solver.hllcMethod then ?>	
		real sStar = (WR.rho * WR.v.x * (sR - WR.v.x) - WL.rho * WL.v.x * (sL - WL.v.x) + WL.P - WR.P) 
						/ (WR.rho * (sR - WR.v.x) - WL.rho * (sL - WL.v.x));
<? end ?>

		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		if (0 <= sL) {
			<?=eqn.cons_t?> FL = fluxFromCons_0(solver, UL, xL);
			*flux = FL;
		
<? if solver.hllcMethod == 0 then ?>
		
		} else if (sL <= 0. && 0. <= sStar) {
			<?=eqn.cons_t?> FL = fluxFromCons_0(solver, UL, xL);
			<?=eqn.cons_t?> ULStar;
			ULStar.rho = UL.rho * (sL - WL.v.x) / (sL - sStar);
			ULStar.m.x = ULStar.rho * sStar;
			ULStar.m.y = ULStar.rho * WL.v.y;
			ULStar.m.z = ULStar.rho * WL.v.z;
			ULStar.ETotal = ULStar.rho * (
				UL.ETotal / UL.rho
				+ (sStar - WL.v.x) 
					* (sStar + WL.P / (UL.rho * (sL - WL.v.x)))
			);
			for (int i = 0; i < numStates; ++i) {
				flux->ptr[i] = FL.ptr[i] + sL * (ULStar.ptr[i] - UL.ptr[i]);
			}
		} else if (sStar <= 0. && 0. <= sR) {
			<?=eqn.cons_t?> FR = fluxFromCons_0(solver, UR, xR);
			<?=eqn.cons_t?> URStar;
			URStar.rho = UR.rho * (sR - WR.v.x) / (sR - sStar);
			URStar.m.x = URStar.rho * sStar;
			URStar.m.y = URStar.rho * WR.v.y;
			URStar.m.z = URStar.rho * WR.v.z;
			URStar.ETotal = URStar.rho * (
				UR.ETotal / UR.rho
				+ (sStar - WR.v.x) 
					* (sStar + WR.P / (UR.rho * (sR - WR.v.x))));
			for (int i = 0; i < numStates; ++i) {
				flux->ptr[i] = FR.ptr[i] + sR * (URStar.ptr[i] - UR.ptr[i]);
			}

<? elseif solver.hllcMethod == 1 then ?>
			
		} else if (sL <= 0. && 0. <= sStar) {
			<?=eqn.cons_t?> FL = fluxFromCons_0(solver, UL, xL);
			flux->rho = (sStar * (sL * UL.rho - FL.rho)) / (sL - sStar);
			flux->m.x = (sStar * (sL * UL.m.x - FL.m.x) + sL * (WL.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x))) / (sL - sStar);
			flux->m.y = (sStar * (sL * UL.m.y - FL.m.y)) / (sL - sStar);
			flux->m.z = (sStar * (sL * UL.m.z - FL.m.z)) / (sL - sStar);
			flux->ETotal = (sStar * (sL * UL.ETotal - FL.ETotal) + sL * (WL.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x)) * sStar) / (sL - sStar);
		} else if (sStar <= 0. && 0. <= sR) {
			<?=eqn.cons_t?> FR = fluxFromCons_0(solver, UR, xR);
			flux->rho = (sStar * (sR * UR.rho - FR.rho)) / (sR - sStar);
			flux->m.x = (sStar * (sR * UR.m.x - FR.m.x) + sR * (WR.P + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x))) / (sR - sStar);
			flux->m.y = (sStar * (sR * UR.m.y - FR.m.y)) / (sR - sStar);
			flux->m.z = (sStar * (sR * UR.m.z - FR.m.z)) / (sR - sStar);
			flux->ETotal = (sStar * (sR * UR.ETotal - FR.ETotal) + sR * (WR.P + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x)) * sStar) / (sR - sStar);

<? elseif solver.hllcMethod == 2 then ?>
	
		} else if (sL <= 0. && 0. <= sStar) {
			<?=eqn.cons_t?> FL = fluxFromCons_0(solver, UL, xL);
			real PLR = .5 * (
				WL.P
				+ WR.P
				+ WL.rho * (sL - WL.v.x) * (sStar - WL.v.x)
				+ WR.rho * (sR - WR.v.x) * (sStar - WR.v.x)
			);
			flux->rho = (sL * UL.rho - FL.rho) * sStar / (sL - sStar);
			flux->m.x = ((sL * UL.m.x - FL.m.x) * sStar + sL * PLR) / (sL - sStar);
			flux->m.y = sStar * (sL * UL.m.y - FL.m.y) / (sL - sStar);
			flux->m.z = sStar * (sL * UL.m.z - FL.m.z) / (sL - sStar);
			flux->ETotal = (sStar * (sL * UL.ETotal - FL.ETotal) + sL * PLR * sStar) / (sL - sStar);
		} else if (sStar <= 0. && 0. <= sR) {
			<?=eqn.cons_t?> FR = fluxFromCons_0(solver, UR, xR);
			real PLR = .5 * (WL.P + WR.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x) + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x));
			flux->rho = sStar * (sR * UR.rho - FR.rho) / (sR - sStar);
			flux->m.x = (sStar * (sR * UR.m.x - FR.m.x) + sR * PLR) / (sR - sStar);
			flux->m.y = sStar * (sR * UR.m.y - FR.m.y) / (sR - sStar);
			flux->m.z = sStar * (sR * UR.m.z - FR.m.z) / (sR - sStar);
			flux->ETotal = (sStar * (sR * UR.ETotal - FR.ETotal) + sR * PLR * sStar) / (sR - sStar);

<? end	--solver.hllcMethod ?>
		
		} else if (sR <= 0) {
			<?=eqn.cons_t?> FR = fluxFromCons_0(solver, UR, xR);
			*flux = FR;
#if 1	//why is this here? for when sStar is not between sL and sR
		} else if (sL <= 0 && 0 <= sR) {
			<?=eqn.cons_t?> FL = fluxFromCons_0(solver, UL, xL);
			<?=eqn.cons_t?> FR = fluxFromCons_0(solver, UR, xR);
			for (int j = 0; j < numIntStates; ++j) {
				flux->ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * (UR.ptr[j] - UL.ptr[j])) / (sR - sL);
			}
#endif		
		}
		flux->m = real3_rotTo<?=side?>(flux->m);
	
	}<? end ?>
}
