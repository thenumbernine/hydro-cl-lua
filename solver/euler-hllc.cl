//HLLC based on
//http://math.lanl.gov/~shenli/publications/hllc_mhd.pdf
//http://marian.fsik.cvut.cz/~bodnar/PragueSum_2012/Toro_2-HLLC-RiemannSolver.pdf

kernel void calcFlux(
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - stepsize.s<?=side?>;
		real3 xL = xR;
		xL.s<?=side?> -= grid_dx<?=side?>;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		int indexInt = side + dim * index;
	
		<?=solver:getULRCode()?>
	
		<?=eqn.prim_t?> WL = primFromCons(*UL, xInt);
		<?=eqn.prim_t?> WR = primFromCons(*UR, xInt);

		// get min/max lambdas of UL, UR, and interface U (based on Roe averaging)
		// TODO this in a more computationally efficient way
		<?=eqn.eigen_t?> eigInt = eigen_forInterface(*UL, *UR, xInt, normalForSide<?=side?>());
	
		real lambdaIntMin, lambdaIntMax;
		{
			<?=eqn:eigenWaveCodePrefix(side, '&eigInt', 'xInt')?>
			lambdaIntMin = <?=eqn:eigenMinWaveCode(side, '&eigInt', 'xInt')?>;
			lambdaIntMax = <?=eqn:eigenMaxWaveCode(side, '&eigInt', 'xInt')?>;
		}

#if 0	//Davis direct
		real sL = lambdaIntMin;
		real sR = lambdaIntMax;
#endif
#if 1	//David direct bounded

#if 1	//inline ... no function calls, but they're replaced with eigen_forCell calls ...
		//goes a bit faster.  getting rid of eigen_forCell would make it go faster.
		real lambdaLMin;
		{
			<?=eqn.eigen_t?> eigL = eigen_forCell_<?=side?>(*UL, xL);
			<?=eqn:eigenWaveCodePrefix(side, '&eigL', 'xL')?>
			lambdaLMin = <?=eqn:eigenMinWaveCode(side, '&eigL', 'xL')?>;
		}

		real lambdaRMax;
		{
			<?=eqn.eigen_t?> eigR = eigen_forCell_<?=side?>(*UR, xR);
			<?=eqn:eigenWaveCodePrefix(side, '&eigR', 'xR')?>
			lambdaRMax = <?=eqn:eigenMaxWaveCode(side, '&eigR', 'xR')?>;
		}
		real sL = min(lambdaLMin, lambdaIntMin);
		real sR = max(lambdaRMax, lambdaIntMax);
#else	//call
		range_t lambdaL = calcCellMinMaxEigenvalues_<?=side?>(UL, xL);
		range_t lambdaR = calcCellMinMaxEigenvalues_<?=side?>(UR, xR);
		real sL = min(lambdaL.min, lambdaIntMin);
		real sR = max(lambdaR.max, lambdaIntMax);
#endif

#endif

<?
--local hllcMethod = 0
--local hllcMethod = 1
local hllcMethod = 2
?>

		//TODO align to x axis

<? if hllcMethod then ?>
		real sStar = (UR->rho * WR.v.x * (sR - WR.v.x) - UL->rho * WL.v.x * (sL - WL.v.x) + WL.P - WR.P) 
						/ (UR->rho * (sR - WR.v.x) - UL->rho * (sL - WL.v.x));
<? end ?>

		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		if (0 <= sL) {
			<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
			*flux = FL;
		
<? if hllcMethod == 0 then ?>
	
	} else if (sL <= 0. && 0. <= sStar) {
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
		<?=eqn.cons_t?> ULStar;
		ULStar.rho = UL->rho * (sL - WL.v.x) / (sL - sStar);
		ULStar.m.x = UL->rho * sStar;
		ULStar.m.y = ULStar.rho * WL.v.y;
		ULStar.m.z = ULStar.rho * WL.v.z;
		ULStar.ETotal = ULStar.rho * (UL->ETotal + (sStar - WL.v.x) * (sStar + WL.P / (UL->rho * (sL - WL.v.x))));
		for (int i = 0; i < numStates; ++i) {
			flux->ptr[i] = FL.ptr[i] + sL * (ULStar.ptr[i] - UL->ptr[i]);
		}
	} else if (sStar <= 0. && 0. <= sR) {
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
		<?=eqn.cons_t?> URStar;
		URStar.rho = UR->rho * (sR - WR.v.x) / (sR - sStar);
		URStar.m.x = URStar.rho * sStar;
		URStar.m.y = URStar.rho * WR.v.y;
		URStar.m.z = URStar.rho * WR.v.z;
		URStar.ETotal = UR->rho * (UR->ETotal + (sStar - WR.v.x) * (sStar + WR.P / (UR->rho * (sR - WR.v.x))));
		for (int i = 0; i < numStates; ++i) {
			flux->ptr[i] = FR.ptr[i] + sR * (URStar.ptr[i] - UR->ptr[i]);
		}

<? elseif hllcMethod == 1 then ?>
		
	} else if (sL <= 0. && 0. <= sStar) {
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
		flux->rho = (sStar * (sL * UL->rho - FL.rho)) / (sL - sStar);
		flux->m.x = (sStar * (sL * UL->m.x - FL.m.x) + sL * (WL.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x))) / (sL - sStar);
		flux->m.y = (sStar * (sL * UL->m.y - FL.m.y)) / (sL - sStar);
		flux->m.z = (sStar * (sL * UL->m.z - FL.m.z)) / (sL - sStar);
		flux->ETotal = (sStar * (sL * UL->ETotal - FL.ETotal) + sL * (WL.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x)) * sStar) / (sL - sStar);
	} else if (sStar <= 0. && 0. <= sR) {
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
		flux->rho = (sStar * (sR * UR->rho - FR.rho)) / (sR - sStar);
		flux->m.x = (sStar * (sR * UR->m.x - FR.m.x) + sR * (WR.P + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x))) / (sR - sStar);
		flux->m.y = (sStar * (sR * UR->m.y - FR.m.y)) / (sR - sStar);
		flux->m.z = (sStar * (sR * UR->m.z - FR.m.z)) / (sR - sStar);
		flux->ETotal = (sStar * (sR * UR->ETotal - FR.ETotal) + sR * (WR.P + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x)) * sStar) / (sR - sStar);

<? elseif hllcMethod == 2 then ?>
	
	} else if (sL <= 0. && 0. <= sStar) {
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
		real PLR = .5 * (WL.P + WR.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x) + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x));
		flux->rho = sStar * (sL * UL->rho - FL.rho) / (sL - sStar);
		flux->m.x = (sStar * (sL * UL->m.x - FL.m.x) + sL * PLR) / (sL - sStar);
		flux->m.y = sStar * (sL * UL->m.y - FL.m.y) / (sL - sStar);
		flux->m.z = sStar * (sL * UL->m.z - FL.m.z) / (sL - sStar);
		flux->ETotal = (sStar * (sL * UL->ETotal - FL.ETotal) + sL * PLR * sStar) / (sL - sStar);
	} else if (sStar <= 0. && 0. <= sR) {
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
		real PLR = .5 * (WL.P + WR.P + WL.rho * (sL - WL.v.x) * (sStar - WL.v.x) + WR.rho * (sR - WR.v.x) * (sStar - WR.v.x));
		flux->rho = sStar * (sR * UR->rho - FR.rho) / (sR - sStar);
		flux->m.x = (sStar * (sR * UR->m.x - FR.m.x) + sR * PLR) / (sR - sStar);
		flux->m.y = sStar * (sR * UR->m.y - FR.m.y) / (sR - sStar);
		flux->m.z = sStar * (sR * UR->m.z - FR.m.z) / (sR - sStar);
		flux->ETotal = (sStar * (sR * UR->ETotal - FR.ETotal) + sR * PLR * sStar) / (sR - sStar);

<? end	--hllcMethod ?>
		
		} else if (sR <= 0) {
			<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
			*flux = FR;
		} else if (sL <= 0 && 0 <= sR) {
			<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
			<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
			for (int j = 0; j < numIntStates; ++j) {
				flux->ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * (UR->ptr[j] - UL->ptr[j])) / (sR - sL);
			}
		}
	}<? end ?>
}
