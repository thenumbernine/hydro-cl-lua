//// MODULE_NAME: weno_l/r
//// MODULE_DEPENDS: cell_x sqr fluxFromCons normal_t eigen_forInterface eigen_left/rightTransform eqn.waveCode

<? 
local clnumber = require "cl.obj.number"

local stencilSize = solver.stencilSize
local coeffs = solver.coeffs

local coeff = coeffs[stencilSize]
local a = coeff.a
local d = coeff.d
local betaCoeffs = coeff.betaCoeffs

for side=0,solver.dim-1 do
	for _,l_or_r in ipairs{"l", "r"} do
		local ak0 = l_or_r == "l" and stencilSize or 1
		local akd = l_or_r == "l" and -1 or 1
		local al0 = l_or_r == "l" and stencilSize or 1
		local ald = l_or_r == "l" and -1 or 1
		local d0 = l_or_r == "l" and stencilSize or 1
		local dd = l_or_r == "l" and -1 or 1

-- how should I generalize this?
-- base pointer
-- number of cells (stencilSize)
-- number of coefficients at each cell ... or what C type it is, and how many variables to cycle across in the ptr field
-- ... or we can just always perform this interpolation in cons/char space
?>
<?=eqn.waves_t?> weno_<?=l_or_r?>_<?=side?>(
	const <?=eqn.waves_t?>* v
) {
	<?=eqn.waves_t?> result;
	for (int k = 0; k < numWaves; ++k) {
<? 		for j=0,stencilSize-1 do 
?>		real beta<?=j?> = 0.<?
			for m=0,stencilSize-1 do
				for n=0,m do ?> 
			+ <?=clnumber(betaCoeffs[j+1][m+1][n+1])?> * v[<?=j+m?>].ptr[k] * v[<?=j+n?>].ptr[k]<?
				end
			end?>;
<?		end
	
		if solver.wenoMethod == "1996 Jiang Shu" then 		-- WENO-JS
			--local epsilon = clnumber(1e-14)
			local epsilon = clnumber(1e-6)
?>
<? 			for i=0,stencilSize-1 do
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> / sqr(<?=epsilon?> + beta<?=i?>);
<? 			end 

		elseif solver.wenoMethod == "2008 Borges" then 	-- WENO-Z
			local epsilon = clnumber(1e-6)
			
			-- TODO different coeffs of betas?
			if false then 	--if stencilSize == 4 then -- for 2016 Rathan, it suggests these:
?>		real tau = fabs(beta0 + 3. * beta1 - 3. * beta2 - beta3);
<?			else	-- for weno7, 2018 Zeytinoglu suggests this:
?>		real tau = fabs(beta0 - beta<?=stencilSize-1?>);
<?			end
 			for i=0,stencilSize-1 do 
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> * (1. + (tau / (beta<?=i?> + <?=epsilon?>)));
<? 			end

		elseif solver.wenoMethod == "2010 Shen Zha" then -- WENO-BS?
			local epsilon = clnumber(1e-10)
		
			-- TODO find these for other order WENO's
			local shen_zha_A = clnumber(50)	-- 0-100
?>		
		real minB = beta0;
<?			for j=1,stencilSize-1 do
?>		minB = min(minB, beta<?=j?>);
<?			end
?>		
		real maxB = beta0;
<?			for j=1,stencilSize-1 do
?>		maxB = max(maxB, beta<?=j?>);
<?			end
?>		
		real R0 = minB / (maxB + <?=epsilon?>);

<? 			for i=0,stencilSize-1 do 
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> / sqr(<?=epsilon?> + beta<?=i?> + R0 * <?=shen_zha_A?> * minB);
<?	 		end 
		else
			error("unknown wenoMethod "..tostring(solver.wenoMethod))
		end
?>
<?		for i=0,stencilSize-1 do 
?>		real vs<?=i?> = 0.<?
			for j=0,stencilSize-1 do ?>
			+ <?=clnumber(a[ak0+akd*i][al0+ald*j])?> * v[<?=i+j?>].ptr[k]<?
			end ?>;
<? 		end 
?>		
		real alphasum = 0.<? for i=0,stencilSize-1 do ?> + alpha<?=i?><? end ?>;
		
		result.ptr[k] = (0.<?
		for i=0,stencilSize-1 do ?> 
			+ alpha<?=i?> * vs<?=i?><?
		end ?>
		) / alphasum;
	}
	return result;
}
<? 	end 
end
?>

//// MODULE_NAME: calcCellFlux
//// MODULE_DEPENDS: cell_x fluxFromCons normal_t

kernel void calcCellFlux(
	constant <?=solver.solver_t?>* solver,
	const global <?=solver.coord.cell_t?>* cellBuf,
	global <?=eqn.cons_t?>* fluxCellBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		real3 xInt = cell_x(i);
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		normal_t n = normal_forSide<?=side?>(xInt);
		global <?=eqn.cons_t?>* F = fluxCellBuf + <?=side?> + dim * index;
		*F = fluxFromCons(solver, *U, xInt, n);
	}<? end ?>
}

//// MODULE_NAME: calcFlux
//// MODULE_DEPENDS: cell_x fluxFromCons normal_t eigen_forInterface eigen_left/rightTransform eqn.waveCode weno_l/r

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf/*<?=solver.getULRArg?>*/,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?
for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - solver->stepsize.s<?=side?>;
		int indexR = index;

#if 0
		<?=solver:getULRCode():gsub('\n', '\n\t\t')?>
#else
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;
		const global <?=eqn.cons_t?>* UR = UBuf + indexR;
#endif

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		normal_t n = normal_forSide<?=side?>(xInt);

		int fluxIndexInt = side + dim * index;
		global <?=eqn.cons_t?>* flux = fluxBuf + fluxIndexInt;

<?
	
	if solver.fluxMethod == 'Lax-Friedrichs' then

?>
		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);
		real maxAbsLambda = 0.;
		for (int j = 0; j < <?=2*stencilSize?>; ++j) {
			const global <?=eqn.cons_t?>* Uj = U + (j - <?=stencilSize?>) * solver->stepsize.s<?=side?>;
			
			<?=eqn:consWaveCodePrefix('n', '*Uj', 'xInt'):gsub('\n', '\n\t\t')?>
			
			real lambdaMin = <?=eqn:consMinWaveCode('n', '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMin));
			
			real lambdaMax = <?=eqn:consMaxWaveCode('n', '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMax));
		}

<? 	
		for j=0,2*stencilSize-1 do
?>		
		const global <?=eqn.cons_t?>* U<?=j?> = U + <?=j - stencilSize?> * solver->stepsize.s<?=side?>;
		//<?=eqn.cons_t?> F<?=j?> = fluxCellBuf[<?=side?> + dim * (index + <?=j - stencilSize?> * solver->stepsize.s<?=side?>)];
		<?=eqn.cons_t?> F<?=j?> = fluxFromCons(solver, *U<?=j?>, xInt, n);
<?
			if j < 2*stencilSize-1 then
?>		<?=eqn.cons_t?> fp<?=j?>;
<?				for k=0,eqn.numStates-1 do 
?>		fp<?=j?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] + maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
<?				end
			end
			
			if j > 0 then
?>		<?=eqn.cons_t?> fm<?=j-1?>;
<?				for k=0,eqn.numStates-1 do 
?>		fm<?=j-1?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] - maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
<? 				end
			end
		end	
?>	
		<?=eqn.waves_t?> afp[<?=2*stencilSize-1?>];
<? 		for j=0,2*stencilSize-2 do
?>		afp[<?=j?>] = eigen_leftTransform(solver, eig, fp<?=j?>, xInt, n);
<?		end
?>
		<?=eqn.waves_t?> afm[<?=2*stencilSize-1?>];
<?		for j=0,2*stencilSize-2 do
?>		afm[<?=j?>] = eigen_leftTransform(solver, eig, fm<?=j?>, xInt, n);
<? 		end 
?>
		<?=eqn.waves_t?> wafp = weno_r_<?=side?>(afp);
		<?=eqn.waves_t?> wafm = weno_l_<?=side?>(afm);
		
		<?=eqn.waves_t?> waf;
		for (int j = 0; j < numWaves; ++j) {
			waf.ptr[j] = wafp.ptr[j] + wafm.ptr[j];
		}
		
		*flux = eigen_rightTransform(solver, eig, waf, xInt, n);
<?
	
	elseif solver.fluxMethod == 'Marquina' then

?>
//// MODULE_DEPENDS: eigen_forCell
		<?=eqn.eigen_t?> eigL = eigen_forCell(solver, *UL, xInt, n);
		<?=eqn.eigen_t?> eigR = eigen_forCell(solver, *UR, xInt, n);

		real lambdaL[<?=eqn.numWaves?>];
		real lambdaR[<?=eqn.numWaves?>];
<? 		for _,lr in ipairs{'L', 'R'} do
?>		{
		<?=eqn:eigenWaveCodePrefix('n', 'eig'..lr, 'xInt'):gsub('\n', '\n\t\t')?>
<? 		
			for k=0,eqn.numWaves-1 do 
?>			lambda<?=lr?>[<?=k?>] = <?=eqn:eigenWaveCode('n', 'U'..lr, 'xInt', k)?>;
<?			end
?>		}
<?		end
?>

<? 		for j=0,2*stencilSize-1 do
?>		const global <?=eqn.cons_t?>* U<?=j?> = U + (<?=j - stencilSize?>) * solver->stepsize.s<?=side?>;
		<?=eqn.cons_t?> F<?=j?> = fluxFromCons(solver, *U<?=j?>, xInt, n);
		<?=eqn.waves_t?> al<?=j?> = eigen_leftTransform(solver, eigL, *U<?=j?>, xInt, n);
		<?=eqn.waves_t?> ar<?=j?> = eigen_leftTransform(solver, eigR, *U<?=j?>, xInt, n);
		<?=eqn.waves_t?> afl<?=j?> = eigen_leftTransform(solver, eigL, F<?=j?>, xInt, n);
		<?=eqn.waves_t?> afr<?=j?> = eigen_leftTransform(solver, eigR, F<?=j?>, xInt, n);

<?		end
?>

		<?=eqn.waves_t?> afp[<?=2*stencilSize-1?>], afm[<?=2*stencilSize-1?>];
<? 
		for k=0,eqn.numWaves-1 do 
?>		if (lambdaL[<?=k?>] > 0.0 && lambdaR[<?=k?>] > 0.0) {
<?			for j=0,2*stencilSize-1 do
?>			afp[<?=j?>].ptr[<?=k?>] = afl<?=j?>.ptr[<?=k?>];
			afm[<?=j?>].ptr[<?=k?>] = 0;
<?			end
?>		} else if (lambdaL[<?=k?>] < 0.0 && lambdaR[<?=k?>] < 0.0) {
<?			for j=0,2*stencilSize-1 do
?>			afp[<?=j?>].ptr[<?=k?>] = 0;
			afm[<?=j?>].ptr[<?=k?>] = afr<?=j?>.ptr[<?=k?>];
<?			end
?>		} else {
			real absLambdaL = fabs(lambdaL[<?=k?>]);
			real absLambdaR = fabs(lambdaR[<?=k?>]);
			real a, absa;
			if (absLambdaL > absLambdaR) {
				a = lambdaL[<?=k?>];
				absa = absLambdaL;
			} else {
				a = lambdaR[<?=k?>];
				absa = absLambdaR;
			}
<?			for j=0,2*stencilSize-1 do
?>			afp[<?=j?>].ptr[<?=k?>] = .5 * (afl<?=j?>.ptr[<?=k?>] + absa * al<?=j?>.ptr[<?=k?>]);
			afm[<?=j?>].ptr[<?=k?>] = .5 * (afr<?=j?>.ptr[<?=k?>] - absa * ar<?=j?>.ptr[<?=k?>]);
<?			end
?>		}
<?		end
?>	

		<?=eqn.waves_t?> wafp = weno_r_<?=side?>(afp);
		<?=eqn.waves_t?> wafm = weno_l_<?=side?>(afm);
		
		<?=eqn.cons_t?> fluxP = eigen_rightTransform(solver, eigL, wafp, xInt, n);
		<?=eqn.cons_t?> fluxM = eigen_rightTransform(solver, eigR, wafm, xInt, n);
		
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] = fluxP.ptr[j] + fluxM.ptr[j];
		}
<?
	
	elseif solver.fluxMethod == 'Roe' then
		
		-- TODO make the following its own function, maybe in solver/roe.cl or solver/roe.lua
		-- so we can call a function instead of copy/paste
		-- but this means making it a function of the lua parameters, like flux limiter.
		-- and a flux limiter means +1 to the numGhost

?>		<?=eqn.waves_t?> afp[<?=2*stencilSize-1?>], afm[<?=2*stencilSize-1?>];
<?		for j=0,2*stencilSize-1 do
?>
		{
			const global <?=eqn.cons_t?>* UL<?=j?> = U + <?=j - stencilSize?> * solver->stepsize.s<?=side?>;
			const global <?=eqn.cons_t?>* UR<?=j?> = U + <?=j+1 - stencilSize?> * solver->stepsize.s<?=side?>;
			<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);

			<?=eqn:eigenWaveCodePrefix('n', 'eig', 'xInt')?>

			<?=eqn.cons_t?> UAvg;
			for (int k = 0; k < numIntStates; ++k) {
				UAvg.ptr[k] = .5 * (UL->ptr[k] + UR->ptr[k]);
			}
			afp[<?=j?>] = eigen_leftTransform(solver, eig, UAvg, xInt, n);
			afm[<?=j?>] = afp[<?=j?>];

			<?=eqn.cons_t?> deltaU;
			for (int k = 0; k < numStates; ++k) {
				deltaU.ptr[k] = UR->ptr[k] - UL->ptr[k];
			}
			
			<?=eqn.waves_t?> deltaUEig = eigen_leftTransform(solver, eig, deltaU, xInt, n);

			<? for k=0,eqn.numWaves-1 do ?>{
				const int k = <?=k?>;
				real lambda = <?=eqn:eigenWaveCode('n', 'eig', 'xInt', k)?>;
				real lambdaPlus = max(lambda, 0.);
				real lambdaMinus = min(lambda, 0.);

				afp[<?=j?>].ptr[k] *= lambdaPlus;
				afm[<?=j?>].ptr[k] *= lambdaMinus;
				real sgnLambda = lambda >= 0 ? 1 : -1;
				afp[<?=j?>].ptr[k] -= .5 * lambdaPlus * deltaUEig.ptr[k] * sgnLambda;
				afm[<?=j?>].ptr[k] -= .5 * lambdaMinus * deltaUEig.ptr[k] * sgnLambda;
			}<? end ?>
		}
<?
		end 
?>
		<?=eqn.waves_t?> waf;
		<?=eqn.waves_t?> wafp = weno_r_<?=side?>(afp);
		<?=eqn.waves_t?> wafm = weno_l_<?=side?>(afm);
		for (int j = 0; j < numWaves; ++j) {
			waf.ptr[j] = wafp.ptr[j] + wafm.ptr[j];
		}
	
		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);
		*flux = eigen_rightTransform(solver, eig, waf, xInt, n);
<?
	
	else
		error("unknown fluxMethod "..tostring(solver.fluxMethod))
	end
?>	}<? 
end ?>
}
