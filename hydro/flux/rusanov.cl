//// MODULE_NAME: calcFluxForInterface
//// MODULE_DEPENDS: solver.macros math cons_parallelPropagate
//// MODULE_DEPENDS: eigen_forInterface eqn.waveCode fluxFromCons

<?=eqn.cons_t?> calcFluxForInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> pUL,
	<?=eqn.cons_t?> pUR,
	real3 xInt,
	normal_t n
) {
	<?=eqn.eigen_t?> eigInt = eigen_forInterface(solver, pUL, pUR, xInt, n);
	<?=eqn:eigenWaveCodePrefix('n', 'eigInt', 'xInt')?>

	real lambdaMax;
	<? for _,U in ipairs{"pUL", "pUR"} do 
		for _,minmax in ipairs{"Min", "Max"} do ?>{
		<?=eqn:consWaveCodePrefix("n", U, "xInt")?>
		real absLambda = fabs(<?=eqn["cons"..minmax.."WaveCode"](eqn, "n", U, "xInt")?>);
		lambdaMax = max(lambdaMax, absLambda);
	}<? end
	end ?>

	<?=eqn.cons_t?> FL = fluxFromCons(solver, pUL, xInt, n);
	<?=eqn.cons_t?> FR = fluxFromCons(solver, pUR, xInt, n);
	<?=eqn.cons_t?> F;
	for (int j = 0; j < numIntStates; ++j) {
		F.ptr[j] = .5 * (FL.ptr[j] + FR.ptr[j] - lambdaMax * (pUR.ptr[j] - pUL.ptr[j]));
	}
	return F;
}
