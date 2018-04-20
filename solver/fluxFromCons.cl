<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(<?=eqn.cons_t?> U, real3 x) {
	return eigen_fluxTransform_<?=side?>(eigen_forCell_<?=side?>(U, x), U, x);
}
<? end ?>
