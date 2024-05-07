//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=coordMap?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);

	real solid = 0;
	real rho = 0;
	real F[<?=#solver.offsets?>];
	for (int i = 0; i < <?=#solver.offsets?>; ++i) {
		F[i] = 0;
	}

<?=initCode()?>

	U->rho = rho;
	U->solid = solid;
	U->v = real3_zero;

<? for i=0,#solver.offsets-1 do
?>	U->F<?=i?> = F[<?=i?>];
<? end
?>
}

//// MODULE_NAME: <?=advect?>

kernel void <?=advect?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UNextBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const UNext = UNextBuf + index;

	UNext->rho = U->rho;
	UNext->v = U->v;
	UNext->solid = U->solid;

<?
	 for i,ofs in ipairs(solver.offsets) do 
		local c = ofs.c
		local ofsindex = i-1
?>	UNext->F<?=ofsindex?> = U[<?=-c.x?> + solver->gridSize.x * (<?=-c.y?> + solver->gridSize.y * <?=-c.z?>)].F<?=ofsindex?>;
<? 	end
?>
}

//// MODULE_NAME: <?=calcPrims?>

// maybe I can combine this with applyCollision ...
kernel void <?=calcPrims?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const U = UBuf + index;
	
	real rho = 0;
	real3 m = real3_zero;

<?
	for i,ofs in ipairs(solver.offsets) do
		local c = ofs.c
		local ofsindex = i-1
	?>{
		real const Fnbhd = U->F<?=ofsindex?>;
		rho += Fnbhd;
		m.x += Fnbhd * <?=c.x?>;
		m.y += Fnbhd * <?=c.y?>;
		m.z += Fnbhd * <?=c.z?>;
	}<?
	end
?>

	U->rho = rho;
	if (U->solid) {
		U->v = real3_zero;
	} else {
		U->v = real3_real_mul(m, 1./rho);
	}
}

//// MODULE_NAME: <?=applyCollision?>

kernel void <?=applyCollision?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	realparam const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const U = UBuf + index;

	real const invdt = 1. / dt;

	if (!U->solid) {
		real3 const v = U->v;
		real const vSq = real3_lenSq(v);
<?
	for i,ofs in ipairs(solver.offsets) do
		local c = ofs.c
		local ofsindex = i-1
?>		{
			real const velDotOfs = v.x * <?=c.x?> + v.y * <?=c.y?> + v.z * <?=c.z?>;
			real const Feq = U->rho * <?=ofs.weight?> * (1. + 3. * velDotOfs + 4.5 * velDotOfs * velDotOfs - 1.5 * vSq);
			U->F<?=ofsindex?> *= 1. - invdt;
			U->F<?=ofsindex?> += invdt *  Feq;
		}
<?	end
?>
	} else {	
		real tmpF[<?=#solver.offsets?>];
<?
	for i,ofs in ipairs(solver.offsets) do
		local c = ofs.c
		local ofsindex = i-1
?>		tmpF[<?=ofsindex?>] = U->F<?=ofs.oppositeOffsetIndex?>;
<?	end
	for ofsindex=0,#solver.offsets-1 do
?>		U->F<?=ofsindex?> = tmpF[<?=ofsindex?>];
<?	end
?>
	}
}
