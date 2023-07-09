//// MODULE_NAME: <?=macros?>

#define maxofsvol (3*3*3)

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=coordMap?>
//// MODULE_DEPENDS: <?=macros?>

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
	real F[maxofsvol];
	for (int i = 0; i < solver->ofsvol; ++i) {
		F[i] = 0;
	}

<?=initCode()?>

	U->rho = rho;
	U->solid = solid;
	U->v = real3_zero;
	for (int i = 0; i < solver->ofsvol; ++i) {
		U->F[i] = F[i];
	}
}

//// MODULE_NAME: <?=advect?>
//// MODULE_DEPENDS: <?=macros?>

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

	int ofsindex = 0;
	int4 c = (int4)(0,0,0,0);		
	int4 ofs = (int4)(0,0,0,0);
	for (ofs.z = 0; ofs.z < solver->ofsmax.z; ++ofs.z) {
		c.z = ofs.z - (solver->ofsmax.z-1)/2;
		for (ofs.y = 0; ofs.y < solver->ofsmax.y; ++ofs.y) {
			c.y = ofs.y - (solver->ofsmax.y-1)/2;
			for (ofs.x = 0; ofs.x < solver->ofsmax.x; ++ofs.x) {
				c.x = ofs.x - (solver->ofsmax.x-1)/2;
				//int Fofs = -dot(c, solver->stepsize);
				int Fofs = 
					- c.x * solver->stepsize.x
					- c.y * solver->stepsize.y
					- c.z * solver->stepsize.z;
				UNext->F[ofsindex] = U[Fofs].F[ofsindex];
				++ofsindex;
			}
		}
	}
}

//// MODULE_NAME: <?=calcPrims?>
//// MODULE_DEPENDS: <?=macros?>

// maybe I can combine this with applyCollision ...
kernel void <?=calcPrims?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const U = UBuf + index;

	int ofsindex = 0;
	real rho = 0;
	real3 m = real3_zero;
	real3 c;
	int4 ofs = (int4)(0,0,0,0);
	for (ofs.z = 0; ofs.z < solver->ofsmax.z; ++ofs.z) {
		c.z = (real)(ofs.z - (solver->ofsmax.z-1)/2);
		for (ofs.y = 0; ofs.y < solver->ofsmax.y; ++ofs.y) {
			c.y = (real)(ofs.y - (solver->ofsmax.y-1)/2);
			for (ofs.x = 0; ofs.x < solver->ofsmax.x; ++ofs.x) {
				c.x = (real)(ofs.x - (solver->ofsmax.x-1)/2);
				real const Fnbhd = U->F[ofsindex];
				++ofsindex;
				rho += Fnbhd;
				m = real3_add(m, real3_real_mul(c, Fnbhd));
			}
		}
	}
	U->rho = rho;
	if (U->solid) {
		U->v = real3_zero;
	} else {
		U->v = real3_real_mul(m, 1./rho);
	}
}

//// MODULE_NAME: <?=applyCollision?>
//// MODULE_DEPENDS: <?=macros?>

kernel void <?=applyCollision?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	realparam const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const U = UBuf + index;

	real const invdt = 1. / dt;

	real tmpF[maxofsvol];

	if (U->solid) {
		int ofsindex = 0;
		int4 ofs = (int4)(0,0,0,0);
		int4 nofs = (int4)(0,0,0,0);
		for (ofs.z = 0; ofs.z < solver->ofsmax.z; ++ofs.z) {
			nofs.z = solver->ofsmax.z - 1 - ofs.z;
			for (ofs.y = 0; ofs.y < solver->ofsmax.y; ++ofs.y) {
				nofs.y = solver->ofsmax.y - 1 - ofs.y;
				for (ofs.x = 0; ofs.x < solver->ofsmax.x; ++ofs.x) {
					nofs.x = solver->ofsmax.x - 1 - ofs.x;
					tmpF[ofsindex] = U->F[nofs.x + 3 * (nofs.y + 3 * nofs.z)];
					++ofsindex;
				}
			}
		}
		for (int ofsindex = 0; ofsindex < solver->ofsvol; ++ofsindex) {
			U->F[ofsindex] = tmpF[ofsindex];
		}
	} else {
		int4 ofs = (int4)(0,0,0,0);
		real3 const v = U->v;
		real const vSq = real3_lenSq(v);
		real3 c;
		int ofsindex = 0;
		for (ofs.z = 0; ofs.z < solver->ofsmax.z; ++ofs.z) {
			c.z = (real)(ofs.z - (solver->ofsmax.z-1)/2);
			for (ofs.y = 0; ofs.y < solver->ofsmax.y; ++ofs.y) {
				c.y = (real)(ofs.y - (solver->ofsmax.y-1)/2);
				for (ofs.x = 0; ofs.x < solver->ofsmax.x; ++ofs.x) {
					c.x = (real)(ofs.x - (solver->ofsmax.x-1)/2);
					real const w = solver->weights[ofs.x + 3 * (ofs.y + 3 * ofs.z)];
					real const velDotOfs = real3_dot(c, v);
					real const Feq = U->rho * w * (1. + 3. * velDotOfs + 4.5 * velDotOfs * velDotOfs - 1.5 * vSq);
					U->F[ofsindex] *= 1. - invdt;
					U->F[ofsindex] += invdt *  Feq;
					++ofsindex;
				}
			}
		}
	}
}
