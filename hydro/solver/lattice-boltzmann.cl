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
	real nbhd[maxofsvol];
	for (int i = 0; i < solver->ofsvol; ++i) {
		nbhd[i] = 0;
	}

<?=initCode()?>

	U->rho = rho;
	U->solid = solid;
	U->v = real3_zero;
	for (int i = 0; i < solver->ofsvol; ++i) {
		U->nbhd[i] = nbhd[i];
	}
}

//// MODULE_NAME: <?=advect?>
//// MODULE_DEPENDS: <?=macros?>

kernel void <?=advect?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UNextBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost-1, solver->numGhost-2);
	
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const UNext = UNextBuf + index;

	UNext->rho = U->rho;
	UNext->v = U->v;
	UNext->solid = U->solid;

	int ofs = 0;
	int4 c = (int4)(0,0,0,0);		
	for (int ofz = 0; ofz < solver->ofsmax.z; ++ofz) {
		c.z = ofz - (solver->ofsmax.z-1)/2;
		for (int ofy = 0; ofy < solver->ofsmax.y; ++ofy) {
			c.y = ofy - (solver->ofsmax.y-1)/2;
			for (int ofx = 0; ofx < solver->ofsmax.x; ++ofx) {
				c.x = ofx - (solver->ofsmax.x-1)/2;
				//int indexofs = -dot(c, solver->stepsize);
				int indexofs = 
					-c.x * solver->stepsize.x
					-c.y * solver->stepsize.y
					-c.z * solver->stepsize.z;
				UNext->nbhd[ofs] = U[indexofs].nbhd[ofs];
				++ofs;
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
	<?=SETBOUNDS?>(solver->numGhost-1, solver->numGhost-2);
	
	global <?=cons_t?> * const U = UBuf + index;

	int ofs = 0;
	real rho = 0;
	real3 m = real3_zero;
	real3 c;
	for (int ofz = 0; ofz < solver->ofsmax.z; ++ofz) {
		c.z = (real)(ofz - 1);
		for (int ofy = 0; ofy < solver->ofsmax.y; ++ofy) {
			c.y = (real)(ofy - 1);
			for (int ofx = 0; ofx < solver->ofsmax.x; ++ofx) {
				c.x = (real)(ofx - 1);
				real rhoNbhd = U->nbhd[ofs];
				rho += rhoNbhd;
				m = real3_add(m, real3_real_mul(c, rhoNbhd));
				++ofs;
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
	global <?=cons_t?> * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost-1, solver->numGhost-2);
	
	global <?=cons_t?> * const U = UBuf + index;

	real const invtau = 1. / .6;	// 1./dt?

	if (U->solid) {
		real tmp[maxofsvol];
		int ofs = 0;
		for (int ofz = 0; ofz < solver->ofsmax.z; ++ofz) {
			int const nofz = solver->ofsmax.z-1-ofz;
			for (int ofy = 0; ofy < solver->ofsmax.y; ++ofy) {
				int const nofy = solver->ofsmax.y-1-ofy;
				for (int ofx = 0; ofx < solver->ofsmax.x; ++ofx) {
					int const nofx = solver->ofsmax.x-1-ofx;	
					tmp[ofs] = U->nbhd[nofx + 3 * (nofy + 3 * nofz)];
					++ofs;
				}
			}
		}
		for (int ofs = 0; ofs < solver->ofsvol; ++ofs) {
			U->nbhd[ofs] = tmp[ofs];
		}
	} else {
		real3 const v = U->v;
		real const vSq = real3_lenSq(v);
		real3 c;
		int ofs = 0;
		for (int ofz = 0; ofz < solver->ofsmax.z; ++ofz) {
			c.z = (real)(ofz - 1);
			for (int ofy = 0; ofy < solver->ofsmax.y; ++ofy) {
				c.y = (real)(ofy - 1);
				for (int ofx = 0; ofx < solver->ofsmax.x; ++ofx) {
					c.x = (real)(ofx - 1);
					real const w = solver->weights[ofx + 3 * (ofy + 3 * ofz)];
					real const velDotOfs = real3_dot(c, v);
					real const Feq = U->rho * w * (1. + 3. * velDotOfs + 4.5 * velDotOfs * velDotOfs - 1.5 * vSq);
					U->nbhd[ofs] *= 1. - invtau;
					U->nbhd[ofs] += invtau *  Feq;
					++ofs;
				}
			}
		}
	}
}
