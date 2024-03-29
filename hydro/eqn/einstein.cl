<? if false then ?>
// TODO all these From/To names should be distinct per function if they will ever be used with composite eqn
//  specifically if they will ever be used with more than one einstein solver
//  but then again, if you already have one solver for your spacetime, why would you need two?
<? end ?>

//// MODULE_NAME: <?=rescaleFromCoord_rescaleToCoord?>
//// MODULE_DEPENDS: <?=coord_dx_i?>
// rescaling from/to diagonalization from the grid metric
// used especially by the bssn solvers
// but anything that wants to use their initial conditions will also need this.
#if 1

//rescaling, used for bssn finite-difference, but I am tempted to try it with other coordinate systems with singularities
//TODO for the initial conditions do this symbolically instead of numerically

//apply this to lower indexes to convert from coordinate metric to better metric
//apply this to upper indexes to convert from better metric to coordinate metric
real3 real3_rescaleFromCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x / coord_dx0(x),
		.y = v.y / coord_dx1(x),
		.z = v.z / coord_dx2(x),
	};
}
#define real3_rescaleToCoord_U real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_L(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y * coord_dx1(x),
		.z = v.z * coord_dx2(x),
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_L

real3s3 real3s3_rescaleFromCoord_ll(real3s3 a, real3 x) {
	return (real3s3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define real3s3_rescaleToCoord_UU real3s3_rescaleFromCoord_ll

real3s3 real3s3_rescaleToCoord_LL(real3s3 a, real3 x) {
	return (real3s3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define real3s3_rescaleFromCoord_uu real3s3_rescaleToCoord_LL

#else	//debugging -- turning it off

#define real3_rescaleFromCoord_l(a,x) a
#define real3_rescaleToCoord_U(a,x) a
#define real3_rescaleToCoord_L(a,x) a
#define real3_rescaleFromCoord_u(a,x) a
#define real3s3_rescaleFromCoord_ll(a,x) a
#define real3s3_rescaleToCoord_UU(a,x) a
#define real3s3_rescaleToCoord_LL(a,x) a
#define real3s3_rescaleFromCoord_uu(a,x) a

#endif

//// MODULE_NAME: <?=cplx3_rescaleFromCoord_cplx3_rescaleToCoord?>
//// MODULE_DEPENDS: cplx3 <?=rescaleFromCoord_rescaleToCoord?>	
// only used by bssnok-fd-num with scalar field
cplx3 cplx3_rescaleFromCoord_l(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		real3_rescaleFromCoord_l(cplx3_re(v), x),
		real3_rescaleFromCoord_l(cplx3_im(v), x));
}

//// MODULE_NAME: <?=real3x3s3_rescaleFromCoord_real3x3s3_rescaleToCoord?>
//// MODULE_DEPENDS: real3x3s3
// <?=rescaleFromCoord_rescaleToCoord?>	// I could use this for sub-member rescaling
// ... but I just did it manually
//// MODULE_DEPENDS: <?=coord_dx_i?>
#if 1

real3x3s3 real3x3s3_rescaleFromCoord_lll(real3x3s3 a, real3 x) {
	return (real3x3s3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3s3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define real3x3s3_rescaleToCoord_UUU real3x3s3_rescaleFromCoord_lll

real3x3s3 real3x3s3_rescaleToCoord_LLL(real3x3s3 a, real3 x) {
	return (real3x3s3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3s3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define real3x3s3_rescaleFromCoord_uuu real3x3s3_rescaleToCoord_LLL

real3x3s3 real3x3s3_rescaleFromCoord_ull(real3x3s3 a, real3 x) {
	return (real3x3s3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3s3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * coord_dx<?=i-1?>(x) / (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}

real3x3s3 real3x3s3_rescaleToCoord_ULL(real3x3s3 a, real3 x) {
	return (real3x3s3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3s3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)) / coord_dx<?=i-1?>(x),
<?	end
?>		},
<? end
?>	};
}

#else

#define real3x3s3_rescaleFromCoord_lll(a,x) a
#define real3x3s3_rescaleToCoord_UUU(a,x) a
#define real3x3s3_rescaleToCoord_LLL(a,x) a
#define real3x3s3_rescaleFromCoord_uuu(a,x) a

#endif

//// MODULE_NAME: <?=real3s3x3s3_rescaleFromCoord_real3s3x3s3_rescaleToCoord?>
//// MODULE_DEPENDS: real3s3x3s3
// <?=rescaleFromCoord_rescaleToCoord?>	// I could use this for sub-member rescaling
// ... but I just did it manually
//// MODULE_DEPENDS: <?=coord_dx_i?>
// only used by bssnok-fd-sym
#if 1

real3s3x3s3 real3s3x3s3_rescaleFromCoord_llll(real3s3x3s3 a, real3 x) {
	return (real3s3x3s3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (real3s3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define real3s3x3s3_rescaleToCoord_UUUU real3s3x3s3_rescaleFromCoord_llll

real3s3x3s3 real3s3x3s3_rescaleToCoord_LLLL(real3s3x3s3 a, real3 x) {
	return (real3s3x3s3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (real3s3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define real3s3x3s3_rescaleFromCoord_uuuu real3s3x3s3_rescaleToCoord_LLLL

#else

#define real3s3x3s3_rescaleFromCoord_llll(a,x) a
#define real3s3x3s3_rescaleToCoord_UUUU(a,x) a
#define real3s3x3s3_rescaleToCoord_LLLL(a,x) a
#define real3s3x3s3_rescaleFromCoord_uuuu (a,x) a

#endif
