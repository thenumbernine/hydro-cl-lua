//#include "math.types.h"

/*
unit conversion variables.
my current convention is this:
- initial conditions should provide variables to eqn pre-converted to unitless if they so desire.
- (therefore) state variables should be considered unitless.
- solver variables are not yet converted to unitless.  every time they are referenced, factor out the units.
- right now I'm dividing by units to convert to unitless, and multiplying by units to convert out.  This might be the inverse of what I should be doing.  The plus side is when inputting units to the conversion, you don't have to invert so often. 'meter = 6.3716' gives you 1 distance unit = 6.3716 meters.
*/

#define unit_m					solver->meter
#define unit_s					solver->second
#define unit_kg					solver->kilogram
#define unit_C					solver->coulomb
#define unit_K					solver->kelvin

#define unit_m2					(unit_m * unit_m)
#define unit_m3					(unit_m * unit_m * unit_m)
#define unit_s2					(unit_s * unit_s)
#define unit_C2					(unit_C * unit_C)
#define unit_m_per_s			(unit_m / unit_s)
#define unit_m2_per_s2			(unit_m2 / unit_s2)
#define unit_m3_per_kg_s2		(unit_m3 / (unit_kg * unit_s2))
#define unit_kg_per_m3			(unit_kg / unit_m3)
#define unit_kg_per_m2_s		(unit_kg / (unit_m2 * unit_s))
#define unit_kg_per_m_s2		(unit_kg / (unit_m * unit_s2))
#define unit_C_per_kg			(unit_C / unit_kg)
#define unit_C_per_m2			(unit_C / unit_m2)
#define unit_kg_per_C_s			(unit_kg / (unit_C * unit_s))
#define unit_kg_m_per_C2		(unit_kg * unit_m / unit_C2)
#define unit_C2_s_per_kg_m3		((unit_C2 * unit_s) / (unit_kg * unit_m3))
#define unit_C2_s2_per_kg_m3	((unit_C2 * unit_s2) / (unit_kg * unit_m3))


#define real_conj(x)		(x)
#define real_from_real(x)	(x)
#define real_from_cplx(x)	((x).re)
#define real_zero			0
#define real_neg(x)			(-(x))
#define real_inv(x)			(1./(x))
#define real_add(a,b)		((a) + (b))
#define real_sub(a,b)		((a) - (b))
#define real_mul(a,b)		((a) * (b))
#define real_real_mul(a,b)	((a) * (b))
#define real_div(a,b)		((a) / (b))
#define real_lenSq(x)		((x) * (x))
#define real_abs			fabs

#define real_sqrt			sqrt

#define _cplx(a,b) 			(cplx){.s={a,b}}
#define cplx_from_real(x)	_cplx(x,0)
#define cplx_from_cplx(x)	(x)
#define cplx_zero 			cplx_from_real(0)


cplx cplx_conj(cplx a);
cplx cplx_neg(cplx a);
real cplx_lenSq(cplx a);
real cplx_abs(cplx a);
real cplx_arg(cplx a);
cplx cplx_add(cplx a, cplx b);
cplx cplx_sub(cplx a, cplx b);
cplx cplx_mul(cplx a, cplx b);
cplx cplx_real_mul(cplx a, real b);
cplx cplx_inv(cplx a);
cplx cplx_div(cplx a, cplx b);
cplx cplx_exp(cplx a);
cplx cplx_log(cplx a);
cplx cplx_pow(cplx a, cplx b);
cplx cplx_sqrt(cplx a);
real4 quatUnitConj(real4 q);
real4 quatMul(real4 q, real4 r);

<?
function makevec3(vec, scalar)
	local add = scalar..'_add'
	local mul = scalar..'_mul'
?>
#define <?=scalar?>_add3(a,b,c)		(<?=add?>(<?=add?>(a,b),c))
#define <?=scalar?>_mul3(a,b,c)		(<?=mul?>(<?=mul?>(a,b),c))
#define _<?=vec?>(a,b,c) 			(<?=vec?>){.s={a,b,c}}
#define <?=vec?>_zero				_<?=vec?>(<?=scalar?>_zero,<?=scalar?>_zero,<?=scalar?>_zero)

<?=scalar?> <?=vec?>_dot(<?=vec?> a, <?=vec?> b);
<?=vec?> <?=vec?>_<?=scalar?>_mul(<?=vec?> a, <?=scalar?> s);
<?=vec?> <?=scalar?>_<?=vec?>_mul(<?=scalar?> a, <?=vec?> b);
<? if scalar ~= 'real' then ?>
<?=vec?> <?=vec?>_real_mul(<?=vec?> a, real b);
<?=vec?> real_<?=vec?>_mul(real a, <?=vec?> b);
<? end ?>
<?=vec?> <?=vec?>_add(<?=vec?> a, <?=vec?> b);
<?=vec?> <?=vec?>_sub(<?=vec?> a, <?=vec?> b);
<?=vec?> <?=vec?>_cross(<?=vec?> a, <?=vec?> b);
<?=vec?> <?=vec?>_neg(<?=vec?> a);

<?
end
makevec3('real3', 'real')
makevec3('cplx3', 'cplx')
?>

#define real3_from_real3(x)	x

real real3_lenSq(real3 a);
real real3_len(real3 a);
real real3_lenSq(real3 a);
real real3_len(real3 a);
real3 real3_swap0(real3 v);
real3 real3_swap1(real3 v);
real3 real3_swap2(real3 v);
real3 real3_rotFrom0(real3 v);
real3 real3_rotFrom1(real3 v);
real3 real3_rotFrom2(real3 v);
real3 real3_rotTo0(real3 v);
real3 real3_rotTo1(real3 v);
real3 real3_rotTo2(real3 v);
real3 real3_rotateFrom(real3 v, real3 n);
real3 real3_rotateTo(real3 v, real3 n);

#define real3_from_cplx3		cplx3_re

cplx3 cplx3_from_real3(real3 re);
cplx3 cplx3_from_real3_real3(real3 re, real3 im);
real3 cplx3_re(cplx3 v);
real3 cplx3_im(cplx3 v);
real cplx3_lenSq(cplx3 v);
real cplx3_len(cplx3 v);

#define _sym3(a,b,c,d,e,f) (sym3){.s={a,b,c,d,e,f}}
#define sym3_zero	_sym3(0,0,0,0,0,0)
#define sym3_ident	_sym3(1,0,0,1,0,1)

sym3 real3_outer(real3 a);
real sym3_det(sym3 m);
sym3 sym3_inv(sym3 m, real det);
real3 sym3_real3_mul(sym3 m, real3 v);
sym3 sym3_add(sym3 a, sym3 b);
sym3 sym3_sub(sym3 a, sym3 b);
sym3 sym3_real_mul(sym3 a, real s);
real sym3_dot(sym3 a, sym3 b);
real3 sym3_x(sym3 m);
real3 sym3_y(sym3 m);
real3 sym3_z(sym3 m);
real sym3_trace(sym3 m);
sym3 sym3_swap0(sym3 m);
sym3 sym3_swap1(sym3 m);
sym3 sym3_swap2(sym3 m);
real real3_weightedDot(real3 a, real3 b, sym3 m);
real real3_weightedLenSq(real3 a, sym3 m);
real real3_weightedLen(real3 a, sym3 m);

#define real3x3_zero (real3x3){.v={real3_zero, real3_zero, real3_zero}}

real3x3 real3_real3_outer(real3 a, real3 b);
real sym3_real3x3_dot(sym3 a, real3x3 b);
real3x3 sym3_sym3_mul(sym3 a, sym3 b);
real3x3 real3x3_sym3_mul(real3x3 a, sym3 b);
real3x3 sym3_real3x3_mul(sym3 a, real3x3 b);
sym3 real3x3_sym3_to_sym3_mul(real3x3 a, sym3 b);
sym3 sym3_real3x3_to_sym3_mul(sym3 a, real3x3 b);
sym3 sym3_from_real3x3(real3x3 a);
real3x3 real3x3_addT(real3x3 a, real3x3 b);
real3x3 real3x3_real3x3_mul(real3x3 a, real3x3 b);
real3 real3x3_real3_mul(real3x3 a, real3 b);
real real3x3_trace(real3x3 m);
real3x3 real3x3_from_real(real x);
real real3x3_det(real3x3 m);
real3x3 real3x3_inv(real3x3 m);

#define _3sym3_zero (_3sym3){.s={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}

_3sym3 _3sym3_add(_3sym3 a, _3sym3 b);
_3sym3 _3sym3_sub(_3sym3 a, _3sym3 b);
_3sym3 _3sym3_real_mul(_3sym3 a, real b);
_3sym3 sym3_3sym3_mul(sym3 a, _3sym3 b);
real3 _3sym3_sym3_dot23(_3sym3 a, sym3 b);
real3 sym3_3sym3_dot12(sym3 a, _3sym3 b);
sym3 real3_3sym3_dot1(real3 a, _3sym3 b);
real3 _3sym3_tr12(_3sym3 a);
real3x3 real3_3sym3_dot2(real3 a, _3sym3 b);

real3x3x3 _3sym3_sym3_mul(_3sym3 a, sym3 b);
real3 sym3_real3x3x3_dot23(sym3 a, real3x3x3 b);
real3x3 _3sym3_real3x3x3_dot12_23(_3sym3 a, real3x3x3 b);
sym3 _3sym3_real3x3x3_dot13_to_sym3(_3sym3 a, real3x3x3 b);

sym3sym3 sym3sym3_add(sym3sym3 a, sym3sym3 b);

real3 normalForSide0();
real3 normalForSide1();
real3 normalForSide2();
float crand();
