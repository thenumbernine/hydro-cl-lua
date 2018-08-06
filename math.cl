<?
local common = require 'common'()	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym
?>

#define real_conj(x)		(x)
#define real_from_real(x)	(x)
#define real_to_real(x)		(x)
#define real_zero			0
#define real_neg(x)			(-(x))
#define real_inv(x)			(1./(x))
#define real_add(a,b)		((a) + (b))
#define real_sub(a,b)		((a) - (b))
#define real_mul(a,b)		((a) * (b))
#define real_real_mul(a,b)	((a) * (b))
#define real_div(a,b)		((a) / (b))
#define real_lenSq(x)		((x) * (x))

#define real_sqrt			sqrt

#define _cplx(a,b) 			(cplx){.s={a,b}}
#define cplx_from_real(x)	_cplx(x,0)
#define cplx_to_real(x)		((x).re)
#define cplx_zero 			cplx_from_real(0)


#define _sym3(a,b,c,d,e,f) (sym3){.s={a,b,c,d,e,f}}

static inline cplx cplx_conj(cplx a) {
	return _cplx(a.re, -a.im);
}

static inline cplx cplx_neg(cplx a) { 
	return _cplx(-a.re, -a.im); 
}

static inline real cplx_lenSq(cplx a) {
	return a.re * a.re + a.im * a.im;
}

static inline real cplx_abs(cplx a) {
	return sqrt(cplx_lenSq(a));
}

static inline real cplx_arg(cplx a) {
	return atan2(a.im, a.re);
}

static inline cplx cplx_add(cplx a, cplx b) {
	return _cplx(
		a.re + b.re,
		a.im + b.im
	);
}

static inline cplx cplx_sub(cplx a, cplx b) {
	return _cplx(
		a.re - b.re,
		a.im - b.im
	);
}

static inline cplx cplx_mul(cplx a, cplx b) {
	return _cplx(
		a.re * b.re - a.im * b.im,
		a.re * b.im + a.im * b.re);
}

static inline cplx cplx_real_mul(cplx a, real b) {
	return _cplx(a.re * b, a.im * b);
}

static inline cplx cplx_inv(cplx a) {
	return cplx_real_mul(cplx_conj(a), 1. / cplx_lenSq(a));
}

static inline cplx cplx_div(cplx a, cplx b) {
	return cplx_mul(a, cplx_inv(b));
}

static inline cplx cplx_exp(cplx a) {
	real expre = exp(a.re);
	return _cplx(
		expre * cos(a.im),
		expre * sin(a.im)
	);
}

static inline cplx cplx_log(cplx a) {
	return _cplx(
		log(cplx_abs(a)),
		cplx_arg(a)
	);
}

static inline cplx cplx_pow(cplx a, cplx b) {
	return cplx_exp(cplx_mul(b, cplx_log(a)));
}

static inline cplx cplx_sqrt(cplx a) {
	return cplx_pow(a, cplx_from_real(.5));
}


//assumes q is unit
//returns the conjugate
real4 quatUnitConj(real4 q) {
	return (real4)(-q.x, -q.y, -q.z, q.x);
}

real4 quatMul(real4 q, real4 r) {
	real a = (q.w + q.x) * (r.w + r.x);
	real b = (q.z - q.y) * (r.y - r.z);
	real c = (q.x - q.w) * (r.y + r.z);
	real d = (q.y + q.z) * (r.x - r.w);
	real e = (q.x + q.z) * (r.x + r.y);
	real f = (q.x - q.z) * (r.x - r.y);
	real g = (q.w + q.y) * (r.w - r.z);
	real h = (q.w - q.y) * (r.w + r.z);

	return (real4)(
		 a - .5 * ( e + f + g + h), //x
		-c + .5 * ( e - f + g - h), //y
		-d + .5 * ( e - f - g + h), //z
		 b + .5 * (-e - f + g + h)); //w
}


<?
function makevec3(vec,scalar)
	local add = scalar..'_add'
	local sub = scalar..'_sub'
	local mul = scalar..'_mul'
	local real_mul = scalar..'_real_mul'
	local conj = scalar..'_conj'
	local add3 = scalar..'_add3'
?>

#define <?=scalar?>_add3(a,b,c)	(<?=add?>(<?=add?>(a,b),c))
#define _<?=vec?>(a,b,c) 		(<?=vec?>){.s={a,b,c}}
#define <?=vec?>_zero			_<?=vec?>(<?=scalar?>_zero,<?=scalar?>_zero,<?=scalar?>_zero)

static inline <?=scalar?> <?=vec?>_dot(<?=vec?> a, <?=vec?> b) {
	return <?=add3?>(
		<?=mul?>(a.x, <?=conj?>(b.x)), 
		<?=mul?>(a.y, <?=conj?>(b.y)),
		<?=mul?>(a.z, <?=conj?>(b.z))
	);
}

//maybe I should just use <?=vec?>_<?=scale?>_mul ?
static inline <?=vec?> <?=vec?>_<?=scalar?>_mul(<?=vec?> a, <?=scalar?> s) {
	return _<?=vec?>(
		<?=mul?>(a.x, s),
		<?=mul?>(a.y, s),
		<?=mul?>(a.z, s));
}

static inline <?=vec?> <?=scalar?>_<?=vec?>_mul(<?=scalar?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=mul?>(a, b.x),
		<?=mul?>(a, b.y),
		<?=mul?>(a, b.z));
}

<? if scalar ~= 'real' then ?>
static inline <?=vec?> <?=vec?>_real_mul(<?=vec?> a, real b) {
	return _<?=vec?>(
		<?=real_mul?>(a.x, b),
		<?=real_mul?>(a.y, b),
		<?=real_mul?>(a.z, b));
}

static inline <?=vec?> real_<?=vec?>_mul(real a, <?=vec?> b) {
	return _<?=vec?>(
		<?=real_mul?>(b.x, a),
		<?=real_mul?>(b.y, a),
		<?=real_mul?>(b.z, a));
}
<? end ?>

static inline <?=vec?> <?=vec?>_add(<?=vec?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=add?>(a.x, b.x), 
		<?=add?>(a.y, b.y), 
		<?=add?>(a.z, b.z));
}

static inline <?=vec?> <?=vec?>_sub(<?=vec?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=sub?>(a.x, b.x),
		<?=sub?>(a.y, b.y),
		<?=sub?>(a.z, b.z));
}

static inline <?=vec?> <?=vec?>_cross(<?=vec?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=sub?>(<?=mul?>(a.y, b.z), <?=mul?>(a.z, b.y)),
		<?=sub?>(<?=mul?>(a.z, b.x), <?=mul?>(a.x, b.z)),
		<?=sub?>(<?=mul?>(a.x, b.y), <?=mul?>(a.y, b.x)));
}

<?
end
makevec3('real3', 'real')
makevec3('cplx3', 'cplx')
?>

#define real3_to_real3(x)	x
#define real3_from_real3(x)	x

static inline real real3_lenSq(real3 a) {
	return real3_dot(a,a);
}

static inline real real3_len(real3 a) {
	return sqrt(real3_lenSq(a));
}

static inline sym3 real3_outer(real3 a, real3 b) {
	return _sym3(
		a.x * b.x, a.x * b.y, a.x * b.z,
		a.y * b.y, a.y * b.z, a.z * b.z);
}

//for swapping dimensions between x and 012
real3 real3_swap0(real3 v) { return v; }
real3 real3_swap1(real3 v) { return _real3(v.y, v.x, v.z); }
real3 real3_swap2(real3 v) { return _real3(v.z, v.y, v.x); }

//rotate from a particular side xyz to put x forward
real3 real3_rotFrom0(real3 v) { return v; }
real3 real3_rotFrom1(real3 v) { return _real3(v.y, -v.x, v.z); }
real3 real3_rotFrom2(real3 v) { return _real3(v.z, v.y, -v.x); }

//rotate to put x back to the side 
real3 real3_rotTo0(real3 v) { return v; }
real3 real3_rotTo1(real3 v) { return _real3(-v.y, v.x, v.z); }
real3 real3_rotTo2(real3 v) { return _real3(-v.z, v.y, v.x); }

//rotate 'n' to x-axis
//assumes 'n' is unit
real3 real3_rotateFrom(real3 v, real3 n) {
#if dim == 1
	return v;
#elif dim == 2
	return _real3(
		v.x * n.x + v.y * n.y,
		-v.x * n.y + v.y * n.x,
		v.z);
#elif dim == 3
	/*
	axis is n cross x-axis
	[ 1  0  0] x [nx ny nz] = [0, -nz, ny] / (ny^2 + nz^2)	
	angle = acos(n.x)
	cos angle = n.x
	sin angle = sqrt(1 - nx^2)
	cos (angle/2) = 

	cos^2 theta + sin^2 theta = 1
	sin^2 theta = 1 - cos^2
	cos^2 theta - sin^2 theta = cos(2 theta) <-> 
	cos(theta/2) = sqrt((1 + cos(theta))/2) <-> 
	2 sin theta cos theta = sin(2 theta)
	*/
	real cosTheta = n.x;
	real cosHalfTheta = sqrt(.5 * (1. + cosTheta));
	real sinHalfTheta = sqrt(1. - cosHalfTheta * cosHalfTheta);
	real n2 = sqrt(n.y * n.y + n.z * n.z);
	real ax = 0.;
	real ay = -n.z / n2;
	real az = n.y / n2;
	real4 q = (real4)(ax * sinHalfTheta, ay * sinHalfTheta, az * sinHalfTheta, cosHalfTheta);
	real4 qInv = quatUnitConj(q);
	real4 _v = (real4)(v.x, v.y, v.z, 0);
	real4 vres = quatMul(quatMul(q, _v), qInv);
	return _real3(vres.x, vres.y, vres.z);
#endif
}

//rotate x-axis to 'n'
real3 real3_rotateTo(real3 v, real3 n) {
#if dim == 1
	return v;
#elif dim == 2
	return _real3(
		v.x * n.x - v.y * n.y,
		v.x * n.y + v.y * n.x,
		v.z);
#elif dim == 3
	//same as above but with negative axis	
	real cosTheta = n.x;
	real cosHalfTheta = sqrt(.5 * (1. + cosTheta));
	real sinHalfTheta = sqrt(1. - cosHalfTheta * cosHalfTheta);
	real n2 = sqrt(n.y * n.y + n.z * n.z);
	real ax = 0.;
	real ay = n.z / n2;
	real az = -n.y / n2;
	real4 q = (real4)(ax * sinHalfTheta, ay * sinHalfTheta, az * sinHalfTheta, cosHalfTheta);
	real4 qInv = quatUnitConj(q);
	real4 _v = (real4)(v.x, v.y, v.z, 0);
	real4 vres = quatMul(quatMul(q, _v), qInv);
	return _real3(vres.x, vres.y, vres.z);
#endif
}

#define cplx3_to_real3		cplx3_re

cplx3 cplx3_from_real3(real3 re) {
	return (cplx3){
		.x = cplx_from_real(re.x),
		.y = cplx_from_real(re.y),
		.z = cplx_from_real(re.z),
	};
}

cplx3 cplx3_from_real3_real3(real3 re, real3 im) {
	return (cplx3){
		.x = {.re = re.x, .im = im.x},
		.y = {.re = re.y, .im = im.y},
		.z = {.re = re.z, .im = im.z},
	};
}

real3 cplx3_re(cplx3 v) { return _real3(v.x.re, v.y.re, v.z.re); }
real3 cplx3_im(cplx3 v) { return _real3(v.x.im, v.y.im, v.z.im); }

real cplx3_lenSq(cplx3 v) {
	return real3_lenSq(cplx3_re(v)) + real3_lenSq(cplx3_im(v));
}

real cplx3_len(cplx3 v) {
	return sqrt(cplx3_lenSq(v));
}

#define sym3_zero	_sym3(0,0,0,0,0,0)
#define sym3_ident	_sym3(1,0,0,1,0,1)

static inline real sym3_det(sym3 m) {
	return m.xx * (m.yy * m.zz - m.yz * m.yz)
		- m.xy * (m.xy * m.zz - m.xz * m.yz)
		+ m.xz * (m.xy * m.yz - m.yy * m.xz);
}

static inline sym3 sym3_inv(sym3 m, real det) {
	real invDet = 1. / det;
	return (sym3){
		.xx = (m.yy * m.zz - m.yz * m.yz) * invDet,
		.xy = (m.xz * m.yz - m.xy * m.zz) * invDet,
		.xz = (m.xy * m.yz - m.xz * m.yy) * invDet,
		.yy = (m.xx * m.zz - m.xz * m.xz) * invDet,
		.yz = (m.xz * m.xy - m.xx * m.yz) * invDet,
		.zz = (m.xx * m.yy - m.xy * m.xy) * invDet,
	};
}

static inline real3 sym3_real3_mul(sym3 m, real3 v) {
	return _real3(
		m.xx * v.x + m.xy * v.y + m.xz * v.z,
		m.xy * v.y + m.yy * v.y + m.yz * v.z,
		m.xz * v.z + m.yz * v.y + m.zz * v.z);
}

static inline sym3 sym3_add(sym3 a, sym3 b) {
	return (sym3){
		.xx = a.xx + b.xx,
		.xy = a.xy + b.xy,
		.xz = a.xz + b.xz,
		.yy = a.yy + b.yy,
		.yz = a.yz + b.yz,
		.zz = a.zz + b.zz,
	};
}

static inline sym3 sym3_sub(sym3 a, sym3 b) {
	return (sym3){
		.xx = a.xx - b.xx,
		.xy = a.xy - b.xy,
		.xz = a.xz - b.xz,
		.yy = a.yy - b.yy,
		.yz = a.yz - b.yz,
		.zz = a.zz - b.zz,
	};
}


static inline sym3 sym3_real_mul(sym3 a, real s) {
	return (sym3){
		.xx = a.xx * s,
		.xy = a.xy * s,
		.xz = a.xz * s,
		.yy = a.yy * s,
		.yz = a.yz * s,
		.zz = a.zz * s,
	};
}

//computes a^ij b_ij
static inline real sym3_dot(sym3 a, sym3 b) {
	return a.xx * b.xx + a.yy * b.yy + a.zz * b.zz
		+ 2. * (a.xy * b.xy + a.xz * b.xz + a.yz * b.yz);
}

static inline real3x3 sym3_sym3_mul(sym3 a, sym3 b) {
	real3x3 m;
<? for i=0,2 do
	for j=0,2 do
?>	m.v<?=i?>.s<?=j?> = 0.
<?		for k=0,2 do
?>		+ a.s<?=i<=k and i..k or k..i?> * b.s<?=k<=j and k..j or j..k?>
<?		end
?>	;
<?	end
end 
?>	return m;
}

static inline real3x3 real3x3_sym3_mul(real3x3 a, sym3 b) {
	real3x3 m;
<? for i=0,2 do
	for j=0,2 do
?>	m.v<?=i?>.s<?=j?> = 0.<?
		for k=0,2 do
?> + a.v<?=i?>.s<?=k?> * b.s<?=k<=j and k..j or j..k?><?
		end
?>;
<?	end
end 
?>	return m;
}

static inline sym3 real3x3_sym3_to_sym3_mul(real3x3 a, sym3 b) {
	sym3 m;
<? for i=0,2 do
	for j=i,2 do
?>	m.s<?=i?><?=j?> = 0.<?
		for k=0,2 do
?> + a.v<?=i?>.s<?=k?> * b.s<?=k<=j and k..j or j..k?><?
		end
?>;
<?	end
end 
?>	return m;
}

//c_ik = a_ij b_jk when you know c_ik is going to be symmetric
static inline sym3 sym3_real3x3_to_sym3_mul(sym3 a, real3x3 b) {
	sym3 m;
<? for i=0,2 do
	for j=i,2 do
?>	m.s<?=i?><?=j?> = 0.<?
		for k=0,2 do
?> + a.s<?=i<=k and i..k or k..i?> * b.v<?=k?>.s<?=j?><?
		end
?>;
<?	end
end 
?>	return m;
}

static inline real3 sym3_x(sym3 m) { return _real3(m.xx, m.xy, m.xz); }
static inline real3 sym3_y(sym3 m) { return _real3(m.xy, m.yy, m.yz); }
static inline real3 sym3_z(sym3 m) { return _real3(m.xz, m.yz, m.zz); }

static inline real sym3_trace(sym3 m) {
	return m.xx + m.yy + m.zz;
}

//for swapping dimensions between x and 012
sym3 sym3_swap0(sym3 m) { return m; }
sym3 sym3_swap1(sym3 m) { return _sym3(m.yy, m.xy, m.yz, m.xx, m.xz, m.zz); }
sym3 sym3_swap2(sym3 m) { return _sym3(m.zz, m.yz, m.xz, m.yy, m.xy, m.xx); }

#define real3x3_zero (real3x3){.v={real3_zero, real3_zero, real3_zero}}

static inline real3x3 real3x3_real3x3_mul(real3x3 a, real3x3 b) {
	real3x3 c;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			real sum = 0.;
			for (int k = 0; k < 3; ++k) {
				sum += a.v[i].s[k] * b.v[k].s[j];
			}
			c.v[i].s[j] = sum;
		}
	}
	return c;
}

static inline real3 real3x3_real3_mul(real3x3 a, real3 b) {
	real3 c;
	for (int i = 0; i < 3; ++i) {
		real sum = 0.;
		for (int j = 0; j < 3; ++j) {
			sum += a.v[i].s[j] * b.s[j];
		}
		c.s[i] = sum;
	}
	return c;
}

static inline real real3x3_trace(real3x3 m) {
	return m.x.x + m.y.y + m.z.z;
}

static inline real3x3 real3x3_from_real(real x) {
	real3x3 m;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			m.v[i].s[j] = (i == j) ? x : 0;
		}
	}
	return m;
}

static inline real real3x3_det(real3x3 m) {
	return m.x.x * (m.y.y * m.z.z - m.y.z * m.z.y)
		- m.x.y * (m.y.x * m.z.z - m.y.z * m.z.x)
		+ m.x.z * (m.y.x * m.z.y - m.y.y * m.z.x);
}

static inline real3x3 real3x3_inv(real3x3 m) {
	real det = real3x3_det(m);
	real3x3 n;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			int i1 = (i+1)%3;
			int i2 = (i+2)%3;
			int j1 = (j+1)%3;
			int j2 = (j+2)%3;
			n.v[i].s[j] = (
				m.v[j1].s[i1] * m.v[j2].s[i2] 
				- m.v[j2].s[i1] * m.v[j1].s[i2]
			) / det;
		}
	}
	//n.x.x = (m.y.y * m.z.z - m.z.y * m.y.z) / det;
	//n.x.y = (m.z.y * m.x.z - m.x.y * m.z.z) / det;
	//n.x.z = (m.x.y * m.y.z - m.y.y * m.x.z) / det;
	//...
	return n;
}


//weighted inner product using 'm'
static inline real real3_weightedDot(real3 a, real3 b, sym3 m) {
	return real3_dot(a, sym3_real3_mul(m, b));
}

static inline real real3_weightedLenSq(real3 a, sym3 m) {
	return real3_weightedDot(a, a, m);
}

static inline real real3_weightedLen(real3 a, sym3 m) {
	return sqrt(real3_weightedLenSq(a, m));
}

#define _3sym3_zero (_3sym3){.s={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}

<? for name,symbol in pairs{add='+', sub='-'} do ?>
static inline _3sym3 _3sym3_<?=name?>(_3sym3 a, _3sym3 b) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for jk,xjk in ipairs(symNames) do
?>			.<?=xjk?> = a.<?=xi?>.<?=xjk?> <?=symbol?> b.<?=xi?>.<?=xjk?>,
<?	end
?>		},
<? end
?>	};
}
<? end ?>

static inline _3sym3 _3sym3_real_mul(_3sym3 a, real b) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for jk,xjk in ipairs(symNames) do
?>			.<?=xjk?> = a.<?=xi?>.<?=xjk?> * b,
<?	end
?>		},
<? end
?>	};
}

//c^i_jk = a^il b_ljk
static inline _3sym3 sym3_3sym3_mul(sym3 a, _3sym3 b) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
?>			.<?=xjk?> = 0.<?
		for l,xl in ipairs(xNames) do
				?> + a.<?=sym(i,l)?> * b.<?=xl?>.<?=xjk?><? 
		end ?>,
<?	end
?>		},
<? end
?>	};
}

//c^i = a^i_jk b^jk
static inline real3 _3sym3_sym3_dot23(_3sym3 a, sym3 b) {
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + a.<?=xi?>.<?=sym(j,k)?> * b.<?=sym(j,k)?><?
		end
	end ?>,
<? end
?>	};
}

real3 normalForSide0() { return _real3(1,0,0); }
real3 normalForSide1() { return _real3(0,1,0); }
real3 normalForSide2() { return _real3(0,0,1); }
