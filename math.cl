<?
local common = require 'common'()	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym
?>

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

#define _real3(a,b,c) (real3){.s={a,b,c}}
#define _sym3(a,b,c,d,e,f) (sym3){.s={a,b,c,d,e,f}}

static inline real real3_dot(real3 a, real3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline real3 real3_cross(real3 a, real3 b) {
	return _real3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

static inline real real3_lenSq(real3 a) {
	return real3_dot(a,a);
}

static inline real real3_len(real3 a) {
	return sqrt(real3_lenSq(a));
}

static inline real3 real3_scale(real3 a, real s) {
	return _real3(a.x * s, a.y * s, a.z * s);
}

static inline real3 real3_add(real3 a, real3 b) {
	return _real3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline real3 real3_sub(real3 a, real3 b) {
	return _real3(a.x - b.x, a.y - b.y, a.z - b.z);
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

static inline sym3 sym3_ident() {
	return _sym3(1,0,0,1,0,1);
}

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


static inline sym3 sym3_scale(sym3 a, real s) {
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

static inline mat3 sym3_sym3_mul(sym3 a, sym3 b) {
	mat3 m;
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

static inline mat3 mat3_sym3_mul(mat3 a, sym3 b) {
	mat3 m;
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

static inline sym3 mat3_sym3_to_sym3_mul(mat3 a, sym3 b) {
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
static inline sym3 sym3_mat3_to_sym3_mul(sym3 a, mat3 b) {
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


static inline mat3 mat3_mat3_mul(mat3 a, mat3 b) {
	mat3 c;
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

static inline real mat3_trace(mat3 m) {
	return m.x.x + m.y.y + m.z.z;
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

static inline _3sym3 _3sym3_scale(_3sym3 a, real b) {
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
