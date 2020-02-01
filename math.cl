// This file needs to be templated, so it must be processed and manually inserted
// I could write out the templated math.h file into the cache-cl folder, and then add that to the CL compiler include folder, and include it from there ...
//#include "math.h"

<?
local common = require 'common'	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local sym = common.sym
local from6to3x3 = common.from6to3x3
?>

////////////////////////// quat //////////////////////////

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

////////////////////////// real3 //////////////////////////

<?
function makevec3(vec,scalar)
	local add = scalar..'_add'
	local sub = scalar..'_sub'
	local mul = scalar..'_mul'
	local real_mul = scalar..'_real_mul'
	local neg = scalar..'_neg'
	local conj = scalar..'_conj'
	local add3 = scalar..'_add3'
	local mul3 = scalar..'_mul3'
?>

<?=vec?> <?=vec?>_<?=scalar?>_mul(<?=vec?> a, <?=scalar?> s) {
	return _<?=vec?>(
		<?=mul?>(a.x, s),
		<?=mul?>(a.y, s),
		<?=mul?>(a.z, s));
}

<?=vec?> <?=scalar?>_<?=vec?>_mul(<?=scalar?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=mul?>(a, b.x),
		<?=mul?>(a, b.y),
		<?=mul?>(a, b.z));
}

<? if scalar ~= 'real' then ?>
<?=vec?> <?=vec?>_real_mul(<?=vec?> a, real b) {
	return _<?=vec?>(
		<?=real_mul?>(a.x, b),
		<?=real_mul?>(a.y, b),
		<?=real_mul?>(a.z, b));
}

<?=vec?> real_<?=vec?>_mul(real a, <?=vec?> b) {
	return _<?=vec?>(
		<?=real_mul?>(b.x, a),
		<?=real_mul?>(b.y, a),
		<?=real_mul?>(b.z, a));
}
<? end ?>

<?=vec?> <?=vec?>_add(<?=vec?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=add?>(a.x, b.x), 
		<?=add?>(a.y, b.y), 
		<?=add?>(a.z, b.z));
}

<?=vec?> <?=vec?>_sub(<?=vec?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=sub?>(a.x, b.x),
		<?=sub?>(a.y, b.y),
		<?=sub?>(a.z, b.z));
}

<?=vec?> <?=vec?>_neg(<?=vec?> a) {
	return _<?=vec?>(<?=neg?>(a.x), <?=neg?>(a.y), <?=neg?>(a.z));
}

<?=vec?> <?=vec?>_cross(<?=vec?> a, <?=vec?> b) {
	return _<?=vec?>(
		<?=sub?>(<?=mul?>(a.y, b.z), <?=mul?>(a.z, b.y)),
		<?=sub?>(<?=mul?>(a.z, b.x), <?=mul?>(a.x, b.z)),
		<?=sub?>(<?=mul?>(a.x, b.y), <?=mul?>(a.y, b.x)));
}

<?
end
makevec3('real3', 'real')
?>

<?
local function resultType(atype, btype)
	if atype == 'real' and btype == 'real' then return 'real' end
	if atype == 'cplx' and btype == 'real' then return 'cplx' end
	if atype == 'real' and btype == 'cplx' then return 'cplx' end
	if atype == 'cplx' and btype == 'cplx' then return 'cplx' end
	error("no resultType for operations of "..atype.." and "..btype)
end

local function make_dot(atype, btype)
	local ctype = resultType(atype, btype)
?>
<?=ctype?> <?=atype?>3_<?=btype?>3_dot(<?=atype?>3 a, <?=btype?>3 b) {
	return <?=ctype?>_add3(
		<?=atype?>_<?=btype?>_mul(a.x, <?=btype?>_conj(b.x)), 
		<?=atype?>_<?=btype?>_mul(a.y, <?=btype?>_conj(b.y)),
		<?=atype?>_<?=btype?>_mul(a.z, <?=btype?>_conj(b.z))
	);
}
<?
end
make_dot('real', 'real')
?>

////////////////////////// real3 //////////////////////////

real real3_lenSq(real3 a) {
	return real3_dot(a,a);
}

real real3_len(real3 a) {
	return sqrt(real3_lenSq(a));
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

////////////////////////// sym3 //////////////////////////

//outer with yourself
sym3 real3_outer(real3 a) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = a.<?=xi?> * a.<?=xj?>,
<? end
?>	};
}

real sym3_det(sym3 m) {
	return m.xx * (m.yy * m.zz - m.yz * m.yz)
		- m.xy * (m.xy * m.zz - m.xz * m.yz)
		+ m.xz * (m.xy * m.yz - m.yy * m.xz);
}

sym3 sym3_inv(sym3 m, real det) {
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

real3 sym3_real3_mul(sym3 m, real3 v) {
	return _real3(
		m.xx * v.x + m.xy * v.y + m.xz * v.z,
		m.xy * v.y + m.yy * v.y + m.yz * v.z,
		m.xz * v.z + m.yz * v.y + m.zz * v.z);
}

sym3 sym3_neg(sym3 a) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = -a.<?=xij?>,
<? end
?>	};
}

sym3 sym3_add(sym3 a, sym3 b) {
	return (sym3){
		.xx = a.xx + b.xx,
		.xy = a.xy + b.xy,
		.xz = a.xz + b.xz,
		.yy = a.yy + b.yy,
		.yz = a.yz + b.yz,
		.zz = a.zz + b.zz,
	};
}

sym3 sym3_sub(sym3 a, sym3 b) {
	return (sym3){
		.xx = a.xx - b.xx,
		.xy = a.xy - b.xy,
		.xz = a.xz - b.xz,
		.yy = a.yy - b.yy,
		.yz = a.yz - b.yz,
		.zz = a.zz - b.zz,
	};
}


sym3 sym3_real_mul(sym3 a, real s) {
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
real sym3_dot(sym3 a, sym3 b) {
	return a.xx * b.xx + a.yy * b.yy + a.zz * b.zz
		+ 2. * (a.xy * b.xy + a.xz * b.xz + a.yz * b.yz);
}

real3 sym3_x(sym3 m) { return _real3(m.xx, m.xy, m.xz); }
real3 sym3_y(sym3 m) { return _real3(m.xy, m.yy, m.yz); }
real3 sym3_z(sym3 m) { return _real3(m.xz, m.yz, m.zz); }

real sym3_trace(sym3 m) {
	return m.xx + m.yy + m.zz;
}

//for swapping dimensions between x and 012
sym3 sym3_swap0(sym3 m) { return m; }
sym3 sym3_swap1(sym3 m) { return _sym3(m.yy, m.xy, m.yz, m.xx, m.xz, m.zz); }
sym3 sym3_swap2(sym3 m) { return _sym3(m.zz, m.yz, m.xz, m.yy, m.xy, m.xx); }

//weighted inner product using 'm'
real real3_weightedDot(real3 a, real3 b, sym3 m) {
	return real3_dot(a, sym3_real3_mul(m, b));
}

real real3_weightedLenSq(real3 a, sym3 m) {
	return real3_weightedDot(a, a, m);
}

real real3_weightedLen(real3 a, sym3 m) {
	return sqrt(real3_weightedLenSq(a, m));
}

////////////////////////// real3x3 //////////////////////////

real3x3 real3x3_add(real3x3 a, real3x3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3_add(a.<?=xi?>, b.<?=xi?>),
<? end
?>	};
}

/*
a^T * b
= a^i b_j
   [a.x b.x, a.x b.y, a.x b.z]
 = [a.y b.x, a.y b.y, a.y b.z]
   [a.z b.x, a.z b.y, a.z b.z]
*/
real3x3 real3_real3_outer(real3 a, real3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = a.<?=xi?> * b.<?=xj?>,
<?	end
?>		},
<? end
?>	};
}

real real3x3_dot(real3x3 a, real3x3 b) {
	return 0.
<? for i,xi in ipairs(xNames) do
?>		<?
	for j,xj in ipairs(xNames) do
?> + a.<?=xi?>.<?=xj?> * b.<?=xi?>.<?=xj?><?
	end ?>
<? end ?>;
}

real real3x3_sym3_dot(real3x3 a, sym3 b) {
	return 0.
<? for i,xi in ipairs(xNames) do
?>		<?
	for j,xj in ipairs(xNames) do
?> + a.<?=xi?>.<?=xj?> * b.<?=sym(i,j)?><?
	end ?>
<? end ?>;
}

real3x3 sym3_sym3_mul(sym3 a, sym3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
?> + a.<?=sym(i,k)?> * b.<?=sym(k,j)?><?
		end ?>,
<?	end
?>		},
<? end 
?>	};
}

real3x3 real3x3_sym3_mul(real3x3 a, sym3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
?> + a.<?=xi?>.<?=xk?> * b.<?=sym(k,j)?><?
		end ?>,
<?	end
?>		},
<? end 
?>	};
}

real3x3 sym3_real3x3_mul(sym3 a, real3x3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
?> + a.<?=sym(i,k)?> * b.<?=xk?>.<?=xj?><?
		end ?>,
<?	end
?>		},
<? end 
?>	};
}

sym3 real3x3_sym3_to_sym3_mul(real3x3 a, sym3 b) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.<?
		for k,xk in ipairs(xNames) do
?> + a.<?=xi?>.<?=xk?> * b.<?=sym(k,j)?><?
		end ?>,
<?	end
?>	};
}

//c_ik = a_ij b_jk when you know c_ik is going to be symmetric
sym3 sym3_real3x3_to_sym3_mul(sym3 a, real3x3 b) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.<?
	for k,xk in ipairs(xNames) do
?> + a.<?=sym(i,k)?> * b.<?=xk?>.<?=xj?><?
	end
?>,
<? end
?>	};
}

//a_ij = b_ij + b_ji
//doesn't do the .5 ... that is left to you
sym3 sym3_from_real3x3(real3x3 a) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = a.<?=xi?>.<?=xj?> + a.<?=xj?>.<?=xi?>,
<? end
?>	};
}


//c_ij = a_ij + b_ji
real3x3 real3x3_addT(real3x3 a, real3x3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = a.<?=xi?>.<?=xj?> + b.<?=xj?>.<?=xi?>,
<?	end
?>		},
<? end
?>	};
}

real3x3 real3x3_real3x3_mul(real3x3 a, real3x3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
?> + a.<?=xi?>.<?=xk?> * b.<?=xk?>.<?=xj?><? 
		end ?>,
<? 	end
?>		},
<? end
?>	};
}

<?
local function make_3x3_3_mul(atype, btype)
	local ctype = resultType(atype, btype)
?>
//c_i = a_ij * b^j
<?=ctype?>3 <?=atype?>3x3_<?=btype?>3_mul(<?=atype?>3x3 a, <?=btype?>3 b) {
	return (<?=ctype?>3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = <?=atype?>3_<?=btype?>3_dot(a.<?=xi?>, b),
<? end
?>	};
}
<?
end
make_3x3_3_mul('real', 'real')	-- real3x3_real3_mul
?>

<?
local function make_3_3x3_mul(atype, btype)
	local ctype = resultType(atype, btype)
?>
//c_j = a^i * b_ij = b_ji * a^i
//so this is the same behavior as transpose(A) * b
<?=ctype?>3 <?=atype?>3_<?=btype?>3x3_mul(<?=atype?>3 a, <?=btype?>3x3 b) {
	return (<?=ctype?>3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = <?=ctype?>_add3(
			<?=atype?>_<?=btype?>_mul(a.x, b.x.<?=xi?>),
			<?=atype?>_<?=btype?>_mul(a.y, b.y.<?=xi?>),
			<?=atype?>_<?=btype?>_mul(a.z, b.z.<?=xi?>)),
<? end
?>	};
}

<?
end
make_3_3x3_mul('real', 'real')	-- real3_real3x3_mul
?>

real real3x3_trace(real3x3 m) {
	return m.x.x + m.y.y + m.z.z;
}

real3x3 real3x3_from_real(real x) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = _real3(<?
	for j,xj in ipairs(xNames) do
		if i~=j then ?>0.<? else ?>x<? end ?><?= j < 3 and ', ' or ''?><?
	end ?>),
<? end
?>	};
}

real real3x3_det(real3x3 m) {
	return m.x.x * (m.y.y * m.z.z - m.y.z * m.z.y)
		- m.x.y * (m.y.x * m.z.z - m.y.z * m.z.x)
		+ m.x.z * (m.y.x * m.z.y - m.y.y * m.z.x);
}

real3x3 real3x3_inv(real3x3 m) {
	real invDet = 1. / real3x3_det(m);
	return (real3x3){
<? 
--n.x.x = (m.y.y * m.z.z - m.z.y * m.y.z) / det;
--n.x.y = (m.z.y * m.x.z - m.x.y * m.z.z) / det;
--n.x.z = (m.x.y * m.y.z - m.y.y * m.x.z) / det;
--...
for i=0,2 do
	local i1 = (i+1)%3;
	local i2 = (i+2)%3;
?>		.v<?=i?> = {
<? 	for j=0,2 do
		local j1 = (j+1)%3;
		local j2 = (j+2)%3;
?>			.s<?=j?> = (m.v<?=j1?>.s<?=i1?> * m.v<?=j2?>.s<?=i2?> - m.v<?=j2?>.s<?=i1?> * m.v<?=j1?>.s<?=i2?>) * invDet,
<? end
?>		},
<? end	
?>	};
}

////////////////////////// _3sym3 //////////////////////////

<? for name,symbol in pairs{add='+', sub='-'} do ?>
_3sym3 _3sym3_<?=name?>(_3sym3 a, _3sym3 b) {
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

_3sym3 _3sym3_real_mul(_3sym3 a, real b) {
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
_3sym3 sym3_3sym3_mul(sym3 a, _3sym3 b) {
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
real3 _3sym3_sym3_dot23(_3sym3 a, sym3 b) {
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(a.<?=xi?>, b),
<? end
?>	};
}

//c_i = a^jk b_jki
real3 sym3_3sym3_dot12(sym3 a, _3sym3 b) {
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + a.<?=sym(j,k)?> * b.<?=xj?>.<?=sym(k,i)?><?
		end
	end ?>,
<? end
?>	};
}

//c_ij = a_k b^k_ij
sym3 real3_3sym3_dot1(real3 a, _3sym3 b) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = 0.<?
	for k,xk in ipairs(xNames) do
?> + a.<?=xk?> * b.<?=xk?>.<?=xij?><?
	end ?>,
<? end
?>	};
}

//c_i = a^j_ji
real3 _3sym3_tr12(_3sym3 a) {
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
?> + a.<?=xj?>.<?=sym(j,i)?><?
	end ?>,
<? end
?>	};
}

//c_ij = a^k b_ikj
real3x3 real3_3sym3_dot2(real3 a, _3sym3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
?> + a.<?=xk?> * b.<?=xi?>.<?=sym(k,j)?><?
		end ?>,
<?	end
?>		},
<? end
?>	};
}


////////////////////////// real3x3x3 //////////////////////////

//c_ij^k = a_ijl b^lk
real3x3x3 _3sym3_sym3_mul(_3sym3 a, sym3 b) {
	return (real3x3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = {
<?		for k,xk in ipairs(xNames) do
?>				.<?=xk?> = 0.<?
			for l,xl in ipairs(xNames) do
?> + a.<?=xi?>.<?=sym(j,l)?> * b.<?=sym(l,k)?><?
			end ?>,
<? end
?>			},
<? end
?>		},
<? end
?>	};
}

real3 real3x3x3_sym3_dot23(real3x3x3 a, sym3 b) {
	return (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_sym3_dot(a.<?=xi?>, b),
<? end
?>	};
}

//c_ij = a^k_li b_jk^l
real3x3 _3sym3_real3x3x3_dot12_23(_3sym3 a, real3x3x3 b) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
?>				+ a.<?=xk?>.<?=sym(l,i)?> * b.<?=xj?>.<?=xk?>.<?=xl?>
<?			end
		end
?>			,
<? 	end
?>		},
<? end
?>	};
}

//c_ij = a^k_il b_kj^l
//assuming the result is symmetric
sym3 _3sym3_real3x3x3_dot13_to_sym3(_3sym3 a, real3x3x3 b) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>			+ a.<?=xk?>.<?=sym(i,l)?> * b.<?=xk?>.<?=xj?>.<?=xl?>
<?		end
	end
?>		,
<? end
?>	};
}

////////////////////////// sym3sym3 //////////////////////////

sym3sym3 sym3sym3_add(sym3sym3 a, sym3sym3 b) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = sym3_add(a.<?=xij?>, b.<?=xij?>),
<? end
?>	};
}

////////////////////////// cplx //////////////////////////

cplx cplx_conj(cplx a) { return _cplx(a.re, -a.im); }
cplx cplx_neg(cplx a) { return _cplx(-a.re, -a.im); }
real cplx_lenSq(cplx a) { return a.re * a.re + a.im * a.im; }
real cplx_abs(cplx a) { return sqrt(cplx_lenSq(a)); }
real cplx_arg(cplx a) { return atan2(a.im, a.re); }

//what do I call this function? complex dot? real component of complex dot?  real component of complex mul of a number with another's conjugate?
real cplx_dot(cplx a, cplx b) {
	return a.re * b.re + a.im * b.im;
}

cplx cplx_add(cplx a, cplx b) {
	return _cplx(
		a.re + b.re,
		a.im + b.im
	);
}

cplx cplx_sub(cplx a, cplx b) {
	return _cplx(
		a.re - b.re,
		a.im - b.im
	);
}

cplx cplx_mul(cplx a, cplx b) {
	return _cplx(
		a.re * b.re - a.im * b.im,
		a.re * b.im + a.im * b.re);
}

cplx cplx_real_add(cplx a, real b) { 
	return _cplx(a.re + b, a.im); 
}

cplx real_cplx_add(real a, cplx b) { 
	return _cplx(b.re + a, b.im); 
}

cplx cplx_real_mul(cplx a, real b) { 
	return _cplx(a.re * b, a.im * b); 
}

cplx cplx_inv(cplx a) { 
	return cplx_real_mul(cplx_conj(a), 1. / cplx_lenSq(a)); 
}

cplx cplx_div(cplx a, cplx b) { 
	return cplx_mul(a, cplx_inv(b)); 
}

cplx cplx_exp(cplx a) {
	real expre = exp(a.re);
	return _cplx(
		expre * cos(a.im),
		expre * sin(a.im)
	);
}

cplx cplx_log(cplx a) {
	return _cplx(
		log(cplx_abs(a)),
		cplx_arg(a)
	);
}

cplx cplx_pow(cplx a, cplx b) { return cplx_exp(cplx_mul(b, cplx_log(a))); }
cplx cplx_sqrt(cplx a) { return cplx_pow(a, cplx_from_real(.5)); }

////////////////////////// cplx3 //////////////////////////

<?
makevec3('cplx3', 'cplx')
?>

cplx3 cplx3_from_real3(real3 re) {
	return (cplx3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = cplx_from_real(re.<?=xi?>),
<? end
?>	};
}

cplx3 cplx3_from_real3_real3(real3 re, real3 im) {
	return (cplx3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = _cplx(re.<?=xi?>, im.<?=xi?>),
<? end
?>	};
}

real3 cplx3_re(cplx3 v) { 
	return _real3(v.x.re, v.y.re, v.z.re); 
}

real3 cplx3_im(cplx3 v) { 
	return _real3(v.x.im, v.y.im, v.z.im); 
}

real cplx3_lenSq(cplx3 v) {
	return cplx_lenSq(v.x) + cplx_lenSq(v.y) + cplx_lenSq(v.z);
}

real cplx3_len(cplx3 v) {
	return sqrt(cplx3_lenSq(v));
}

cplx3 real3_cplx_mul(real3 a, cplx b) {
	return _cplx3(
		real_cplx_mul(a.x, b),
		real_cplx_mul(a.y, b),
		real_cplx_mul(a.z, b));
}

<?
make_dot('cplx', 'real')	-- cplx3_real3_dot
?>

/*
assumes the weights are real and only returns the real component ...
what to properly name this
a^i b^j* g_ij
= (re(a^i) re(b^j) + im(a^i) im(a^j)) g_ij
= re(a^i) re(b^j) g_ij + im(a^i) im(a^j) g_ij
*/
real cplx3_re_weightedDot(cplx3 a, cplx3 b, sym3 m) {
	return real3_weightedDot(cplx3_re(a), cplx3_re(b), m)
		+ real3_weightedDot(cplx3_im(a), cplx3_im(b), m);
}

real cplx3_weightedLenSq(cplx3 a, sym3 m) {
	return real3_weightedLenSq(cplx3_re(a), m)
		+ real3_weightedLenSq(cplx3_im(a), m);
}

////////////////////////// cplx3x3 //////////////////////////

cplx3x3 cplx3x3_from_real3x3_real3x3(real3x3 re, real3x3 im) {
	return (cplx3x3) {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = cplx3_from_real3_real3(re.<?=xi?>, im.<?=xi?>),
<? end
?>	};
}

real3x3 cplx3x3_re(cplx3x3 v) { 
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = cplx3_re(v.<?=xi?>),
<? end
?>	};
}

real3x3 cplx3x3_im(cplx3x3 v) { 
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = cplx3_im(v.<?=xi?>),
<? end
?>	};
}

<?
make_3x3_3_mul('cplx', 'real')	-- cplx3x3_real3_mul
make_3x3_3_mul('real', 'cplx')	-- real3x3_cplx3_mul
make_3_3x3_mul('cplx', 'real')	-- cplx3_real3x3_mul
?>

cplx cplx3x3_sym3_dot(cplx3x3 a, sym3 b) {
	cplx sum = cplx_zero;
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?> 	sum = cplx_add(sum, cplx_real_mul(a.<?=xi?>.<?=xj?>, b.<?=sym(i,j)?>));
<?	end
end 
?>	return sum;
}

////////////////////////// extra //////////////////////////


real3 normalForSide0() { return _real3(1,0,0); }
real3 normalForSide1() { return _real3(0,1,0); }
real3 normalForSide2() { return _real3(0,0,1); }

// https://community.amd.com/thread/169701
// meh, not so great
float crand() {
	unsigned int seed = get_global_id(0)
		+ 13 * get_global_id(1)
		+ 87 * get_global_id(2);
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	seed = ((seed) * 16807 ) % 2147483647;
	return (float)(seed) * 4.6566129e-10;
}

real sech(real x) {
	return 2. / (exp(x) + exp(-x));
}
