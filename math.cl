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

//c^i_jk = a^il b_ljk
<? for _,addr in ipairs{'', 'global'} do
?>static void sym3_3sym3_mul_<?=addr?>(_3sym3 c, const sym3 a, <?=addr?> const _3sym3 b) {
<? 
for i,xi in ipairs(xNames) do
?>	c[<?=i-1?>] = (sym3){
<?	for jk,xjk in ipairs(symNames) do
?>		.<?=xjk?> = 0.<?
		for l,xl in ipairs(xNames) do
			?> + a.<?=sym(i,l)?> * b[<?=l-1?>].<?=xjk?><? 
		end ?>,
<?	end
?>	};
<? 
end
?>}
<? end
?>
