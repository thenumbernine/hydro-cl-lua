#define _real3(a,b,c) (real3){.s={a,b,c}}

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

#define _sym3(a,b,c,d,e,f) (sym3){.s={a,b,c,d,e,f}}

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


sym3 sym3_scale(sym3 a, real s) {
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

mat3 sym3_sym3_mul(sym3 a, sym3 b) {
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

mat3 mat3_sym3_mul(mat3 a, sym3 b) {
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

sym3 mat3_sym3_to_sym3_mul(mat3 a, sym3 b) {
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
sym3 sym3_mat3_to_sym3_mul(sym3 a, mat3 b) {
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

real3 sym3_x(sym3 m) { return _real3(m.xx, m.xy, m.xz); }
real3 sym3_y(sym3 m) { return _real3(m.xy, m.yy, m.yz); }
real3 sym3_z(sym3 m) { return _real3(m.xz, m.yz, m.zz); }

real sym3_trace(sym3 m) {
	return m.xx + m.yy + m.zz;
}

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

mat3 mat3_mat3_mul(mat3 a, mat3 b) {
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

real mat3_trace(mat3 m) {
	return m.x.x + m.y.y + m.z.z;
}
