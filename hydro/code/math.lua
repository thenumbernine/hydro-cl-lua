-- This file needs to be templated, so it must be processed and manually inserted
-- I could write out the templated hydro/code/math.h file into the cache-bin folder, and then add that to the CL compiler include folder, and include it from there ...

local template = require 'template'

local common = require 'hydro.common'	-- xNames, symNames
local xNames = common.xNames
local sym = common.sym
local symNames = common.symNames
local from6to3x3 = common.from6to3x3
local from3x3to6 = common.from3x3to6

return function(modules, app)
	-- used for passing parameters into real, distinct from real in the case that real=half which can't be passed as kernel arguments
	modules:add{
		name = 'realparam',
		typecode = template([[
typedef <?=app.realparam?> realparam;
typedef <?=app.realparam?>2 realparam2;
typedef <?=app.realparam?>4 realparam4;
typedef <?=app.realparam?>8 realparam8;
]], {app=app}),
	}

	modules:add{
		name = 'quat',
		headercode = [[
real4 quatUnitConj(real4 q);
real4 quatMul(real4 q, real4 r);
]],
		code = [[

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

]],
	}

	local function makevec3type(name, scalar)
		if scalar == 'real' 
		and (	
			app.real == 'float'
			or app.real == 'double'
		) then
			-- use vec-ffi when we can 
			assert(name == 'real3')
			
			local vecType
			if app.real == 'float' then
				vecType = require 'vec-ffi.vec3f'
			elseif app.real == 'double' then
				vecType = require 'vec-ffi.vec3d'
			end
			-- use the vec-ffi type code
			-- granted if we ffi.cdef this, it will have already been ffi.cdef'd from the require 'vec-ffi.vec3x'
			return template([[
<?=vecType.typeCode?>
typedef <?=vecType.type?> <?=name?>;
]], {name=name, vecType=vecType})		
		else 
			-- normal
			return template([[
typedef union {
	<?=scalar?> s[3];
	struct { <?=scalar?> s0, s1, s2; };
	struct { <?=scalar?> x, y, z; };
} <?=app.real=='half' and '__attribute__ ((packed))' or ''
-- __attribute__ ((packed)) seems to need to be here with real=half
?> <?=name?>;
]], {name=name, scalar=scalar, app=app})
		end
	end 

	local function makevec3header(vec, scalar)
		return template([[

<? -- TODO move this to makescalar ?>
#define <?=scalar?>_add3(a,b,c)		(<?=add?>(a, <?=add?>(b,c)))
#define <?=scalar?>_add4(a,b,c,d)	(<?=add?>(a, <?=scalar?>_add3(b,c,d)))
#define <?=scalar?>_add5(a,b,c,d,e) (<?=add?>(a, <?=scalar?>_add4(b,c,d,e)))

#define <?=scalar?>_mul3(a,b,c)		(<?=mul?>(<?=mul?>(a,b),c))

//is buggy with doubles on intel opencl ubuntu compiler
//#define _<?=vec?>(a,b,c) 			((<?=vec?>){.s={a,b,c}})
//so we do this instead and are safe:
#define _<?=vec?>(a,b,c) 			((<?=vec?>){.x=a, .y=b, .z=c})

#define <?=vec?>_zero				_<?=vec?>(<?=scalar?>_zero,<?=scalar?>_zero,<?=scalar?>_zero)

<?=scalar?> <?=vec?>_<?=vec?>_dot(<?=vec?> a, <?=vec?> b);
#define	<?=vec?>_dot <?=vec?>_<?=vec?>_dot
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

#define <?=vec?>_add3(a,b,c)	<?=vec?>_add(<?=vec?>_add(a,b),c)
#define <?=vec?>_add4(a,b,c,d)	<?=vec?>_add(<?=vec?>_add(a,b),<?=vec?>_add(c,d))

]], 	{
			vec = vec,
			scalar = scalar,
			add = scalar..'_add',
			mul = scalar..'_mul',
		})
	end

	local function makevec3code(vec,scalar)
		return template([[

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

<?=vec?> <?=vec?>_unit(<?=vec?> a) {
	return <?=vec?>_real_mul(a, 1. / <?=vec?>_len(a));
}
]], {
			vec = vec,
			scalar = scalar,
			add = scalar..'_add',
			sub = scalar..'_sub',
			mul = scalar..'_mul',
			real_mul = scalar..'_real_mul',
			neg = scalar..'_neg',
			conj = scalar..'_conj',
			add3 = scalar..'_add3',
			mul3 = scalar..'_mul3',
		})
	end


	local function make3x3type(scalar)
		return template([[
typedef union {
	<?=scalar?> s[9];
	<?=vec3?> v[3];
	struct {
		<?=vec3?> v0,v1,v2;
	};
	struct {
		<?=vec3?> x,y,z;
	};
} <?=name?>;
]], 	{
			scalar = scalar,
			vec3 = scalar..'3',
			name = scalar..'3x3',
		})
	end

	local function makecplx(name, real) 
		return template([[
typedef union {
	<?=real?> s[2];
	struct { <?=real?> s0, s1; };
	struct { <?=real?> re, im; };
} <?=name?>;
]], 	{
			name = name,
			real = real,
		})
	end 

	local function resultType(atype, btype)
		if atype == 'real' and btype == 'real' then return 'real' end
		if atype == 'cplx' and btype == 'real' then return 'cplx' end
		if atype == 'real' and btype == 'cplx' then return 'cplx' end
		if atype == 'cplx' and btype == 'cplx' then return 'cplx' end
		error("no resultType for operations of "..atype.." and "..btype)
	end

	local function make_dot(atype, btype)
		return template([[
<?=ctype?> <?=atype?>3_<?=btype?>3_dot(<?=atype?>3 a, <?=btype?>3 b) {
	return <?=ctype?>_add3(
		<?=atype?>_<?=btype?>_mul(a.x, <?=btype?>_conj(b.x)), 
		<?=atype?>_<?=btype?>_mul(a.y, <?=btype?>_conj(b.y)),
		<?=atype?>_<?=btype?>_mul(a.z, <?=btype?>_conj(b.z))
	);
}
]],		{
			atype = atype,
			btype = btype,
			ctype = resultType(atype, btype),
		})
	end


	local function make_3x3_3_mul(atype, btype)
		return template([[
//c_i = a_ij * b^j
<?=ctype?>3 <?=atype?>3x3_<?=btype?>3_mul(<?=atype?>3x3 a, <?=btype?>3 b) {
	return (<?=ctype?>3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = <?=atype?>3_<?=btype?>3_dot(a.<?=xi?>, b),
<? end
?>	};
}
]], 	{
			atype = atype,
			btype = btype,
			ctype = resultType(atype, btype),
			xNames = xNames,
		})
	end

	local function make_3_3x3_mul(atype, btype)
		return template([[
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
]], 	{
			atype = atype,
			btype = btype,
			ctype = resultType(atype, btype),
			xNames = xNames,
		})
	end

	modules:add{
		name = 'units',
		headercode = [[
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
]],
	}

	modules:add{
		name = 'real',
		headercode = [[

#define real_conj(x)		(x)
#define real_from_real(x)	(x)
#define real_zero			0
#define real_neg(x)			(-(x))
#define real_inv(x)			(1./(x))
#define real_add(a,b)		((a) + (b))
#define real_sub(a,b)		((a) - (b))
#define real_mul(a,b)		((a) * (b))
#define real_div(a,b)		((a) / (b))
#define real_real_add(a,b)	((a) + (b))
#define real_real_sub(a,b)	((a) - (b))
#define real_real_mul(a,b)	((a) * (b))
#define real_real_div(a,b)	((a) / (b))
#define real_lenSq(x)		((x) * (x))
#define real_abs			fabs

#define real_sqrt			sqrt

]],
	}

	modules:add{
		name = 'real3',
		depends = {'real'},	
		typecode = makevec3type('real3', 'real'),
		headercode = template([[

<?=makevec3header('real3', 'real')?>

#define real3_from_real3(x)	x

real real3_lenSq(real3 a);
real real3_len(real3 a);

]], 	{
			app = app,
			makevec3header = makevec3header,
		}),
		
		code = template([[

<?=makevec3code('real3', 'real')?>
<?=make_dot('real', 'real')?>

real real3_lenSq(real3 a) {
	return real3_dot(a,a);
}

real real3_len(real3 a) {
	return sqrt(real3_lenSq(a));
}

]], 	{
			app = app,
			makevec3code = makevec3code,
			make_dot = make_dot,
		}),
	}

	modules:add{
		name = 'rotate',
		depends = {'real3'},
		headercode = [[

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

]],
		code = [[

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

real3 real3_rotateX(real3 v, real theta) {
	real cosTheta = cos(theta);
	real sinTheta = sin(theta);
	real vy = v.y;
	real vz = v.z;
	v.y = vy * cosTheta - vz * sinTheta;
	v.z = vy * sinTheta + vz * cosTheta;
	return v;
}

real3 real3_rotateY(real3 v, real theta) {
	real cosTheta = cos(theta);
	real sinTheta = sin(theta);
	real vx = v.x;
	real vz = v.z;
	v.x = vx * cosTheta + vz * sinTheta;
	v.z = -vx * sinTheta + vz * cosTheta;
	return v;
}

real3 real3_rotateZ(real3 v, real theta) {
	real cosTheta = cos(theta);
	real sinTheta = sin(theta);
	real vx = v.x;
	real vy = v.y;
	v.x = vx * cosTheta - vy * sinTheta;
	v.y = vx * sinTheta + vy * cosTheta;
	return v;
}



]],
	}

	modules:add{
		name = 'sym3',
		typecode = [[
typedef union {
	real s[6];
	struct {
		real xx, xy, xz, yy, yz, zz;
	};
	struct {
		real s00, s01, s02, s11, s12, s22;
	};
} sym3;
]],
		headercode = [[

//buggy on intel opencl ubuntu
//#define _sym3(a,b,c,d,e,f) ((sym3){.s={a,b,c,d,e,f}})
//fix
#define _sym3(a,b,c,d,e,f) ((sym3){.xx=a, .xy=b, .xz=c, .yy=d, .yz=e, .zz=f})

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

]],
		code = template([[

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

]],		{
			xNames = xNames,
			symNames = symNames,
			from6to3x3 = from6to3x3,
		}),
	}

	modules:add{
		name = 'real3x3',
		depends = {'real3','sym3'},
		typecode = '// row vectors, so a.i.j = a_ij\n'
					..make3x3type'real',
		headercode = [[

//specified in row-order, like you would a C array
#define _real3x3(xx,xy,xz,yx,yy,yz,zx,zy,zz) ((real3x3){.x=_real3(xx,xy,xz), .y=_real3(yx,yy,yz), .z=_real3(zx,zy,zz)})

#define real3x3_zero ((real3x3){.v={real3_zero, real3_zero, real3_zero}})

real3x3 real3x3_add(real3x3 a, real3x3 b);
real3x3 real3_real3_outer(real3 a, real3 b);
real real3x3_dot(real3x3 a, real3x3 b);
real real3x3_sym3_dot(real3x3 a, sym3 b);
real3x3 sym3_sym3_mul(sym3 a, sym3 b);
real3x3 real3x3_sym3_mul(real3x3 a, sym3 b);
real3x3 sym3_real3x3_mul(sym3 a, real3x3 b);
sym3 real3x3_sym3_to_sym3_mul(real3x3 a, sym3 b);
sym3 sym3_real3x3_to_sym3_mul(sym3 a, real3x3 b);
sym3 sym3_from_real3x3(real3x3 a);
real3x3 real3x3_from_sym3(sym3 a);
real3x3 real3x3_addT(real3x3 a, real3x3 b);
real3x3 real3x3_real3x3_mul(real3x3 a, real3x3 b);
real3 real3x3_real3_mul(real3x3 a, real3 b);
real real3x3_trace(real3x3 m);
real3x3 real3x3_from_real(real x);
real real3x3_det(real3x3 m);
real3x3 real3x3_inv(real3x3 m);

]],
		code = template([[


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

real3x3 real3x3_from_sym3(sym3 a) {
	return (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
		local ij,xij = from3x3to6(i,j)
?>			.<?=xj?> = a.<?=xij?>,
<?	end
?>		},
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

<?=
make_3x3_3_mul('real', 'real')
-- real3x3_real3_mul
?>

<?=
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

]], 	{
			xNames = xNames,
			sym = sym,
			symNames = symNames,
			from6to3x3 = from6to3x3,
			from3x3to6 = from3x3to6,
			make_3x3_3_mul = make_3x3_3_mul,
			make_3_3x3_mul = make_3_3x3_mul,
		}),
	}

	modules:add{
		name = '_3sym3',
		depends = {'sym3', 'real3x3'},
		typecode = [[
typedef union {
	real s[18];
	struct {
		sym3 v0,v1,v2;	//why not s0,s1,s2?
	};
	struct {
		sym3 x,y,z;
	};
} _3sym3;
]],
		headercode = [[
//buggy
//#define _3sym3_zero ((_3sym3){.s={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}})
//fixed
#define _3sym3_zero ((_3sym3){.x=sym3_zero, .y=sym3_zero, .z=sym3_zero})

_3sym3 _3sym3_add(_3sym3 a, _3sym3 b);
_3sym3 _3sym3_sub(_3sym3 a, _3sym3 b);
_3sym3 _3sym3_real_mul(_3sym3 a, real b);
_3sym3 sym3_3sym3_mul(sym3 a, _3sym3 b);
real3 _3sym3_sym3_dot23(_3sym3 a, sym3 b);
real3 sym3_3sym3_dot12(sym3 a, _3sym3 b);
sym3 real3_3sym3_dot1(real3 a, _3sym3 b);
real3 _3sym3_tr12(_3sym3 a);
real3x3 real3_3sym3_dot2(real3 a, _3sym3 b);
]],
		code = template([[

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



]], 	{
			xNames = xNames,
			sym = sym,
			symNames = symNames,
		}),
	}

	modules:add{
		name = 'real3x3x3',
		depends = {'real3', 'sym3', 'real3x3', '_3sym3'},
		typecode = [[

typedef union {
	real s[27];
	real3x3 v[3];
	struct {
		real3x3 v0,v1,v2;
	};
	struct {
		real3x3 x,y,z;
	};
} real3x3x3;

]],
		headercode = [[

real3x3x3 _3sym3_sym3_mul(_3sym3 a, sym3 b);
real3 real3x3x3_sym3_dot23(real3x3x3 a, sym3 b);
real3x3 _3sym3_real3x3x3_dot12_23(_3sym3 a, real3x3x3 b);
sym3 _3sym3_real3x3x3_dot13_to_sym3(_3sym3 a, real3x3x3 b);

]],
		code = template([[

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


]],		{
			xNames = xNames,
			sym = sym,
			symNames = symNames,
			from6to3x3 = from6to3x3,
		}),
	}

	modules:add{
		name = 'sym3sym3',
		depends = {'sym3'},
		typecode = [[

typedef union {
	real s[36];
	sym3 v[6];
	struct {
		sym3 v0, v1, v2, v3, v4, v5;
	};
	struct {
		sym3 xx, xy, xz, yy, yz, zz;
	};
} sym3sym3;

]],
		headercode = [[

sym3sym3 sym3sym3_add(sym3sym3 a, sym3sym3 b);

]],
		code = template([[

sym3sym3 sym3sym3_add(sym3sym3 a, sym3sym3 b) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = sym3_add(a.<?=xij?>, b.<?=xij?>),
<? end
?>	};
}

]],		{
			symNames = symNames,
		}),
	}

	modules:add{
		name = 'cplx',
		typecode = makecplx('cplx', 'real'),
		headercode = [[

#define real_from_cplx(x)	((x).re)

//buggy
//#define _cplx(a,b) 			((cplx){.s={a,b}})
//fixed
#define _cplx(a,b) 			((cplx){.re=a, .im=b})

#define cplx_from_real(x)	_cplx(x,0)
#define cplx_from_cplx(x)	(x)
#define cplx_1				_cplx(1,0)
#define cplx_i	 			_cplx(0,1)
#define cplx_zero 			cplx_from_real(0)


cplx cplx_conj(cplx a);
cplx cplx_neg(cplx a);
real cplx_lenSq(cplx a);
real cplx_abs(cplx a);
real cplx_arg(cplx a);
real cplx_dot(cplx a, cplx b);
cplx cplx_add(cplx a, cplx b);
cplx cplx_sub(cplx a, cplx b);
cplx cplx_mul(cplx a, cplx b);
cplx cplx_real_mul(cplx a, real b);
#define real_cplx_mul(a,b)	cplx_real_mul(b,a)
cplx cplx_inv(cplx a);
cplx cplx_div(cplx a, cplx b);
cplx cplx_exp(cplx a);
cplx cplx_log(cplx a);
cplx cplx_pow(cplx a, cplx b);
cplx cplx_sqrt(cplx a);

]],
		code = [[


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



]],
	}

	modules:add{
		name = 'cplx3',
		depends = {'cplx', 'sym3'},
		typecode = makevec3type('cplx3', 'cplx'),
		headercode = makevec3header('cplx3', 'cplx')..[[

#define real3_from_cplx3		cplx3_re

cplx3 cplx3_from_real3(real3 re);
cplx3 cplx3_from_real3_real3(real3 re, real3 im);
real3 cplx3_re(cplx3 v);
real3 cplx3_im(cplx3 v);
real cplx3_lenSq(cplx3 v);
real cplx3_len(cplx3 v);
cplx3 real3_cplx_mul(real3 a, cplx b);
cplx cplx3_real3_dot(cplx3 a, real3 b);
#define real3_cplx3_dot(a,b) cplx3_real3_dot(b,a)

//TODO instead of sym3 as a depends, put this into a cplx+sym3 module? 
real cplx3_re_weightedDot(cplx3 a, cplx3 b, sym3 m);
real cplx3_weightedLenSq(cplx3 a, sym3 m);

]],
		code = template([[

<?=makevec3code('cplx3', 'cplx')
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

<?=make_dot('cplx', 'real')	-- cplx3_real3_dot
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


]],		{
			app = app,
			xNames = xNames,
			makevec3code = makevec3code,
			make_dot = make_dot,
		}),
	}

	modules:add{
		name = 'cplx3x3',
		depends = {'sym3', 'real3x3', 'cplx3'},
		typecode = make3x3type'cplx'..[[

cplx3x3 cplx3x3_from_real3x3_real3x3(real3x3 re, real3x3 im);
cplx3 cplx3x3_real3_mul(cplx3x3 a, real3 b);
cplx3 cplx3_real3x3_mul(cplx3 a, real3x3 b);
real3x3 cplx3x3_re(cplx3x3 v);
real3x3 cplx3x3_im(cplx3x3 v);
]],
		code = template([[
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

<?=
make_3x3_3_mul('cplx', 'real')		-- cplx3x3_real3_mul
..make_3x3_3_mul('real', 'cplx')	-- real3x3_cplx3_mul
..make_3_3x3_mul('cplx', 'real')	-- cplx3_real3x3_mul
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

]], 	{
			xNames = xNames,
			sym = sym,
			make_3x3_3_mul = make_3x3_3_mul,
			make_3_3x3_mul = make_3_3x3_mul,
		}),
	}

	modules:add{
		name = 'sech',
		code = [[
real sech(real x) {
	return 2. / (exp(x) + exp(-x));
}
]],
	}

	modules:add{
		name = 'Bessel',
		code = [[

// Reference: From Numath Library By Tuan Dang Trong in Fortran 77.
// C++ Release 1.0 By J-P Moreau, Paris.
// (www.jpmoreau.fr)
real BESSJ0(real X) {
	/***********************************************************************
	This subroutine calculates the First Kind Bessel Function of
	order 0, for any real number X. The polynomial approximation by
	series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	REFERENCES:
	M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	VOL.5, 1962.
	************************************************************************/
	if (X==0.0) return 1.0;
	real AX = fabs(X);
	if (AX < 8.0) {
		const real
			R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
			R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
			S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
			S4= 59272.64853, S5=267.8532712, S6=1.0;
		real Y = X*X;
		real FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
		real FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
		return FR/FS;
	} else {
		const real
			P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
			P4=-0.2073370639E-5, P5= 0.2093887211E-6,
			Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
			Q4= 0.7621095161E-6, Q5=-0.9349451520E-7;	
		real Z = 8./AX;
		real Y = Z*Z;
		real XX = AX-0.785398164;
		real FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
		real FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
		return sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
	}
}


real Sign(real X, real Y) {
  if (Y<0.0) return (-fabs(X));
  else return (fabs(X));
}

// Reference: From Numath Library By Tuan Dang Trong in Fortran 77.
// C++ Release 1.0 By J-P Moreau, Paris.
// (www.jpmoreau.fr)
real BESSJ1(real X) {
	/**********************************************************************
	This subroutine calculates the First Kind Bessel Function of
	order 1, for any real number X. The polynomial approximation by
	series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	REFERENCES:
	M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	VOL.5, 1962.
	***********************************************************************/
	const real
		P1=1.0, P2=0.183105E-2, P3=-0.3516396496E-4, P4=0.2457520174E-5,
		P5=-0.240337019E-6,  P6=0.636619772,
		Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
		Q4=-0.88228987E-6, Q5= 0.105787412E-6,
		R1= 72362614232.0, R2=-7895059235.0, R3=242396853.1,
		R4=-2972611.439,   R5=15704.48260,  R6=-30.16036606,
		S1=144725228442.0, S2=2300535178.0, S3=18583304.74,
		S4=99447.43394,    S5=376.9991397,  S6=1.0;

	real AX = fabs(X);
	if (AX < 8.0) {
		real Y = X*X;
		real FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
		real FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
		return X*(FR/FS);
	} else {
		real Z = 8.0/AX;
		real Y = Z*Z;
		real XX = AX-2.35619491;
		real FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
		real FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
		return sqrt(P6/AX)*(cos(XX)*FP-Z*sin(XX)*FQ)*Sign(S6,X);
	}
}


]],
	}

	modules:add{
		name = 'getPerpendicularBasis',
		depends = {'real3x3'},
		code = [[

//TODO include metric weight
//TODO I'm using this for n^u and n_u
void getPerpendicularBasis(real3 n, real3* n2, real3* n3) {
	real3 n_x_x = real3_cross(n, _real3(1,0,0));
	real3 n_x_y = real3_cross(n, _real3(0,1,0));
	real3 n_x_z = real3_cross(n, _real3(0,0,1));
	real n_x_xSq = real3_lenSq(n_x_x);
	real n_x_ySq = real3_lenSq(n_x_y);
	real n_x_zSq = real3_lenSq(n_x_z);
	if (n_x_xSq > n_x_ySq) {
		if (n_x_xSq > n_x_zSq) {
			*n2 = n_x_x;	//use x
		} else {
			*n2 = n_x_z;	//use z
		}
	} else {
		if (n_x_ySq > n_x_zSq) {
			*n2 = n_x_y;	//use y
		} else {
			*n2 = n_x_z;	//use z
		}
	}
	*n2 = real3_unit(*n2);
	*n3 = real3_cross(n, *n2);
}

//based on n->x, calculate n->y and n->z
//same as above, but for the row vectors of a 3x3 matrix
void getPerpendicularBasis3x3(real3x3* n) {
	getPerpendicularBasis(n->x, &n->y, &n->z);
}

]],
	}

	modules:add{
		name = 'normalForSide',
		headercode = [[
		
real3 normalForSide0();
real3 normalForSide1();
real3 normalForSide2();

]],
		code = [[

#define normalForSide0 _real3(1,0,0)
#define normalForSide1 _real3(0,1,0)
#define normalForSide2 _real3(0,0,1)

#define normalBasisForSide0 _real3x3(1,0,0, 0,1,0, 0,0,1)
#define normalBasisForSide1 _real3x3(0,1,0, 0,0,1, 1,0,0)
#define normalBasisForSide2 _real3x3(0,0,1, 1,0,0, 0,1,0)

]],
	}

	-- this is for GLSL esp which doesn't have these defs
	modules:add{
		name = 'sinh',
		code = [[
real sinh(real x) { 
	real ex = exp(x);
	return .5 * (ex - 1. / ex);
}
]],
	}
	
	modules:add{
		name = 'cosh',
		code = [[
real cosh(real x) { 
	real ex = exp(x);
	return .5 * (ex + 1. / ex);
}
]],
	}

	modules:add{
		name = 'asinh',
		code = [[
real asinh(real x) { 
	return log(x + sqrt(x*x + 1.));
}
]],
	}


--[[
TODO don't just use 'math', 
put a global depends list somewhere, build on it as we add eqn and solver
but app has its own modules which it needs typecode for first (which other solvers use ffi instances of)
and solvers modules depend on those app modules...
--]]	
	modules:add{
		name = 'math',
		depends = {
			'realparam',
			'units',
			'real',
			'real3',	-- used by coord
			
			-- used by coord, but not needed by {flux=roe, eqn=euler, coord=cartesian}
			--'sym3',		-- used by bssnok-fd 
			--'real3x3',-- used by calcCellMinMaxEigenvalues, which is only used by the calcDT of mhd and glm-mhd eqns
			--'_3sym3',	-- used by conn
			--'real3x3x3',
			--'quat',
		
			-- not needed
			--'sym3sym3',
			--'cplx',
			--'cplx3',
			--'cplx3x3',
			--'rotate',
			--'sech',
			--'Bessel',
			--'getPerpendicularBasis',
			--'normalForSide',
		},
	}
end
