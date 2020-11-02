-- This file needs to be templated, so it must be processed and manually inserted
-- I could write out the templated hydro/code/math.h file into the cache-bin folder, and then add that to the CL compiler include folder, and include it from there ...

local file = require 'ext.file'
local template = require 'template'

local common = require 'hydro.common'	-- xNames, symNames
local xNames = common.xNames
local sym = common.sym
local symNames = common.symNames
local from6to3x3 = common.from6to3x3
local from3x3to6 = common.from3x3to6

return function(modules, app)
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

	local function makescalarheader(scalar)
		return template([[
#define <?=scalar?>_add3(a,b,c)		(<?=add?>(a, <?=add?>(b,c)))
#define <?=scalar?>_add4(a,b,c,d)	(<?=add?>(a, <?=scalar?>_add3(b,c,d)))
#define <?=scalar?>_add5(a,b,c,d,e) (<?=add?>(a, <?=scalar?>_add4(b,c,d,e)))

#define <?=scalar?>_mul3(a,b,c)		(<?=mul?>(<?=mul?>(a,b),c))
]], 	{
			scalar = scalar,
			add = scalar..'_add',
			mul = scalar..'_mul',
		})
	end

	modules:addFromMarkup(template(file['hydro/code/math.cl'], {
		app = app,
		
		xNames = xNames,
		symNames = symNames,
		sym = sym,
		from6to3x3 = from6to3x3,
		from3x3to6 = from3x3to6,
		
		make_3x3_3_mul = make_3x3_3_mul,
		make_3_3x3_mul = make_3_3x3_mul,
		make3x3type = make3x3type,
		makevec3code = makevec3code,
		make_dot = make_dot,
		makevec3type = makevec3type,
		makevec3header = makevec3header,
		makecplx = makecplx,
		makescalarheader = makescalarheader,
	}))
end
