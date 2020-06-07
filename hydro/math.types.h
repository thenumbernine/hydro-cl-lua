typedef <?=app.realparam?> realparam;
typedef <?=app.realparam?>2 realparam2;
typedef <?=app.realparam?>4 realparam4;
typedef <?=app.realparam?>8 realparam8;

<? 
local function makevec3(name, scalar)
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
?><?=vecType.typeCode?>
typedef <?=vecType.type?> <?=name?>;
<?
	else 
		-- normal
?>
typedef union {
	<?=scalar?> s[3];
	struct { <?=scalar?> s0, s1, s2; };
	struct { <?=scalar?> x, y, z; };
} <?=app.real=='half' and '__attribute__ ((packed))' or ''
-- __attribute__ ((packed)) seems to need to be here with real=half
?> <?=name?>;
<? 
	end
end 
makevec3('real3', 'real')
?> 

typedef union {
	real s[6];
	struct {
		real xx, xy, xz, yy, yz, zz;
	};
	struct {
		real s00, s01, s02, s11, s12, s22;
	};
} sym3;

typedef union {
	real s[18];
	struct {
		sym3 v0,v1,v2;	//why not s0,s1,s2?
	};
	struct {
		sym3 x,y,z;
	};
} _3sym3;

//row vectors, so a.i.j = a_ij
<?
local function make3x3(scalar)
	local vec3 = scalar..'3'
	local name = scalar..'3x3'
?>
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
<?
end
make3x3'real'	-- real3x3
?>

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


// cplx

<? 
local function makecplx(name, real) 
?>
typedef union {
	<?=real?> s[2];
	struct { <?=real?> s0, s1; };
	struct { <?=real?> re, im; };
} <?=name?>;
<? 
end 
makecplx('cplx', 'real')
makevec3('cplx3', 'cplx')
make3x3'cplx'
?>


// used by the display code

typedef union {
	real	ptr[9];
	real	vreal;
	sym3	vsym3;
	cplx	vcplx;
	real3	vreal3;
	cplx3	vcplx3;
	real3x3	vreal3x3;
} displayValue_t;
