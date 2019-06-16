typedef <?=app.realparam?> realparam;

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
?> 

<? 
local function makevec3(name, real) 
?>
typedef union {
	<?=real?> s[3];
	struct { <?=real?> s0, s1, s2; };
	struct { <?=real?> x, y, z; };
} <?=name?>;
<? 
end 
makevec3('real3', 'real')
makevec3('cplx3', 'cplx')
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
typedef union {
	real s[9];
	real3 v[3];
	struct {
		real3 v0,v1,v2;
	};
	struct {
		real3 x,y,z;
	};
} real3x3;

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
