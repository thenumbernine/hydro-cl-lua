typedef union {
	real s[3];
	struct { real s0, s1, s2; };
	struct { real x, y, z; };
} real3;

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
		sym3 v0,v1,v2;
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
} mat3;
