typedef struct {
	real rho;
	real v;
	real P;
} prim_t;

typedef struct {
	real rho;
	real m;
	real ETotal;
} cons_t;

cons_t consFromPrim(prim_t prim) {
	cons_t cons;
	cons.rho = prim.rho;
	cons.m = prim.rho * prim.v;
	cons.ETotal = .5 * prim.rho * prim.v * prim.v + prim.P / (gamma - 1);
	return cons;
}

__kernel void initState(
	__global cons_t* us
) {
	BEGIN_KERNEL();
	real x = (real)(i.x + .5) / (real)SIZE_X * (xmax - xmin) + xmin;
	prim_t prim;
	prim.rho = x < 0 ? 1 : .125;
	prim.v = 0;
	prim.P = x < 0 ? 1 : .1;
	us[index] = consFromPrim(prim); 
}

__kernel void calcFlux(
	__global cons_t* fs,
	const __global cons_t* us
) {
	BEGIN_KERNEL();
}

__kernel void convertToTex(
	__write_only IMAGETYPE tex,
	const __global cons_t* us
) {
	BEGIN_KERNEL();
	real value = us[index].rho;
	write_imagef(tex, WRITEIMAGEARGS, (float4)(value, 0., 0., 0.));
}
