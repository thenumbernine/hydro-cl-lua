//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
kernel void calcDT(
	global real* dtBuf,
	const global cons_t* UBuf,
	const global real* ePotBuf
) {
	int4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int index = INDEXV(i);
	if (OOB(2,2)) {
		dtBuf[index] = INFINITY;
		return;
	}
	
	cons_t U = UBuf[index];
	real ePot = ePotBuf[index];
	prim_t W = primFromCons(U, ePot);
	real Cs = calc_Cs(W);

	real dt = INFINITY;

	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		real lambdaMin = min((real)0., W.v.s[side] - Cs);
		real lambdaMax = max((real)0., W.v.s[side] + Cs);
		dt = min(dt, dx<?=side?>_at(i) / (fabs(lambdaMax - lambdaMin) + (real)1e-9));
	}<? end ?>
	dtBuf[index] = dt;
}

//used by PLM
<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	eigen_t* eig,
	const global cons_t* U
) {
	real ePot = 0; //TODO need ePot...
	prim_t W = primFromCons(*U, ePot);
	real vSq = coordLenSq(W.v);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U->ETotal);
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	eig->rho = W.rho;
	eig->v = W.v;
	eig->hTotal = hTotal;
	eig->vSq = vSq;
	eig->Cs = Cs;
}
<? end ?>

//used by PLM
<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr0?> real* wave,
	<?=addr1?> const eigen_t* eig
) {
	real v_n = eig->v.s[<?=side?>];
	wave[0] = v_n - eig->Cs;
	wave[1] = v_n;
	wave[2] = v_n;
	wave[3] = v_n;
	wave[4] = v_n + eig->Cs;
}
<?		end
	end
end ?>

//used for interface eigen basis
void eigen_forSide(
	global eigen_t* eig,
	global const cons_t* UL, 
	real ePotL, 
	global const cons_t* UR, 
	real ePotR
) {
	prim_t WL = primFromCons(*UL, ePotL);
	real sqrtRhoL = sqrt(WL.rho);
	real3 vL = WL.v;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL->ETotal);
	
	prim_t WR = primFromCons(*UR, ePotR);
	real sqrtRhoR = sqrt(UR->rho);
	real3 vR = WR.v;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR->ETotal);
	
	real invDenom = 1./(sqrtRhoL + sqrtRhoR);
	
	//Roe-averaged
	real rho = sqrtRhoL * sqrtRhoR;
	real3 v = real3_add(
			real3_scale(vL, sqrtRhoL * invDenom),
			real3_scale(vR, sqrtRhoR * invDenom));
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);
	
	//derived:
	real vSq = coordLenSq(v);
	real eKin = .5 * vSq;
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	
	eig->rho = rho; 
	eig->v = v;
	eig->hTotal = hTotal;
	eig->vSq = vSq;
	eig->Cs = Cs;
}

kernel void calcEigenBasis(
	global real* waveBuf,			//[volume][dim][numWaves]
	global eigen_t* eigenBuf,		//[volume][dim]
	const global consLR_t* ULRBuf,	//[volume]
	const global real* ePotBuf
) {
	SETBOUNDS(2,1);
	int indexR = index;
	real ePotR = ePotBuf[indexR];
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize[side];
		real ePotL = ePotBuf[indexL];
		
		global const cons_t* UL = &ULRBuf[side + dim * indexL].R;
		global const cons_t* UR = &ULRBuf[side + dim * indexR].L;
		
		int intindex = side + dim * index;	

		global eigen_t* eig = eigenBuf + intindex;
		eigen_forSide(eig, UL, ePotL, UR, ePotR);
		
		global real* wave = waveBuf + numWaves * intindex;
		eigen_calcWaves_<?=side?>_global_global(wave, eig);
	}<? end ?>
}

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,solver.dim-1 do 
				local prefix
				if side == 0 then
					prefix = [[
	const real nx = 1, ny = 0, nz = 0;
	const real n1x = 0, n1y = 1, n1z = 0;
	const real n2x = 0, n2y = 0, n2z = 1;
	real v_n = v.x, v_n1 = v.y, v_n2 = v.z;
]] 
				elseif side == 1 then
					prefix = [[
	const real nx = 0, ny = 1, nz = 0;
	const real n1x = 0, n1y = 0, n1z = 1;
	const real n2x = 1, n2y = 0, n2z = 0;
	real v_n = v.y, v_n1 = v.z, v_n2 = v.x;
]] 
				elseif side == 2 then
					prefix = [[
	const real nx = 0, ny = 0, nz = 1;
	const real n1x = 1, n1y = 0, n1z = 0;
	const real n2x = 0, n2y = 1, n2z = 0;
	real v_n = v.z, v_n1 = v.x, v_n2 = v.y;
]]
				end
				prefix = [[
	real3 v = eig->v;
	real hTotal = eig->hTotal;
	real vSq = eig->vSq;
	real Cs = eig->Cs;
]] .. prefix
?>

void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) { 
	<?=prefix?>

	real invDenom = .5 / (Cs * Cs);
	y[0] = (x[0] * ((heatCapacityRatio - 1.) * .5 * vSq + Cs * v_n)
		+ x[1] * -(nx * Cs + (heatCapacityRatio - 1.) * v.x) 
		+ x[2] * -(ny * Cs + (heatCapacityRatio - 1.) * v.y)
		+ x[3] * -(nz * Cs + (heatCapacityRatio - 1.) * v.z)
		+ x[4] * (heatCapacityRatio - 1.)
	) * invDenom;
	y[1] = (x[0] * (2.*Cs*Cs - (heatCapacityRatio - 1.) * vSq)
		+ x[1] * (heatCapacityRatio - 1.) * v.x * 2
		+ x[2] * (heatCapacityRatio - 1.) * v.y * 2
		+ x[3] * (heatCapacityRatio - 1.) * v.z * 2
		+ x[4] * -(heatCapacityRatio - 1.) * 2
	) * invDenom;
	y[2] = x[0] * -v_n1 + x[1] * n1x + x[2] * n1y + x[3] * n1z;
	y[3] = x[0] * -v_n2 + x[1] * n2x + x[2] * n2y + x[3] * n2z;
	y[4] = (x[0] * ((heatCapacityRatio - 1.) * .5 * vSq - Cs * v_n) 
		+ x[1] * (nx * Cs - (heatCapacityRatio - 1.) * v.x) 
		+ x[2] * (ny * Cs - (heatCapacityRatio - 1.) * v.y) 
		+ x[3] * (nz * Cs - (heatCapacityRatio - 1.) * v.z) 
		+ x[4] * (heatCapacityRatio - 1.)
	) * invDenom;
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	<?=prefix?>

	y[0] = x[0] + x[1] + x[4];
	y[1] = x[0] * (v.x - nx * Cs) + x[1] * v.x + x[2] * n1x + x[3] * n2x + x[4] * (v.x + nx * Cs);
	y[2] = x[0] * (v.y - ny * Cs) + x[1] * v.y + x[2] * n1y + x[3] * n2y + x[4] * (v.y + ny * Cs);
	y[3] = x[0] * (v.z - nz * Cs) + x[1] * v.z + x[2] * n1z + x[3] * n2z + x[4] * (v.z + nz * Cs);
	y[4] = x[0] * (hTotal - v_n * Cs) + x[1] * .5 * vSq + x[2] * v_n1 + x[3] * v_n2 + x[4] * (hTotal + v_n * Cs);
}

<?	if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	<?=prefix?>
	y[0] = x[1] * nx + x[2] * ny + x[3] * nz;
	y[1] = x[0] * (-v_n * v.x + (heatCapacityRatio - 1.) * .5 * vSq * nx)
		+ x[1] * (v.x * nx - (heatCapacityRatio - 1.) * nx * v.x + v_n)
		+ x[2] * (v.x * ny - (heatCapacityRatio - 1.) * nx * v.y)
		+ x[3] * (v.x * nz - (heatCapacityRatio - 1.) * nx * v.z)
		+ x[4] * (heatCapacityRatio - 1.) * nx;
	y[2] = x[0] * (-v_n * v.y + (heatCapacityRatio - 1.) * .5 * vSq * ny)
		+ x[1] * (v.y * nx - (heatCapacityRatio - 1.) * ny * v.x)
		+ x[2] * (v.y * ny - (heatCapacityRatio - 1.) * ny * v.y + v_n)
		+ x[3] * (v.y * nz - (heatCapacityRatio - 1.) * ny * v.z)
		+ x[4] * (heatCapacityRatio - 1.) * ny;
	y[3] = x[0] * (-v_n * v.z + (heatCapacityRatio - 1.) * .5 * vSq * nz)
		+ x[1] * (v.z * nx - (heatCapacityRatio - 1.) * nz * v.x)
		+ x[2] * (v.z * ny - (heatCapacityRatio - 1.) * nz * v.y)
		+ x[3] * (v.z * nz - (heatCapacityRatio - 1.) * nz * v.z + v_n)
		+ x[4] * (heatCapacityRatio - 1.) * nz;
	y[4] = x[0] * v_n * ((heatCapacityRatio - 1.) * .5 * vSq - hTotal)
		+ x[1] * (-(heatCapacityRatio - 1.) * v_n * v.x + nx * hTotal)
		+ x[2] * (-(heatCapacityRatio - 1.) * v_n * v.y + ny * hTotal)
		+ x[3] * (-(heatCapacityRatio - 1.) * v_n * v.z + nz * hTotal)
		+ x[4] * heatCapacityRatio * v_n;
}
<?
				end
			end
		end
	end
end
?>
