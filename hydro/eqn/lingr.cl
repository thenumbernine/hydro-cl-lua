//// MODULE_NAME: applyInitCond
//// MODULE_DEPENDS: cartesianToCoord

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;
	
	<?=initCode()?>

	UBuf[index] = (<?=eqn.cons_t?>){
<? if eqn.usePressure then
?>		.Pi_tt = P,
<? else		
?>		.Pi_tt = rho,
<? end		
?>		.Pi_ti = real3_zero,
		.Pi_ij = sym3_zero,
		
		.Psi_ttk = cartesianToCoord(v, x),
		.Psi_tik = real3x3_zero,
		.Psi_ijk = _3sym3_zero,
	};
}

//// MODULE_NAME: fluxFromCons
//// MODULE_DEPENDS: solver_t normal_t cons_t

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	real c = solver->wavespeed;
	real3 nL = normal_l1(n);
	real3 nU = normal_u1(n);
	
	<?=eqn.cons_t?> F;

<?
for _,fields in ipairs{
	{"Pi_tt", "Psi_ttk", {""}},
	{"Pi_ti", "Psi_tik", xNames:mapi(function(xi) return "."..xi end)},
	{"Pi_ij", "Psi_ijk", symNames:mapi(function(xij) return "."..xij end)},
} do
	local Pi, Psi, xabs = table.unpack(fields)
	for _,xab in ipairs(xabs) do
?>
	//F^Pi_ab = -c Psi_ab,k n^k
	F.<?=Pi?><?=xab?> = -c * (U.<?=Psi?>.x<?=xab?> * nU.x + U.<?=Psi?>.y<?=xab?> * nU.y + U.<?=Psi?>.z<?=xab?> * nU.z);
	
	//F^{Psi_ab,k} = -c (Pi_ab n_k)
<? 		for k,xk in ipairs(xNames) do
?>	F.<?=Psi?>.<?=xk?><?=xab?> = -c * (U.<?=Pi?><?=xab?> * nL.<?=xk?>);
<? 		end
	end
end
?>

	return F;
}

//// MODULE_NAME: eigen_forInterface
//// MODULE_DEPENDS: eigen_t

#define eigen_forInterface(solver, UL, UR, x, n) ((eigen_t){})

//// MODULE_NAME: eigen_forCell
//// MODULE_DEPENDS: eigen_t

#define eigen_forCell(solver, U, x, n) ((eigen_t){})

//// MODULE_NAME: eigen_left/rightTransform
//// MODULE_DEPENDS: eigen_t waves_t

waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x,
	normal_t n
) { 
	waves_t Y;
	real* Yp = Y.ptr;

	real nLen = normal_len(n);
	real invDenom = 1. / (nLen * nLen);

<?
for _,fields in ipairs{
	{"Pi_tt", "Psi_ttk", {""}},
	{"Pi_ti", "Psi_tik", xNames:mapi(function(xi) return "."..xi end)},
	{"Pi_ij", "Psi_ijk", symNames:mapi(function(xij) return "."..xij end)},
} do
	local Pi, Psi, xabs = table.unpack(fields)
	for _,xab in ipairs(xabs) do
?>
	Yp[0] = .5 * invDenom * (
		X.<?=Pi?><?=xab?> * nLen
		+ X.<?=Psi?>.x<?=xab?> * normal_u1x(n)
		+ X.<?=Psi?>.y<?=xab?> * normal_u1y(n)
		+ X.<?=Psi?>.z<?=xab?> * normal_u1z(n)
	);
	Yp[1] = invDenom * (
		X.<?=Psi?>.x<?=xab?> * normal_u2x(n)
		+ X.<?=Psi?>.y<?=xab?> * normal_u2y(n)
		+ X.<?=Psi?>.z<?=xab?> * normal_u2z(n)
	);
	Yp[2] = invDenom * (
		X.<?=Psi?>.x<?=xab?> * normal_u3x(n)
		+ X.<?=Psi?>.y<?=xab?> * normal_u3y(n)
		+ X.<?=Psi?>.z<?=xab?> * normal_u3z(n)
	);
	Yp[3] = .5 * invDenom * (
		-X.<?=Pi?><?=xab?> * nLen
		+ X.<?=Psi?>.x<?=xab?> * normal_u1x(n)
		+ X.<?=Psi?>.y<?=xab?> * normal_u1y(n)
		+ X.<?=Psi?>.z<?=xab?> * normal_u1z(n)
	);
	Yp += 4;
<?
	end
end
?>
	return Y;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x,
	normal_t n
) {
	real* Xp = X.ptr;
	real nLen = normal_len(n);

	cons_t Y;
<?
for _,fields in ipairs{
	{"Pi_tt", "Psi_ttk", {""}},
	{"Pi_ti", "Psi_tik", xNames:mapi(function(xi) return "."..xi end)},
	{"Pi_ij", "Psi_ijk", symNames:mapi(function(xij) return "."..xij end)},
} do
	local Pi, Psi, xabs = table.unpack(fields)
	for _,xab in ipairs(xabs) do
?>
	Y.<?=Pi?><?=xab?> = (Xp[0] - Xp[3]) * nLen;
	Y.<?=Psi?>.x<?=xab?> = Xp[0] * normal_l1x(n) 
			+ Xp[1] * normal_l2x(n) 
			+ Xp[2] * normal_l3x(n) 
			+ Xp[3] * normal_l1x(n);
	Y.<?=Psi?>.y<?=xab?> = Xp[0] * normal_l1y(n) 
			+ Xp[1] * normal_l2y(n) 
			+ Xp[2] * normal_l3y(n) 
			+ Xp[3] * normal_l1y(n);
	Y.<?=Psi?>.z<?=xab?> = Xp[0] * normal_l1z(n) 
			+ Xp[1] * normal_l2z(n) 
			+ Xp[2] * normal_l3z(n) 
			+ Xp[3] * normal_l1z(n);
	Xp += 4;
<?
	end
end
?>
	return Y;
}

//// MODULE_NAME: eigen_fluxTransform

// by default in hydro/eqn/eqn.lua, fluxFromCons is defined by eigen_fluxTransform
// but since eig is empty, we can define eigen_fluxTransform with fluxFromCons
#define eigen_fluxTransform(solver, eig, X, x, n) fluxFromCons(solver, X, x, n)

//// MODULE_NAME: addSource

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cellBuf[index].pos;
	
	// TODO add stress-energy influence here
	// TODO TODO coupled solvers, once again, but it does seem messy
}
