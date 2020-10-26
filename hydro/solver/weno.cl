//// MODULE_NAME: calcFlux
//// MODULE_DEPENDS: cell_x fluxFromCons normal_t eigen_forInterface eigen_left/rightTransform eqn.waveCode

/*
WENO solver:
sources:
1996 Jiang, Shu, "Efficient Implementation of Weighted ENO Schemes"
1998 Shu "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws"
"A hybrid approach for the regularized long wave-Burgers equation"
2016 Rathan, Raju "An improved Non-linear Weights for Seventh-Order WENO Scheme"
https://github.com/jzrake/Mara for weno5 examples
https://github.com/wme7/WENO7-Z/blob/master/WENO7ZresAdv1d.m for weno7 examples
https://github.com/python-hydro/hydro_examples/blob/master/compressible/weno_coefficients.py likewise

TODO incorporate parallel propagators
*/

static inline real sqr(real x) { return x * x; }

<? 
local clnumber = require 'cl.obj.number'

local stencilSize = solver.stencilSize

local coeffs = {
	-- [r] of table 1 of 1996 Jiang, Shu
	[2] = {
		d = {1/3, 2/3},
		-- coeffs[r].a[k][l] == a^r_k,l of table 1 of 1996 Jiang, Shu
		a = {
			{-1/2,	3/2,	},
			{1/2,	1/2,	},
		},
		betaCoeffs = {
			{
				{1},
				{-2, 1},
			}, {
				{1},
				{-2, 1},
			},
		},
	},
	
	[3] = {
-- [[		
		d = {1/10, 6/10, 3/10},
		a = {
			{ 2/6,	-7/6, 11/6, },
			{-1/6,	 5/6,  2/6, },
			{ 2/6,	 5/6, -1/6, },
		},
--]]
--[[ which paper is this from?
		d = {1/16,	10/16,	5/16},
		a = {
			{3/8,	-10/8,	15/8},
			{-1/8,	6/8,	3/8},
			{3/8,	6/8,	-1/8},
		},
--]]	
		betaCoeffs = {
			{
				{4/3},
				{-19/3,	25/3},
				{11/3,	-31/3,	10/3},
			}, {
				{4/3},
				{-13/3,	13/3},
				{5/3,	-13/3,	4/3},
			}, {
				{10/3},
				{-31/3,	25/3},
				{11/3,	-19/3,	4/3},
			},
		},
	},
	

	
	[4] = {
		d = {1/35, 12/35, 18/35, 4/35},
		a = {	
			{-3/12, 13/12, -23/12, 25/12, },
			{ 1/12, -5/12,  13/12,  3/12, },
			{-1/12,  7/12,   7/12, -1/12, },
			{ 3/12, 13/12,  -5/12,  1/12, },
		},
		betaCoeffs = {

--[[
			{
				{ 2107/240},
				{-9402/240,  11003/240},
				{ 7042/240, -17246/240,   7043/240},
				{-1854/240,   4642/240,  -3882/240,    547/240},
			},
			{
				{  547/240},
				{-2522/240,   3443/240},
				{ 1922/240,  -5966/240,   2843/240},
				{ -494/240,   1602/240,  -1642/240,    267/240},
			},
			{
				{  267/240},
				{-1642/240,   2843/240},
				{ 1602/240,  -5966/240,   3443/240},
				{ -494/240,   1922/240,  -2522/240,    547/240},
			},
			{
				{  547/240},
				{-3882/240,   7043/240},
				{ 4642/240, -17246/240,  11003/240},
				{-1854/240,   7042/240,  -9402/240,   2107/240},
			},
--]]		
-- [[ 2016 Rathan, Raju after eqn 12 ... 
-- https://github.com/wme7/WENO7-Z/blob/master/WENO7ZresAdv1d.m
			{
				{  6649},
				{-47214,	  85641},
				{ 56694,	-210282,	134241},
				{-22778,	  86214,	-114894,	25729},
			},
			{
				{  3169},
				{-19374,	 33441},
				{ 19014,	-70602,	 41001},
				{ -5978,	 23094,	-30414,	6649},
			},
			{
				{  6649},
				{-30414,	 41001},
				{ 23094,	-70602,	 33441},
				{ -5978,	 19014,	-19374,	3169},
			},
			{
				{  25729},
				{-114894,	134241},
				{  86214,	-210282,	85641},
				{ -22778,	  56694,	-47214,	6649},
			},
--]]	
		},
	},
	[5] = {
		d = {  1/126,  20/126,  60/126,  40/126,   5/126 },
		a = {
			{12/60, -63/60, 137/60, -163/60, 137/60, },
			{-3/60,  17/60, -43/60,   77/60,  12/60, },
			{ 2/60, -13/60,  47/60,   27/60,  -3/60, },
			{-3/60,  27/60,  47/60,  -13/60,   2/60, },
			{12/60,  77/60, -43/60,   17/60,  -3/60, },
		},
		betaCoeffs = {
			{
				{  107918/5040},
				{ -649501/5040,  1020563/5040},
				{  758823/5040, -2462076/5040,  1521393/5040},
				{ -411487/5040,  1358458/5040, -1704396/5040,   482963/5040},
				{   86329/5040,  -288007/5040,   364863/5040,  -208501/5040,    22658/5040},
			},                                                                                 
			{                                                                                         
				{   22658/5040},
				{ -140251/5040,   242723/5040},
				{  165153/5040,  -611976/5040,   406293/5040},
				{  -88297/5040,   337018/5040,  -464976/5040,   138563/5040},
				{   18079/5040,   -70237/5040,    99213/5040,   -60871/5040,     6908/5040},
			},                                                                                 
			{                                                                                         
				{    6908/5040},
				{  -51001/5040,   104963/5040},
				{   67923/5040,  -299076/5040,   231153/5040},
				{  -38947/5040,   179098/5040,  -299076/5040,   104963/5040},
				{    8209/5040,   -38947/5040,    67923/5040,   -51001/5040,     6908/5040},
			},                                                                                 
			{                                                                                         
				{    6908/5040},
				{  -60871/5040,   138563/5040},
				{   99213/5040,  -464976/5040,   406293/5040},
				{  -70237/5040,   337018/5040,  -611976/5040,   242723/5040},
				{   18079/5040,   -88297/5040,   165153/5040,  -140251/5040,    22658/5040},
			},                                                                                 
			{                                                                                         
				{   22658/5040},
				{ -208501/5040,   482963/5040},
				{  364863/5040, -1704396/5040,  1521393/5040},
				{ -288007/5040,  1358458/5040, -2462076/5040,  1020563/5040},
				{   86329/5040,  -411487/5040,   758823/5040,  -649501/5040,   107918/5040},
			},
		},
	},

	[6] = {
		d = { 1/462,   30/462,  150/462,  200/462,   75/462,    6/462 },
		a = {
			{-10/60,  62/60, -163/60, 237/60, -213/60, 147/60, },
			{  2/60, -13/60,   37/60, -63/60,   87/60,  10/60, },
			{ -1/60,   7/60,  -23/60,  57/60,   22/60,  -2/60, },
			{  1/60,  -8/60,   37/60,  37/60,   -8/60,   1/60, },
			{ -2/60,  22/60,   57/60, -23/60,    7/60,  -1/60, },
			{ 10/60,  87/60,  -63/60,  37/60,  -13/60,   2/60, },
		},
		betaCoeffs = {
			{
				{   6150211/120960,},
				{ -47460464/120960,   94851237/120960,},
				{  76206736/120960, -311771244/120960,  260445372/120960,},
				{ -63394124/120960,  262901672/120960, -444003904/120960,  190757572/120960,},
				{  27060170/120960, -113206788/120960,  192596472/120960, -166461044/120960,   36480687/120960,},
				{  -4712740/120960,   19834350/120960,  -33918804/120960,   29442256/120960,  -12950184/120960,    1152561/120960,},
			}, {
				{   1152561/120960,},
				{  -9117992/120960,   19365967/120960,},
				{  14742480/120960,  -65224244/120960,   56662212/120960,},
				{ -12183636/120960,   55053752/120960,  -97838784/120960,   43093692/120960,},
				{   5134574/120960,  -23510468/120960,   42405032/120960,  -37913324/120960,    8449957/120960,},
				{   -880548/120960,    4067018/120960,   -7408908/120960,    6694608/120960,   -3015728/120960,     271779/120960,},
			}, {
				{    271779/120960,},
				{  -2380800/120960,    5653317/120960,},
				{   4086352/120960,  -20427884/120960,   19510972/120960,},
				{  -3462252/120960,   17905032/120960,  -35817664/120960,   17195652/120960,},
				{   1458762/120960,   -7727988/120960,   15929912/120960,  -15880404/120960,    3824847/120960,},
				{   -245620/120960,    1325006/120960,   -2792660/120960,    2863984/120960,   -1429976/120960,     139633/120960,},
			}, {
				{    139633/120960,},
				{  -1429976/120960,    3824847/120960,},
				{   2863984/120960,  -15880404/120960,   17195652/120960,},
				{  -2792660/120960,   15929912/120960,  -35817664/120960,   19510972/120960,},
				{   1325006/120960,   -7727988/120960,   17905032/120960,  -20427884/120960,    5653317/120960,},
				{   -245620/120960,    1458762/120960,   -3462252/120960,    4086352/120960,   -2380800/120960,     271779/120960,},
			}, {
				{    271779/120960,},
				{  -3015728/120960,    8449957/120960,},
				{   6694608/120960,  -37913324/120960,   43093692/120960,},
				{  -7408908/120960,   42405032/120960,  -97838784/120960,   56662212/120960,},
				{   4067018/120960,  -23510468/120960,   55053752/120960,  -65224244/120960,   19365967/120960,},
				{   -880548/120960,    5134574/120960,  -12183636/120960,   14742480/120960,   -9117992/120960,    1152561/120960,},
			}, {
				{   1152561/120960,},
				{ -12950184/120960,   36480687/120960,},
				{  29442256/120960, -166461044/120960,  190757572/120960,},
				{ -33918804/120960,  192596472/120960, -444003904/120960,  260445372/120960,},
				{  19834350/120960, -113206788/120960,  262901672/120960, -311771244/120960,   94851237/120960,},
				{  -4712740/120960,   27060170/120960,  -63394124/120960,   76206736/120960,  -47460464/120960,    6150211/120960,},
			},
		},
	},

	[7] = {
		d = {  1/1716,   42/1716,  315/1716,  700/1716,  525/1716,  126/1716,    7/1716 },
		a = {
			{  60/420, -430/420, 1334/420, -2341/420, 2559/420, -1851/420, 1089/420, },
			{ -10/420,   74/420, -241/420,   459/420, -591/420,   669/420,   60/420, },
			{   4/420,  -31/420,  109/420,  -241/420,  459/420,   130/420,  -10/420, },
			{  -3/420,   25/420, -101/420,   319/420,  214/420,   -38/420,    4/420, },
			{   4/420,  -38/420,  214/420,   319/420, -101/420,    25/420,   -3/420, },
			{ -10/420,  130/420,  459/420,  -241/420,  109/420,   -31/420,    4/420, },
			{  60/420,  669/420, -591/420,   459/420, -241/420,    74/420,  -10/420, },
		},
		betaCoeffs = {
			{
				{    7177657304/59875200},
				{  -68289277071/59875200,   166930543737/59875200},
				{  140425750893/59875200,  -698497961463/59875200,   739478564460/59875200},
				{ -158581758572/59875200,   797280592452/59875200, -1701893556420/59875200,   985137198380/59875200},
				{  102951716988/59875200,  -521329653333/59875200,  1119254208255/59875200, -1301580166020/59875200,   431418789360/59875200},
				{  -36253275645/59875200,   184521097818/59875200,  -397822832973/59875200,   464200620612/59875200,  -308564463663/59875200,    55294430841/59875200},
				{    5391528799/59875200,   -27545885877/59875200,    59577262788/59875200,   -69700128812/59875200,    46430779053/59875200,   -16670007831/59875200,     1258225940/59875200},
			}, {
				{    1258225940/59875200},
				{  -12223634361/59875200,    31090026771/59875200},
				{   25299603603/59875200,  -132164397513/59875200,   143344579860/59875200},
				{  -28498553012/59875200,   151212114012/59875200,  -332861569020/59875200,   195601143380/59875200},
				{   18375686988/59875200,   -98508059523/59875200,   219064013505/59875200,  -259838403420/59875200,    86959466460/59875200},
				{   -6414710427/59875200,    34632585198/59875200,   -77574968883/59875200,    92646554652/59875200,   -62392325913/59875200,    11250068787/59875200},
				{     945155329/59875200,    -5128661355/59875200,    11548158588/59875200,   -13862429972/59875200,     9380155443/59875200,    -3397272201/59875200,      257447084/59875200},
			}, {
				{     257447084/59875200},
				{   -2659103847/59875200,     7257045753/59875200},
				{    5684116173/59875200,   -32164185663/59875200,    36922302360/59875200},
				{   -6473137292/59875200,    37531128132/59875200,   -88597133220/59875200,    54531707180/59875200},
				{    4158865908/59875200,   -24530177853/59875200,    59045150655/59875200,   -74236325220/59875200,    25788772260/59875200},
				{   -1432622085/59875200,     8555779674/59875200,   -20891234853/59875200,    26694456132/59875200,   -18869146983/59875200,     3510366201/59875200},
				{     206986975/59875200,    -1247531949/59875200,     3078682188/59875200,    -3982402892/59875200,     2854088973/59875200,    -1077964287/59875200,       84070496/59875200},
			}, {
				{      84070496/59875200},
				{    -969999969/59875200,     2927992563/59875200},
				{    2283428883/59875200,   -14296379553/59875200,    18133963560/59875200},
				{   -2806252532/59875200,    18083339772/59875200,   -47431870620/59875200,    32154783380/59875200},
				{    1902531828/59875200,   -12546315963/59875200,    33820678305/59875200,   -47431870620/59875200,    18133963560/59875200},
				{    -676871859/59875200,     4550242446/59875200,   -12546315963/59875200,    18083339772/59875200,   -14296379553/59875200,     2927992563/59875200},
				{      99022657/59875200,     -676871859/59875200,     1902531828/59875200,    -2806252532/59875200,     2283428883/59875200,     -969999969/59875200,       84070496/59875200},
			}, {
				{      84070496/59875200},
				{   -1077964287/59875200,     3510366201/59875200},
				{    2854088973/59875200,   -18869146983/59875200,    25788772260/59875200},
				{   -3982402892/59875200,    26694456132/59875200,   -74236325220/59875200,    54531707180/59875200},
				{    3078682188/59875200,   -20891234853/59875200,    59045150655/59875200,   -88597133220/59875200,    36922302360/59875200},
				{   -1247531949/59875200,     8555779674/59875200,   -24530177853/59875200,    37531128132/59875200,   -32164185663/59875200,     7257045753/59875200},
				{     206986975/59875200,    -1432622085/59875200,     4158865908/59875200,    -6473137292/59875200,     5684116173/59875200,    -2659103847/59875200,      257447084/59875200},
			}, {
				{     257447084/59875200},
				{   -3397272201/59875200,    11250068787/59875200},
				{    9380155443/59875200,   -62392325913/59875200,    86959466460/59875200},
				{  -13862429972/59875200,    92646554652/59875200,  -259838403420/59875200,   195601143380/59875200},
				{   11548158588/59875200,   -77574968883/59875200,   219064013505/59875200,  -332861569020/59875200,   143344579860/59875200},
				{   -5128661355/59875200,    34632585198/59875200,   -98508059523/59875200,   151212114012/59875200,  -132164397513/59875200,    31090026771/59875200},
				{     945155329/59875200,    -6414710427/59875200,    18375686988/59875200,   -28498553012/59875200,    25299603603/59875200,   -12223634361/59875200,     1258225940/59875200},
			}, {
				{    1258225940/59875200},
				{  -16670007831/59875200,    55294430841/59875200},
				{   46430779053/59875200,  -308564463663/59875200,   431418789360/59875200},
				{  -69700128812/59875200,   464200620612/59875200, -1301580166020/59875200,   985137198380/59875200},
				{   59577262788/59875200,  -397822832973/59875200,  1119254208255/59875200, -1701893556420/59875200,   739478564460/59875200},
				{  -27545885877/59875200,   184521097818/59875200,  -521329653333/59875200,   797280592452/59875200,  -698497961463/59875200,   166930543737/59875200},
				{    5391528799/59875200,   -36253275645/59875200,   102951716988/59875200,  -158581758572/59875200,   140425750893/59875200,   -68289277071/59875200,     7177657304/59875200},
			},
		},
	},
}

local coeff = coeffs[stencilSize]
local a = coeff.a
local d = coeff.d
local betaCoeffs = coeff.betaCoeffs

for side=0,solver.dim-1 do
	for _,l_or_r in ipairs{'l', 'r'} do
		local ak0 = l_or_r == 'l' and stencilSize or 1
		local akd = l_or_r == 'l' and -1 or 1
		local al0 = l_or_r == 'l' and stencilSize or 1
		local ald = l_or_r == 'l' and -1 or 1
		local d0 = l_or_r == 'l' and stencilSize or 1
		local dd = l_or_r == 'l' and -1 or 1

-- how should I generalize this?
-- base pointer
-- number of cells (stencilSize)
-- number of coefficients at each cell ... or what C type it is, and how many variables to cycle across in the ptr field
-- ... or we can just always perform this interpolation in cons/char space
?>
<?=eqn.waves_t?> weno_<?=l_or_r?>_<?=side?>(
	const <?=eqn.waves_t?>* v
) {
	<?=eqn.waves_t?> result;
	for (int k = 0; k < numWaves; ++k) {
<? 		for j=0,stencilSize-1 do 
?>		real beta<?=j?> = 0.<?
			for m=0,stencilSize-1 do
				for n=0,m do
			?> + <?=clnumber(betaCoeffs[j+1][m+1][n+1])?> * v[<?=j+m?>].ptr[k] * v[<?=j+n?>].ptr[k]<?
				end
			end?>;
<?		end
	
		if solver.wenoMethod == '1996 Jiang Shu' then 		-- WENO-JS
			--local epsilon = clnumber(1e-14)
			local epsilon = clnumber(1e-6)
?>
<? 			for i=0,stencilSize-1 do
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> / sqr(<?=epsilon?> + beta<?=i?>);
<? 			end 

		elseif solver.wenoMethod == '2008 Borges' then 	-- WENO-Z
			local epsilon = clnumber(1e-6)
			
			-- TODO different coeffs of betas?
			if false then 	--if stencilSize == 4 then -- for 2016 Rathan, it suggests these:
?>		real tau = fabs(beta0 + 3. * beta1 - 3. * beta2 - beta3);
<?			else	-- for weno7, 2018 Zeytinoglu suggests this:
?>		real tau = fabs(beta0 - beta<?=stencilSize-1?>);
<?			end
 			for i=0,stencilSize-1 do 
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> * (1. + (tau / (beta<?=i?> + <?=epsilon?>)));
<? 			end

		elseif solver.wenoMethod == '2010 Shen Zha' then -- WENO-BS?
			local epsilon = clnumber(1e-10)
		
			-- TODO find these for other order WENO's
			local shen_zha_A = clnumber(50)	-- 0-100
?>		
		real minB = beta0, maxB = beta1;
<?			for j=1,stencilSize-1 do
?>		minB = min(minB, beta<?=j?>); maxB = max(maxB, beta<?=j?>);
<?			end
?>		real R0 = minB / (maxB + <?=epsilon?>);
<? 			for i=0,stencilSize-1 do 
?>		beta<?=i?> += R0 * <?=shen_zha_A?> * minB;
<? 			end 
 			for i=0,stencilSize-1 do 
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> / sqr(<?=epsilon?> + beta<?=i?>);
<?	 		end 
 
		else
			error("unknown wenoMethod "..tostring(solver.wenoMethod))
		end

		for i=0,stencilSize-1 do 
?>		real vs<?=i?> = 0.<?
			for j=0,stencilSize-1 do
		?> + <?=clnumber(a[ak0+akd*i][al0+ald*j])?> * v[<?=i+j?>].ptr[k]<?
			end ?>;
<? 
		end 
?>		real alphasum= 0.<? for i=0,stencilSize-1 do ?> + alpha<?=i?><? end ?>;
		result.ptr[k] = (0.<?
		for i=0,stencilSize-1 do
		?> + alpha<?=i?> * vs<?=i?><?
		end ?>) / alphasum;
	}
	return result;
}
<? 	end 
end
?>

kernel void calcCellFlux(
	constant <?=solver.solver_t?>* solver,
	const global <?=solver.coord.cell_t?>* cellBuf,
	global <?=eqn.cons_t?>* fluxCellBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		real3 xInt = cell_x(i);
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		normal_t n = normal_forSide<?=side?>(xInt);
		global <?=eqn.cons_t?>* F = fluxCellBuf + <?=side?> + dim * index;
		*F = fluxFromCons(solver, *U, xInt, n);
	}<? end ?>
}

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	const global <?=solver.coord.cell_t?>* cellBuf,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf
//	,const global <?=eqn.cons_t?>* fluxCellBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?
for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - solver->stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;
	
		int indexR = index;
		const global <?=eqn.cons_t?>* UR = UBuf + indexR;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		normal_t n = normal_forSide<?=side?>(xInt);

		int fluxIndexInt = side + dim * index;
		global <?=eqn.cons_t?>* flux = fluxBuf + fluxIndexInt;

<?
	
	if solver.fluxMethod == 'Lax-Friedrichs' then

?>
		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);
		real maxAbsLambda = 0.;
		for (int j = 0; j < <?=2*stencilSize?>; ++j) {
			const global <?=eqn.cons_t?>* Uj = U + (j - <?=stencilSize?>) * solver->stepsize.s<?=side?>;
			
			<?=eqn:consWaveCodePrefix('n', '*Uj', 'xInt'):gsub('\n', '\n\t\t')?>
			
			real lambdaMin = <?=eqn:consMinWaveCode('n', '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMin));
			
			real lambdaMax = <?=eqn:consMaxWaveCode('n', '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMax));
		}

<? 	
		for j=0,2*stencilSize-1 do
?>		const global <?=eqn.cons_t?>* U<?=j?> = U + <?=j - stencilSize?> * solver->stepsize.s<?=side?>;
		//<?=eqn.cons_t?> F<?=j?> = fluxCellBuf[<?=side?> + dim * (index + <?=j - stencilSize?> * solver->stepsize.s<?=side?>)];
		<?=eqn.cons_t?> F<?=j?> = fluxFromCons(solver, *U<?=j?>, xInt, n);
<?
			if j < 2*stencilSize-1 then
?>		<?=eqn.cons_t?> fp<?=j?>;
<?				for k=0,eqn.numStates-1 do 
?>		fp<?=j?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] + maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
<?				end
			end
			
			if j > 0 then
?>		<?=eqn.cons_t?> fm<?=j-1?>;
<?				for k=0,eqn.numStates-1 do 
?>		fm<?=j-1?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] - maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
<? 				end
			end
		end	
?>	
		
		<?=eqn.waves_t?> afp[<?=2*stencilSize-1?>], afm[<?=2*stencilSize-1?>];
<? 		for j=0,2*stencilSize-2 do
?>		afp[<?=j?>] = eigen_leftTransform(solver, eig, fp<?=j?>, xInt, n);
		afm[<?=j?>] = eigen_leftTransform(solver, eig, fm<?=j?>, xInt, n);
<? 		end 
?>
		<?=eqn.waves_t?> waf;
		<?=eqn.waves_t?> wafp = weno_r_<?=side?>(afp);
		<?=eqn.waves_t?> wafm = weno_l_<?=side?>(afm);
		for (int j = 0; j < numWaves; ++j) {
			waf.ptr[j] = wafp.ptr[j] + wafm.ptr[j];
		}
		
		*flux = eigen_rightTransform(solver, eig, waf, xInt, n);
<?
	
	elseif solver.fluxMethod == 'Marquina' then

?>
//// MODULE_DEPENDS: eigen_forCell
		<?=eqn.eigen_t?> eigL = eigen_forCell(solver, *UL, xInt, n);
		<?=eqn.eigen_t?> eigR = eigen_forCell(solver, *UR, xInt, n);

		real lambdaL[<?=eqn.numWaves?>];
		real lambdaR[<?=eqn.numWaves?>];
<? 		for _,lr in ipairs{'L', 'R'} do
?>		{
		<?=eqn:eigenWaveCodePrefix('n', 'eig'..lr, 'xInt'):gsub('\n', '\n\t\t')?>
<? 		
			for k=0,eqn.numWaves-1 do 
?>			lambda<?=lr?>[<?=k?>] = <?=eqn:eigenWaveCode('n', 'U'..lr, 'xInt', k)?>;
<?			end
?>		}
<?		end
?>

<? 		for j=0,2*stencilSize-1 do
?>		const global <?=eqn.cons_t?>* U<?=j?> = U + (<?=j - stencilSize?>) * solver->stepsize.s<?=side?>;
		<?=eqn.cons_t?> F<?=j?> = fluxFromCons(solver, *U<?=j?>, xInt, n);
		<?=eqn.waves_t?> al<?=j?> = eigen_leftTransform(solver, eigL, *U<?=j?>, xInt, n);
		<?=eqn.waves_t?> ar<?=j?> = eigen_leftTransform(solver, eigR, *U<?=j?>, xInt, n);
		<?=eqn.waves_t?> afl<?=j?> = eigen_leftTransform(solver, eigL, F<?=j?>, xInt, n);
		<?=eqn.waves_t?> afr<?=j?> = eigen_leftTransform(solver, eigR, F<?=j?>, xInt, n);

<?		end
?>

		<?=eqn.waves_t?> afp[<?=2*stencilSize-1?>], afm[<?=2*stencilSize-1?>];
<? 
		for k=0,eqn.numWaves-1 do 
?>		if (lambdaL[<?=k?>] > 0.0 && lambdaR[<?=k?>] > 0.0) {
<?			for j=0,2*stencilSize-1 do
?>			afp[<?=j?>].ptr[<?=k?>] = afl<?=j?>.ptr[<?=k?>];
			afm[<?=j?>].ptr[<?=k?>] = 0;
<?			end
?>		} else if (lambdaL[<?=k?>] < 0.0 && lambdaR[<?=k?>] < 0.0) {
<?			for j=0,2*stencilSize-1 do
?>			afp[<?=j?>].ptr[<?=k?>] = 0;
			afm[<?=j?>].ptr[<?=k?>] = afr<?=j?>.ptr[<?=k?>];
<?			end
?>		} else {
			real absLambdaL = fabs(lambdaL[<?=k?>]);
			real absLambdaR = fabs(lambdaR[<?=k?>]);
			real a, absa;
			if (absLambdaL > absLambdaR) {
				a = lambdaL[<?=k?>];
				absa = absLambdaL;
			} else {
				a = lambdaR[<?=k?>];
				absa = absLambdaR;
			}
<?			for j=0,2*stencilSize-1 do
?>			afp[<?=j?>].ptr[<?=k?>] = .5 * (afl<?=j?>.ptr[<?=k?>] + absa * al<?=j?>.ptr[<?=k?>]);
			afm[<?=j?>].ptr[<?=k?>] = .5 * (afr<?=j?>.ptr[<?=k?>] - absa * ar<?=j?>.ptr[<?=k?>]);
<?			end
?>		}
<?		end
?>	

		<?=eqn.waves_t?> wafp = weno_r_<?=side?>(afp);
		<?=eqn.waves_t?> wafm = weno_l_<?=side?>(afm);
		
		<?=eqn.cons_t?> fluxP = eigen_rightTransform(solver, eigL, wafp, xInt, n);
		<?=eqn.cons_t?> fluxM = eigen_rightTransform(solver, eigR, wafm, xInt, n);
		
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] = fluxP.ptr[j] + fluxM.ptr[j];
		}
<?
	
	elseif solver.fluxMethod == 'Roe' then
		
		-- TODO make the following its own function, maybe in solver/roe.cl or solver/roe.lua
		-- so we can call a function instead of copy/paste
		-- but this means making it a function of the lua parameters, like flux limiter.
		-- and a flux limiter means +1 to the numGhost

?>		<?=eqn.waves_t?> afp[<?=2*stencilSize-1?>], afm[<?=2*stencilSize-1?>];
<?		for j=0,2*stencilSize-1 do
?>
		{
			const global <?=eqn.cons_t?>* UL<?=j?> = U + <?=j - stencilSize?> * solver->stepsize.s<?=side?>;
			const global <?=eqn.cons_t?>* UR<?=j?> = U + <?=j+1 - stencilSize?> * solver->stepsize.s<?=side?>;
			<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);

			<?=eqn:eigenWaveCodePrefix('n', 'eig', 'xInt')?>

			<?=eqn.cons_t?> UAvg;
			for (int k = 0; k < numIntStates; ++k) {
				UAvg.ptr[k] = .5 * (UL->ptr[k] + UR->ptr[k]);
			}
			afp[<?=j?>] = eigen_leftTransform(solver, eig, UAvg, xInt, n);
			afm[<?=j?>] = afp[<?=j?>];

			<?=eqn.cons_t?> deltaU;
			for (int k = 0; k < numStates; ++k) {
				deltaU.ptr[k] = UR->ptr[k] - UL->ptr[k];
			}
			
			<?=eqn.waves_t?> deltaUEig = eigen_leftTransform(solver, eig, deltaU, xInt, n);

			<? for k=0,eqn.numWaves-1 do ?>{
				const int k = <?=k?>;
				real lambda = <?=eqn:eigenWaveCode('n', 'eig', 'xInt', k)?>;
				real lambdaPlus = max(lambda, 0.);
				real lambdaMinus = min(lambda, 0.);

				afp[<?=j?>].ptr[k] *= lambdaPlus;
				afm[<?=j?>].ptr[k] *= lambdaMinus;
				real sgnLambda = lambda >= 0 ? 1 : -1;
				afp[<?=j?>].ptr[k] -= .5 * lambdaPlus * deltaUEig.ptr[k] * sgnLambda;
				afm[<?=j?>].ptr[k] -= .5 * lambdaMinus * deltaUEig.ptr[k] * sgnLambda;
			}<? end ?>
		}
<?
		end 
?>
		<?=eqn.waves_t?> waf;
		<?=eqn.waves_t?> wafp = weno_r_<?=side?>(afp);
		<?=eqn.waves_t?> wafm = weno_l_<?=side?>(afm);
		for (int j = 0; j < numWaves; ++j) {
			waf.ptr[j] = wafp.ptr[j] + wafm.ptr[j];
		}
	
		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, n);
		*flux = eigen_rightTransform(solver, eig, waf, xInt, n);
<?
	
	else
		error("unknown fluxMethod "..tostring(solver.fluxMethod))
	end
?>	}<? 
end ?>
}
