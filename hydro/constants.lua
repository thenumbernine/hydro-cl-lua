-- TODO use symmath.vars for units?
return {
	speedOfLight_in_m_per_s = 299792458,
	gravitationalConstant_in_m3_per_kg_s2 = 6.6740831e-11,
	CoulombConstant_in_kg_m3_per_C2_s2 = 8.9875517873681764e+9,
	BoltzmannConstant_in_kg_m2_per_K_s2 = 1.3806485279e-23,

	-- TODO calculate from eps0 = 1 / (4 * pi * k_e), mu0 = 1 / (eps0 * c^2)
	vacuumPermeability_in_kg_m_per_C2 = 1.2566370621219e-6,		-- 2018 CODATA
	vacuumPermittivity_in_C2_s2_per_kg_m3 = 8.854187812813e-12,	-- 2018 CODATA

	protonMass_in_kg = 1.6726219236951e-27,						-- 2018 CODATA
	protonCharge_in_C = 1.602176634e-19,						-- 2018 CODATA

	electronMass_in_kg = 9.109383701528e-31,					-- 2018 CODATA
	electronCharge_in_C = -1.602176634e-19,						-- 2018 CODATA

	EarthRadius_in_m = 6.371e+6,
	EarthMass_in_kg = 5.9722e+24,

	SolarRadius_in_m = 6.960e+8, 
	SolarMass_in_kg = 1.9891e+30,
}
