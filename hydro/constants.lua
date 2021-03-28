-- TODO use symmath.vars for units?
local units = {}

units.speedOfLight_in_m_per_s = 299792458
units.gravitationalConstant_in_m3_per_kg_s2 = 6.6740831e-11
units.CoulombConstant_in_kg_m3_per_C2_s2 = 8.9875517873681764e+9
units.BoltzmannConstant_in_kg_m2_per_K_s2 = 1.3806485279e-23

	-- TODO calculate from eps0 = 1 / (4 * pi * k_e), mu0 = 1 / (eps0 * c^2)
units.vacuumPermeability_in_kg_m_per_C2 = 1.2566370621219e-6		-- 2018 CODATA
units.vacuumPermittivity_in_C2_s2_per_kg_m3 = 8.854187812813e-12	-- 2018 CODATA

units.protonMass_in_kg = 1.6726219236951e-27						-- 2018 CODATA
units.protonCharge_in_C = 1.602176634e-19						-- 2018 CODATA

units.electronMass_in_kg = 9.109383701528e-31					-- 2018 CODATA
units.electronCharge_in_C = -1.602176634e-19						-- 2018 CODATA

units.EarthRadius_in_m = 6.371e+6
units.EarthMass_in_kg = 5.9722e+24

units.SolarRadius_in_m = 6.960e+8
units.SolarMass_in_kg = 1.9891e+30


units.Mpc_in_m = 648000/math.pi*149597870700*1000000	--megaparsec
units.Kpc_in_m = 648000/math.pi*149597870700*1000	--kiloparsec
units.pc_in_m = 648000/math.pi*149597870700	--parsec
units.lyr_in_m = 9460730472580800	--light year
units.AU_in_m = 149597870700	--astronomical unit
units.ls_in_m = 299792458	--light second
units.km_in_m = 1000	--kilometer
units.m_in_m = 1		--meter

	--1 = G m^3 / (kg s^2) <=> G m^3 / (c m)^2 = kg <=> G/c^2 = kg/m
units.kg_in_m = units.gravitationalConstant_in_m3_per_kg_s2 / (units.speedOfLight_in_m_per_s * units.speedOfLight_in_m_per_s)

	--1 = c m/s  <-> c m = 1 s
units.m_in_s = units.speedOfLight_in_m_per_s

return units
