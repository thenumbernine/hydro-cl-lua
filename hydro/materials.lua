-- https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity
-- resistivity is in ohm-meters, conductivity is in siemens/meters, temperature coefficient is in 1/Kelvin
local materials = {
	['Carbon (graphene)'] = {resistivity=1.00e-8, conductivity=1.00e+8, tempCoeff=-0.0002},
	['Silver'] = {resistivity=1.59e-8, conductivity=6.30e+7, tempCoeff=0.0038},
	['Copper'] = {resistivity=1.68e-8, conductivity=5.96e+7, tempCoeff=0.003862},
	['Annealed copper'] = {resistivity=1.72e-8, conductivity=5.80e+7, tempCoeff=0.00393},
	['Gold'] = {resistivity=2.44e-8, conductivity=4.10e+7, tempCoeff=0.0034},
	['Aluminium'] = {resistivity=2.82e-8, conductivity=3.50e+7, tempCoeff=0.0039},
	['Calcium'] = {resistivity=3.36e-8, conductivity=2.98e+7, tempCoeff=0.0041},
	['Tungsten'] = {resistivity=5.60e-8, conductivity=1.79e+7, tempCoeff=0.0045},
	['Zinc'] = {resistivity=5.90e-8, conductivity=1.69e+7, tempCoeff=0.0037},
	['Nickel'] = {resistivity=6.99e-8, conductivity=1.43e+7, tempCoeff=0.006},
	['Lithium'] = {resistivity=9.28e-8, conductivity=1.08e+7, tempCoeff=0.006},
	['Iron'] = {resistivity=9.71e-8, conductivity=1.00e+7, tempCoeff=0.005},
	['Platinum'] = {resistivity=1.06e-7, conductivity=9.43e+6, tempCoeff=0.00392},
	['Tin'] = {resistivity=1.09e-7, conductivity=9.17e+6, tempCoeff=0.0045},
	['Carbon steel (1010)'] = {resistivity=1.43e-7, conductivity=6.99e+6},
	['Lead'] = {resistivity=2.20e-7, conductivity=4.55e+6, tempCoeff=0.0039},
	['Titanium'] = {resistivity=4.20e-7, conductivity=2.38e+6, tempCoeff=0.0038},
	['Grain oriented electrical steel'] = {resistivity=4.60e-7, conductivity=2.17e+6},
	['Manganin'] = {resistivity=4.82e-7, conductivity=2.07e+6, tempCoeff=0.000002},
	['Constantan'] = {resistivity=4.90e-7, conductivity=2.04e+6, tempCoeff=0.000008},
	['Stainless steel'] = {resistivity=6.90e-7, conductivity=1.45e+6, tempCoeff=0.00094},
	['Mercury'] = {resistivity=9.80e-7, conductivity=1.02e+6, tempCoeff=0.0009},
	['Nichrome'] = {resistivity=1.10e-6, conductivity=6.7e+5, tempCoeff=0.0004},
	['GaAs'] = {resistivity=1.00e-3, conductivity=1.00e-8},
	['Carbon (amorphous)'] = {resistivity=5.00e-4, conductivity=1.25e+3, tempCoeff=-0.0005},
	['Carbon (graphite)'] = {resistivity=2.50e-6, conductivity=2.00e+5},
	['PEDOT:PSS'] = {resistivity=2e-6, conductivity=1e+1},
	['Germanium'] = {resistivity=4.60e-1, conductivity=2.17, tempCoeff=-0.048},
	['Sea water'] = {resistivity=2.00e-1, conductivity=4.80},
	['Swimming pool water'] = {resistivity=3.33e-1, conductivity=0.25},
	['Drinking water'] = {resistivity=2.00e+1, conductivity=5.00e-4},
	['Silicon'] = {resistivity=6.40e+2, conductivity=1.56e-3, tempCoeff=-0.075},
	['Wood (damp)'] = {resistivity=1.00e+3, conductivity=1e-4},
	['Deionized water'] = {resistivity=1.80e+5, conductivity=5.50e-6},
	['Glass'] = {resistivity=1.00e+11, conductivity=1e-15},
	['Hard rubber'] = {resistivity=1.00e+13, conductivity=1e-14},
	['Wood (oven dry)'] = {resistivity=1.00e+14, conductivity=1e-16},
	['Sulfur'] = {resistivity=1.00e+15, conductivity=1e-16},
	['Air'] = {resistivity=1.30e+14, conductivity=3e-15},
	['Carbon (diamond)'] = {resistivity=1.00e+12, conductivity=1e-13},
	['Fused quartz'] = {resistivity=7.50e+17, conductivity=1.30e-18},
	['PET'] = {resistivity=1.00e+21, conductivity=1e-21},
	['Teflon'] = {resistivity=1.00e+23, conductivity=1e-25},
}
-- convert everything to meters
local ohm = 0.033356409519815	-- m^0
local S = 29.9792458	-- m^0
for _,material in pairs(materials) do
	material.resistivity = material.resistivity * ohm	-- ohm-m = m^1
	material.conductivity = material.conductivity * S	-- S/m = m^-1
end

-- https://www.engineeringtoolbox.com/air-properties-d_156.html
-- measurements at 0' C at 1 bar
-- C_p = specific heat at constant pressure, in J / (kg K) = m^2 / (K s^2)
-- C_v = specific heat at constant volume, in J / (kg K)
materials.Air.C_p = 1006
materials.Air.C_v = 717.1
materials.Air.heatCapacityRatio = 1006 / 717.1 -- ~ 1.4

-- sea level air pressure
-- = 101.325 kPa
-- = 101325 Pa = N / m^2 = kg / (m s^2)
materials.Air.seaLevelPressure = 101325

-- https://en.wikipedia.org/wiki/Density_of_air
-- sea level 0' C, density is 1.2754 kg / m^3
materials.Air.seaLevelDensity = 1.2754

-- sqrt(gamma P / rho) = 333.8445024974 m/s
materials.Air.speedOfSound = math.sqrt(materials.Air.heatCapacityRatio * materials.Air.seaLevelPressure / materials.Air.seaLevelDensity)

-- from https://en.wikipedia.org/wiki/Viscosity
-- at 25 C and 1 bar of pressure
materials.Air.shearViscosity = 18.5e-6	-- Pa s = kg / (m s)

-- from https://en.wikipedia.org/wiki/List_of_thermal_conductivities
materials.Air.heatConductivity = 0.0235	-- W / (m K) = kg m / (s^3 K)

-- from https://en.wikipedia.org/wiki/Viscosity
materials['Drinking water'].shearViscosity = 8.9e-4	-- Pa s = kg / (m s) ... at 25 C

-- https://github.com/Bowserinator/Periodic-Table-JSON
materials.Lead.seaLevelDensity = 11.34e+3	-- kg/m^3 ... as a solid.  what about as a gas / plasma? at higher temps?
materials.Lead.boilingPoint = 2022			-- K
materials.Mercury.seaLevelDensity = 11.34e+3-- kg/m^3
materials.Mercury.boilingPoint = 629.88		-- K

return materials
