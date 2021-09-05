# Feynman output sheet 12

"""
General engine parameters
"""

# Burn time 36 seconds
# Total mass flow 12 kg/s
# Chamber pressure 900 kPa
# Cooling channel pressure drop 450 kPa
# Injector pressure drop 180 kPa

# Chamber dimensions. All units are meters and degrees.
exit_diameter = 0.28197
throat_diameter = 0.1685
diverging_half_angle = 15
converging_half_angle = 30
chamber_length = 0.3434
chamber_diameter = 0.27788
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.02

fuel_flow_rate = 3.871  # kg/s
lox_flow_rate = 8.129  # kg/s

# Design chamber pressure (Feynman)
# Used as initial estimate
chamber_pressure = 900e3  # Pa

"""
Regenerative cooling parameters
"""
helical_channels = False
helix_pitch = 7.6e-2  # m
helix_number = 1

inner_wall_thickness = 0.0015875  # m
channel_height = 1e-3  # m
wall_thermal_conductivity = 50  # W / (m K)
surface_roughness = 7.5e-6  # m
fuel_input_temperature = 300  # K

# Name of file containing discrete data for fuel properties, which will be interpolated.
# Pressure should be constant and temperature varied.
# A CSV file with columns titled Temperature, Density, Specific Heat, Thermal Conductivity, Viscosity, Enthalpy.
# Standard SI units with no prefixes should be used.
fuel_properties_name = "Kerosene_12E5Pa.csv"

"""
Inputs for stress calculations
"""
coolant_pressure = 1.12e6  # Pa
elastic_modulus = 200e9  # Pa
thermal_expansion_coeff = 13e-6  # m / m / K
poissons_ratio = 0.29  # unitless

"""
These inputs are only used for the CEARUN option
Read manual, GAMMAs is not gas_gamma
"""
chamber_temperature = 3304.4  # K
gas_gamma = 1.1755  # unitless
gas_viscosity = 1.0224e-4  # Pa S
gas_specific_heat = 2.0755e3  # J / (kg K)
gas_prandtl_number = 0.6190  # unitless