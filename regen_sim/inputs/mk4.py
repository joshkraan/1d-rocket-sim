# MK4 engine inputs

"""
General engine parameters
"""

# Chamber dimensions. All units are meters and degrees.
exit_diameter = 0.1447
throat_diameter = 0.0934
diverging_half_angle = 15
converging_half_angle = 30
chamber_length = 0.3315
chamber_diameter = 0.1615
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.02

fuel_flow_rate = 0.961  # kg/s
lox_flow_rate = 2.20  # kg/s

# Design chamber pressure (Feynman)
# Used as initial estimate
chamber_pressure = 770000  # Pa

"""
Regenerative cooling parameters
"""
inner_wall_thickness = 0.0015875  # m
channel_height = 0.001  # m
wall_thermal_conductivity = 50  # W / (m K)
surface_roughness = 0.001  # mm
fuel_input_temperature = 300  # K

# Name of file containing discrete data for fuel properties, which will be interpolated.
# Pressure should be constant and temperature varied.
# A CSV file with columns titled Temperature, Density, Specific Heat, Thermal Conductivity, Viscosity, Enthalpy.
# Standard SI units with no prefixes should be used.
fuel_properties_file = "Kerosene_11E5Pa.csv"

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
