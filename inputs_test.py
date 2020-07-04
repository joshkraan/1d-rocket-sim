"""
High level inputs
"""
num_stations = 1000

"""
Path to file containing discrete data for fuel properties, which will be interpolated.
Pressure should be constant and temperature varied.
Should be a CSV file with columns titled Temperature, Density, Specific Heat, Thermal Conductivity, Viscosity, Enthalpy.
Standard SI units with no prefixes should be used.
"""
fuel_properties_file = "kerosene_17E5Pa_fit.csv"

"""
Chamber dimensions. All units are meters and degrees.
"""
exit_diameter = 0.5
throat_diameter = 0.28448
diverging_half_angle = 15
converging_half_angle = 30
chamber_length = 0.34272
chamber_diameter = 0.45
throat_bevel_radius = 0.133858
chamber_bevel_radius = 0.02

"""
Double-wall cooling parameters
"""
inner_wall_thickness = 0.0015875  # m
channel_height = 0.001  # m
wall_thermal_conductivity = 50  # W / (m K)
surface_roughness = 0.001  # mm
fuel_flow_rate = 3.8  # kg/s
lox_flow_rate = 8.7  # kg/s

# fuel_input_pressure = 1820000  # Pa
fuel_input_temperature = 300  # K

"""
CEARUN Results
"""
chamber_pressure = 5.516e6  # Pa
chamber_temperature = 3188.9  # K
gas_gamma = 1.213  # unitless
gas_viscosity = 5.215e-5  # Pa S
gas_specific_heat = 3948.15  # J / (kg K)
gas_prandtl_number = 0.82  # unitless

"""
Stress calcs
"""
coolant_pressure = 1.3e6  # Pa
#chamber_pressure = 9.5e5  # Pa
elastic_modulus = 200e9  # Pa
thermal_expansion_coeff = 13e-6  # m / m / K
poissons_ratio = 0.29  # unitless



