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
# TODO update fuel properties to correct pressure
"""
Chamber dimensions. All units are meters and degrees.
"""
exit_diameter = 0.0981366
throat_diameter = 0.05356
diverging_half_angle = 15
converging_half_angle = 30
chamber_length = 0.23679
chamber_diameter = 0.10940
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.02

"""
Double-wall cooling parameters
"""
inner_wall_thickness = 0.0015875  # m
channel_height = 0.0001  # m
wall_thermal_conductivity = 385  # W / (m K)
surface_roughness = 0.001  # mm
fuel_flow_rate = 0.7707  # kg/s
lox_flow_rate = 1.76498  # kg/s

# fuel_input_pressure = 1820000  # Pa
fuel_input_temperature = 300  # K

"""
Heat flux calc info
For CEARUN Results use Frozen comp (DOLPRE page 70)
"""
# 2.29 OF
chamber_pressure = 1877777.8  # Pa
chamber_temperature = 3411.57  # K
gas_gamma = 1.139  # unitless
gas_viscosity = 1.0469e-4  # Pa S
gas_specific_heat = 2082.2  # J / (kg K)
gas_prandtl_number = 0.6223  # unitless
cstar = 1786.0  # m/s

# 1.8 OF
# chamber_pressure = 1877777.8  # Pa
# chamber_temperature = 3014.6  # K
# gas_gamma = 1.1888  # unitless
# gas_viscosity = 0.91417e-4  # Pa S
# gas_specific_heat = 2180.8  # J / (kg K)
# gas_prandtl_number = 0.5689  # unitless
# cstar = 1747.8  # m/s

# Thermal resistance of carbon deposit (DOLPRE 87)
# Use low value for all to be conservative, should be experimentally determined
# 3.4e-4 corresponds to ~1000 in ^2 s F / BTU
#carbon_resistance = 7e-4  # m^2 K / W
carbon_resistance = 0  # m^2 K / W

"""
Stress calcs
"""
coolant_pressure = 2.54e6  # Pa

elastic_modulus = 117e9  # Pa
thermal_expansion_coeff = 17e-6  # m / m / K
poissons_ratio = 0.33  # unitless


