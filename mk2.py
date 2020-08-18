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
exit_diameter = 0.147
throat_diameter = 0.0766
diverging_half_angle = 15
converging_half_angle = 30
chamber_length = 0.305
chamber_diameter = 0.146
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.02

"""
Double-wall cooling parameters
"""
inner_wall_thickness = 0.0015875  # m
channel_height = 0.001  # m
wall_thermal_conductivity = 385  # W / (m K)
surface_roughness = 0.001  # mm
# TODO incorrect, this should be divided by two as there were 2 engines
fuel_flow_rate = 2.16  # kg/s
lox_flow_rate = 4.95  # kg/s

# fuel_input_pressure = 1820000  # Pa
fuel_input_temperature = 300  # K

"""
Heat flux calc info
For CEARUN Results use Frozen comp (DOLPRE page 70)
"""
chamber_pressure = 1.29e6  # Pa
chamber_temperature = 3365.31  # K
gas_gamma = 1.1366  # unitless
gas_viscosity = 1.0363e-4  # Pa S
gas_specific_heat = 2.0794e3  # J / (kg K)
gas_prandtl_number = 0.6209  # unitless
cstar = 1779.2  # m/s

# Thermal resistance of carbon deposit (DOLPRE 87)
# Use low value for all to be conservative, should be experimentally determined
# 3.4e-4 corresponds to ~1000 in ^2 s F / BTU
#carbon_resistance = 7e-4  # m^2 K / W
carbon_resistance = 0e-4  # m^2 K / W

"""
Stress calcs
"""
coolant_pressure = 1.74e6  # Pa

elastic_modulus = 117e9  # Pa
thermal_expansion_coeff = 17e-6  # m / m / K
poissons_ratio = 0.33  # unitless


