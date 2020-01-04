"""
High level inputs
"""
num_stations = 1000
guess_for_exit_wall_temp = 400

"""
Chamber dimensions. All units are meters and degrees.
"""
exit_diameter = 0.211
throat_diameter = 0.0944
diverging_half_angle = 20
converging_half_angle = 30
chamber_length = 0.33
chamber_diameter = 0.172
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.01

"""
Double-wall cooling parameters
"""
inner_wall_thickness = 0.0015875  # m
channel_height = 0.003  # m
wall_thermal_conductivity = 401  # W / (m K)
surface_roughness = 0.001  # mm
fuel_input_pressure = 1820000  # Pa
fuel_input_temperature = 300  # k

"""
CEARUN Results
"""
chamber_pressure = 16.5  # Bar
chamber_temperature = 3111.90  # kelvin
gamma = 1.123  # unitless
viscosity = 0.97784  # millipoise
specific_heat = 9.1653 # kJ / (kg K)
prandtl_number = 0.4133  # unitless

"""
Miscellaneous inputs
"""
fuel_flow_rate = 2.11  # kg/s
ox_flow_rate = 4.82  # kg/s
input_enthalpy = -408900 #J/s
