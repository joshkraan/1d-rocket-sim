"""
High level inputs
"""
num_stations = 1000
guess_for_exit_heat_flux = 5
guess_for_exit_wall_temp = 500

"""
Chamber dimensions. All units are meters and degrees.
"""
exit_diameter = 0.214
throat_diameter = 0.114
diverging_half_angle = 15
converging_half_angle = 30
chamber_length = 0.354
chamber_diameter = 0.201
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.05

"""
Double-wall cooling parameters
"""
inner_wall_thickness = 0.0016  # m
channel_percentage = 1  # channel height as a % of chamber diameter
channel_height = 0.006350  # m
wall_thermal_conductivity = 401  # W / (m K)
surface_roughness = 0.001  # mm
fuel_input_pressure = 1820000  # Pa
fuel_input_temperature = 300  # k

"""
CEARUN Results
"""
chamber_pressure = 13.2  # Bar
chamber_temperature = 3368.15  # kelvin
gamma = 1.13  # unitless
viscosity = 1.037  # millipoise
specific_heat = 6.5462  # kJ / (kg K)
prandtl_number = 0.4364  # unitless

"""
Miscellaneous inputs
"""
fuel_flow_rate = 2.56  # kg/s
ox_flow_rate = 5.65  # kg/s
input_enthalpy = -408900 #J/s
