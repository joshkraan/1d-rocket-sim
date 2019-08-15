"""
High level inputs
"""
num_stations = 1000

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
Double-wall cooling parameters (meters)
"""
inner_wall_thickness = 0.01
channel_thickness = 0.01
wall_thermal_conductivity = 401  # W / (m K)

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
mass_flow_rate = 8.12  # kg/s
