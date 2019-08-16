import math

import inputs

"""
Chamber dimension processing for easier use
All units are in meters and radians
"""
exit_radius = inputs.exit_diameter / 2
throat_radius = inputs.throat_diameter / 2
diverging_half_angle = math.radians(inputs.diverging_half_angle)
converging_half_angle = math.radians(inputs.converging_half_angle)
chamber_length = inputs.chamber_length
chamber_radius = inputs.chamber_diameter / 2
throat_bevel_radius = inputs.throat_bevel_radius
chamber_bevel_radius = inputs.chamber_bevel_radius

"""
Cooling properties
"""
inner_wall_thickness = inputs.inner_wall_thickness
channel_thickness = inputs.channel_thickness
wall_thermal_conductivity = inputs.wall_thermal_conductivity

"""
NASA CEARUN unit conversions & c-star
"""
chamber_pressure = inputs.chamber_pressure * 100000  # pascal
chamber_temperature = inputs.chamber_temperature  # kelvin
gamma = inputs.gamma  # unitless
viscosity = inputs.viscosity / 10000  # Pa s
specific_heat = inputs.specific_heat  # J / (g K)
prandtl_number = inputs.prandtl_number  # unitless

characteristic_velocity = chamber_pressure * math.pi * throat_radius ** 2 / inputs.mass_flow_rate  # m/s

"""
Bounds separating sections of chamber
"""
bounds = [chamber_length]
bounds.append(bounds[0] + chamber_bevel_radius * math.sin(converging_half_angle))
bounds.append(bounds[1] + (chamber_radius + chamber_bevel_radius * math.cos(converging_half_angle)
                           - chamber_bevel_radius - throat_radius - throat_bevel_radius + throat_bevel_radius
                           * math.cos(converging_half_angle)) / math.tan(converging_half_angle))
bounds.append(bounds[2] + throat_bevel_radius * (math.sin(diverging_half_angle) + math.sin(converging_half_angle)))
bounds.append(bounds[3] + (exit_radius - (throat_radius + throat_bevel_radius - throat_bevel_radius * math.cos(diverging_half_angle))) / math.tan(diverging_half_angle) )

"""
High level processing
"""
num_stations = inputs.num_stations
station_length = bounds[4] / num_stations

"""
Calculate throat position for convenience
"""
throat_position = chamber_length + chamber_bevel_radius * math.sin(converging_half_angle) + (chamber_radius + chamber_bevel_radius * math.cos(converging_half_angle) - chamber_bevel_radius - throat_radius - throat_bevel_radius + throat_bevel_radius * math.cos(converging_half_angle)) / math.tan(converging_half_angle) + throat_bevel_radius * math.sin(converging_half_angle)
