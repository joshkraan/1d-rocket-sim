import math
import numpy
import matplotlib.pyplot

from input_processing import (
    # Chamber dimensions
    throat_radius, diverging_half_angle, converging_half_angle,
    chamber_length, chamber_radius, throat_bevel_radius, chamber_bevel_radius,
    # Cooling
    inner_wall_thickness, channel_thickness, wall_thermal_conductivity,
    # Combustion gas properties
    gamma, viscosity, specific_heat, prandtl_number,
    # Engine properties
    chamber_pressure, chamber_temperature, characteristic_velocity,
    # Bounds
    bounds, throat_position,
    # High level
    num_stations, station_length
)


def inner_radius(x):
    # x=0 is defined at the injector plate
    if x <= bounds[0]:
        return chamber_radius
    elif bounds[0] < x <= bounds[1]:
        return math.sqrt(
            throat_bevel_radius ** 2 - (x - chamber_length) ** 2) + chamber_radius - chamber_bevel_radius
    elif bounds[1] < x <= bounds[2]:
        return - math.tan(converging_half_angle) * (x - chamber_length - chamber_bevel_radius * math.sin(converging_half_angle)) + chamber_radius + chamber_bevel_radius * math.cos(converging_half_angle) - chamber_bevel_radius
    elif bounds[2] < x <= bounds[3]:
        return - math.sqrt(throat_bevel_radius ** 2 - (x - bounds[2] - throat_bevel_radius * math.sin(converging_half_angle)) ** 2) + throat_radius + throat_bevel_radius
    elif bounds[3] < x:
        return math.tan(diverging_half_angle) * (x - bounds[3]) + throat_radius + throat_bevel_radius - throat_bevel_radius * math.cos(diverging_half_angle)


def area_ratio(x):
    local_area = math.pi * inner_radius(x) ** 2
    throat_area = math.pi * throat_radius**2
    return local_area / throat_area


def mach_number(x):
    if x < throat_position:
        mach_current = 0.3
        mach_last = 0
        while abs(mach_current - mach_last) > 0.00001:
            mach_last = mach_current
            mach_current = 1 / area_ratio(x) * ((2 + (gamma - 1) * mach_last ** 2) / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
        return mach_current
    elif x == throat_position:
        return 1
    elif x > throat_position:
        mach_current = 1.2
        mach_last = 0
        while abs(mach_current - mach_last) > 0.00001:
            mach_last = mach_current
            mach_current = math.sqrt(((gamma + 1) * (area_ratio(x) * mach_last) ** ((2 * (gamma - 1)) / (gamma + 1)) - 2) / (gamma - 1))
        return mach_current


def adiabatic_wall_temp(x):
    recovery_factor = prandtl_number ** 0.33
    local_mach = mach_number(x)

    a = 1 + recovery_factor * (gamma - 1) / 2 * local_mach ** 2
    b = 1 + (gamma - 1) / 2 * local_mach ** 2

    return chamber_temperature * (a / b)


def gas_transfer_coefficient(x, gas_wall_temp):
    local_mach = mach_number(x)
    a = 1 + (gamma - 1) / 2 * local_mach ** 2
    sigma = 1/((0.5 * (gas_wall_temp / chamber_temperature) * a + 0.5) ** 0.68 * a ** 0.12)

    return 0.026 / (throat_radius * 2) ** 0.2 * viscosity ** 0.2 * specific_heat / prandtl_number ** 0.6\
           * (chamber_pressure * 9.8 / characteristic_velocity) ** 0.8 * (throat_radius * 2 / throat_bevel_radius) ** 0.1 * (1 / area_ratio(x)) ** 0.9 * sigma


def coolant_transfer_coefficient(x, coolant_wall_temp, coolant_properties):
    hydraulic_diameter = channel_thickness
    coolant_wall_viscosity = 1
    re = coolant_density * coolant_velocity * hydraulic_diameter / coolant_viscosity
    pr = coolant_viscosity * coolant_capacity / coolant_conductivity
    nu = 0.027 * re ** 0.8 * pr ** (1/3) * (coolant_viscosity / coolant_wall_viscosity) ** 0.14
    return nu * coolant_conductivity / hydraulic_diameter


def calculate_heat_balance(x, coolant_properties):
    # Initial guesses for gas wall temp and coolant temp, well below actual
    gas_wall_temp = 300
    coolant_temp = 250

    step_size = 1
    while coolant_properties[1] - coolant_temp > 0.00001:
        q = gas_transfer_coefficient(x, gas_wall_temp) * (adiabatic_wall_temp(x) - gas_wall_temp)
        conduction_resistance = math.log((inner_radius(x) + inner_wall_thickness) / inner_radius(x)) / (2 * math.pi * station_length * wall_thermal_conductivity)
        coolant_wall_temp = gas_wall_temp - q * conduction_resistance
        convection_resistance = 1 / (coolant_transfer_coefficient(x, coolant_wall_temp, coolant_properties) * 2 * math.pi * station_length * (inner_radius(x) + inner_wall_thickness))
        coolant_temp = coolant_wall_temp - q * convection_resistance
        if coolant_temp > coolant_properties[1]:
            gas_wall_temp -= step_size
            step_size /= 10
        else:
            gas_wall_temp += step_size

def update_coolant_props():
    #TODO

# Initialize coolant properties here:

for i in range(1, num_stations + 1):
    # Starting from nozzle end
    position = bounds[4] - i * station_length
    calculate_heat_balance(position, #coolant properties)
    update_coolant_props()


# test = numpy.arange(0, 0.7, 0.001)
# matplotlib.pyplot.scatter(test, list(map(heat_flux, test)))
# matplotlib.pyplot.axvline(x=throat_position)
# matplotlib.pyplot.show()

