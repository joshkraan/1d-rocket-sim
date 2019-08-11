import math
import numpy
import matplotlib.pyplot

from input_processing import (
    # Chamber dimensions
    exit_radius, throat_radius, diverging_half_angle, converging_half_angle,
    chamber_length, chamber_radius, throat_bevel_radius, chamber_bevel_radius,
    # Combustion gas properties
    gamma, viscosity, specific_heat, prandtl_number,
    # Engine properties
    chamber_pressure, chamber_temperature, characteristic_velocity,
    # Bounds
    bounds, throat_position
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


def gas_wall_temp():
    # Need to fix this
    return 300


def gas_transfer_coefficient(x):
    local_mach = mach_number(x)
    a = 1 + (gamma - 1) / 2 * local_mach ** 2
    sigma = 1/((0.5 * (gas_wall_temp() / chamber_temperature) * a + 0.5) ** 0.68 * a ** 0.12)

    return 0.026 / (throat_radius * 2) ** 0.2 * viscosity ** 0.2 * specific_heat / prandtl_number ** 0.6\
           * (chamber_pressure * 9.8 / characteristic_velocity) ** 0.8 * (throat_radius * 2 / throat_bevel_radius) ** 0.1 * (1 / area_ratio(x)) ** 0.9 * sigma


def heat_flux(x):
    return gas_transfer_coefficient(x) * (adiabatic_wall_temp(x) - gas_wall_temp())

test = numpy.arange(0, 0.7, 0.001)
matplotlib.pyplot.scatter(test, list(map(heat_flux, test)))
matplotlib.pyplot.axvline(x=throat_position)
matplotlib.pyplot.show()

