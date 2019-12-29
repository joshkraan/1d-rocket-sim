import math
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from scipy.optimize import fixed_point
import matplotlib.pyplot as plt

from input_processing import (
    # Chamber dimensions
    throat_radius, diverging_half_angle, converging_half_angle,
    chamber_length, chamber_radius, throat_bevel_radius, chamber_bevel_radius,
    # Cooling
    inner_wall_thickness, channel_percent, wall_thermal_conductivity,
    coolant_flow_rate, surface_roughness, fuel_input_pressure, fuel_input_temperature,
    # Combustion gas properties
    gamma, viscosity, specific_heat, prandtl_number,
    # Engine properties
    chamber_pressure, chamber_temperature, characteristic_velocity,
    # Bounds
    bounds, throat_position,
    # High level
    num_stations, station_length
)

FLUID = "Water"


def test(station, coolant_properties, data):

    # Calculate the radius of the chamber at the current position
    # x = 0 is defined at the injector plate

    x = bounds[4] - station * station_length

    def inner_radius1(position):
        if x <= bounds[0]:
            return chamber_radius
        elif bounds[0] < x <= bounds[1]:
            return math.sqrt(throat_bevel_radius ** 2 - (x - chamber_length) ** 2) + chamber_radius - chamber_bevel_radius
        elif bounds[1] < x <= bounds[2]:
            return - math.tan(converging_half_angle) * (x - chamber_length - chamber_bevel_radius * math.sin(converging_half_angle)) + chamber_radius + chamber_bevel_radius * math.cos(converging_half_angle) - chamber_bevel_radius
        elif bounds[2] < x <= bounds[3]:
            return - math.sqrt(throat_bevel_radius ** 2 - (x - bounds[2] - throat_bevel_radius * math.sin(converging_half_angle)) ** 2) + throat_radius + throat_bevel_radius
        elif bounds[3] < x:
            return math.tan(diverging_half_angle) * (x - bounds[3]) + throat_radius + throat_bevel_radius - throat_bevel_radius * math.cos(diverging_half_angle)

    # Calculate the area ratio

    local_area = math.pi * inner_radius1(x) ** 2
    throat_area = math.pi * throat_radius ** 2
    area_ratio1 = local_area / throat_area

    # Calculate the mach number
    # Solved iteratively

    if x < throat_position:
        mach_current = 0.3
        mach_last = 0
        while abs(mach_current - mach_last) > 0.00001:
            mach_last = mach_current
            mach_current = 1 / area_ratio(x) * ((2 + (gamma - 1) * mach_last ** 2) / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
        mach_number1 = mach_current
    elif x == throat_position:
        mach_number1 = 1
    elif x > throat_position:
        mach_current = 1.2
        mach_last = 0
        while abs(mach_current - mach_last) > 0.00001:
            mach_last = mach_current
            mach_current = math.sqrt(((gamma + 1) * (area_ratio(x) * mach_last) ** ((2 * (gamma - 1)) / (gamma + 1)) - 2) / (gamma - 1))
        mach_number1 = mach_current

    # Calculate the adiabatic wall temperature

    recovery_factor = prandtl_number ** 0.33
    a = 1 + recovery_factor * (gamma - 1) / 2 * mach_number1 ** 2
    b = 1 + (gamma - 1) / 2 * mach_number1 ** 2
    adiabatic_wall_temp1 = chamber_temperature * (a / b)

    # Calculate the gas transfer coefficient

    def gas_transfer_coefficient1(gas_wall_temp):
        g = 1 + (gamma - 1) / 2 * mach_number1 ** 2
        sigma = 1 / ((0.5 * (gas_wall_temp / chamber_temperature) * g + 0.5) ** 0.68 * g ** 0.12)

        h_0 = 0.026 * prandtl_number ** 0.4 * (4 / math.pi) ** 0.8 * 1.555 * viscosity ** -0.8
        h = h_0 * 8.12 ** 0.8 * (inner_radius(x) * 2) ** -1.8

        bartz = 0.026 / (throat_radius * 2) ** 0.2 * viscosity ** 0.2 * specific_heat / prandtl_number ** 0.6 \
                * (chamber_pressure / characteristic_velocity) ** 0.8 * (
                            throat_radius * 2 / throat_bevel_radius) ** 0.1 * (1 / area_ratio(x)) ** 0.9 * sigma

        # print("Bartz = ", bartz)
        # print("Dittus = ", h)

        return bartz

    # Calculate the channel height
    # Currently 1% of the local chamber diameter

    def channel_height1(position):
        return 2 * inner_radius1(position) * channel_percent / 100

    # Calculate the coolant transfer coefficient

    def coolant_transfer_coefficient1(wall_temp):
        hydraulic_diameter = channel_height1(x)
        coolant_wall_viscosity = PropsSI('V', 'T', wall_temp, 'P', coolant_properties[station][1], FLUID)
        re = coolant_properties[station][3] * coolant_properties[station][4] * hydraulic_diameter / coolant_properties[station][5]
        pr = coolant_properties[station][5] * coolant_properties[station][6] / coolant_properties[station][7]
        nu = 0.027 * re ** 0.8 * pr ** (1 / 3) * (coolant_properties[station][5] / coolant_wall_viscosity) ** 0.14
        return nu * coolant_properties[station][7] / hydraulic_diameter

    # Calculate the heat balance

    # Initial guesses for gas wall temp and coolant temp, well below actual
    gas_wall_temp = 800
    coolant_temp = 250

    q = 5

    step_size = 100
    while abs(coolant_properties[station][2] - coolant_temp) > 0.00001:
        q = gas_transfer_coefficient(x, gas_wall_temp) * station_length * 2 * math.pi * inner_radius(x) * (
                    adiabatic_wall_temp(x) - gas_wall_temp)  # TODO update latex document with changes
        conduction_resistance = math.log((inner_radius(x) + inner_wall_thickness) / inner_radius(x)) / (
                    2 * math.pi * station_length * wall_thermal_conductivity)
        coolant_wall_temp = gas_wall_temp - q * conduction_resistance
        convection_resistance = 1 / (coolant_transfer_coefficient(station, coolant_wall_temp,
                                                                  coolant_properties) * 2 * math.pi * station_length * (
                                                 inner_radius(x) + inner_wall_thickness))
        coolant_temp = coolant_wall_temp - q * convection_resistance
        if coolant_temp > coolant_properties[station][2]:
            gas_wall_temp -= step_size
            step_size /= 10
        else:
            gas_wall_temp += step_size

    print("coolant temp = ", coolant_temp)

    data[station][0] = gas_wall_temp
    data[station][1] = q

    # Calculate the coolant pressure at the next station
    # Based off of Naraghi and Foulon, need to double check

#TODO fix this mess
    def coolant_area1(position):
        outer_area = math.pi * (inner_radius1(position) + inner_wall_thickness + channel_height(position)) ** 2
        inner_area = math.pi * (inner_radius1(position) + inner_wall_thickness) ** 2

        return outer_area - inner_area

    # Pressure drop for expansion and contraction
    k_expansion = ((channel_height1(x) / channel_height1(x - station_length)) ** 2 - 1) ** 2
    k_contraction = 0.5 - 0.167 * (channel_height1(x - station_length) / channel_height1(x)) \
                    - 0.125 * (channel_height1(x - station_length) / channel_height1(x)) ** 2 \
                    - 0.208 * (channel_height1(x - station_length) / channel_height1(x)) ** 3
    k = k_contraction if x > throat_position else k_expansion
    expansion_contraction_drop = k * coolant_properties[station][3] * coolant_properties[station][4] ** 2 / 2

    # Fluid acceleration pressure drop
    density = coolant_properties[station][3]  # TODO same density for both?
    accel_drop = 2 / (coolant_area1(x) + coolant_area1(x - station_length)) * coolant_flow_rate ** 2 \
                 * (1 / (density * coolant_area1(x)) - 1 / (density * coolant_area1(x - station_length)))

    # Viscous pressure drop
    def friction(f):
        re = coolant_properties[station][3] * coolant_properties[station][4] * channel_height1(x) / \
             coolant_properties[station][5]
        c = -2 * math.log((surface_roughness / (3.7065 * channel_height1(x)) + 2.5226 / (re * math.sqrt(f))), 10)
        return 1 / c ** 2

    friction_factor = fixed_point(friction, 0.02)

    viscous_drop = friction_factor * coolant_properties[station][3] * coolant_properties[station][
        4] ** 2 * station_length / (2 * channel_height1(x))

    pressure_drop1 = expansion_contraction_drop + accel_drop + viscous_drop

    enthalpy_change = data[station][1] / coolant_flow_rate
    coolant_properties[station + 1][0] = coolant_properties[station][0] + enthalpy_change
    coolant_properties[station + 1][1] = coolant_properties[station][1] - pressure_drop1
    coolant_properties[station + 1][2] = PropsSI('T', 'H', coolant_properties[station + 1][0], 'P',
                                                 coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][3] = PropsSI('D', 'H', coolant_properties[station + 1][0], 'P',
                                                 coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][4] = coolant_flow_rate / (
                coolant_properties[station][3] * coolant_area(bounds[4] - station * station_length))
    coolant_properties[station + 1][5] = PropsSI('V', 'H', coolant_properties[station + 1][0], 'P',
                                                 coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][6] = PropsSI('C', 'H', coolant_properties[station + 1][0], 'P',
                                                 coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][7] = PropsSI('L', 'H', coolant_properties[station + 1][0], 'P',
                                                 coolant_properties[station + 1][1], FLUID)


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
    throat_area = math.pi * throat_radius ** 2
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

    # h_0 = 0.026 * prandtl_number ** 0.4 * (4 / math.pi) ** 0.8 * 1.555 * viscosity ** -0.8
    # h = h_0 * 8.12 ** 0.8 * (inner_radius(x) * 2) ** -1.8

    bartz = 0.026 / (throat_radius * 2) ** 0.2 * viscosity ** 0.2 * specific_heat / prandtl_number ** 0.6\
           * (chamber_pressure / characteristic_velocity) ** 0.8 * (throat_radius * 2 / throat_bevel_radius) ** 0.1 * (1 / area_ratio(x)) ** 0.9 * sigma

    # print("Bartz = ", bartz)
    # print("Dittus = ", h)

    return bartz


def channel_height(x):
    return 2 * inner_radius(x) * channel_percent / 100


def coolant_transfer_coefficient(station, coolant_wall_temp, coolant_properties):
    x = bounds[4] - station * station_length
    hydraulic_diameter = channel_height(x)  # TODO
    coolant_wall_viscosity = PropsSI('V', 'T', coolant_wall_temp, 'P', coolant_properties[station][1], FLUID)
    re = coolant_properties[station][3] * coolant_properties[station][4] * hydraulic_diameter / coolant_properties[station][5]
    pr = coolant_properties[station][5] * coolant_properties[station][6] / coolant_properties[station][7]
    nu = 0.027 * re ** 0.8 * pr ** (1/3) * (coolant_properties[station][5] / coolant_wall_viscosity) ** 0.14
    print("coolant viscosity = ", coolant_properties[station][5], "wall viscosity = ", coolant_wall_viscosity)
    return nu * coolant_properties[station][7] / hydraulic_diameter


def coolant_transfer_coefficient_test(station, coolant_properties):
    x = bounds[4] - station * station_length
    hydraulic_diameter = channel_height(x)
    re = coolant_properties[station][3] * coolant_properties[station][4] * hydraulic_diameter / coolant_properties[station][5]
    pr = coolant_properties[station][5] * coolant_properties[station][6] / coolant_properties[station][7]
    nu = 0.023 * re ** 0.8 * pr ** 0.4
    return nu * coolant_properties[station][7] / hydraulic_diameter


def calculate_heat_balance(station, coolant_properties, data):
    # Initial guesses for gas wall temp and coolant temp, well below actual
    gas_wall_temp = 800
    coolant_temp = 250

    x = bounds[4] - station * station_length

    q = 5

    step_size = 100
    while abs(coolant_properties[station][2] - coolant_temp) > 0.00001:
        #TODO: q in this case is actually in Watts, not heat flux  q W/m^2, rename variable
        q = gas_transfer_coefficient(x, gas_wall_temp) * station_length * 2 * math.pi * inner_radius(x) * (adiabatic_wall_temp(x) - gas_wall_temp)  # TODO update latex document with changes
        conduction_resistance = math.log((inner_radius(x) + inner_wall_thickness) / inner_radius(x)) / (2 * math.pi * station_length * wall_thermal_conductivity)
        coolant_wall_temp = gas_wall_temp - q * conduction_resistance
        convection_resistance = 1 / (coolant_transfer_coefficient(station, coolant_wall_temp, coolant_properties) * 2 * math.pi * station_length * (inner_radius(x) + inner_wall_thickness))
        coolant_temp = coolant_wall_temp - q * convection_resistance
        if coolant_temp > coolant_properties[station][2]:
            gas_wall_temp -= step_size
            step_size /= 10
        else:
            gas_wall_temp += step_size

    print("coolant temp = ", coolant_temp)

    data[station][0] = gas_wall_temp
    data[station][1] = q
    

def coolant_area(station):
    x = bounds[4] - station * station_length
    outer_area = math.pi * (inner_radius(x) + inner_wall_thickness + channel_height(x)) ** 2
    inner_area = math.pi * (inner_radius(x) + inner_wall_thickness) ** 2

    return outer_area - inner_area


def coolant_pressure_drop(station, coolant_properties):
    # Based off of naraghi and foulon
    x = bounds[4] - station * station_length

    # Pressure drop for expansion and contraction
    k_expansion = ((channel_height(x) / channel_height(x - station_length)) ** 2 - 1) ** 2
    k_contraction = 0.5 - 0.167 * (channel_height(x - station_length) / channel_height(x))\
                    - 0.125 * (channel_height(x - station_length) / channel_height(x)) ** 2\
                    - 0.208 * (channel_height(x - station_length) / channel_height(x)) ** 3
    k = k_contraction if x > throat_position else k_expansion
    expansion_contraction_drop = k * coolant_properties[station][3] * coolant_properties[station][4] ** 2 / 2

    # Fluid acceleration pressure drop
    density = coolant_properties[station][3]  # TODO same density for both?
    accel_drop = 2 / (coolant_area(station) + coolant_area(station + 1)) * coolant_flow_rate ** 2 \
                 * (1 / (density * coolant_area(station)) - 1 / (density * coolant_area(station + 1)))

    # Viscous pressure drop
    def friction(f):
        re = coolant_properties[station][3] * coolant_properties[station][4] * channel_height(x) / \
             coolant_properties[station][5]
        a = -2 * math.log((surface_roughness / (3.7065 * channel_height(x)) + 2.5226 / (re * math.sqrt(f))), 10)
        return 1 / a ** 2

    friction_factor = fixed_point(friction, 0.02)

    viscous_drop = friction_factor * coolant_properties[station][3] * coolant_properties[station][4] ** 2 * station_length / (2 * channel_height(x))

    return expansion_contraction_drop + accel_drop + viscous_drop


# Calculates coolant properties at next station
def update_coolant_props(station, coolant_properties, data):
    enthalpy_change = data[station][1] / coolant_flow_rate
    pressure_drop = coolant_pressure_drop(station, coolant_properties)
    coolant_properties[station + 1][0] = coolant_properties[station][0] + enthalpy_change
    coolant_properties[station + 1][1] = coolant_properties[station][1] - pressure_drop
    coolant_properties[station + 1][2] = PropsSI('T', 'H', coolant_properties[station + 1][0], 'P', coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][3] = PropsSI('D', 'H', coolant_properties[station + 1][0], 'P', coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][4] = coolant_flow_rate / (coolant_properties[station][3] * coolant_area(bounds[4] - station * station_length))
    coolant_properties[station + 1][5] = PropsSI('V', 'H', coolant_properties[station + 1][0], 'P', coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][6] = PropsSI('C', 'H', coolant_properties[station + 1][0], 'P', coolant_properties[station + 1][1], FLUID)
    coolant_properties[station + 1][7] = PropsSI('L', 'H', coolant_properties[station + 1][0], 'P', coolant_properties[station + 1][1], FLUID)


def calculate_stress(x, temp):
    max_stress = (fuel_input_pressure - chamber_pressure) * inner_radius(x) / inner_wall_thickness  # TODO pressures (propagate), thermal etc
    yield_strength = 1570.3 - 14.184 * temp + 5.6410 * 10 ** (-2) * temp ** 2 - 1.0592 * 10 ** (-4) * temp ** 3 + 9.2881 * 10 ** (-8) * temp ** 4 - 3.086 * 10 ** (-11) * temp ** 5
    print("stress = ", max_stress / 1000000)
    print("yield strength = ", yield_strength)

# Initialize coolant properties
input_pressure = fuel_input_pressure
input_temperature = fuel_input_temperature

cool_props = [[0 for x in range(8)] for y in range(num_stations)]
cool_props[0][0] = PropsSI('H', 'P', input_pressure, 'T', input_temperature, FLUID)  # Stagnation enthalpy
cool_props[0][1] = input_pressure  # Pressure
cool_props[0][2] = input_temperature  # Temperature
cool_props[0][3] = PropsSI('D', 'P', input_pressure, 'T', input_temperature, FLUID)  # Density
cool_props[0][4] = coolant_flow_rate / (cool_props[0][3] * coolant_area(bounds[4]))  # Velocity
cool_props[0][5] = PropsSI('V', 'P', input_pressure, 'T', input_temperature, FLUID)  # Viscosity
cool_props[0][6] = PropsSI('C', 'P', input_pressure, 'T', input_temperature, FLUID)  # Heat capacity
cool_props[0][7] = PropsSI('L', 'P', input_pressure, 'T', input_temperature, FLUID)  # Thermal conductivity

heat_data = [[0 for x in range(3)] for y in range(num_stations)]
# gas wall temp is first column
# heat flux is second column

# for i in range(num_stations - 1):
    # test(i, cool_props, heat_data)
    # print("calculating ", i, "of ", num_stations)

print(gas_transfer_coefficient(throat_position, 300) * (adiabatic_wall_temp(throat_position)))
print(adiabatic_wall_temp(throat_position))
