import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve


from input_processing import (
    # Chamber dimensions
    throat_radius, diverging_half_angle, converging_half_angle,
    chamber_length, chamber_radius, throat_bevel_radius, chamber_bevel_radius,
    # Cooling
    inner_wall_thickness, channel_height, wall_thermal_conductivity,
    coolant_flow_rate, surface_roughness, fuel_input_pressure, fuel_input_temperature,
    # Combustion gas properties
    gas_gamma, gas_viscosity, gas_specific_heat, gas_prandtl_number,
    # Engine properties
    chamber_pressure, chamber_temperature, characteristic_velocity,
    # Bounds
    bounds, throat_position,
    # High level
    num_stations, station_length, fuel_input_enthalpy, guess_for_exit_wall_temp
)


class BoundsError(Exception):
    pass


# Polynomial fits for kerosene properties at 1.7Mpa for temperatures between 270K and 650K.
def density(temperature):
    if temperature < 270 or temperature > 650:
        raise BoundsError("Temperature out of range for density: " + str(temperature))
    return -7.35246E-13 * temperature ** 6 + 1.89632E-09 * temperature ** 5 - 2.01450E-06 * temperature ** 4 + 1.12544E-03 * temperature ** 3 - 3.48382E-01 * temperature ** 2 + 5.59004E+01 * temperature - 2.75598E+03


def specific_heat(temperature):
    return 2500
    if temperature < 270 or temperature > 650:
        raise BoundsError("Temperature out of range for specific heat: " + str(temperature))
    return 1.40053E-11 * temperature ** 6 - 3.66091E-08 * temperature ** 5 + 3.93278E-05 * temperature ** 4 - 2.22085E-02 * temperature ** 3 + 6.94910E+00 * temperature ** 2 - 1.13776E+03 * temperature + 7.76898E+04


def viscosity(temperature):
    return 0.003
    if temperature < 270 or temperature > 650:
        raise BoundsError("Temperature out of range for viscosity: " + str(temperature))
    return 2.39318E-17 * temperature ** 6 - 7.07636E-14 * temperature ** 5 + 8.64557E-11 * temperature ** 4 - 5.58968E-08 * temperature ** 3 + 2.01936E-05 * temperature ** 2 - 3.87381E-03 * temperature + 3.09736E-01


def thermal_conductivity(temperature):
    return 0.1
    if temperature < 270 or temperature > 650:
        raise BoundsError("Temperature out of range for thermal conductivity: " + str(temperature))
    return 6.27385E-17 * temperature ** 6 - 1.63841E-13 * temperature ** 5 + 1.76034E-10 * temperature ** 4 - 9.95363E-08 * temperature ** 3 + 3.13833E-05 * temperature ** 2 - 5.41772E-03 * temperature + 5.25984E-01


def kerosene_temperature(enthalpy):
    if enthalpy < -466070 or enthalpy > 547770:
        raise BoundsError("Enthalpy out of range for temperature: " + str(enthalpy))
    return -1.10030E-34 * enthalpy ** 6 + 1.57293E-29 * enthalpy ** 5 - 3.53013E-23 * enthalpy ** 4 + 6.01412E-17 * enthalpy ** 3 - 1.08504E-10 * enthalpy ** 2 + 3.70175E-04 * enthalpy + 4.75311E+02


def inner_radius(x):
    # x=0 is defined at the injector plate
    if x <= bounds[0]:
        return chamber_radius
    elif bounds[0] < x <= bounds[1]:
        return math.sqrt(
            chamber_bevel_radius ** 2 - (x - chamber_length) ** 2) + chamber_radius - chamber_bevel_radius
    elif bounds[1] < x <= bounds[2]:
        return - math.tan(converging_half_angle) * (x - chamber_length - chamber_bevel_radius * math.sin(
            converging_half_angle)) + chamber_radius + chamber_bevel_radius * math.cos(
            converging_half_angle) - chamber_bevel_radius
    elif bounds[2] < x <= bounds[3]:
        return - math.sqrt(throat_bevel_radius ** 2 - (x - bounds[2] - throat_bevel_radius * math.sin(
            converging_half_angle)) ** 2) + throat_radius + throat_bevel_radius
    elif bounds[3] < x:
        return math.tan(diverging_half_angle) * (
                x - bounds[3]) + throat_radius + throat_bevel_radius - throat_bevel_radius * math.cos(
            diverging_half_angle)


def test(station, data):
    # x=0 is defined at the injector plate
    x = bounds[4] - station * station_length

    # Calculate the radius of the chamber at the current position
    radius = inner_radius(x)

    # Calculate the area ratio
    local_area = math.pi * radius ** 2
    throat_area = math.pi * throat_radius ** 2
    area_ratio = local_area / throat_area

    # Calculate the mach number
    # Solved iteratively

    mach_number = 0
    if x < throat_position:
        mach_current = 0.3
        mach_last = 0
        while abs(mach_current - mach_last) > 0.00001:
            mach_last = mach_current
            mach_current = 1 / area_ratio * ((2 + (gas_gamma - 1) * mach_last ** 2) / (gas_gamma + 1)) ** (
                    (gas_gamma + 1) / (2 * (gas_gamma - 1)))
        mach_number = mach_current
    elif x == throat_position:
        mach_number = 1
    elif x > throat_position:
        mach_current = 1.2
        mach_last = 0
        while abs(mach_current - mach_last) > 0.00001:
            mach_last = mach_current
            mach_current = math.sqrt(
                ((gas_gamma + 1) * (area_ratio * mach_last) ** ((2 * (gas_gamma - 1)) / (gas_gamma + 1)) - 2) / (
                        gas_gamma - 1))
        mach_number = mach_current

    # Calculate the adiabatic wall temperature
    recovery_factor = gas_prandtl_number ** 0.33
    a = 1 + recovery_factor * (gas_gamma - 1) / 2 * mach_number ** 2
    b = 1 + (gas_gamma - 1) / 2 * mach_number ** 2
    adiabatic_wall_temp = chamber_temperature * (a / b)

    # Calculate the gas transfer coefficient
    def gas_transfer_coefficient(gas_wall_temp1):
        g = 1 + (gas_gamma - 1) / 2 * mach_number ** 2
        sigma = 1 / ((0.5 * (gas_wall_temp1 / chamber_temperature) * g + 0.5) ** 0.68 * g ** 0.12)

        return 0.026 / (
                throat_radius * 2) ** 0.2 * gas_viscosity ** 0.2 * gas_specific_heat / gas_prandtl_number ** 0.6 * \
               (chamber_pressure / characteristic_velocity) ** 0.8 * (throat_radius * 2 / throat_bevel_radius) ** 0.1 * \
               (1 / area_ratio) ** 0.9 * sigma

    fuel_temperature = kerosene_temperature(data[station][3])
    data[station][0] = fuel_temperature

    # Calculate the coolant transfer coefficient
    def coolant_transfer_coefficient(wall_temp):
        hydraulic_diameter = 2 * channel_height
        coolant_wall_viscosity = viscosity(wall_temp)
        velocity = coolant_flow_rate / ((math.pi * (radius + inner_wall_thickness + channel_height) ** 2 - math.pi * (
                radius + inner_wall_thickness) ** 2) * density(fuel_temperature))
        re = density(fuel_temperature) * velocity * hydraulic_diameter / viscosity(fuel_temperature)
        pr = viscosity(fuel_temperature) * specific_heat(fuel_temperature) / thermal_conductivity(fuel_temperature)
        nu = 0.027 * re ** 0.8 * pr ** (1 / 3) * (viscosity(fuel_temperature) / coolant_wall_viscosity) ** 0.14
        return nu * thermal_conductivity(fuel_temperature) / hydraulic_diameter

    # Calculate the heat balance

    conduction_resistance = math.log((radius + inner_wall_thickness) / radius) / (2 * math.pi * station_length * wall_thermal_conductivity)

    def heat_flux(gas_wall_temp1):
        return gas_transfer_coefficient(gas_wall_temp1) * (adiabatic_wall_temp - gas_wall_temp1)

    def heat_flow(gas_wall_temp1):
        return heat_flux(gas_wall_temp1) * station_length * 2 * math.pi * radius

    def coolant_wall_temp(gas_wall_temp1):
        return gas_wall_temp1 - heat_flow(gas_wall_temp1) * conduction_resistance

    def convection_resistance(gas_wall_temp1):
        return 1 / (coolant_transfer_coefficient(coolant_wall_temp(gas_wall_temp1)) * 2 * math.pi * station_length * (
                    radius + inner_wall_thickness))

    def coolant_temp(gas_wall_temp1):
        return coolant_wall_temp(gas_wall_temp1) - heat_flow(gas_wall_temp1) * convection_resistance(gas_wall_temp1)

    def funct2solve(gas_wall_temp1):
        try:
            return fuel_temperature - coolant_temp(gas_wall_temp1)
        except BoundsError:
            return 10000

    # Initial guesses for gas wall temp, well below actual
    gas_wall_temp_guess = data[station - 1][1] if station != 0 else guess_for_exit_wall_temp

    # print("gas transfer: " + str(gas_transfer_coefficient(300)))
    # print("convection resistance: " + str(convection_resistance(400)))
    # print("conduction resistance: " + str(conduction_resistance))
    # print("heat flow " + str(heat_flow(400)))
    # print("heat flux " + str(heat_flux(400)))
    # print("area " + str(station_length * 2 * math.pi * radius))
    #
    # plt.scatter([temp + 270 for temp in range(800)], [funct2solve(temp + 270) for temp in range(800)])
    # plt.show()

    gas_wall_temp, _, flag, message = fsolve(funct2solve, gas_wall_temp_guess, full_output=True)

    if flag != 1:
        raise Exception("Equilibrium not found for station " + str(station) + ", solver error: " + message + str(gas_wall_temp))

    # Storing Calculated Values
    data[station][1] = gas_wall_temp
    data[station][2] = heat_flux(gas_wall_temp)
    data[station][4] = gas_transfer_coefficient(gas_wall_temp) * (adiabatic_wall_temp - gas_wall_temp)
    data[station][5] = coolant_wall_temp(gas_wall_temp)
    # Enthalpy of station + 1 is current station enthalpy plus change in enthalpy
    data[station + 1][3] = data[station][3] + heat_flow(gas_wall_temp) / coolant_flow_rate


# Driving Code
calc_data = [[0 for x in range(8)] for y in range(num_stations + 1)]
calc_data[0][0] = fuel_input_temperature
calc_data[0][3] = fuel_input_enthalpy

# This is the array for storing calculated data at each station
# Column 0 is Temperature of Fluid (k)
# Column 1 is Gas Wall Temp (k)
# Column 2 is Heat Transfer (W)
# Column 3 is Enthalpy (J/Kg)

# for i in range(num_stations):
#     test(i, calc_data)
#
# fig = plt.figure()
# ax1 = fig.add_subplot(211)
# ax1 = plt.scatter([bounds[4] - station_length * station for station in range(num_stations + 1)],
#                   [row[2] for row in calc_data])
# axes = plt.gca()
#
#
# ax2 = fig.add_subplot(212)
# ax2 = plt.scatter([bounds[4] - station_length * station for station in range(num_stations + 1)],
#                   [row[1] for row in calc_data])
# axes = plt.gca()
#
#
# plt.show()

# print(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4])

pos = np.linspace(0, 0.5, 500)
radius = np.vectorize(inner_radius)(pos)

plt.scatter(pos, radius)
plt.show()

# TODO gas wall temp is fucked => anything else could be
