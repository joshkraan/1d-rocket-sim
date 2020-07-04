import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math

import downsized as inp

sns.set()
sns.set_style("ticks")

div_ang_rad = math.radians(inp.diverging_half_angle)
conv_ang_rad = math.radians(inp.converging_half_angle)

chamber_radius = inp.chamber_diameter / 2
throat_radius = inp.throat_diameter / 2

converging_start = inp.chamber_length + inp.chamber_bevel_radius * math.sin(conv_ang_rad)
throat_start = converging_start + \
               (chamber_radius + inp.chamber_bevel_radius * math.cos(conv_ang_rad) - inp.chamber_bevel_radius -
                throat_radius - inp.throat_bevel_radius + inp.throat_bevel_radius * math.cos(conv_ang_rad)) / math.tan(conv_ang_rad)
throat_end = throat_start + inp.throat_bevel_radius * (math.sin(div_ang_rad) + math.sin(conv_ang_rad))
diverging_end = throat_end + \
                (inp.exit_diameter / 2 -
                 (throat_radius + inp.throat_bevel_radius - inp.throat_bevel_radius *
                  math.cos(div_ang_rad))) / math.tan(div_ang_rad)

throat_position = inp.chamber_length + inp.chamber_bevel_radius * math.sin(conv_ang_rad) + \
                  (chamber_radius + inp.chamber_bevel_radius * math.cos(conv_ang_rad) - inp.chamber_bevel_radius -
                   throat_radius - inp.throat_bevel_radius + inp.throat_bevel_radius * math.cos(conv_ang_rad)) / \
                  math.tan(conv_ang_rad) + inp.throat_bevel_radius * math.sin(conv_ang_rad)

station_width = diverging_end / inp.num_stations

data = pd.read_csv(inp.fuel_properties_file)
fuel_density = interp1d(data["Temperature"], data["Density"], kind='cubic', bounds_error=True)
fuel_specific_heat = interp1d(data["Temperature"], data["Specific Heat"], kind='cubic', bounds_error=True)
fuel_conductivity = interp1d(data["Temperature"], data["Thermal Conductivity"], kind='cubic', bounds_error=True)
fuel_viscosity = interp1d(data["Temperature"], data["Viscosity"], kind='cubic', bounds_error=True)
fuel_enthalpy = interp1d(data["Temperature"], data["Enthalpy"], kind='cubic', bounds_error=True)
fuel_temperature = interp1d(data["Enthalpy"], data["Temperature"], kind='cubic', bounds_error=True)


def inner_radius(pos):
    conds = [pos <= inp.chamber_length, (pos > inp.chamber_length) & (pos <= converging_start),
             (pos > converging_start) & (pos <= throat_start), (pos > throat_start) & (pos <= throat_end),
             (pos > throat_end) & (pos <= diverging_end)]
    funcs = [chamber_radius,
             lambda x: np.sqrt(inp.chamber_bevel_radius ** 2 - (x - inp.chamber_length) ** 2) + chamber_radius - inp.chamber_bevel_radius,
             lambda x: - math.tan(conv_ang_rad) * (x - inp.chamber_length - inp.chamber_bevel_radius * math.sin(conv_ang_rad)) + chamber_radius + inp.chamber_bevel_radius * math.cos(conv_ang_rad) - inp.chamber_bevel_radius,
             lambda x: - np.sqrt(inp.throat_bevel_radius ** 2 - (x - throat_start - inp.throat_bevel_radius * math.sin(conv_ang_rad)) ** 2) + throat_radius + inp.throat_bevel_radius,
             lambda x: math.tan(div_ang_rad) * (x - throat_end) + throat_radius + inp.throat_bevel_radius - inp.throat_bevel_radius * math.cos(div_ang_rad)]

    return np.piecewise(pos, conds, funcs)


# def wall_temp(pos, fuel_temp):
#     """
#     Non-vectorized
#     :param pos:
#     :param fuel_temp:
#     :return:
#     """
#     radius = inner_radius(pos)
#
#     hydraulic_diameter = 2 * inp.channel_height
#     density = fuel_density(fuel_temp)
#     viscosity = fuel_viscosity(fuel_temp)
#
#     velocity = inp.fuel_flow_rate / ((np.pi * (radius + inp.inner_wall_thickness + inp.channel_height) ** 2 - np.pi * (
#             radius + inp.inner_wall_thickness) ** 2) * density)
#     re = density * velocity * hydraulic_diameter / viscosity
#     pr = viscosity * fuel_specific_heat(fuel_temp) / fuel_conductivity(fuel_temp)
#
#     heat_flow = heat_flux(pos) * station_width * 2 * np.pi * radius
#
#     def coolant_difference(coolant_wall_temp):
#         coolant_wall_viscosity = fuel_viscosity(coolant_wall_temp)
#         nu = 0.027 * re ** 0.8 * pr ** (1 / 3) * (viscosity / coolant_wall_viscosity) ** 0.14
#         coolant_transfer_coefficient = nu * fuel_conductivity(fuel_temp) / hydraulic_diameter
#         convection_resistance = 1 / (coolant_transfer_coefficient * 2 * np.pi * station_width * (
#                     radius + inp.inner_wall_thickness))
#         temp = coolant_wall_temp - heat_flow * convection_resistance
#         return fuel_temp - temp
#
#     coolant_wall_temperature = scipy.optimize.brentq(coolant_difference, inp.temp_lower_bound, inp.temp_upper_bound)
#     conduction_resistance = math.log((radius + inp.inner_wall_thickness) / radius) / (2 * math.pi * station_width * inp.wall_thermal_conductivity)
#     gas_wall_temperature = coolant_wall_temperature + heat_flow * conduction_resistance
#     return coolant_wall_temperature, gas_wall_temperature


def wall_temp1(radius, fuel_temp, heatflux):
    density, viscosity, conductivity, velocity, re, pr = calc_fuel_props(radius, fuel_temp)

    heat_flow = heatflux * station_width * 2 * np.pi * radius

    base_gnielinski = gnielinkski(re, pr)
    max_temp = data["Temperature"].iloc[-1]

    coolant_wall_temperature = np.zeros(radius.size)

    for i in range(radius.size):
        def coolant_difference(coolant_wall_temp):
            hydraulic_diameter = 2 * inp.channel_height
            pr_w = fuel_viscosity(coolant_wall_temp) * fuel_specific_heat(coolant_wall_temp) / fuel_conductivity(coolant_wall_temp)
            coolant_transfer_coeff = base_gnielinski[i] * (pr[i] / pr_w) ** 0.11 * conductivity[i] / hydraulic_diameter
            convection_resistance = 1 / (coolant_transfer_coeff * 2 * np.pi * station_width * (radius[i] + inp.inner_wall_thickness))
            return fuel_temp[i] - coolant_wall_temp + heat_flow[i] * convection_resistance

        coolant_wall_temperature[i] = scipy.optimize.brentq(coolant_difference, fuel_temp[i], max_temp)

    conduction_resistance = np.log((radius + inp.inner_wall_thickness) / radius) / (2 * np.pi * station_width * inp.wall_thermal_conductivity)
    gas_wall_temperature = coolant_wall_temperature + heat_flow * conduction_resistance

    return coolant_wall_temperature, gas_wall_temperature


def gnielinkski(re, pr):
    f = (0.79 * np.log(re) - 1.64) ** (-2)
    nu = (f / 8) * (re - 1000) * pr / (1 + 12.7 * (f / 8) ** 0.5 * (pr ** (2 / 3) - 1))
    return nu


def wall_temp2(radius, fuel_temp, heatflux):
    density, viscosity, conductivity, velocity, re, pr = calc_fuel_props(radius, fuel_temp)

    heat_flow = heatflux * station_width * 2 * np.pi * radius

    max_temp = data["Temperature"].iloc[-1]

    coolant_wall_temperature = np.zeros(radius.size)

    for i in range(radius.size):
        def coolant_difference(coolant_wall_temp):
            hydraulic_diameter = 2 * inp.channel_height
            pr_w = fuel_viscosity(coolant_wall_temp) * fuel_specific_heat(coolant_wall_temp) / fuel_conductivity(coolant_wall_temp)

            # Annular form from gnielinski paper
            a = 2 * radius[i] / (2 * radius[i] + 2 * inp.channel_height)
            k1 = 1.07 + 900 / re[i] - 0.63 / (1 + 10 * pr[i])
            modified_re = re[i] * ((1 + a ** 2) * np.log(a) + (1 - a ** 2)) / ((1 - a) ** 2 * np.log(a))
            f = (1.8 * np.log10(modified_re) - 1.5) ** -2
            f_ann = 0.75 * a ** -0.17
            k = (pr[i] / pr_w) ** 0.11
            nu = (f / 8) * re[i] * pr[i] / (k1 + 12.7 * np.sqrt(f/8) * (pr[i] ** (2/3) - 1)) * f_ann * k

            coolant_transfer_coeff = nu * conductivity[i] / hydraulic_diameter
            convection_resistance = 1 / (coolant_transfer_coeff * 2 * np.pi * station_width * (radius[i] + inp.inner_wall_thickness))
            return fuel_temp[i] - coolant_wall_temp + heat_flow[i] * convection_resistance

        coolant_wall_temperature[i] = scipy.optimize.brentq(coolant_difference, fuel_temp[i], max_temp)

    conduction_resistance = np.log((radius + inp.inner_wall_thickness) / radius) / (2 * np.pi * station_width * inp.wall_thermal_conductivity)
    gas_wall_temperature = coolant_wall_temperature + heat_flow * conduction_resistance

    return coolant_wall_temperature, gas_wall_temperature


def calc_fuel_temp(radius, heatflux):
    heat_flow = heatflux * station_width * 2 * math.pi * radius
    enthalpy_change = heat_flow / inp.fuel_flow_rate
    enthalpy = np.cumsum(enthalpy_change[::-1])[::-1] + fuel_enthalpy(inp.fuel_input_temperature)
    return fuel_temperature(enthalpy)


def calc_fuel_props(radius, fuel_temp):
    hydraulic_diameter = 2 * inp.channel_height
    density = fuel_density(fuel_temp)
    viscosity = fuel_viscosity(fuel_temp)
    conductivity = fuel_conductivity(fuel_temp)

    velocity = inp.fuel_flow_rate / ((np.pi * (radius + inp.inner_wall_thickness + inp.channel_height) ** 2 - np.pi * (
            radius + inp.inner_wall_thickness) ** 2) * density)
    re = density * velocity * hydraulic_diameter / viscosity
    pr = viscosity * fuel_specific_heat(fuel_temp) / conductivity
    return density, viscosity, conductivity, velocity, re, pr


def bartz_heat_flux(pos):
    gas_wall_temp = 800  # Conservative estimate

    radius = inner_radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * throat_radius ** 2

    mach_number = calc_mach_number(pos)
    recovery_factor = inp.gas_prandtl_number ** 0.33
    a = 1 + recovery_factor * (inp.gas_gamma - 1) / 2 * mach_number ** 2
    b = 1 + (inp.gas_gamma - 1) / 2 * mach_number ** 2
    adiabatic_wall_temp = inp.chamber_temperature * (a / b)

    g = 1 + (inp.gas_gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / inp.chamber_temperature) * g + 0.5) ** 0.68 * g ** 0.12)
    htc = 0.026 / inp.throat_diameter ** 0.2 * (inp.gas_viscosity ** 0.2 * inp.gas_specific_heat / inp.gas_prandtl_number ** 0.6) * \
          (inp.chamber_pressure / inp.cstar) ** 0.8 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma

    overall_htc = 1 / (1 / htc + inp.carbon_resistance)

    return overall_htc * (adiabatic_wall_temp - gas_wall_temp)


def calc_mach_number(pos):
    radius = inner_radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * throat_radius ** 2
    area_ratio = local_area / throat_area

    gamma = inp.gas_gamma

    mach_number = np.zeros(pos.size)
    for i in range(pos.size):
        def mach_equation(m):
            return 1 / m ** 2 * (2 / (gamma + 1) * (1 + 0.5 * (gamma - 1) * m ** 2)) ** ((gamma + 1) / (gamma - 1)) - area_ratio[i] ** 2
        if pos[i] < throat_position:
            mach_number[i] = brentq(mach_equation, 1e-6, 1)
        else:
            mach_number[i] = brentq(mach_equation, 1, 50)

    return mach_number


def stress(radius, heat_flux):
    hoop_stress = (inp.coolant_pressure - inp.chamber_pressure) * radius / inp.inner_wall_thickness
    thermal_stress = inp.elastic_modulus * inp.thermal_expansion_coeff * heat_flux * \
                     inp.inner_wall_thickness / (2 * (1 - inp.poissons_ratio) * inp.wall_thermal_conductivity)
    return hoop_stress + thermal_stress


def pressure_drop(radius, fuel_temp):
    density, viscosity, conductivity, velocity, re, pr = calc_fuel_props(radius, fuel_temp)
    hydraulic_diameter = 2 * inp.channel_height
    a = -2 * np.log10(inp.surface_roughness / hydraulic_diameter / 3.7 + 12 / re)
    b = -2 * np.log10(inp.surface_roughness / hydraulic_diameter / 3.7 + 2.51 * a / re)
    c = -2 * np.log10(inp.surface_roughness / hydraulic_diameter / 3.7 + 2.51 * b / re)
    f = (a - (b - a) ** 2 / (c - 2 * b + a)) ** (-2)
    drop = f * density * station_width * velocity ** 2 / (2 * hydraulic_diameter)
    return np.sum(drop)


def main():
    df = pd.DataFrame()
    df['Position'] = np.linspace(0, diverging_end, inp.num_stations)
    df['Radius'] = inner_radius(df['Position'].to_numpy())
    df['Heat Flux'] = bartz_heat_flux(df['Position'].to_numpy())
    #df['Fuel Temperature'] = calc_fuel_temp(df['Radius'].to_numpy(), df['Heat Flux'].to_numpy())
    #df['Coolant Wall Temperature'], df['Gas Wall Temperature'] = wall_temp2(df['Radius'].to_numpy(), df['Fuel Temperature'].to_numpy(), df['Heat Flux'].to_numpy())

    df['Mach'] = calc_mach_number(df['Position'].to_numpy())
    df['Stress'] = stress(df['Radius'].to_numpy(), df['Heat Flux'].to_numpy())

    df.plot(x='Position', y='Radius')
    sns.despine()

    #print("Pressure Drop: ", pressure_drop(df['Radius'].to_numpy(), df['Fuel Temperature'].to_numpy()))

    plt.show()


if __name__ == "__main__":
    main()
