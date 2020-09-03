import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import cantera as ct
import math
import time

import inputs as inp
import engine_geometry as geom
from gas_properties import calc_gas_properties

sns.set()
sns.set_style("ticks")

data = pd.read_csv(inp.fuel_properties_file)
fuel_density = interp1d(data["Temperature"], data["Density"], kind='cubic', bounds_error=True)
fuel_specific_heat = interp1d(data["Temperature"], data["Specific Heat"], kind='cubic', bounds_error=True)
fuel_conductivity = interp1d(data["Temperature"], data["Thermal Conductivity"], kind='cubic', bounds_error=True)
fuel_viscosity = interp1d(data["Temperature"], data["Viscosity"], kind='cubic', bounds_error=True)
fuel_enthalpy = interp1d(data["Temperature"], data["Enthalpy"], kind='cubic', bounds_error=True)
fuel_temperature = interp1d(data["Enthalpy"], data["Temperature"], kind='cubic', bounds_error=True)


def wall_temp(radius, fuel_temp, heatflux):
    density, viscosity, conductivity, velocity, re, pr = calc_fuel_props(radius, fuel_temp)

    heat_flow = heatflux * geom.station_width * 2 * np.pi * radius

    max_temp = data["Temperature"].iloc[-1]

    coolant_wall_temperature = np.zeros(radius.size)

    for i in range(radius.size):
        def coolant_difference(coolant_wall_temp):
            hydraulic_diameter = 2 * inp.channel_height
            pr_w = fuel_viscosity(coolant_wall_temp) * fuel_specific_heat(coolant_wall_temp) / fuel_conductivity(coolant_wall_temp)

            if re[i] < 3e3:
                print("Warning, Re out of range:", re[i])

            #nu = gnielinski_annular(re[i], pr[i], pr_w, 2 * radius[i] / (2 * radius[i] + 2 * inp.channel_height))

            nu = gnielinski(re[i], pr[i], pr_w)

            #nu = sieder_tate(re[i], pr[i], viscosity[i], fuel_viscosity(coolant_wall_temp))

            coolant_transfer_coeff = nu * conductivity[i] / hydraulic_diameter
            convection_resistance = 1 / (coolant_transfer_coeff * 2 * np.pi * geom.station_width * (radius[i] + inp.inner_wall_thickness))
            return fuel_temp[i] - coolant_wall_temp + heat_flow[i] * convection_resistance

        # coolant_wall_temperature[i] = scipy.optimize.brentq(coolant_difference, 300, max_temp)
        try:
            coolant_wall_temperature[i] = scipy.optimize.brentq(coolant_difference, inp.fuel_input_temperature, max_temp)
        except ValueError:
            print(radius[i])
            print(coolant_difference(300))
            print(coolant_difference(max_temp))
            raise Exception("test")

    conduction_resistance = np.log((radius + inp.inner_wall_thickness) / radius) / (2 * np.pi * geom.station_width * inp.wall_thermal_conductivity)
    gas_wall_temperature = coolant_wall_temperature + heat_flow * conduction_resistance

    return coolant_wall_temperature, gas_wall_temperature


def sieder_tate(re, pr, viscosity, wall_temp_viscosity):
    return 0.027 * re ** 0.8 * pr ** (1 / 3) * (viscosity / wall_temp_viscosity) ** 0.14


def dittus_boelter(re, pr):
    return 0.023 * re ** 0.8 * pr ** 0.3


def gnielinski(re, pr, pr_w):
    f = (0.79 * np.log(re) - 1.64) ** -2
    return (f / 8) * (re - 1000) * pr / (1 + 12.7 * (f / 8) ** 0.5 * (pr ** (2 / 3) - 1)) * (pr / pr_w) ** 0.11


def gnielinski_annular(re, pr, pr_w, diameter_ratio):
    a = diameter_ratio
    k1 = 1.07 + 900 / re - 0.63 / (1 + 10 * pr)
    modified_re = re * ((1 + a ** 2) * np.log(a) + (1 - a ** 2)) / ((1 - a) ** 2 * np.log(a))
    f = (1.8 * np.log10(modified_re) - 1.5) ** -2
    f_ann = 0.75 * a ** -0.17
    k = (pr / pr_w) ** 0.11
    return (f / 8) * re * pr / (k1 + 12.7 * np.sqrt(f / 8) * (pr ** (2 / 3) - 1)) * f_ann * k


def calc_fuel_temp(radius, heatflux):
    heat_flow = heatflux * geom.station_width * 2 * math.pi * radius
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
    drop = f * density * geom.station_width * velocity ** 2 / (2 * hydraulic_diameter)
    return np.sum(drop)


def main():
    # TODO remove station_width, directly calculate in functions
    position = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)
    chamber_only = position[position < inp.chamber_length]

    # bartz = bartz_heat_flux(position)
    # bartz1 = bartz_heat_flux1(position)
    # sns.lineplot(position, bartz)
    # sns.lineplot(position, bartz1)
    # pressure, density, viscosity, velocity, gamma, thermal_cond, prandtl = calc_gas_properties(position)
    # sns.lineplot(position, viscosity)
    gas, states, mach = calc_gas_properties(position)

    t_wall = 500

    # temperature, pressure = states.TP
    # density = states.density
    # specific_heat = states.cp
    # thermal_conductivity = states.thermal_conductivity
    # gamma = specific_heat / states.cv
    # viscosity = states.viscosity
    # velocity = mach * np.sqrt(gamma * pressure / density)
    # re = density * velocity * (2 * geom.radius(position)) / viscosity
    # pr = specific_heat * viscosity / thermal_conductivity
    #
    # states.TP = t_wall, pressure
    # viscosity_wall = states.viscosity
    #
    # states.TP = (t_wall + temperature) / 2, pressure
    # viscosity_am = states.viscosity
    # density_am = states.density
    # re_am = density_am * velocity * (2 * geom.radius(position)) / viscosity_am
    # pr_am = specific_heat * viscosity_am / thermal_conductivity
    #
    # sieder = sieder_tate(re, pr, viscosity, viscosity_wall)
    # dittus = dittus_boelter(re_am, pr_am)
    #
    # plt.plot(position, sieder)
    # plt.plot(position, dittus)
    # plt.show()

    # radius = geom.radius(position)
    # heat_flux = scaled_heat_flux(position)
    # fuel_temp = calc_fuel_temp(radius, heat_flux)
    # coolant_wall_temp, gas_wall_temp = wall_temp(radius, fuel_temp, heat_flux)
    # mach = calc_mach_number(position)
    # wall_stress = stress(radius, heat_flux)
    # density, viscosity, conductivity, velocity, re, pr = calc_fuel_props(radius, fuel_temp)
    #
    # df = pd.DataFrame({'Position': position, 'Radius': radius, 'Heat Flux': heat_flux, 'Fuel Temperature': fuel_temp,
    #                    'Coolant Wall Temperature': coolant_wall_temp, 'Gas Wall Temperature': gas_wall_temp,
    #                    'Mach': mach, 'Stress': wall_stress, 'Re': re})
    # df.plot(x='Position', y='Re')
    # sns.despine()
    #
    # plt.show()


if __name__ == "__main__":
    main()
