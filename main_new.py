import inputs_new as inp
import numpy as np
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math

sns.set()

# TODO convert to math?
converging_start = inp.chamber_length + inp.chamber_bevel_radius * np.sin(np.radians(inp.converging_half_angle))
throat_start = converging_start + \
               (inp.chamber_diameter / 2 + inp.chamber_bevel_radius * np.cos(np.radians(inp.converging_half_angle)) -
                inp.chamber_bevel_radius - inp.throat_diameter / 2 - inp.throat_bevel_radius + inp.throat_bevel_radius *
                np.cos(np.radians(inp.converging_half_angle))) / np.tan(np.radians(inp.converging_half_angle))
throat_end = throat_start + inp.throat_bevel_radius * (np.sin(np.radians(inp.diverging_half_angle)) +
                                                       np.sin(np.radians(inp.converging_half_angle)))
diverging_end = throat_end + \
                (inp.exit_diameter / 2 -
                 (inp.throat_diameter / 2 + inp.throat_bevel_radius - inp.throat_bevel_radius *
                  np.cos(np.radians(inp.diverging_half_angle)))) / np.tan(np.radians(inp.diverging_half_angle))

throat_position = inp.chamber_length + inp.chamber_bevel_radius * np.sin(np.radians(inp.converging_half_angle)) + (inp.chamber_diameter / 2 + inp.chamber_bevel_radius * np.cos(np.radians(inp.converging_half_angle)) - inp.chamber_bevel_radius - inp.throat_diameter / 2 - inp.throat_bevel_radius + inp.throat_bevel_radius * np.cos(np.radians(inp.converging_half_angle))) / np.tan(np.radians(inp.converging_half_angle)) + inp.throat_bevel_radius * np.sin(np.radians(inp.converging_half_angle))


class BoundsError(Exception):
    pass


def props(property_name, temp):
    coeffs = {'dens': inp.density_coefficients, 'cp': inp.density_coefficients,
              'visc': inp.viscosity_coefficients, 'cond': inp.thermal_conductivity_coefficients}
    if temp < inp.temp_lower_bound or temp > inp.temp_upper_bound:
        raise BoundsError("Temperature out of range for " + property_name + ": " + str(temp))
    return np.poly1d(coeffs[property_name])(temp)


def fuel_temperature(enthalpy):
    if np.any(enthalpy < inp.enthalpy_lower_bound) or np.any(enthalpy > inp.enthalpy_upper_bound):
        raise BoundsError("Enthalpy out of range for temperature: " + str(enthalpy))
    return np.poly1d(inp.temperature_coefficients)(enthalpy)


def inner_radius(pos):
    conds = [pos <= inp.chamber_length, (pos > inp.chamber_length) & (pos <= converging_start),
             (pos > converging_start) & (pos <= throat_start), (pos > throat_start) & (pos <= throat_end),
             (pos > throat_end) & (pos <= diverging_end)]
    funcs = [inp.chamber_diameter / 2,
             lambda x: np.sqrt(inp.chamber_bevel_radius ** 2 - (x - inp.chamber_length) ** 2) + inp.chamber_diameter / 2 - inp.chamber_bevel_radius,
             lambda x: - np.tan(np.radians(inp.converging_half_angle)) * (x - inp.chamber_length - inp.chamber_bevel_radius * np.sin(np.radians(inp.converging_half_angle))) + inp.chamber_diameter / 2 + inp.chamber_bevel_radius * np.cos(np.radians(inp.converging_half_angle)) - inp.chamber_bevel_radius,
             lambda x: - np.sqrt(inp.throat_bevel_radius ** 2 - (x - throat_start - inp.throat_bevel_radius * np.sin(np.radians(inp.converging_half_angle))) ** 2) + inp.throat_diameter / 2 + inp.throat_bevel_radius,
             lambda x: np.tan(np.radians(inp.converging_half_angle)) * (x - throat_end) + inp.throat_diameter / 2 + inp.throat_bevel_radius - inp.throat_bevel_radius * np.cos(np.radians(inp.diverging_half_angle))]

    return np.piecewise(pos, conds, funcs)


def wall_temp(pos, fuel_temp):
    """
    Non-vectorized (just need to change optimization method)
    :param pos:
    :param fuel_temp:
    :return:
    """
    radius = inner_radius(pos)

    hydraulic_diameter = 2 * inp.channel_height
    fuel_density = props('dens', fuel_temp)
    fuel_viscosity = props('visc', fuel_temp)

    velocity = inp.fuel_flow_rate / ((np.pi * (radius + inp.inner_wall_thickness + inp.channel_height) ** 2 - np.pi * (
            radius + inp.inner_wall_thickness) ** 2) * fuel_density)
    re = fuel_density * velocity * hydraulic_diameter / fuel_viscosity
    pr = fuel_viscosity * props('cp', fuel_temp) / props('cond', fuel_temp)

    heat_flow = heat_flux(pos) * (diverging_end / inp.num_stations) * 2 * np.pi * radius

    def coolant_difference(coolant_wall_temp):
        coolant_wall_viscosity = props('visc', coolant_wall_temp)
        nu = 0.027 * re ** 0.8 * pr ** (1 / 3) * (fuel_viscosity / coolant_wall_viscosity) ** 0.14
        coolant_transfer_coefficient = nu * props('cond', fuel_temp) / hydraulic_diameter
        convection_resistance = 1 / (coolant_transfer_coefficient * 2 * np.pi * (diverging_end / inp.num_stations) * (
                    radius + inp.inner_wall_thickness))
        temp = coolant_wall_temp - heat_flow * convection_resistance
        return fuel_temp - temp

    coolant_wall_temperature = scipy.optimize.brentq(coolant_difference, inp.temp_lower_bound, inp.temp_upper_bound)
    conduction_resistance = math.log((radius + inp.inner_wall_thickness) / radius) / (2 * math.pi * (diverging_end / inp.num_stations) * inp.wall_thermal_conductivity)
    gas_wall_temperature = coolant_wall_temperature + heat_flow * conduction_resistance
    return coolant_wall_temperature, gas_wall_temperature


def calc_fuel_temp(pos):
    """
    Vectorized
    :param pos:
    :return:
    """
    heat_flow = heat_flux(pos) * (diverging_end / inp.num_stations) * 2 * math.pi * inner_radius(pos)
    enthalpy_change = heat_flow / inp.fuel_flow_rate
    enthalpy = np.cumsum(enthalpy_change[::-1])[::-1] + inp.input_enthalpy
    return fuel_temperature(enthalpy)


def heat_flux(pos):
    # TODO
    return 5E4


def main():
    df = pd.DataFrame()
    df['Position'] = np.linspace(0, diverging_end, inp.num_stations)
    df['Radius'] = inner_radius(df['Position'].to_numpy())
    df['Heat Flux'] = heat_flux(df['Position'].to_numpy())
    df['Fuel Temperature'] = calc_fuel_temp(df['Position'].to_numpy())
    df['Coolant Wall Temperature'], df['Gas Wall Temperature'] = \
        np.vectorize(wall_temp)(df['Position'].to_numpy(), df['Fuel Temperature'].to_numpy())
    df.plot(x='Position', y='Gas Wall Temperature')
    plt.show()


if __name__ == "__main__":
    main()
