import inputs_new as inp
import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math

sns.set()

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


def wall_temp(pos, fuel_temp):
    """
    Non-vectorized
    :param pos:
    :param fuel_temp:
    :return:
    """
    radius = inner_radius(pos)

    hydraulic_diameter = 2 * inp.channel_height
    density = fuel_density(fuel_temp)
    viscosity = fuel_viscosity(fuel_temp)

    velocity = inp.fuel_flow_rate / ((np.pi * (radius + inp.inner_wall_thickness + inp.channel_height) ** 2 - np.pi * (
            radius + inp.inner_wall_thickness) ** 2) * density)
    re = density * velocity * hydraulic_diameter / viscosity
    pr = viscosity * fuel_specific_heat(fuel_temp) / fuel_conductivity(fuel_temp)

    heat_flow = heat_flux(pos) * station_width * 2 * np.pi * radius

    def coolant_difference(coolant_wall_temp):
        coolant_wall_viscosity = fuel_viscosity(coolant_wall_temp)
        nu = 0.027 * re ** 0.8 * pr ** (1 / 3) * (viscosity / coolant_wall_viscosity) ** 0.14
        coolant_transfer_coefficient = nu * fuel_conductivity(fuel_temp) / hydraulic_diameter
        convection_resistance = 1 / (coolant_transfer_coefficient * 2 * np.pi * station_width * (
                    radius + inp.inner_wall_thickness))
        temp = coolant_wall_temp - heat_flow * convection_resistance
        return fuel_temp - temp

    coolant_wall_temperature = scipy.optimize.brentq(coolant_difference, inp.temp_lower_bound, inp.temp_upper_bound)
    conduction_resistance = math.log((radius + inp.inner_wall_thickness) / radius) / (2 * math.pi * station_width * inp.wall_thermal_conductivity)
    gas_wall_temperature = coolant_wall_temperature + heat_flow * conduction_resistance
    return coolant_wall_temperature, gas_wall_temperature


def calc_fuel_temp(pos):
    """
    Vectorized
    :param pos:
    :return:
    """
    heat_flow = heat_flux(pos) * station_width * 2 * math.pi * inner_radius(pos)
    enthalpy_change = heat_flow / inp.fuel_flow_rate
    enthalpy = np.cumsum(enthalpy_change[::-1])[::-1] + fuel_enthalpy(inp.fuel_input_temperature)
    return fuel_temperature(enthalpy)


def heat_flux(pos):
    # TODO
    return 500000 * (np.pi * throat_radius ** 2) / (np.pi * inner_radius(pos) ** 2)


def main():
    df = pd.DataFrame()
    df['Position'] = np.linspace(0, diverging_end, inp.num_stations)
    df['Radius'] = inner_radius(df['Position'].to_numpy())
    df['Heat Flux'] = heat_flux(df['Position'].to_numpy())
    df['Fuel Temperature'] = calc_fuel_temp(df['Position'].to_numpy())
    df['Coolant Wall Temperature'], df['Gas Wall Temperature'] = \
        np.vectorize(wall_temp)(df['Position'].to_numpy(), df['Fuel Temperature'].to_numpy())
    df.plot(x='Position')
    plt.show()


if __name__ == "__main__":
    main()
