import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math

import inputs as inp
import engine_geometry as geom
import heat_flux as hf

sns.set()
sns.set_style("ticks")

data = pd.read_csv(inp.fuel_properties_file)
fuel_density = interp1d(data["Temperature"], data["Density"], kind='cubic', bounds_error=True)
fuel_specific_heat = interp1d(data["Temperature"], data["Specific Heat"], kind='cubic', bounds_error=True)
fuel_conductivity = interp1d(data["Temperature"], data["Thermal Conductivity"], kind='cubic', bounds_error=True)
fuel_viscosity = interp1d(data["Temperature"], data["Viscosity"], kind='cubic', bounds_error=True)
fuel_enthalpy = interp1d(data["Temperature"], data["Enthalpy"], kind='cubic', bounds_error=True)
fuel_temperature = interp1d(data["Enthalpy"], data["Temperature"], kind='cubic', bounds_error=True)


def wall_temp(radius, fuel_props, heatflux):
    fuel_temp, density, viscosity, conductivity, velocity, re, pr, hydraulic_diameter = fuel_props
    heat_flow = heatflux * geom.station_width * 2 * np.pi * radius
    max_temp = data["Temperature"].iloc[-1]
    coolant_wall_temperature = np.zeros(radius.size)

    for i in range(radius.size):
        def coolant_difference(coolant_wall_temp):
            if re[i] < 3e3:
                print("Warning, Re out of range:", re[i])

            pr_w = fuel_viscosity(coolant_wall_temp) * fuel_specific_heat(coolant_wall_temp) / fuel_conductivity(coolant_wall_temp)
            nu = gnielinski(re[i], pr[i], pr_w)
            # nu = sieder_tate(re[i], pr[i], viscosity[i], fuel_viscosity(coolant_wall_temp))
            # nu = gnielinski_annular(re[i], pr[i], pr_w, radius[i] / (radius[i] + inp.channel_height))

            coolant_transfer_coeff = nu * conductivity[i] / hydraulic_diameter
            convection_resistance = 1 / (coolant_transfer_coeff * 2 * np.pi * geom.station_width * (radius[i] + inp.inner_wall_thickness))
            return fuel_temp[i] - coolant_wall_temp + heat_flow[i] * convection_resistance

        try:
            coolant_wall_temperature[i] = scipy.optimize.root_scalar(coolant_difference,
                                                                     bracket=[inp.fuel_input_temperature, max_temp],
                                                                     method='brentq').root
        except ValueError:
            raise Exception("Coolant wall temperature exceeded boiling point of fuel.")

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


def calc_fuel_props(radius, heat_flux):
    heat_flow = heat_flux * geom.station_width * 2 * math.pi * radius
    enthalpy_change = heat_flow / inp.fuel_flow_rate
    enthalpy = np.cumsum(enthalpy_change[::-1])[::-1] + fuel_enthalpy(inp.fuel_input_temperature)
    fuel_temp = fuel_temperature(enthalpy)

    if inp.helical_channels is True:
        hydraulic_diameter = 2 * inp.helix_pitch * inp.channel_height / (inp.helix_pitch + inp.channel_height)
        area = inp.helix_number * inp.helix_pitch * inp.channel_height
    else:
        hydraulic_diameter = 2 * inp.channel_height
        area = np.pi * (radius + inp.inner_wall_thickness + inp.channel_height) ** 2 - np.pi * (
            radius + inp.inner_wall_thickness) ** 2

    density = fuel_density(fuel_temp)
    viscosity = fuel_viscosity(fuel_temp)
    conductivity = fuel_conductivity(fuel_temp)

    velocity = inp.fuel_flow_rate / (area * density)
    re = density * velocity * hydraulic_diameter / viscosity
    pr = viscosity * fuel_specific_heat(fuel_temp) / conductivity
    return fuel_temp, density, viscosity, conductivity, velocity, re, pr, hydraulic_diameter


def stress(radius, heat_flux):
    # TODO shell buckling
    hoop_stress = (inp.coolant_pressure - inp.chamber_pressure) * radius / inp.inner_wall_thickness
    thermal_stress = inp.elastic_modulus * inp.thermal_expansion_coeff * heat_flux * \
                     inp.inner_wall_thickness / (2 * (1 - inp.poissons_ratio) * inp.wall_thermal_conductivity)
    return hoop_stress + thermal_stress


def pressure_drop(radius, fuel_props):
    fuel_temp, density, viscosity, conductivity, velocity, re, pr, hydraulic_diameter = fuel_props

    if inp.helical_channels is True:
        delta_l = geom.station_width * \
                  (1 + 2 * np.pi * (radius + inp.inner_wall_thickness + inp.channel_height / 2) / (inp.helix_number * inp.helix_pitch))
    else:
        delta_l = geom.station_width

    a = -2 * np.log10(inp.surface_roughness / hydraulic_diameter / 3.7 + 12 / re)
    b = -2 * np.log10(inp.surface_roughness / hydraulic_diameter / 3.7 + 2.51 * a / re)
    c = -2 * np.log10(inp.surface_roughness / hydraulic_diameter / 3.7 + 2.51 * b / re)
    f = (a - (b - a) ** 2 / (c - 2 * b + a)) ** (-2)
    drop = f * density * delta_l * velocity ** 2 / (2 * hydraulic_diameter)
    return np.sum(drop)


def main():

    position = np.arange(0, inp.chamber_length, geom.station_width, dtype=np.double)

    radius = geom.radius(position)
    heat_flux = 0.6 * hf.heat_flux(position, 700)

    fuel_props = calc_fuel_props(radius, heat_flux)
    coolant_wall_temp, gas_wall_temp = wall_temp(radius, fuel_props, heat_flux)
    #wall_stress = stress(radius, heat_flux)

    plt.plot(position, coolant_wall_temp, color='black', label='Coolant Wall')
    plt.plot(position, gas_wall_temp, color='black', label='Gas Wall', linestyle=':')
    plt.plot(position, fuel_props[0], color='black', label='Fuel', linestyle='--')

    plt.scatter([0, 0.26, 0.11], [611, 654, 600], color='red', label='High ANSYS Coolant Wall')
    plt.scatter([0, 0.26, 0.11], [554, 625, 560.5], color='blue', label='Low ANSYS Coolant Wall')

    plt.xlabel('Axial Position (m)')
    plt.ylabel('Temperature (K)')
    plt.tight_layout(pad=0.5)
    plt.legend()
    sns.despine()
    #plt.savefig('Temperatures_3mm.png')

    print("Pressure drop = " + str(pressure_drop(radius, fuel_props) / 1e3) + " kPa")

    plt.show()

    plt.plot(position, fuel_props[4])
    plt.ylabel('Velocity (m/s)')
    plt.show()


if __name__ == "__main__":
    main()
