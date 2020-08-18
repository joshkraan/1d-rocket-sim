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

import downsized2 as inp
import engine_geometry as geom
from gas_properties import calc_gas_properties


def bartz_heat_flux(pos, mach_number):
    gas_wall_temp = 300  # Conservative estimate

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2

    recovery_factor = inp.gas_prandtl_number ** 0.33
    a = 1 + recovery_factor * (inp.gas_gamma - 1) / 2 * mach_number ** 2
    b = 1 + (inp.gas_gamma - 1) / 2 * mach_number ** 2
    adiabatic_wall_temp = inp.chamber_temperature * (a / b)

    g = 1 + (inp.gas_gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / inp.chamber_temperature) * g + 0.5) ** 0.68 * g ** 0.12)
    htc = 0.026 / inp.throat_diameter ** 0.2 * (inp.gas_viscosity ** 0.2 * inp.gas_specific_heat / inp.gas_prandtl_number ** 0.6) * \
          (inp.chamber_pressure / inp.cstar) ** 0.8 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma

    return htc * (adiabatic_wall_temp - gas_wall_temp)


def bartz_heat_flux1(pos, mach_number):
    gas1 = ct.Solution('lox_kero.cti')
    of = inp.lox_flow_rate / inp.fuel_flow_rate
    fuel_fraction = 1 / (1 + of)
    ox_fraction = 1 - fuel_fraction
    gas1.TPY = 195, inp.chamber_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
    gas1.basis = "mass"
    gas1.equilibrate('HP')

    specific_heat = gas1.cp
    viscosity = gas1.viscosity
    gamma = gas1.cp / gas1.cv
    temperature = gas1.T
    prandtl = specific_heat * viscosity / gas1.thermal_conductivity

    gas_wall_temp = 300  # Conservative estimate

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2

    recovery_factor = prandtl ** 0.33
    a = 1 + recovery_factor * (gamma - 1) / 2 * mach_number ** 2
    b = 1 + (gamma - 1) / 2 * mach_number ** 2
    adiabatic_wall_temp = temperature * (a / b)

    g = 1 + (gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / temperature) * g + 0.5) ** 0.68 * g ** 0.12)
    htc = 0.026 / inp.throat_diameter ** 0.2 * (viscosity ** 0.2 * specific_heat / prandtl ** 0.6) * \
          ((inp.lox_flow_rate + inp.fuel_flow_rate) / throat_area) ** 0.8 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma

    #TODO
    return htc


def bartz_heat_flux2(pos, mach_number):
    gas1 = ct.Solution('lox_kero.cti')
    of = inp.lox_flow_rate / inp.fuel_flow_rate
    fuel_fraction = 1 / (1 + of)
    ox_fraction = 1 - fuel_fraction
    gas1.TPY = 195, inp.chamber_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
    gas1.basis = "mass"
    gas1.equilibrate('HP')

    specific_heat = gas1.cp
    viscosity = gas1.viscosity
    gamma = gas1.cp / gas1.cv
    temperature = gas1.T
    pr = specific_heat * viscosity / gas1.thermal_conductivity

    gas_wall_temp = 500  # Conservative estimate

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2

    recovery_factor = pr ** 0.33
    a = 1 + recovery_factor * (gamma - 1) / 2 * mach_number ** 2
    b = 1 + (gamma - 1) / 2 * mach_number ** 2
    adiabatic_wall_temp = temperature * (a / b)

    density = gas1.density
    mass_flow_rate = inp.lox_flow_rate + inp.fuel_flow_rate
    throat_velocity = mass_flow_rate / (throat_area * density)
    re = density * throat_velocity * inp.throat_diameter / viscosity

    g = 1 + (gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / temperature) * g + 0.5) ** 0.68 * g ** 0.12)

    nu = 0.026 * re ** 0.8 * pr ** 0.4 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma
    htc = nu * gas1.thermal_conductivity / inp.throat_diameter
    #TODO
    return nu


def colburn(pos):
    gas, state, mach = calc_gas_properties(pos)

    gas_wall_temp = 500  # Conservative estimate
    t_mean = (gas_wall_temp + state.T) / 2

    diameter = 2 * geom.radius(pos)
    x_eff = 3.53 * diameter * (1 + (pos / (3.53 * diameter)) ** (-1.2)) ** (1/1.2)

def scaled_heat_flux(pos):
    # Experimentally measured parameters for similar engine
    exp_pressure = 2e6  # Pa
    exp_diameter = 5e-2  # m
    exp_heat_flux = 3954e3  # W / m^2
    ind_heat_flux_param = exp_heat_flux * exp_diameter ** 0.1 / exp_pressure ** 0.8
    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2
    throat_flux = ind_heat_flux_param * inp.chamber_pressure ** 0.8 / inp.throat_diameter ** 0.1
    return throat_flux * (throat_area / local_area) ** 0.9


def main():
    position = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)
    gas, states, mach = calc_gas_properties(position)
    bartz1 = bartz_heat_flux1(position, mach)
    bartz2 = bartz_heat_flux2(position, mach)
    plt.plot(position, bartz1, linestyle=':')
    plt.plot(position, bartz2)
    plt.show()


if __name__ == "__main__":
    main()
