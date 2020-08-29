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
import gas_properties
from gas_properties import calc_gas_properties


def aw_temp(gas, states):
    pr = states.cp * states.viscosity / states.thermal_conductivity
    r = pr ** (1/3)
    return states.T + r * (gas.T - states.T)


def cea_aw_temp(pos):
    r = inp.gas_prandtl_number ** 0.33

    mach_number = np.zeros(pos.size)
    for i in range(pos.size):
        mach_number[i] = gas_properties.mach_number(inp.gas_gamma, pos[i])

    a = 1 + r * (inp.gas_gamma - 1) / 2 * mach_number ** 2
    b = 1 + (inp.gas_gamma - 1) / 2 * mach_number ** 2
    return inp.chamber_temperature * (a / b)


def cea_bartz(pos, gas_wall_temp):

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2

    mach_number = np.zeros(pos.size)
    for i in range(pos.size):
        mach_number[i] = gas_properties.mach_number(inp.gas_gamma, pos[i])

    g = 1 + (inp.gas_gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / inp.chamber_temperature) * g + 0.5) ** 0.68 * g ** 0.12)
    htc = 0.026 / inp.throat_diameter ** 0.2 * (inp.gas_viscosity ** 0.2 * inp.gas_specific_heat / inp.gas_prandtl_number ** 0.6) * \
          ((inp.lox_flow_rate + inp.fuel_flow_rate) / throat_area) ** 0.8 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma

    return htc


def bartz(pos, gas, mach_number, gas_wall_temp):

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2

    mass_flow_rate = inp.lox_flow_rate + inp.fuel_flow_rate
    throat_velocity = mass_flow_rate / (throat_area * gas.density)
    re = gas.density * throat_velocity * inp.throat_diameter / gas.viscosity
    pr = gas.cp * gas.viscosity / gas.thermal_conductivity

    m = 0.75
    gamma = gas.cp / gas.cv
    g = 1 + (gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / gas.T) * g + 0.5) ** (0.8 - m / 5) * g ** (m / 5))

    nu = 0.026 * re ** 0.8 * pr ** 0.4 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma
    htc = nu * gas.thermal_conductivity / inp.throat_diameter
    return htc


def bartz_free_stream(pos, gas, states, mach, gas_wall_temp):

    am_gas = ct.Solution("lox_kero.cti")
    am_gas.basis = "mass"
    am_states = ct.SolutionArray(am_gas, pos.size)
    am_states.TPX = states.TPX
    am_states.TP = (states.T + gas_wall_temp) / 2, states.P

    velocity = mach * np.sqrt((states.cp / states.cv) * states.P / states.density)

    re = gas.density * velocity * 2 * geom.radius(pos) / gas.viscosity
    pr = gas.cp * gas.viscosity / gas.thermal_conductivity

    sigma = (am_states.density / gas.density) ** 0.8 * (am_states.viscosity / gas.viscosity) ** 0.2

    nu = 0.026 * re ** 0.8 * pr ** 0.4 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * sigma
    htc = nu * gas.thermal_conductivity / (2 * geom.radius(pos))
    return htc


def bartz_free_stream1(pos, gas, states, velocity, gas_wall_temp):

    am_gas = ct.Solution("lox_kero.cti")
    am_gas.basis = "mass"
    am_states = ct.SolutionArray(am_gas, pos.size)
    am_states.TPX = (states.T + gas_wall_temp) / 2, states.P, states.X

    re = gas.density * velocity * 2 * geom.radius(pos) / gas.viscosity
    pr = gas.cp * gas.viscosity / gas.thermal_conductivity

    sigma = (am_states.density / gas.density) ** 0.8 * (am_states.viscosity / gas.viscosity) ** 0.2

    nu = 0.026 * re ** 0.8 * pr ** 0.4 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * sigma
    htc = nu * gas.thermal_conductivity / (2 * geom.radius(pos))
    return htc


def dittus_boelter(pos, states, mach, gas_wall_temp, prandtl_correction, stagnation_temp):
    pr_free_stream = states.cp * states.viscosity / states.thermal_conductivity

    if prandtl_correction:
        reference_temp = 0.5 * (states.T + gas_wall_temp) + 0.22 * pr_free_stream ** (1/3) * (stagnation_temp - states.T)
    else:
        reference_temp = 0.5 * (states.T + gas_wall_temp)

    ref_gas = ct.Solution("lox_kero.cti")
    ref_gas.basis = "mass"
    ref_states = ct.SolutionArray(ref_gas, pos.size)
    ref_states.TPX = states.TPX
    ref_states.TP = reference_temp, states.P

    velocity = mach * np.sqrt((states.cp / states.cv) * states.P / states.density)

    re = ref_states.density * velocity * 2 * geom.radius(pos) / ref_states.viscosity
    pr = ref_states.cp * ref_states.viscosity / ref_states.thermal_conductivity

    nu = 0.023 * re ** 0.8 * pr ** 0.4
    htc = nu * ref_states.thermal_conductivity / (2 * geom.radius(pos))
    return htc


def sieder_tate(pos, states, mach, gas_wall_temp):

    wall_gas = ct.Solution("lox_kero.cti")
    wall_gas.basis = "mass"
    wall_states = ct.SolutionArray(wall_gas, pos.size)
    wall_states.TPX = states.TPX
    wall_states.TP = gas_wall_temp, states.P

    diameter = 2 * geom.radius(pos)

    velocity = mach * np.sqrt((states.cp / states.cv) * states.P / states.density)

    re = states.density * velocity * diameter / states.viscosity
    pr = states.cp * states.viscosity / states.thermal_conductivity

    nu = 0.027 * re ** 0.8 * pr ** (1/3) * (states.viscosity / wall_states.viscosity) ** 0.14

    htc = nu * states.thermal_conductivity / diameter

    return htc


def sieder_tate1(pos, states, velocity, gas_wall_temp):

    wall_gas = ct.Solution("lox_kero.cti")
    wall_gas.basis = "mass"
    wall_states = ct.SolutionArray(wall_gas, pos.size)
    wall_states.TPX = gas_wall_temp, states.P, states.X

    diameter = 2 * geom.radius(pos)

    re = states.density * velocity * diameter / states.viscosity
    pr = states.cp * states.viscosity / states.thermal_conductivity

    nu = 0.027 * re ** 0.8 * pr ** (1/3) * (states.viscosity / wall_states.viscosity) ** 0.14

    htc = nu * states.thermal_conductivity / diameter

    return htc


def colburn(pos, states, mach, gas_wall_temp, effective_position):

    am_gas = ct.Solution("lox_kero.cti")
    am_gas.basis = "mass"
    am_states = ct.SolutionArray(am_gas, pos.size)
    am_states.TPX = states.TPX
    am_states.TP = (states.T + gas_wall_temp) / 2, states.P

    if effective_position:
        diameter = 2 * geom.radius(pos)
        x_eff = 3.53 * diameter * (1 + (pos / (3.53 * diameter)) ** (-1.2)) ** (-1 / 1.2)
    else:
        x_eff = pos

    velocity = mach * np.sqrt((states.cp / states.cv) * states.P / states.density)

    re = am_states.density * velocity * x_eff / am_states.viscosity
    pr = am_states.cp * am_states.viscosity / am_states.thermal_conductivity

    nu = 0.0296 * re ** 0.8 * pr ** (1/3)

    htc = nu * am_states.thermal_conductivity / x_eff

    return htc


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

    aw_cea = cea_aw_temp(position)
    aw = aw_temp(gas, states)

    gas_wall_temp = 500

    col = colburn(position, states, mach, gas_wall_temp, False) * aw
    col_eff = colburn(position, states, mach, gas_wall_temp, True) * aw
    dittus = dittus_boelter(position, states, mach, gas_wall_temp, False, gas.T) * aw
    dittus_pr = dittus_boelter(position, states, mach, gas_wall_temp, True, gas.T) * aw
    sieder = sieder_tate(position, states, mach, gas_wall_temp) * aw
    bartz_cea = cea_bartz(position, gas_wall_temp) * aw_cea
    bartz_free = bartz_free_stream(position, gas, states, mach, gas_wall_temp) * aw
    bartz0 = bartz(position, gas, mach, gas_wall_temp) * aw

    gas1, states1, velocity = gas_properties.calc_gas_properties1(position)
    sieder_new = sieder_tate1(position, states1, velocity, gas_wall_temp) * aw_temp(gas1, states1)
    bartz_new = bartz_free_stream1(position, gas1, states1, velocity, gas_wall_temp) * aw_temp(gas1, states1)

    plt.plot(position, col, color='black', alpha=0.2)
    plt.plot(position, col_eff, color='black', alpha=0.2)
    plt.plot(position, dittus, color='black', alpha=0.2)
    plt.plot(position, dittus_pr, color='black', alpha=0.2)
    plt.plot(position, sieder, color='black', alpha=0.2)
    plt.plot(position, bartz_cea, color='black', alpha=1, linestyle=':')
    plt.plot(position, bartz_free, color='black', alpha=1, linestyle='--')
    plt.plot(position, bartz0, color='black', alpha=1)

    plt.plot(position, sieder_new, color='red')
    plt.plot(position, bartz_new, color='red')
    plt.vlines(geom.throat_position, 0.4e7, 1.4e7)

    plt.ylim(3e6, 1.4e7)
    sns.despine()
    plt.show()

    plt.fill_between(position, bartz0, dittus, color='black', alpha=0.2)
    plt.show()


if __name__ == "__main__":
    main()
