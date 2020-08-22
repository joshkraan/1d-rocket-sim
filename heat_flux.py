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


def cea_bartz(pos):
    gas_wall_temp = 500  # Conservative estimate

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


def bartz(pos, gas, mach_number):

    gas_wall_temp = 500  # Conservative estimate

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


def bartz_free_stream(pos, gas, states, mach):

    gas_wall_temp = 500  # Conservative estimate
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


def dittus_boelter(pos, states, mach):
    gas_wall_temp = 500  # Conservative estimate
    am_gas = ct.Solution("lox_kero.cti")
    am_gas.basis = "mass"
    am_states = ct.SolutionArray(am_gas, pos.size)
    am_states.TPX = states.TPX
    am_states.TP = (states.T + gas_wall_temp) / 2, states.P

    velocity = mach * np.sqrt((states.cp / states.cv) * states.P / states.density)

    re = am_states.density * velocity * 2 * geom.radius(pos) / am_states.viscosity
    pr = am_states.cp * am_states.viscosity / am_states.thermal_conductivity

    nu = 0.023 * re ** 0.8 * pr ** 0.4
    htc = nu * am_states.thermal_conductivity / (2 * geom.radius(pos))
    return htc


def sieder_tate(pos, states, mach):
    #gas, states, mach = calc_gas_properties(pos)

    gas_wall_temp = 500  # Conservative estimate

    temperature, pressure = states.TP

    diameter = 2 * geom.radius(pos)

    gamma = states.cp / states.cv
    density = states.density
    viscosity = states.viscosity
    thermal_cond = states.thermal_conductivity

    velocity = mach * np.sqrt(gamma * pressure / density)

    re = density * velocity * diameter / states.viscosity
    pr = states.cp * states.viscosity / thermal_cond

    states.TP = gas_wall_temp, pressure

    nu = 0.027 * re ** 0.8 * pr ** (1/3) * (viscosity / states.viscosity) ** 0.14

    htc = nu * thermal_cond / diameter

    states.TP = temperature, pressure

    return htc


def colburn(pos, states, mach):
    #gas, states, mach = calc_gas_properties(pos)

    gas_wall_temp = 500  # Conservative estimate

    temperature, pressure = states.TP
    t_mean = (gas_wall_temp + temperature) / 2

    diameter = 2 * geom.radius(pos)
    x_eff = 3.53 * diameter * (1 + (pos / (3.53 * diameter)) ** (-1.2)) ** (-1/1.2)
    #x_eff = pos

    gamma = states.cp / states.cv
    density = states.density

    velocity = mach * np.sqrt(gamma * pressure / density)

    states.TP = t_mean, pressure

    re = states.density * velocity * x_eff / states.viscosity
    pr = states.cp * states.viscosity / states.thermal_conductivity

    nu = 0.0296 * re ** 0.8 * pr ** (1/3)

    htc = nu * states.thermal_conductivity / x_eff

    states.TP = temperature, pressure

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
    pr = states.cp * states.viscosity / states.thermal_conductivity
    gamma = states.cp / states.cv

    # aw_cea = cea_aw_temp(position)
    # aw = aw_temp(gas, states)
    #
    # plt.plot(position, aw_cea)
    # plt.plot(position, aw, linestyle=':')
    # plt.show()

    # col = colburn(position, states, mach)
    dittus = dittus_boelter(position, states, mach)
    # sieder = sieder_tate(position, states, mach)
    bartz0 = cea_bartz(position)
    bartz1 = bartz_free_stream(position, gas, states, mach)
    bartz2 = bartz(position, gas, mach)
    # plt.plot(position, col, linestyle=':')
    plt.plot(position, bartz2, linestyle='-.')
    plt.plot(position, bartz1, linestyle=':')
    plt.plot(position, bartz0)
    plt.plot(position, dittus)
    # plt.plot(position, sieder, linestyle='--')
    plt.show()


if __name__ == "__main__":
    main()
