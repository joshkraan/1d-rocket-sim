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

sns.set()
sns.set_style("ticks")


def calc_gas_properties(pos):
    # TODO calculate chamber and exit pressure

    stagnation_gas = ct.Solution('lox_kero.cti')
    stagnation_gas.basis = "mass"
    of = inp.lox_flow_rate / inp.fuel_flow_rate
    fuel_fraction = 1 / (1 + of)
    ox_fraction = 1 - fuel_fraction
    stagnation_gas.TPY = 195, inp.chamber_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
    stagnation_gas.equilibrate('HP')

    p_0 = inp.chamber_pressure

    entropy = stagnation_gas.entropy_mass

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    gas.state = stagnation_gas.state
    states = ct.SolutionArray(gas)

    ma = np.zeros(pos.size, dtype=np.float)

    #TODO
    exit_pressure = 67730

    for i in range(pos.size):

        def residual(p):
            gas.SP = entropy, p
            gas.equilibrate('SP')
            gamma = gas.cp / gas.cv
            mach = mach_number(gamma, pos[i])
            return (1 + (gamma - 1) / 2 * mach ** 2) ** (- gamma / (gamma - 1)) - p / p_0

        pressure = brentq(residual, p_0, exit_pressure)
        gas.SP = entropy, pressure
        gas.equilibrate('SP')

        states.append(gas.state)
        ma[i] = mach_number(gas.cp / gas.cv, pos[i])

    return stagnation_gas, states, ma


#TODO chamber mass fractions
def calc_mole_fractions(pos, of):
    stagnation_gas = ct.Solution('lox_kero.cti')
    stagnation_gas.basis = "mass"
    stagnation_states = ct.SolutionArray(stagnation_gas)

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    states = ct.SolutionArray(gas)
    #TODO
    exit_pressure = 67730

    for i in range(of.size):
        fuel_fraction = 1 / (1 + of[i])
        ox_fraction = 1 - fuel_fraction
        stagnation_gas.TPY = 195, inp.chamber_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
        stagnation_gas.equilibrate('HP')
        stagnation_states.append(stagnation_gas.state)

        entropy = stagnation_gas.entropy_mass

        p_0 = inp.chamber_pressure

        gas.state = stagnation_gas.state

        def residual(p):
            gas.SP = entropy, p
            gas.equilibrate('SP')
            gamma = gas.cp / gas.cv
            mach = mach_number(gamma, pos)
            return (1 + (gamma - 1) / 2 * mach ** 2) ** (- gamma / (gamma - 1)) - p / p_0

        pressure = brentq(residual, p_0, exit_pressure)
        gas.SP = entropy, pressure
        gas.equilibrate('SP')

        states.append(gas.state)

    return stagnation_states, states


def mach_number(gamma, position):
    radius = geom.radius(position)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2
    ratio = local_area / throat_area

    def mach_equation(m):
        return 1 / m ** 2 * (2 / (gamma + 1) * (1 + 0.5 * (gamma - 1) * m ** 2)) ** ((gamma + 1) / (gamma - 1)) - ratio ** 2

    if position < geom.throat_position:
        return brentq(mach_equation, 1e-6, 1)
    else:
        return brentq(mach_equation, 1, 50)


def main():
    position = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)
    gas, states, mach = calc_gas_properties(position)
    print(states[781].cp / states[781].cv)

    # gas = calc_gas_properties(position)[0]
    # temperature = np.linspace(500, 4000, 1000)
    # pressure = gas.P
    # viscosity = np.zeros(temperature.size)
    # for i in range(temperature.size):
    #     gas.TP = temperature[i], pressure
    #     viscosity[i] = gas.viscosity
    #
    # plt.plot(temperature, viscosity)
    # m = 0.75
    # plt.plot(temperature, (viscosity[0] / temperature[0] ** m) * temperature ** m)
    # plt.show()
    #
    # of = np.linspace(1.5, 2.5, 100)
    # stagnation_states, states = calc_mole_fractions(geom.diverging_end, of)
    #
    # mole_fractions = stagnation_states.X
    # plt.plot(of, mole_fractions[:, states.species_index('H')], color='black', label='H')
    # plt.plot(of, mole_fractions[:, states.species_index('O')], color='black', linestyle='--', label='O')
    # plt.plot(of, mole_fractions[:, states.species_index('OH')], color='black', linestyle='-.', label='OH')
    # plt.xlabel('O/F Ratio')
    # plt.ylabel('Mole Fraction')
    # plt.rcParams["figure.figsize"] = (1, 1)
    # sns.despine()
    # plt.legend()
    # plt.show()
    #
    # flame_temp = states.T
    # plt.plot(of, flame_temp)
    # plt.show()
    #
    # gamma = states.cp / states.cv
    # plt.plot(of, gamma)
    # plt.show()


if __name__ == "__main__":
    main()
