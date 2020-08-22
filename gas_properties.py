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


def stagnation_properties():
    #TODO maybe return stagnation properties from here?
    #TODO bounds on root finder
    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    of = inp.lox_flow_rate / inp.fuel_flow_rate
    fuel_fraction = 1 / (1 + of)
    ox_fraction = 1 - fuel_fraction

    throat_area = np.pi * geom.throat_radius ** 2
    mass_flow_rate = inp.lox_flow_rate + inp.fuel_flow_rate

    def residual(p):
        gas.TPY = 195, p, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
        gas.equilibrate('HP')
        gas_constant = ct.gas_constant / gas.mean_molecular_weight
        g = gas.cp / gas.cv
        return throat_area * p * np.sqrt(g * (2 / (g + 1)) ** ((g + 1) / (g - 1)) / (gas_constant * gas.T)) - mass_flow_rate

    stagnation_pressure = scipy.optimize.fsolve(residual, inp.chamber_pressure)
    gas.TPY = 195, stagnation_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
    gas.equilibrate('HP')

    return gas, stagnation_pressure


def exit_pressure_estimate(p_0, gamma):
    exit_mach = mach_number(gamma, geom.diverging_end)
    return p_0 * (1 + (gamma - 1) / 2 * exit_mach ** 2) ** (- gamma / (gamma - 1))


def calc_gas_properties(pos):
    stagnation_gas, p_0 = stagnation_properties()

    p_0 = inp.chamber_pressure

    entropy = stagnation_gas.entropy_mass

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    gas.state = stagnation_gas.state
    states = ct.SolutionArray(gas)

    ma = np.zeros(pos.size, dtype=np.float)

    estimated_exit_pressure = exit_pressure_estimate(p_0, stagnation_gas.cp / stagnation_gas.cv)

    for i in range(pos.size):

        def residual(p):
            gas.SP = entropy, p
            gas.equilibrate('SP')
            gamma = gas.cp / gas.cv
            mach = mach_number(gamma, pos[i])
            return (1 + (gamma - 1) / 2 * mach ** 2) ** (- gamma / (gamma - 1)) - p / p_0

        pressure = brentq(residual, p_0, 0.9 * estimated_exit_pressure)
        gas.SP = entropy, pressure
        gas.equilibrate('SP')

        states.append(gas.state)
        ma[i] = mach_number(gas.cp / gas.cv, pos[i])

    return stagnation_gas, states, ma


def calc_gas_properties1(pos):
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

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2
    ratio = local_area / throat_area

    #TODO
    exit_pressure = 67730

    def throat_residual(p):
        gas.SP = entropy, p
        gas.equilibrate('SP')
        gamma = gas.cp / gas.cv
        return (2 / (gamma + 1)) ** (gamma / (gamma - 1)) - p / p_0

    throat_pressure = brentq(throat_residual, p_0, exit_pressure)

    for i in range(pos.size):

        def area_ratio(g, p):
            if pos[i] < geom.throat_position:
                return (2 / (g + 1) * (p_0 / p) ** ((g-1)/g)) ** ((g + 1) / (2 * g - 2)) / np.sqrt(2 / (g - 1) * ((p_0 / p) ** ((g - 1) / g) - 1))
            else:
                return (2 / (g + 1)) ** (1 / (g - 1)) * (p_0 / p) ** (1 / g) / np.sqrt((g + 1) / (g - 1) * (1 - (p / p_0) ** ((g - 1) / g)))

        def residual(p):
            gas.SP = entropy, p
            gas.equilibrate('SP')
            gamma = gas.cp / gas.cv
            return area_ratio(gamma, p) - ratio[i]

        upper_bound = p_0 if pos[i] < geom.throat_position else throat_pressure
        lower_bound = throat_pressure if pos[i] < geom.throat_position else exit_pressure
        pressure = brentq(residual, upper_bound, lower_bound)
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
    print(stagnation_properties()[1])
    position = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)
    # t = time.time()
    gas, states, mach = calc_gas_properties(position)
    print(states[999].P)
    # print(time.time() - t)
    # t = time.time()
    # gas, states1, mach = calc_gas_properties1(position)
    # print(time.time() - t)
    #
    # plt.plot(position, states.viscosity)
    # plt.plot(position, states1.viscosity, linestyle=':')
    # plt.show()


    #print(gas.viscosity, states[781].viscosity)


    # pressure = np.linspace(1e6, 30000, 1000)
    # fraction = gas.X
    # gamma = np.zeros(pressure.size)
    # gamma_eq = np.zeros(pressure.size)
    # entropy = gas.entropy_mass
    # for i in range(pressure.size):
    #     gas.SPX = entropy, pressure[i], fraction
    #     gamma[i] = gas.cp / gas.cv
    #     gas.equilibrate('SP')
    #     gamma_eq[i] = gas.cp / gas.cv
    #
    # plt.plot(pressure, gamma)
    # plt.plot(pressure, gamma_eq)
    # plt.show()

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
