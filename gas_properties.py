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


    #TODO remove this
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
            if inp.equilibrium_expansion:
                gas.equilibrate('SP')
            gamma = gas.cp / gas.cv
            mach = mach_number(gamma, pos[i])
            return (1 + (gamma - 1) / 2 * mach ** 2) ** (- gamma / (gamma - 1)) - p / p_0

        pressure = brentq(residual, p_0, 0.9 * estimated_exit_pressure)
        gas.SP = entropy, pressure
        if inp.equilibrium_expansion:
            gas.equilibrate('SP')

        states.append(gas.state)
        ma[i] = mach_number(gas.cp / gas.cv, pos[i])

    return stagnation_gas, states, ma


def calc_gas_properties1(pos):
    stagnation_gas, throat_pressure = solve_chamber_throat()

    p_0 = stagnation_gas.P

    entropy = stagnation_gas.entropy_mass

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    gas.state = stagnation_gas.state
    states = ct.SolutionArray(gas)

    estimated_exit_pressure = exit_pressure_estimate(p_0, stagnation_gas.cp / stagnation_gas.cv)

    average_radius = 0.5 * (geom.radius(pos - geom.station_width / 2) + geom.radius(pos + geom.station_width / 2))

    for i in range(pos.size):

        def residual(p):
            gas.SP = entropy, p
            if inp.equilibrium_expansion:
                gas.equilibrate('SP')
            velocity = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(pos[i]) ** 2 * gas.density)
            h0 = gas.enthalpy_mass + velocity ** 2 / 2
            return h0 - stagnation_gas.enthalpy_mass

        #TODO figure out average radius discontinuity
        upper_bound = p_0 if pos[i] < geom.throat_position else throat_pressure
        lower_bound = throat_pressure if pos[i] < geom.throat_position else 0.9 * estimated_exit_pressure
        pressure = brentq(residual, lower_bound, upper_bound, disp=True)
        gas.SP = entropy, pressure
        if inp.equilibrium_expansion:
            gas.equilibrate('SP')

        states.append(gas.state)

    return stagnation_gas, states


# Function is out of date
def calc_mole_fractions(pos, of):
    stagnation_gas = ct.Solution('lox_kero.cti')
    stagnation_gas.basis = "mass"
    stagnation_states = ct.SolutionArray(stagnation_gas)

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    states = ct.SolutionArray(gas)

    for i in range(of.size):
        fuel_fraction = 1 / (1 + of[i])
        ox_fraction = 1 - fuel_fraction
        stagnation_gas.TPY = 195, inp.chamber_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
        stagnation_gas.equilibrate('HP')
        stagnation_states.append(stagnation_gas.state)

        entropy = stagnation_gas.entropy_mass

        p_0 = inp.chamber_pressure

        gas.state = stagnation_gas.state

        estimated_exit_pressure = exit_pressure_estimate(p_0, gas.cp / gas.cv)

        def residual(p):
            gas.SP = entropy, p
            gas.equilibrate('SP')
            gamma = gas.cp / gas.cv
            mach = mach_number(gamma, pos)
            return (1 + (gamma - 1) / 2 * mach ** 2) ** (- gamma / (gamma - 1)) - p / p_0

        pressure = brentq(residual, p_0, 0.9 * estimated_exit_pressure)
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


def solve_chamber_throat():
    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"

    gas1 = ct.Solution('lox_kero.cti')
    gas1.basis = "mass"

    of = inp.lox_flow_rate / inp.fuel_flow_rate
    fuel_fraction = 1 / (1 + of)
    ox_fraction = 1 - fuel_fraction

    def stagnation_residual(p):
        gas.TPY = 195, p, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
        gas.equilibrate('HP')
        entropy = gas.entropy_mass
        stagnation_enthalpy = gas.enthalpy_mass

        p_throat = brentq(throat_residual, inp.chamber_pressure, 67730, args=(entropy, gas.X))
        gas.SP = entropy, p_throat
        if inp.equilibrium_expansion:
            gas.equilibrate('SP')
        #sound_speed = np.sqrt((gas.cp / gas.cv) * gas.P / gas.density)
        return gas.enthalpy_mass + sound_speed(gas) ** 2 / 2 - stagnation_enthalpy

    def throat_residual(p, entropy, mole_fractions):
        gas1.SPX = entropy, p, mole_fractions
        if inp.equilibrium_expansion:
            gas1.equilibrate('SP')
        velocity = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.throat_radius ** 2 * gas1.density)
        #a = np.sqrt((gas1.cp / gas1.cv) * gas1.P / gas1.density)
        return velocity - sound_speed(gas1)

    stagnation_pressure = scipy.optimize.fsolve(stagnation_residual, inp.chamber_pressure)
    gas.TPY = 195, stagnation_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
    gas.equilibrate('HP')
    #TODO exit pressure (also above)
    throat_pressure = brentq(throat_residual, inp.chamber_pressure, 67730, args=(gas.entropy_mass, gas.X))
    return gas, throat_pressure


def sound_speed(gas):
    perturbed_gas = ct.Solution('lox_kero.cti')
    perturbed_gas.basis = "mass"

    perturbed_gas.state = gas.state

    perturbed_gas.SP = gas.entropy_mass, gas.P * 1.0001

    if inp.equilibrium_expansion:
        perturbed_gas.equilibrate('SP')

    return np.sqrt((perturbed_gas.P - gas.P) / (perturbed_gas.density - gas.density))


def main():
    # print(solve_chamber_throat()[0].P, solve_chamber_throat()[1])
    # gas = solve_chamber_throat()[0]
    # enthalpy = gas.enthalpy_mass
    # gas.SP = gas.entropy_mass, solve_chamber_throat()[1]
    # gas.equilibrate('SP')
    # print(np.sqrt((gas.cp / gas.cv) * gas.P / gas.density))
    # print((inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.throat_radius ** 2 * gas.density))
    # print(np.sqrt(2 * (enthalpy - gas.enthalpy_mass)))
    position = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)
    gas, throat_pressure = solve_chamber_throat()
    print(gas.P, throat_pressure)

    # chamber_pressure = gas.P
    # entropy = gas.entropy_mass
    # enthalpy = gas.enthalpy_mass
    #
    # def residual(p):
    #     gas.SP = entropy, p
    #     if inp.equilibrium_expansion:
    #         gas.equilibrate('SP')
    #     vel = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(position[779]) ** 2 * gas.density)
    #     h0 = gas.enthalpy_mass + vel ** 2 / 2
    #     return (h0 ** 2 - enthalpy ** 2) / h0 ** 2
    #
    # pressure = np.linspace(67730, chamber_pressure, 1000, dtype=np.double)
    # result = np.zeros(pressure.size)
    # for i in range(pressure.size):
    #     result[i] = residual(pressure[i])
    # plt.plot(pressure, result)
    # plt.vlines(throat_pressure, -5e4, 5e4)
    # plt.hlines(0, 6773, chamber_pressure)
    # plt.show()

    # throat_pressure = brentq(throat_residual, inp.chamber_pressure, 67730)
    # gas.SP = entropy, throat_pressure
    # gas.equilibrate('SP')
    # print((inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.throat_radius ** 2 * gas.density))
    # print(np.sqrt((gas.cp / gas.cv) * gas.P / gas.density))
    # print(throat_pressure)
    #
    # def residual(p):
    #     gas.SP = entropy, p
    #     if inp.equilibrium_expansion:
    #         gas.equilibrate('SP')
    #     gamma = gas.cp / gas.cv
    #     velocity = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(position[770]) ** 2 * gas.density)
    #     a = np.sqrt((gas.cp / gas.cv) * gas.P / gas.density)
    #     m = velocity / a
    #     return (1 + (gamma - 1) / 2 * m ** 2) ** (- gamma / (gamma - 1)) - p / inp.chamber_pressure
    #
    # pressure = np.linspace(67730, inp.chamber_pressure, 1000, dtype=np.double)
    # result = np.zeros(pressure.size)
    # for i in range(pressure.size):
    #     result[i] = residual(pressure[i])
    # plt.plot(pressure, result)
    # plt.show()
    #
    # # t = time.time()
    #gas, states, mach = calc_gas_properties(position)

    gas1, states1 = calc_gas_properties1(position)
    velocity = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(position) ** 2 * states1.density)
    re = states1.density * velocity * 2 * geom.radius(position) / states1.viscosity
    plt.vlines(geom.throat_position, 200, 500)
    plt.plot(position, velocity)
    #plt.hlines(throat_pressure, position[0], position[999])
    plt.show()


    # #velocity_mass_1 = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(position) ** 2 * states.density)
    # velocity_mass_2 = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(position) ** 2 * states1.density)
    # #velocity_energy_1 = np.sqrt(2 * (gas.enthalpy_mass - states.enthalpy_mass))
    # velocity_energy_2 = np.sqrt(2 * (gas1.enthalpy_mass - states1.enthalpy_mass))
    # #
    # #plt.plot(position, velocity_mass_1)
    # plt.plot(position, velocity_mass_2)
    # #plt.plot(position, velocity_energy_1, linestyle='--')
    # plt.plot(position, velocity_energy_2, linestyle=':')
    #
    # plt.show()

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
