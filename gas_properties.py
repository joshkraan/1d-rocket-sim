import numpy as np
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
import cantera as ct

import downsized2 as inp
import engine_geometry as geom

sns.set()
sns.set_style("ticks")


def calc_gas_properties(pos, lox_flow_rate=inp.lox_flow_rate, fuel_flow_rate=inp.fuel_flow_rate):
    stagnation_gas, throat_pressure = solve_chamber_throat(lox_flow_rate, fuel_flow_rate)

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    gas.state = stagnation_gas.state
    states = ct.SolutionArray(gas)

    gamma = stagnation_gas.cp / stagnation_gas.cv
    exit_mach = mach_number(gamma, geom.diverging_end)
    estimated_exit_pressure = stagnation_gas.P * (1 + (gamma - 1) / 2 * exit_mach ** 2) ** (- gamma / (gamma - 1))

    for i in range(pos.size):

        def residual(p):
            gas.SP = stagnation_gas.entropy_mass, p
            if inp.equilibrium_expansion:
                gas.equilibrate('SP')
            vel = (fuel_flow_rate + lox_flow_rate) / (np.pi * geom.radius(pos[i]) ** 2 * gas.density)
            h0 = gas.enthalpy_mass + vel ** 2 / 2
            return h0 - stagnation_gas.enthalpy_mass

        upper_bound = stagnation_gas.P if pos[i] < geom.throat_position else throat_pressure
        lower_bound = throat_pressure if pos[i] < geom.throat_position else 0.9 * estimated_exit_pressure
        pressure = float(scipy.optimize.root_scalar(residual, bracket=[lower_bound, upper_bound]).root)

        gas.SP = stagnation_gas.entropy_mass, pressure
        if inp.equilibrium_expansion:
            gas.equilibrate('SP')

        states.append(gas.state)

    velocity = (fuel_flow_rate + lox_flow_rate) / (np.pi * geom.radius(pos) ** 2 * states.density)

    return stagnation_gas, states, velocity


def calc_mole_fractions(pos, of):

    stagnation_gas = ct.Solution('lox_kero.cti')
    stagnation_gas.basis = "mass"
    stagnation_states = ct.SolutionArray(stagnation_gas)

    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"
    states = ct.SolutionArray(gas)

    for i in range(of.size):
        mass_flow = inp.lox_flow_rate + inp.fuel_flow_rate
        fuel_fraction = 1 / (1 + of[i])
        ox_fraction = 1 - fuel_fraction
        gas1, throat_pressure = solve_chamber_throat(mass_flow * ox_fraction, mass_flow * fuel_fraction)
        stagnation_gas.state = gas1.state
        stagnation_states.append(stagnation_gas.state)

        gas.state = stagnation_gas.state

        gamma = stagnation_gas.cp / stagnation_gas.cv
        exit_mach = mach_number(gamma, geom.diverging_end)
        estimated_exit_pressure = stagnation_gas.P * (1 + (gamma - 1) / 2 * exit_mach ** 2) ** (- gamma / (gamma - 1))

        def residual(p):
            gas.SP = stagnation_gas.entropy_mass, p
            if inp.equilibrium_expansion:
                gas.equilibrate('SP')
            vel = (inp.fuel_flow_rate + inp.lox_flow_rate) / (np.pi * geom.radius(pos) ** 2 * gas.density)
            h0 = gas.enthalpy_mass + vel ** 2 / 2
            return h0 - stagnation_gas.enthalpy_mass

        upper_bound = stagnation_gas.P if pos < geom.throat_position else throat_pressure
        lower_bound = throat_pressure if pos < geom.throat_position else 0.9 * estimated_exit_pressure
        pressure = float(scipy.optimize.root_scalar(residual, bracket=[lower_bound, upper_bound]).root)

        gas.SP = stagnation_gas.entropy_mass, pressure
        if inp.equilibrium_expansion:
            gas.equilibrate('SP')

        states.append(gas.state)

    return stagnation_states, states


def solve_chamber_throat(lox_flow_rate, fuel_flow_rate):
    gas = ct.Solution('lox_kero.cti')
    gas.basis = "mass"

    gas_t = ct.Solution('lox_kero.cti')
    gas_t.basis = "mass"

    of = lox_flow_rate / fuel_flow_rate
    fuel_fraction = 1 / (1 + of)
    ox_fraction = 1 - fuel_fraction

    def stagnation_residual(p):
        gas.TPY = 195, p, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
        gas.equilibrate('HP')

        gas_t.SP = gas.entropy_mass, throat_pressure(stagnation_gas=gas)
        if inp.equilibrium_expansion:
            gas_t.equilibrate('SP')
        return gas_t.enthalpy_mass + sound_speed(gas_t) ** 2 / 2 - gas.enthalpy_mass

    def throat_pressure(stagnation_gas):
        gas_t.state = stagnation_gas.state
        gamma = stagnation_gas.cp / stagnation_gas.cv
        throat_estimate = stagnation_gas.P * ((gamma + 1) / 2) ** (- gamma / (gamma - 1))

        def throat_residual(p):
            gas_t.SP = stagnation_gas.entropy_mass, p
            if inp.equilibrium_expansion:
                gas_t.equilibrate('SP')
            velocity = (fuel_flow_rate + lox_flow_rate) / (np.pi * geom.throat_radius ** 2 * gas_t.density)
            return velocity - sound_speed(gas_t)

        return float(scipy.optimize.root(throat_residual, throat_estimate).x)

    stagnation_pressure = float(scipy.optimize.root(stagnation_residual, inp.chamber_pressure).x)
    gas.TPY = 195, stagnation_pressure, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
    gas.equilibrate('HP')

    return gas, throat_pressure(stagnation_gas=gas)


def sound_speed(gas):
    perturbed_gas = ct.Solution('lox_kero.cti')
    perturbed_gas.basis = "mass"

    perturbed_gas.state = gas.state

    perturbed_gas.SP = gas.entropy_mass, gas.P * 1.0001

    if inp.equilibrium_expansion:
        perturbed_gas.equilibrate('SP')

    return np.sqrt((perturbed_gas.P - gas.P) / (perturbed_gas.density - gas.density))


def mach_number(gamma, position):
    radius = geom.radius(position)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2
    ratio = local_area / throat_area

    def mach_equation(m):
        return 1 / m ** 2 * (2 / (gamma + 1) * (1 + 0.5 * (gamma - 1) * m ** 2)) ** ((gamma + 1) / (gamma - 1)) - ratio ** 2

    if position < geom.throat_position:
        return scipy.optimize.root_scalar(mach_equation, bracket=[1e-6, 1], method='brentq').root
    else:
        return scipy.optimize.root_scalar(mach_equation, bracket=[1, 50], method='brentq').root


def main():
    position = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)

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
    of = np.linspace(1.5, 2.5, 100)
    stagnation_states, states = calc_mole_fractions(geom.diverging_end, of)

    mole_fractions = stagnation_states.X
    plt.plot(of, mole_fractions[:, states.species_index('H')], color='black', label='H')
    plt.plot(of, mole_fractions[:, states.species_index('O')], color='black', linestyle='--', label='O')
    plt.plot(of, mole_fractions[:, states.species_index('OH')], color='black', linestyle='-.', label='OH')
    plt.xlabel('O/F Ratio')
    plt.ylabel('Mole Fraction')
    plt.rcParams["figure.figsize"] = (1, 1)
    sns.despine()
    plt.legend()
    plt.show()
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
