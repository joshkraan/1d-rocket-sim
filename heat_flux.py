import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import cantera as ct

import downsized2 as inp
import engine_geometry as geom
import gas_properties
from gas_properties import calc_gas_properties
from main import calc_fuel_temp


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


def bartz(pos, gas, velocity, gas_wall_temp):

    radius = geom.radius(pos)
    local_area = np.pi * radius ** 2
    throat_area = np.pi * geom.throat_radius ** 2

    mass_flow_rate = inp.lox_flow_rate + inp.fuel_flow_rate
    throat_velocity = mass_flow_rate / (throat_area * gas.density)
    re = gas.density * throat_velocity * inp.throat_diameter / gas.viscosity
    pr = gas.cp * gas.viscosity / gas.thermal_conductivity

    mach_number = velocity / gas_properties.sound_speed(gas)

    m = 0.75
    gamma = gas.cp / gas.cv
    g = 1 + (gamma - 1) / 2 * mach_number ** 2
    sigma = 1 / ((0.5 * (gas_wall_temp / gas.T) * g + 0.5) ** (0.8 - m / 5) * g ** (m / 5))

    nu = 0.026 * re ** 0.8 * pr ** 0.4 * (inp.throat_diameter / inp.throat_bevel_radius) ** 0.1 * (throat_area / local_area) ** 0.9 * sigma
    htc = nu * gas.thermal_conductivity / inp.throat_diameter
    return htc


def bartz_free_stream(pos, gas, states, velocity, gas_wall_temp):

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


def dittus_boelter(pos, states, velocity, gas_wall_temp):
    am_gas = ct.Solution("lox_kero.cti")
    am_gas.basis = "mass"
    am_states = ct.SolutionArray(am_gas, pos.size)
    am_states.TPX = (states.T + gas_wall_temp) / 2, states.P, states.X

    re = am_states.density * velocity * 2 * geom.radius(pos) / am_states.viscosity
    pr = am_states.cp * am_states.viscosity / am_states.thermal_conductivity

    nu = 0.023 * re ** 0.8 * pr ** 0.4
    htc = nu * am_states.thermal_conductivity / (2 * geom.radius(pos))
    return htc


def sieder_tate(pos, states, velocity, gas_wall_temp):

    wall_gas = ct.Solution("lox_kero.cti")
    wall_gas.basis = "mass"
    wall_states = ct.SolutionArray(wall_gas, pos.size)
    wall_states.TPX = gas_wall_temp, states.P, states.X

    re = states.density * velocity * 2 * geom.radius(pos) / states.viscosity
    pr = states.cp * states.viscosity / states.thermal_conductivity

    nu = 0.027 * re ** 0.8 * pr ** (1/3) * (states.viscosity / wall_states.viscosity) ** 0.14

    htc = nu * states.thermal_conductivity / (2 * geom.radius(pos))
    return htc


#TODO fix infinity at 0
def colburn(pos, states, velocity, gas_wall_temp, effective_position):

    am_gas = ct.Solution("lox_kero.cti")
    am_gas.basis = "mass"
    am_states = ct.SolutionArray(am_gas, pos.size)
    am_states.TPX = (states.T + gas_wall_temp) / 2, states.P, states.X

    if effective_position:
        diameter = 2 * geom.radius(pos)
        x_eff = 3.53 * diameter * (1 + (pos / (3.53 * diameter)) ** (-1.2)) ** (-1 / 1.2)
    else:
        x_eff = pos

    re = am_states.density * velocity * x_eff / am_states.viscosity
    pr = am_states.cp * am_states.viscosity / am_states.thermal_conductivity

    nu = 0.0296 * re ** 0.8 * pr ** (1/3)

    htc = nu * am_states.thermal_conductivity / x_eff

    return htc


def scaled_heat_flux(pos):
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
    #gas, states, velocity = calc_gas_properties(position)

    # aw_cea = cea_aw_temp(position)
    # aw = aw_temp(gas, states)
    #
    # gas_wall_temp = 300
    #
    # # col = colburn(position, states, velocity, gas_wall_temp, False) * (aw - gas_wall_temp)
    # # col_eff = colburn(position, states, velocity, gas_wall_temp, True) * (aw - gas_wall_temp)
    # dittus = dittus_boelter(position, states, velocity, gas_wall_temp) * (aw - gas_wall_temp)
    # sieder = sieder_tate(position, states, velocity, gas_wall_temp) * (aw - gas_wall_temp)
    # bartz_cea = cea_bartz(position, gas_wall_temp) * (aw_cea - gas_wall_temp)
    # bartz_free = bartz_free_stream(position, gas, states, velocity, gas_wall_temp) * (aw - gas_wall_temp)
    # bartz0 = bartz(position, gas, velocity, gas_wall_temp) * (aw - gas_wall_temp)

    # #plt.plot(position, col, color='black', alpha=0.2)
    # #plt.plot(position, col_eff, color='black', alpha=0.2)
    # plt.plot(position, bartz_cea, color='black', linestyle=':', label='Bartz CEA')
    # plt.plot(position, bartz0, color='black', alpha=1, label='Bartz Cantera')
    # plt.plot(position, bartz_free, color='black', alpha=1, linestyle='--', label='Bartz Free Stream')
    # plt.plot(position, dittus, color='red', label='Dittus-Boelter')
    # plt.plot(position, sieder, color='blue', label='Sieder-Tate')
    #
    # plt.xlabel('Axial Position (m)')
    # plt.ylabel('Heat Flux (W/m$^2$)')
    #
    # plt.legend()
    #
    # plt.ylim(2e6, 1.4e7)
    # plt.xlim(0, geom.diverging_end)
    # sns.despine()
    #
    # plt.rcParams["figure.figsize"] = (8, 1)
    # plt.savefig('Heat_Flux.png')
    # plt.show()

    # plt.plot(position, bartz_free, label='2.29', color='black')
    # mass_flow_rate = inp.fuel_flow_rate + inp.lox_flow_rate
    # fuel_fraction = 1 / (1 + 1.8)
    # ox_fraction = 1 - fuel_fraction
    # gas, states, velocity = calc_gas_properties(position, mass_flow_rate * ox_fraction, mass_flow_rate * fuel_fraction)
    # plt.plot(position, bartz_free_stream(position, gas, states, velocity, gas_wall_temp) * (aw_temp(gas, states) - gas_wall_temp), label='1.80', color='black', linestyle='--')
    # fuel_fraction = 1 / (1 + 1.5)
    # ox_fraction = 1 - fuel_fraction
    # gas, states, velocity = calc_gas_properties(position, ox_fraction * mass_flow_rate, fuel_fraction * mass_flow_rate)
    # plt.plot(position, bartz_free_stream(position, gas, states, velocity, gas_wall_temp) * (aw_temp(gas, states) - gas_wall_temp), label='1.50', color='black', linestyle=':')
    # plt.xlabel('Axial Position (m)')
    # plt.ylabel('Heat Flux (W/m$^2$)')
    # plt.legend()
    # plt.rcParams["figure.figsize"] = (1, 1)
    # sns.despine()
    # plt.savefig('OF_Constant.png')
    # plt.show()

    # gas_wall_temp = 300
    #
    # inp.equilibrium_expansion = True
    #
    # gas, states, velocity = calc_gas_properties(position)
    #
    # aw = aw_temp(gas, states)
    # bartz_free = bartz_free_stream(position, gas, states, velocity, gas_wall_temp) * (aw - gas_wall_temp)
    # plt.plot(position, bartz_free, color='black', label='Equilibrium')
    #
    # inp.equilibrium_expansion = False
    #
    # gas, states, velocity = calc_gas_properties(position)
    #
    # aw = aw_temp(gas, states)
    # bartz_free = bartz_free_stream(position, gas, states, velocity, gas_wall_temp) * (aw - gas_wall_temp)
    # plt.plot(position, bartz_free, color='black', linestyle='--', label='Frozen')
    #
    # plt.xlabel('Axial Position (m)')
    # plt.ylabel('Heat Flux (W/m$^2$)')
    # plt.legend()
    # plt.rcParams["figure.figsize"] = (1, 1)
    # sns.despine()
    # plt.savefig('Equilibrium_Frozen.png')
    # plt.show()

    gas, states, velocity = calc_gas_properties(position)

    aw = aw_temp(gas, states)
    plt.plot(position, bartz_free_stream(position, gas, states, velocity, 300) * (aw - 300), color='black', label='300K')
    plt.plot(position, bartz_free_stream(position, gas, states, velocity, 897) * (aw - 897), color='black', linestyle=':', label='897K')

    plt.xlabel('Axial Position (m)')
    plt.ylabel('Heat Flux (W/m$^2$)')
    plt.legend()
    plt.rcParams["figure.figsize"] = (1, 1)
    sns.despine()
    plt.savefig('Wall_Temp.png')
    plt.show()

    # position1 = np.linspace(0, geom.diverging_end, inp.num_stations, dtype=np.double)
    # position = position1[position1 < inp.chamber_length]
    #
    # gas, states, velocity = calc_gas_properties(position)
    #
    # aw = aw_temp(gas, states)
    # heat_flux = bartz_free_stream(position, gas, states, velocity, 300) * (aw - 300)
    # fuel_temp = calc_fuel_temp(geom.radius(position), 1.75 * heat_flux)
    # plt.plot(position, fuel_temp)
    # plt.show()





if __name__ == "__main__":
    main()
