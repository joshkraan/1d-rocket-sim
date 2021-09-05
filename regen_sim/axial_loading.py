import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

import inputs as inp
import engine_geometry as geom
import gas_properties
from gas_properties import calc_gas_properties


def main():
    position = np.arange(0, geom.diverging_end, geom.station_width, dtype=np.double)
    gas, states, velocity = calc_gas_properties(position)

    infinitesimal_forces = states.P * 2 * np.pi * geom.radius(position) * geom.station_width

    injector_force = gas.P * np.pi * geom.chamber_radius ** 2
    converging_force = np.sum(infinitesimal_forces[int(geom.converging_start / geom.station_width):int(
        geom.throat_start / geom.station_width)]) * np.sin(np.radians(inp.converging_half_angle))
    diverging_force = np.sum(infinitesimal_forces[int(geom.throat_end / geom.station_width):int(
        geom.diverging_end/ geom.station_width)]) * np.sin(np.radians(inp.diverging_half_angle))

    print(injector_force)
    print(converging_force)
    print(diverging_force)
    print(injector_force + diverging_force - converging_force)


if __name__ == "__main__":
    main()
