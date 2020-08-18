import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

import downsized2 as inp

sns.set()
sns.set_style("ticks")

div_ang_rad = math.radians(inp.diverging_half_angle)
conv_ang_rad = math.radians(inp.converging_half_angle)

chamber_radius = inp.chamber_diameter / 2
throat_radius = inp.throat_diameter / 2

converging_start = inp.chamber_length + inp.chamber_bevel_radius * math.sin(conv_ang_rad)
throat_start = converging_start + \
               (chamber_radius + inp.chamber_bevel_radius * math.cos(conv_ang_rad) - inp.chamber_bevel_radius -
                throat_radius - inp.throat_bevel_radius + inp.throat_bevel_radius * math.cos(conv_ang_rad)) / math.tan(conv_ang_rad)
throat_end = throat_start + inp.throat_bevel_radius * (math.sin(div_ang_rad) + math.sin(conv_ang_rad))
diverging_end = throat_end + \
                (inp.exit_diameter / 2 -
                 (throat_radius + inp.throat_bevel_radius - inp.throat_bevel_radius *
                  math.cos(div_ang_rad))) / math.tan(div_ang_rad)

throat_position = inp.chamber_length + inp.chamber_bevel_radius * math.sin(conv_ang_rad) + \
                  (chamber_radius + inp.chamber_bevel_radius * math.cos(conv_ang_rad) - inp.chamber_bevel_radius -
                   throat_radius - inp.throat_bevel_radius + inp.throat_bevel_radius * math.cos(conv_ang_rad)) / \
                  math.tan(conv_ang_rad) + inp.throat_bevel_radius * math.sin(conv_ang_rad)

station_width = diverging_end / inp.num_stations


def radius(pos):
    conds = [pos <= inp.chamber_length, (pos > inp.chamber_length) & (pos <= converging_start),
             (pos > converging_start) & (pos <= throat_start), (pos > throat_start) & (pos <= throat_end),
             (pos > throat_end) & (pos <= diverging_end)]
    funcs = [chamber_radius,
             lambda x: np.sqrt(inp.chamber_bevel_radius ** 2 - (x - inp.chamber_length) ** 2) + chamber_radius - inp.chamber_bevel_radius,
             lambda x: - math.tan(conv_ang_rad) * (x - inp.chamber_length - inp.chamber_bevel_radius * math.sin(conv_ang_rad)) + chamber_radius + inp.chamber_bevel_radius * math.cos(conv_ang_rad) - inp.chamber_bevel_radius,
             lambda x: - np.sqrt(inp.throat_bevel_radius ** 2 - (x - throat_start - inp.throat_bevel_radius * math.sin(conv_ang_rad)) ** 2) + throat_radius + inp.throat_bevel_radius,
             lambda x: math.tan(div_ang_rad) * (x - throat_end) + throat_radius + inp.throat_bevel_radius - inp.throat_bevel_radius * math.cos(div_ang_rad)]

    return np.piecewise(pos, conds, funcs)


def main():
    position = np.linspace(0, diverging_end, inp.num_stations, dtype=np.double)
    plt.plot(position, radius(position))
    plt.show()


if __name__ == "__main__":
    main()
