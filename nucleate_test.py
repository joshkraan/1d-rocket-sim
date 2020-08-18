import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math

density = 800  # kg / m^3
velocity = 10  # m/s
hydraulic_diameter = 1.7e-3  # m
viscosity = 0.0015460  # Pa s
specific_heat = 1967.6  # J / kg / K
thermal_cond = 0.1109  # W / m / K
boiling_point = 564.77  # K
heat_of_vaporization = 246460  # J / kg

surface_tension = 0.025  # N / m very uncertain, need to calculate with REFPROP

fluid_surface_param = 1  # stainless


def main():
    re = density * velocity * hydraulic_diameter / viscosity
    pr = viscosity * specific_heat / thermal_cond
    nu = 0.023 * re ** 0.8 * pr ** 0.4

    mass_flux = density * velocity

    coeff = nu * thermal_cond / hydraulic_diameter

    temp = np.linspace(boiling_point + 1, 600, 100, dtype=np.double)
    tsat = temp - boiling_point

    heat_flux = (1058 * (mass_flux * heat_of_vaporization) ** -0.7 * fluid_surface_param * coeff * tsat) ** (1/0.3)
    heat_flux1 = coeff * (np.linspace(350, 550, 100, dtype=np.double) - 300)

    x = [348.936298803336, 372.0152082898365, 401.2340252519147, 443.8507521730052, 480.9759018933441, 510.8238189949781, 530.3149346724833, 559.5914420403722, 569.2889245409492, 580.2198828610046, 580.1594452930124, 580.0852719141128, 589.7992373877785, 594.6342428271596, 599.4857312396293, 601.8592793644166, 609.1310176587585]
    y = [373564.62718589563, 494103.336152924, 620619.6345531865, 855742.0827902809, 1158698.430127688, 1407259.6305304456, 1601922.883088366, 1804922.580793521, 2108229.532957906, 2412092.6594305555, 2702912.0699141393, 3108176.993969801, 3519508.222041664, 3903440.2382721812, 4196913.451378728, 4801397.336587869, 5408610.69316477]

    kin_vis = viscosity / density
    t_onb = 4 * surface_tension * tsat * kin_vis * coeff / (thermal_cond * heat_of_vaporization) * (1 + np.sqrt(1 + thermal_cond * heat_of_vaporization * (boiling_point - 300) / (2 * surface_tension * tsat * kin_vis * coeff)))
    q_onb = thermal_cond * heat_of_vaporization / (8 * surface_tension * kin_vis * tsat) * t_onb ** 2

    plt.scatter(x, y)
    plt.scatter(temp, q_onb)
    plt.scatter(np.linspace(350, 550, 100, dtype=np.double), heat_flux1)
    # plt.yscale("log")
    # plt.ylim((0.1e6, 10e6))
    # plt.xlim((0, 1400))
    plt.show()


if __name__ == "__main__":
    main()
