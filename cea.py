import cantera as ct
import numpy as np

gas1 = ct.Solution('lox_kero.cti')
gas1.basis = "mass"
of = 2.29
fuel_fraction = 1 / (1 + of)
ox_fraction = 1 - fuel_fraction

gas1.TPY = 195, 7.8888e5, 'RP-1:' + str(fuel_fraction) + ',' + 'O2(L):' + str(ox_fraction)
#gas1.HPY = average_specific_enthalpy, 7.8888e5, 'O:' + str(ox_fraction) + ',' + 'C:' + str(carbon_fraction) + ',' + 'H:' + str(hydrogen_fraction)
print(gas1.T)
gas1.equilibrate('HP')
print(gas1())

