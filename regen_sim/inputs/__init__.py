import importlib_resources
from regen_sim.inputs.output11 import *

cantera_properties_name = "lox_kero.cti"
num_stations = 1000
equilibrium_expansion = True

fuel_properties_file = str(importlib_resources.files('regen_sim.resources').joinpath(fuel_properties_name))
cantera_properties_file = str(importlib_resources.files('regen_sim.resources').joinpath(cantera_properties_name))
