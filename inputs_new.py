"""
High level inputs
"""
num_stations = 1000

"""
Polynomial fits for fuel
"""
temp_lower_bound = 270
temp_upper_bound = 650
density_coefficients = [-7.35246E-13, 1.89632E-09, -2.01450E-06, 1.12544E-03, -3.48382E-01, 5.59004E+01, -2.75598E+03]
specific_heat_coefficients = [1.40053E-11, -3.66091E-08, 3.93278E-05, -2.22085E-02, 6.94910E+00, -1.13776E+03, 7.76898E+04]
viscosity_coefficients = [2.39318E-17, -7.07636E-14, 8.64557E-11, -5.58968E-08, 2.01936E-05, -3.87381E-03, 3.09736E-01]
thermal_conductivity_coefficients = [6.27385E-17, -1.63841E-13, 1.76034E-10, -9.95363E-08, 3.13833E-05, -5.41772E-03, 5.25984E-01]

enthalpy_lower_bound = -466070
enthalpy_upper_bound = 547770
temperature_coefficients = [-1.10030E-34, 1.57293E-29, -3.53013E-23, 6.01412E-17, -1.08504E-10, 3.70175E-04, 4.75311E+02]

"""
Path to file containing discrete data for fuel properties, which will be interpolated.
Pressure should be constant and temperature varied.
Should be a CSV file with columns titled Temperature, Density, Specific Heat, Thermal Conductivity, Viscosity, Enthalpy.
Standard SI units with no prefixes should be used.
"""
fuel_properties_file = "kerosene_17E5Pa_fit.csv"

"""
Chamber dimensions. All units are meters and degrees.
"""
exit_diameter = 0.211
throat_diameter = 0.0944
diverging_half_angle = 20
converging_half_angle = 30
chamber_length = 0.33
chamber_diameter = 0.172
throat_bevel_radius = 0.05
chamber_bevel_radius = 0.01

"""
Double-wall cooling parameters
"""
inner_wall_thickness = 0.0015875  # m
channel_height = 0.003  # m
wall_thermal_conductivity = 401  # W / (m K)
surface_roughness = 0.001  # mm
fuel_flow_rate = 2.11  # kg/s

# TODO calculate enthalpy from input temperature
input_enthalpy = -408900 #J/s
fuel_input_pressure = 1820000  # Pa
fuel_input_temperature = 300  # k

