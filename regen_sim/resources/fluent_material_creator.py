import pandas as pd

csv_file = "Kerosene_12E5Pa.csv"

lines = []


def piecewise_writer(column_name):
    df = pd.read_csv(csv_file)
    temp = df['Temperature']
    data = df[column_name]

    piecewise_list = []
    for i in range(len(temp)):
        piecewise_list.append(" (" + str(temp[i]) + ". . " + str(data[i]) + ")")
    return ''.join(piecewise_list)


lines.append("(")
lines.append("\t(kerosene-model-liquid fluid")
lines.append("\t\t(chemical-formula . c12h22<l>)")
lines.append("\t\t(density (polynomial piecewise-linear" + piecewise_writer("Density") + "))")
lines.append("\t\t(specific-heat (polynomial piecewise-linear" + piecewise_writer("Specific Heat") + "))")
lines.append("\t\t(thermal-conductivity (polynomial piecewise-linear" + piecewise_writer("Thermal Conductivity") + "))")
lines.append("\t\t(viscosity (polynomial piecewise-linear" + piecewise_writer("Viscosity") + "))")
lines.append("\t\t(molecular-weight (constant . 167.31))")
lines.append("\t)")
lines.append(")")


with open(csv_file[:-4] + "_Fluent_Material", 'w') as f:
    f.writelines(s + '\n' for s in lines)
