import numpy as np

generate_data_file = True

plot_data = False


# Read data from Minton and Bessire and create csv file

# Import `os`
import os

# Retrieve current working directory (`cwd`)
cwd = os.getcwd()

# Change directory
os.chdir("C:/Users/Joffrey/Documents/University/References/Minton/TGA")

# List all files and directories in current directory
# (os.listdir('.')


import pandas as pd

# Assign spreadsheet filename to `file`
file = "Molar Yield Table 6.1 C.xlsx"

# Load spreadsheet
xl = pd.ExcelFile(file)

# Print the sheet names
print(xl.sheet_names)

# Load a sheet into a DataFrame by name: df1
df1 = xl.parse("Sheet1")

# List of all indices
print(df1.columns)

# Apply averaged std for H2
meanstd_val = np.mean(df1["0.2"].values)
for i, val in enumerate(df1["0.2"].values):
    if val <= 0.0:
        df1["0.2"][i] = meanstd_val
# df1['0.2'][:] = np.mean(df1['0.2'].values)


# Create the data list to be put in csv
data = {
    "T": df1[0] + 273.4,
    "std_T": df1["0.1"],
    "H2": df1["H2"],
    "std_H2": df1["0.2"],
    "CH4": df1["CH4"],
    "std_CH4": df1["0.3"],
    "CO": df1["CO"],
    "std_CO": df1["0.4"],
    "CO2": df1["CO2"],
    "std_CO2": df1["0.5"],
    "Phenol": df1["Phenol"],
    "std_Phenol": df1["0.6"],
    "H2O": df1["H2O"],
    "std_H2O": df1["0.7"],
    "T_act": 520,
}

print(data["T_act"])
if generate_data_file:
    df = pd.DataFrame(data, columns=data.keys())

    df.to_csv(cwd + "/data_file.csv")


if plot_data:

    import pybit
    import matplotlib.pyplot as plt

    os.chdir(cwd)

    # Colors
    lineColor = [["C0"], ["C1"], ["C2"], ["C3"], ["C4"], ["C5"], ["C6"], ["C7"]]

    data_keys = ["H2", "CO"]
    for i, keys in enumerate(data_keys):
        # plt.figure(1)
        (line,) = plt.plot(
            data["T"].values, data[keys].values, "-o", color=lineColor[i][0], mfc="none", label=keys
        )
        pybit.post_process.error_bar(
            data["T"].values, data[keys].values, data["std_" + keys].values, lineColor[i][0]
        )

    # plt.legend()
    pybit.post_process.saveToTikz("bessire_exp.tex")
    plt.show()
