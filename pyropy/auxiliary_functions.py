import re

import pandas as pd


def data_to_csv(name, pyrolysis_object):
    file_save = pd.DataFrame(data=None, columns=[], index=None)
    file_save["time"] = pyrolysis_object.time
    file_save["temperature"] = pyrolysis_object.temperature
    file_save["rho"] = pyrolysis_object.rho_solid
    file_save["dRho"] = pyrolysis_object.drho_solid
    file_save.to_csv(name)


def write_file_scheme(vector, param_names, filename, symbolleft="(", symbolright=")", folder="./"):
    with open(filename + ".template", "r") as fileread:
        with open(folder + filename, "w+") as filewrite:
            for line in fileread:
                line_new = line
                for val, param in zip(vector, param_names):
                    line_new = line_new.replace(symbolleft + param + symbolright, str(val))
                line = line_new
                filewrite.write(line)


def replace_results(vector, param_names, filename_template, filename_out, symbolleft="(", symbolright=")"):
    with open(filename_template, "r") as fileread:
        with open(filename_out, "w+") as filewrite:
            for line in fileread:
                line_new = line
                for val, param in zip(vector, param_names):
                    line_new = line_new.replace(symbolleft + param + symbolright, str(val))
                line = line_new
                filewrite.write(line)


def get_numbers_from_filename(filename):
    """
    This is used to get the heating rate from the filename

    :param filename: str (with the heating rate)
    :return: str
    """
    return re.findall(r"[-+]?\d*\.\d+|\d+", filename)
