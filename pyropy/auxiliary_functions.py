import re

import pandas as pd


def data_to_csv(name: str, pyrolysis_object) -> None:
    """
    Save pyrolysis data to CSV.

    Parameters
    ----------
    name : str
        Output CSV file name.
    pyrolysis_object : object
        Must have `time`, `temperature`, `rho_solid`, and `drho_solid` attributes.

    Returns
    -------
    None
    """
    file_save = pd.DataFrame(data=None, columns=[], index=None)
    file_save["time"] = pyrolysis_object.time
    file_save["temperature"] = pyrolysis_object.temperature
    file_save["rho"] = pyrolysis_object.rho_solid
    file_save["dRho"] = pyrolysis_object.drho_solid
    file_save.to_csv(name)


def write_file_scheme(
        vector: list[float],
        param_names: list[str],
        filename: str,
        symbolleft: str = "(",
        symbolright: str = ")",
        folder: str = "./",
) -> None:
    """
    Replace placeholders in a template and write the result to file.

    Parameters
    ----------
    vector : list of float
        Values to insert.
    param_names : list of str
        Names of placeholders.
    filename : str
        Base filename (template: `filename.template`; output: `folder/filename`).
    symbolleft : str, optional
        Placeholder opening symbol (default: `'('`).
    symbolright : str, optional
        Placeholder closing symbol (default: `')'`).
    folder : str, optional
        Output directory (default: `'./'`).

    Returns
    -------
    None
    """
    with open(filename + ".template", "r") as fileread:
        with open(folder + filename, "w+") as filewrite:
            for line in fileread:
                line_new = line
                for val, param in zip(vector, param_names):
                    line_new = line_new.replace(
                        symbolleft + param + symbolright, str(val)
                    )
                line = line_new
                filewrite.write(line)


def replace_results(
        vector_vals: list[float],
        param_names: list[str],
        filename_template: str,
        filename_out: str,
        symbolleft: str = "(",
        symbolright: str = ")",
) -> None:
    """
    Replace placeholders in a template file and write the result to a new file.

    Parameters
    ----------
    vector_vals : list of float
        Values to insert.
    param_names : list of str
        Placeholder names.
    filename_template : str
        Path to the template file.
    filename_out : str
        Path for the output file.
    symbolleft : str, optional
        Placeholder start symbol (default: `'('`).
    symbolright : str, optional
        Placeholder end symbol (default: `')'`).

    Returns
    -------
    None
    """
    with open(filename_template, "r") as fin, open(filename_out, "w") as fout:
        for line in fin:
            for val, param in zip(vector_vals, param_names):
                line = line.replace(f"{symbolleft}{param}{symbolright}", str(val))
            fout.write(line)


def get_numbers_from_filename(filename: str) -> list[str]:
    """
    Extract numeric values from a filename.

    Parameters
    ----------
    filename : str
        Filename containing numeric values.

    Returns
    -------
    list of str
        All numbers found, as strings (including integers and floats).
    """
    return re.findall(r"[-+]?\d*\.\d+|\d+", filename)
