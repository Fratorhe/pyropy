import json

from jsmin import jsmin


def replace_results(
    vector, param_names, filename_template, filename_out, symbolleft="(", symbolright=")"
):
    """
    Replaces parameter placeholders in a template file with values from a vector.

    Parameters
    ----------
    vector : list or numpy.ndarray
        A list or array of values to replace the parameters with.
    param_names : list of str
        A list of parameter names corresponding to the values in `vector`.
    filename_template : str
        The path to the template file.
    filename_out : str
        The path to the output file where the replaced content will be written.
    symbolleft : str, optional
        The left symbol surrounding the parameter names in the template. Default is "(".
    symbolright : str, optional
        The right symbol surrounding the parameter names in the template. Default is ")".

    Returns
    -------
    None
    """
    with open(filename_template, "r") as fileread, open(filename_out, "w") as filewrite:
        for line in fileread:
            line_new = line
            for val, param in zip(vector, param_names):
                line_new = line_new.replace(symbolleft + param + symbolright, str(val))
            filewrite.write(line_new)


class ReactManager:
    """
    Manages reaction-related files and folders.
    """

    def __init__(self, filename: str = "", folder: str = "./"):
        """
        Initializes the ReactManager.

        Parameters
        ----------
        filename : str, optional
            The filename associated with the reactions. Default is "".
        folder : str, optional
            The folder path where the reaction files are located. Default is "./".
        """
        self.filename = filename
        self.folder = folder
        self.dict_params = dict()

    def react_reader(self):
        """
        Reads the reaction processes from the input file and processes them.

        This function extracts and processes solid reactants and products,
        and generates a list of unique gases involved in the reactions.
        """
        # Read and parse the JSON data from the input file
        with open(self.folder + self.filename) as f:
            data = json.loads(jsmin(f.read()))

        self.solids = data["solids"]
        self.reactions = data["reactions"]
        self.n_reactions = len(self.reactions)
        self.rhoIni = data["rhoIni"]

        # Initialize lists to store reactants, products, and gases
        reactants, rhsList, solid_product, gases_product = [], [], [], []

        # Process reactions to extract reactants and products
        for reaction in self.reactions:
            for key, rhs in reaction.items():
                if key in self.solids:
                    reactants.append(key)
                    rhsList.append(rhs)

        # Process rhsList to separate solids and gases
        for rhs in rhsList:
            rhsSplit = rhs.strip().replace(" ", "").split("+")
            gases = [product for product in rhsSplit if product not in self.solids]
            solid_product.extend(
                [product for product in rhsSplit if product in self.solids]
            )
            gases_product.append(gases)

        # Store unique gases and other attributes
        self.unique_gases = list(
            set([item for sublist in gases_product for item in sublist])
        )
        self.solid_product = solid_product
        self.solid_reactant = reactants
        self.gas_product = gases_product
        self.n_solids = len(self.solids)
        self.rhs = rhsList

    def param_reader(self):
        """Read the parameter values after the reaction reader"""
        # Read and parse the JSON data from the input file
        with open(self.folder + self.filename) as f:
            data = json.load(f)

        parameters = data["parameters"]
        self.param_names = parameters[0].keys()

        # Populate the parameter dictionary with values for each parameter name
        for react in parameters:
            for key, value in react.items():
                self.dict_params.setdefault(key, []).append(value)

        # If 'g' is a key in the parameter dictionary, compute complementary values
        if "g" in self.dict_params:
            self.g_sol = [1 - sum(g_react) for g_react in self.dict_params["g"]]

    def react_writer(self, filename):
        """
        Writes reaction data, including parameters and reactions, to a JSON file.

        Constructs a dictionary with the reaction data and writes it to the specified file.

        Parameters
        ----------
        filename : str
        The path to the output file where the data will be written.
        """
        # Rebuild the reactions and parameters
        reactions_rebuilt = [
            {self.solid_reactant[idx]: self.rhs[idx]} for idx in range(self.n_reactions)
        ]
        param_rebuilt = [
            {
                parameter: self.dict_params[parameter][idx]
                for parameter in self.param_names
            }
            for idx in range(self.n_reactions)
        ]

        # Create the data dictionary
        data = {
            "rhoIni": self.rhoIni,
            "solids": self.solids,
            "reactions": reactions_rebuilt,
            "parameters": param_rebuilt,
        }

        # Write the data to the file
        with open(filename, "w") as outfile:
            json.dump(data, outfile, indent=2)
