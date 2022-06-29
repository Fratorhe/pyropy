import json

from jsmin import jsmin


class ReactManager:
    def __init__(self, filename: str = None, folder: str = "./"):
        self.filename = filename
        self.folder = folder

    def react_reader(self):
        """Read the reaction processes from input file"""
        self.dict_params = dict()
        with open(self.folder + self.filename) as f:
            minified = jsmin(f.read())
            data = json.loads(minified)
        self.solids = data["solids"]
        # self.gas = self.data["gas"]
        self.reactions = data["reactions"]
        self.n_reactions = len(self.reactions)
        self.rhoIni = data["rhoIni"]
        reactants = []
        rhsList = []  # this is the list of all the products of one reaction (right hand side)

        # Build the list of reactants and products (contained in rhs)
        for reaction in self.reactions:
            keys = reaction.keys()
            for key in keys:
                if key in self.solids:
                    reactants.append(key)
                    rhsList.append(reaction[key])

        self.rhs = rhsList
        products = []
        solid_product = []
        gases_product = []

        for rhs in rhsList:
            rhsSplit = rhs.strip().replace(" ", "").split("+")
            gases = rhsSplit
            products.append(rhsSplit)
            for product in rhsSplit:
                if product in self.solids:
                    solid_product.append(product)
                    gases.remove(product)
            gases_product.append(gases)
        self.unique_gases = list(
            set([item for sublist in gases_product for item in sublist])
        )  # flattens list, and gets unique gases to make gas list

        self.solid_product = solid_product
        self.solid_reactant = reactants
        self.gas_product = gases_product
        self.n_solids = len(self.solids)
        self.g_sol = []
        return 0

    def param_reader(self):
        """Read the parameter values after the reaction reader"""
        with open(self.folder + self.filename) as f:
            data = json.load(f)
        parameters = data["parameters"]
        self.param_names = parameters[0].keys()
        for react in parameters:
            for key in self.param_names:
                if key in self.dict_params:
                    self.dict_params[key].append(react[key])
                else:
                    self.dict_params[key] = []
                    self.dict_params[key].append(react[key])

        if "g" in self.dict_params.keys():
            for g_react in self.dict_params["g"]:
                g_sol_temp = sum(g_react)
                self.g_sol.append(1 - g_sol_temp)
        return 0

    def react_writer(self, filename):
        data = {}
        data["rhoIni"] = self.rhoIni
        data["solids"] = self.solids

        reactions_rebuilt = []
        for idx in range(0, self.n_reactions):
            react_dict = {}
            react_dict[self.solid_reactant[idx]] = self.rhs[idx]
            reactions_rebuilt.append(react_dict)
        data["reactions"] = reactions_rebuilt

        param_rebuilt = []
        for idx in range(0, self.n_reactions):
            param_dict = {}
            for parameter in self.param_names:
                param_dict[parameter] = self.dict_params[parameter][idx]
            param_rebuilt.append(param_dict)
        data["parameters"] = param_rebuilt

        with open(filename, "w") as outfile:
            json.dump(data, outfile, indent=2)

    # def __str__(self):
