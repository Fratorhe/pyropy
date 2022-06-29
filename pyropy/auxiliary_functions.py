import pandas as pd

def data_to_csv(name, pyrolysis_object):
    file_save = pd.DataFrame(data=None, columns=[], index=None)
    file_save['time'] = pyrolysis_object.time
    file_save['temperature'] = pyrolysis_object.temperature
    file_save['rho'] = pyrolysis_object.rho_solid
    file_save['dRho'] = pyrolysis_object.drho_solid
    file_save.to_csv(name)