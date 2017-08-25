# Program name: atomic1D/reference/build_json.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 14 July 2017
# 
# 
# Makes data_dict and copies it into a .json file 'ADAS/0.2.json'

def replace_guards(var):
	"""
	This in-place replaces the points in the guard cells with the points on the boundary
	
	"""
	var[0] = 0.5*(var[0] + var[1])
	var[-1] = 0.5*(var[-1] + var[-2])

filename = 'ADAS/0.2'

from boutdata.collect import collect
import numpy as np

data_dict = {}

dy = collect("dy", path=filename, yguards=True)[-1,1:-1]
n = len(dy)
pos = np.zeros(n)
# position at the centre of the grid cell
pos[0] = -0.5*dy[1]
pos[1] = 0.5*dy[1]
for i in range(2,n):
	pos[i] = pos[i-1] + 0.5*dy[i-1] + 0.5*dy[i]
replace_guards(pos)

data_dict["pos"] = pos

# Normalisation factor for temperature - T * Tnorm returns in eV
data_dict["Tnorm"] = collect("Tnorm", path=filename, tind=-1)
# Normalisation factor for density - N * Nnorm returns in m^-3
data_dict["Nnorm"] = collect("Nnorm", path=filename, tind=-1)
# Plasma pressure (normalised). Pe = 2 Ne Te => P/Ne = Te (and assume Ti=Te)
data_dict["P"] = collect("P", path=filename, tind=-1)
# Electron density (normalised)
data_dict["Ne"] = collect("Ne", path=filename, tind=-1)
print("len(data_dict[Ne])")
print(len(data_dict["Ne"]))
# Neutral density (normalised)
data_dict["Nn"] = collect("Nn", path=filename, tind=-1)

data_dict["nnorm"] = collect("Nnorm", path=filename)  # m^-3
data_dict["tnorm"] = collect("Tnorm", path=filename)  # eV
data_dict["pnorm"] = data_dict["nnorm"]*1.602e-19*data_dict["tnorm"] # Pressure normalisation [Pa]
data_dict["cs0"] = collect("Cs0", path=filename) # m/s
data_dict["timenorm"] = collect("Omega_ci", path=filename)

data_dict["Enorm"] = 1.602e-19*data_dict["tnorm"]*data_dict["nnorm"]*data_dict["timenorm"]

data_dict["Rzrad"] = collect("Rzrad", path=filename, tind=-1) #Impurity radiation
data_dict["Rex"] = collect("Rex",   path=filename, tind=-1) #Hydrogen-excitation radiation

# Help for user
data_dict["help"] = "Contains outputs from Boutprojects/SD1D/MAST-U/ex-C-1.0/area-2.0/nloss-0.0/tn-0.5/0.2. Created with data_dict_export.py - stored in Github.com/TBody/atomic1D/reference"

from copy import deepcopy
import numpy as np
import json

# Need to 'jsonify' the numpy arrays (i.e. convert to nested lists) so that they can be stored in plain-text
# Deep-copy data to a new dictionary and then edit that one (i.e. break the data pointer association - keep data_dict unchanged in case you want to run a copy-verify on it)

data_dict_jsonified = deepcopy(data_dict)

numpy_ndarrays = [];
for key, element in data_dict.items():
    if type(element) == np.ndarray:
        # Store which keys correspond to numpy.ndarray, so that you can de-jsonify the arrays when reading
        numpy_ndarrays.append(key)
        data_dict_jsonified[key] = data_dict_jsonified[key].tolist()

data_dict_jsonified['numpy_ndarrays'] = numpy_ndarrays

# Encode help
# >> data_dict['help'] = 'help string'

# <<Use original filename, except with .json instead of .dat extension>>
with open('{}.json'.format(filename),'w') as fp:
    json.dump(data_dict_jsonified, fp, sort_keys=True, indent=4)