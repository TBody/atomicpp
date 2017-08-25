# Program name: atomic1D/reference/build_json.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 14 July 2017
# 
# 
# Makes data and copies it into a .json file 'ADAS/0.2.json'

def replace_guards(var):
	"""
	This in-place replaces the points in the guard cells with the points on the boundary
	
	"""
	var[0] = 0.5*(var[0] + var[1])
	var[-1] = 0.5*(var[-1] + var[-2])

path = 'ADAS/0.5'

from boutdata.collect import collect
import numpy as np

t = collect("t_array", path=path)
tind = len(t)-1 # Get the last time point

var_list = ["Ne", "P", "Nn",      # Plasma profiles
            # "Srec", "Siz",   # Particle sources / sinks
            # "Frec", "Fiz", "Fcx", "Fel",   # Momentum source / sinks to neutrals
            "Rrec", "Riz", "Rzrad", "Rex", # Radiation, energy loss from system
            # "Erec", "Eiz", "Ecx", "Eel",   # Energy transfer between neutrals and plasma
            "Nnorm", "Tnorm", "Omega_ci", "rho_s0", "Cs0"]  # Normalisations

########################################################
# Position
dy = collect("dy", path=path)[0,:]
n = len(dy)
pos = np.zeros(n)

# position at the centre of the grid cell
pos[0] = 0.5*dy[0]
for i in range(1,n):
    pos[i] = pos[i-1] + 0.5*dy[i-1] + 0.5*dy[i]
    
########################################################
# Read the data into a dictionary


data = {}
data["pos"] = pos
for var in var_list:
    try:
        data[var] = collect(var, tind=tind, path=path)
        
        if len(data[var].shape) == 4:
            # 4D variable
            data[var] = data[var][0,0,:,0] # Make 1D [y]
    except:
        print("Variable '%s' not found" % (var,))
        data[var] = None

from copy import deepcopy
# import numpy as np
import json


# Need to 'jsonify' the numpy arrays (i.e. convert to nested lists) so that they can be stored in plain-text
# Deep-copy data to a new dictionary and then edit that one (i.e. break the data pointer association - keep data unchanged in case you want to run a copy-verify on it)

data_jsonified = deepcopy(data)

numpy_ndarrays = [];
for key, element in data.items():
    if type(element) == np.ndarray:
        # Store which keys correspond to numpy.ndarray, so that you can de-jsonify the arrays when reading
        numpy_ndarrays.append(key)
        data_jsonified[key] = data_jsonified[key].tolist()

data_jsonified['numpy_ndarrays'] = numpy_ndarrays

# Encode help
# >> data['help'] = 'help string'

# <<Use original filename, except with .json instead of .dat extension>>
with open('{}.json'.format(path),'w') as fp:
    json.dump(data_jsonified, fp, sort_keys=False, indent=4)