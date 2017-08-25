import numpy as np
from atomicpp import atomicpy
from scipy.integrate import odeint #ODEPACK, for numerical integration
# from scipy.integrate import simps #Simpson's rule, for definite integrals
import pickle
from atomicpp import atomicpy
import matplotlib.pyplot as plt

def retrieveFromJSON(file_name):
	# Inputs - a JSON file corresponding to an OpenADAS .dat file or SD1D output file
	# file_name can be either relative or absolute path to JSON file
	# Must have .json extension and match keys of creation
	# Not need for the .dat -> .json conversion, but included for reference
	import json
	from warnings import warn
	from copy import deepcopy
	import numpy as np

	file_extension  = file_name.split('.')[-1] #Look at the extension only (last element of split on '.')
	if file_extension != 'json':
		raise NotImplementedError('File extension (.{}) is not .json'.format(file_extension))

	with open(file_name,'r') as fp:
		data_dict = json.load(fp)

	# if  set(data_dict.keys()) != {'numpy_ndarrays', 'charge', 'help', 'log_density', 'number_of_charge_states', 'log_temperature', 'element', 'log_coeff', 'name', 'class'}\
	# and set(data_dict.keys()) != {'numpy_ndarrays', 'help', 'Ne', 'Tnorm', 'P', 'Nnorm', 'Nn'}:
	# 	warn('Imported JSON file {} does not have the expected set of keys - could result in an error'.format(file_name))

	# Convert jsonified numpy.ndarrays back from nested lists
	data_dict_dejsonified = deepcopy(data_dict)

	for key in data_dict['numpy_ndarrays']:
		data_dict_dejsonified[key] = np.array(data_dict_dejsonified[key])

	return data_dict_dejsonified

class SD1DData(object):
	# For storing the data output from SD1D. To create the required JSON run the function
	# data_dict_export.py in an I/O (case) folder in SD1D.
	#
	def __init__(self,input_file):
		# process a input JSON file to extract Te(s,t), ne(s,t), ne/nn (s,t)
		# n.b. s refers to the upstream distance from the strike-point
		#      t is time (will need to add normalisation factor <<TODO>> to convert to real-time)
		#
		# input:    input_file -> JSON file from SD1D run
		# return:   Te, ne, neutral_fraction (from data)
		#           impurity_fraction (fixed fraction impurity density, set programmatically)

		# input_file can be either relative or absolute path to JSON file

		data_dict = retrieveFromJSON(input_file)

		# Retrieve (normalised values)
		Ne = data_dict['Ne']
		Nn = data_dict['Nn']
		# P = 2*Ne*Te => Te = P/(2*Ne)
		# N.b. division between two numpy ndarrays is piecewise
		T  = data_dict['P']/(2*data_dict['Ne'])
		p = data_dict['P']
		nnorm = data_dict['Nnorm']
		tnorm = data_dict['Tnorm']
		pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
		Enorm = 1.602e-19*data_dict["Tnorm"]*data_dict["Nnorm"]*data_dict["Omega_ci"]
		# Neutral fraction affects charge exchange

		# Retrieve normalisation factors
		# Converts N into m^-3, T into eV

		Enorm = Enorm
		Rzrad = data_dict["Rzrad"]
		Rex = data_dict["Rex"]

		# Dimensions are [t, x, y, z]
		#                [0, 1, 2, 3]
		#
		# SD1D only outputs on time (t) and distance from strike-point (currently stored as y)
		#
		# Apply normalisations, and then use np.squeeze to remove single-dimensional entries
		# Should return a 2D numpy array with no length-1 dimensions
		self.position		  = np.squeeze(np.array(data_dict["pos"]))
		self.temperature      = np.squeeze(np.array(T*tnorm))
		self.density          = np.squeeze(np.array(Ne*nnorm))
		self.pressure 		  = np.squeeze(np.array(p*pnorm))
		self.neutral_density  = np.squeeze(np.array(Nn*nnorm))
		self.p_rad_hydrogen   = np.squeeze(np.array(Rex*Enorm))
		self.p_rad_carbon     = np.squeeze(np.array(Rzrad*Enorm))

		# data_shape = self.temperature.shape
		# assert data_shape == self.density.shape
		# assert data_shape == self.neutral_fraction.shape

		# self.data_shape = data_shape

	# def setImpurityFraction(self,impurity_fraction):
	# 	# Set the impurity density (will be scalar multiplication if fixed fraction,
	# 	# and piecewise multiplication (i.e. Hadamard product) if impurity_fraction is
	# 	# an array of shape equal to density, otherwise will result in error)

	# 	self.impurity_fraction = impurity_fraction
	# 	self.setImpurityDensity(self.impurity_fraction * self.density)

	# def setImpurityDensity(self,impurity_density):
	# 	# Set the impurity density, and check that it has the same shape as the other attributes
	# 	self.impurity_density  = impurity_density

	# 	if type(self.data_shape) == int:
	# 		assert self.data_shape == len(self.impurity_density)
	# 	elif type(self.data_shape) == np.ndarray:
	# 		assert self.data_shape == self.impurity_density.shape
	# 	else:
	# 		raise NotImplementedError('Error checking data_shape match')

	# def selectSingleTime(self,t):
	# 	try:
	# 		self.temperature      = self.temperature[t,:]
	# 		self.density          = self.density[t,:]
	# 		self.neutral_fraction = self.neutral_fraction[t,:]
	# 	except IndexError:
	# 		# If using the length of the time array, will get an out-of-bounds error (since python indexes from 0)
	# 		t -= 1
	# 		self.temperature      = self.temperature[t,:]
	# 		self.density          = self.density[t,:]
	# 		self.neutral_fraction = self.neutral_fraction[t,:]

	# 	self.data_shape = self.data_shape[1]


if __name__ == '__main__':

	# input_file = "sd1d-MAST-U-nups_0.2.json"
	# input_file = "0.2.json"
	input_file = "0.5.json"

	fixed_fraction = 1e-2

	test_data = SD1DData(input_file)

	position = -(test_data.position - max(test_data.position))

	impurity_symbol = b'c'
	impurity = atomicpy.PyImpuritySpecies(impurity_symbol)

	impurity_derivatives = atomicpy.PyRateEquations(impurity)
	impurity_derivatives.setThresholdDensity(-1.0) #Don't use a threshold density at first
	impurity_derivatives.setDominantIonMass(1.0)

	Prad = []
	

	for index in range(position):


	if False:
		f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)

		# Top plot - density and temperature

		ax1.semilogy(position, test_data.density, '-', label=r'$N_e$')
		ax1.semilogy(position, test_data.neutral_density, '-', label=r'$N_n$')
		ax1.semilogy(position, test_data.density*fixed_fraction, '-', label=r'$\Sigma$ $N_z$')
		ax1.set_ylim(1e16, 1e23)
		ax1.grid()
		ax1.set_ylabel(r"Density [$m^{-3}$]")

		ax1_twin = ax1.twinx()
		ax1_twin.plot(position, test_data.temperature, 'r--', label=r'T_e')
		ax1_twin.set_ylabel(r"Temperature [$eV$]", color='r')
		ax1_twin.tick_params('y', colors='r')
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax1_twin.get_legend_handles_labels()
		ax1.set_zorder(1)
		ax1.legend(h1+h2, l1+l2, loc=2)

		# Centre plot - power

		ax2.semilogy(position, test_data.p_rad_hydrogen, label=r'H')
		ax2.semilogy(position, test_data.p_rad_carbon, label=r'C')
		ax2.grid()
		ax2.legend()
		ax2.set_ylabel(r"P$_{rad}$ [$W/m^{3}$]")

		# Bottom plot - instantaneous change

		ax3.set_xlabel("Position [$m$]")

		# ax1.set_xbound(lower = min(position), upper = max(position))
		ax1.set_xbound(lower = 0, upper = 5)
		ax1.invert_xaxis()

		plt.show()












