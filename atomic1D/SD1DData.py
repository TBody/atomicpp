import numpy as np
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
		from atomic1D import sharedFunctions

		data_dict = sharedFunctions.retrieveFromJSON(input_file)

		# Retrieve (normalised values)
		Ne = data_dict['Ne']
		Nn = data_dict['Nn']
		# P = 2*Ne*Te => Te = P/(2*Ne)
		# N.b. division between two numpy ndarrays is piecewise
		T  = data_dict['P']/(2*data_dict['Ne'])
		
		# Neutral fraction affects charge exchange
		neutral_fraction = Nn/Ne
		
		# Retrieve normalisation factors
		Nnorm = data_dict['Nnorm']
		Tnorm = data_dict['Tnorm']
		# Converts N into m^-3, T into eV 
		
		# Dimensions are [t, x, y, z]
		#                [0, 1, 2, 3]
		#
		# SD1D only outputs on time (t) and distance from strike-point (currently stored as y)
		# 
		# Apply normalisations, and then use np.squeeze to remove single-dimensional entries
		# Should return a 2D numpy array with no length-1 dimensions
		self.temperature      = np.squeeze(np.array(T*Tnorm))
		self.density          = np.squeeze(np.array(Ne*Nnorm))
		self.neutral_fraction = np.squeeze(np.array(neutral_fraction))

		data_shape = self.temperature.shape
		assert data_shape == self.density.shape
		assert data_shape == self.neutral_fraction.shape
		
		self.data_shape = data_shape

	def setImpurityFraction(self,impurity_fraction):
		# Set the impurity density (will be scalar multiplication if fixed fraction,
		# and piecewise multiplication (i.e. Hadamard product) if impurity_fraction is
		# an array of shape equal to density, otherwise will result in error)
		
		self.impurity_fraction = impurity_fraction
		self.setImpurityDensity(self.impurity_fraction * self.density)

	def setImpurityDensity(self,impurity_density):
		# Set the impurity density, and check that it has the same shape as the other attributes
		self.impurity_density  = impurity_density
		
		if type(self.data_shape) == int:
			assert self.data_shape == len(self.impurity_density)
		elif type(self.data_shape) == np.ndarray:
			assert self.data_shape == self.impurity_density.shape
		else:
			raise NotImplementedError('Error checking data_shape match')

	def selectSingleTime(self,t):
		try:
			self.temperature      = self.temperature[t,:]
			self.density          = self.density[t,:]
			self.neutral_fraction = self.neutral_fraction[t,:]
		except IndexError:
			# If using the length of the time array, will get an out-of-bounds error (since python indexes from 0)
			t -= 1
			self.temperature      = self.temperature[t,:]
			self.density          = self.density[t,:]
			self.neutral_fraction = self.neutral_fraction[t,:]

		self.data_shape = self.data_shape[1]




















