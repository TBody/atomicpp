import numpy as np
from scipy.interpolate import RectBivariateSpline

class RateCoefficient(object):
	# For storing the RateCoefficients encoded in an OpenADAS data file
	# Intended to be called from the .makeRateCoefficients method of an ImpuritySpecies object
	# 
	# Closely based on the cfe316/atomic/atomic_data.py/RateCoefficient class
	# 
	# Interpolation tables for the rate of some physical process.
	# Contains one 2D spline interpolator for each charge state of an element,
	# per  process like 'ionisation', 'recombination'.

	# Attributes:
	#     atomic_number (int) : The element's Z.
	#     element (str)       : Short element name like 'c'
	#     adf11_file (str)    : The /full/filename it came from (link to .json, not .dat)
	#     log_temperature     : np.array of log10 of temperature values
	#     log_density         : np.array of log10 of density values
	#     log_coeff           : a 3D np.array with shape (Z, temp, dens)
	#     splines             : list of scipy.interpolate.fitpack2.RectBivariateSpline
	#         The list has length Z and is interpolations of log_coeff.
	
	def __init__(self,filename):
		# Create an instance of RateCoefficient by reading an OpenADAS JSON file
		from atomic1D import sharedFunctions

		data_dict = sharedFunctions.retrieveFromJSON(filename)

		self.atomic_number   = data_dict['charge']
		self.element         = data_dict['element']
		self.adf11_file      = filename
		self.log_temperature = data_dict['log_temperature']
		self.log_density     = data_dict['log_density']
		self.log_coeff       = data_dict['log_coeff']

		self._compute_interpolating_splines()

	def __str__(self):
		return 'RateCoefficient object with attributes'+\
		'\n'+'{:>25} = {}'.format('atomic_number',		self.atomic_number)+\
		'\n'+'{:>25} = {}'.format('element',			self.element)+\
		'\n'+'{:>25} = {}'.format('adf11_file',			self.adf11_file)+\
		'\n'+'{:>25} = {}'.format('log_temperature',	'{} numpy array'.format(self.log_temperature.shape))+\
		'\n'+'{:>25} = {}'.format('log_density',		'{} numpy array'.format(self.log_density.shape))+\
		'\n'+'{:>25} = {}'.format('log_coeff',			'{} numpy array'.format(self.log_coeff.shape))+\
		'\n'+'{:>25} = {}'.format('splines[k]',			'RectBivariateSpline interpolations')

	def _compute_interpolating_splines(self):
		# Generate the interpolation functions for log_coeff
		
		self.splines = []
		for k in range(self.atomic_number):
			x = self.log_temperature
			y = self.log_density
			z = self.log_coeff[k]
			self.splines.append(RectBivariateSpline(x, y, z))

	def call1D(self, k, Te, ne):
		"""Evaluate the ionisation/recombination coefficients of
		k'th atomic state at a given temperature and density.

		For 1D case - either Te or ne vary while the other is fixed (i.e. parameter scan),
		or both vary but are linked (i.e. analysis of 1D output)
		>>Key difference to call2D is the grid=False argument - otherwise will return
		  interpolations over the len(Te)*len(ne) 2D array

		Args:
			k  (int): Ionising or recombined ion stage,
				between 0 and k=Z-1, where Z is atomic number.
			Te (array_like): Temperature in [eV].
			ne (array_like): Density in [m-3].

		Returns:
			c (array_like): Rate coefficent in [m3/s].
		"""

		# broadcast_arrays ensures that Te and ne are of the same shape
		# if they are originally equal then they remain so, while if len(A) = L
		# and len(B) = 1 then B' = L repeats of B
		Te, ne = np.broadcast_arrays(Te, ne)

		# Need to convert both temp and density to log-scale, since this is what the spline-interpolation is performed for
		log_temperature = np.log10(Te)
		log_density = np.log10(ne)

		# Find the logarithm of the rate-coefficient
		log_coeff = self.splines[k](log_temperature, log_density, grid=False)
		# Raise (piecewise) to the power 10 to return in m3/s
		coeffs = np.power(10,log_coeff)

		return coeffs

	def call2D(self, k, Te, ne):
		"""Evaluate the ionisation/recombination coefficients of
		k'th atomic state at a given temperature and density.

		For 2D case - Te and ne vary independantly - return a 2D array of results across
		parameter space
		>>Key difference to call1D is the implicit grid=True argument - will return
		  interpolations over the len(Te)*len(ne) 2D array

		Args:
			k  (int): Ionising or recombined ion stage,
				between 0 and k=Z-1, where Z is atomic number.
			Te (array_like): Temperature in [eV].
			ne (array_like): Density in [m-3].

		Returns:
			c (array_like): Rate coefficent in [m3/s].
		"""

		# broadcast_arrays ensures that Te and ne are of the same shape
		# if they are originally equal then they remain so, while if len(A) = L
		# and len(B) = 1 then B' = L repeats of B
		Te, ne = np.broadcast_arrays(Te, ne)

		# Need to convert both temp and density to log-scale, since this is what the spline-interpolation is performed for
		log_temperature = np.log10(Te)
		log_density = np.log10(ne)

		# Find the logarithm of the rate-coefficient
		log_coeff = self.splines[k](log_temperature, log_density)
		# Raise (piecewise) to the power 10 to return in m3/s
		coeffs = np.power(10,log_coeff)

		return coeffs

	@property
	def temperature_grid(self):
		#Get a np.array of temperatures in [eV]
		#These are the points for which the rate-coefficients are given in the file
		#
		return np.power(10, self.log_temperature)

	@property
	def density_grid(self):
		#Get an np.array of densities in [m^3]
		#These are the points for which the rate-coefficients are given in the file
		#
		return np.power(10, self.log_density)

	def inspect_with_plot(self, ionisation_stage):
		import matplotlib as mpl
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		# For TeX labelling
		from matplotlib import rc
		plt.rc('text', usetex=True)
		plt.rc('font', family='sans-serif')

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		X, Y = np.meshgrid(self.log_temperature, self.log_density)
		Z = np.transpose(self.log_coeff[ionisation_stage, :, :])

		# Plot a basic wireframe.
		ax.plot_wireframe(X, Y, Z)

		ax.set_xlabel(r'$log(T_e) [eV]$')
		ax.set_ylabel(r'$log(n_e) [m^{-3}]$')
		ax.set_zlabel(r'$log(R) [m^3s^{-1}]$')

		plt.show()

	def inspect_interp_with_plot(self, ionisation_stage):
		import matplotlib as mpl
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		# For TeX labelling
		from matplotlib import rc
		plt.rc('text', usetex=True)
		plt.rc('font', family='sans-serif')

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		X, Y = np.meshgrid(self.log_temperature, self.log_density)
		Z = np.transpose(self.log_coeff[ionisation_stage, :, :])

		# Plot a basic wireframe.
		ax.plot_wireframe(X, Y, Z)

		ax.set_xlabel(r'$log(T_e) [eV]$')
		ax.set_ylabel(r'$log(n_e) [m^{-3}]$')
		ax.set_zlabel(r'$log(R) [m^3s^{-1}]$')

		plt.show()

	def compare_interp_with_plot(self,ionisation_stage):
		import matplotlib as mpl
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		# For TeX labelling
		from matplotlib import rc
		plt.rc('text', usetex=True)
		plt.rc('font', family='sans-serif')

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		X, Y = np.meshgrid(self.log_temperature, self.log_density)
		Z = np.transpose(self.splines[ionisation_stage](self.log_temperature, self.log_density))
		Z_interp = np.transpose(self.splines[ionisation_stage](self.log_temperature, self.log_density))

		# Plot a basic wireframe.
		# ax.plot_wireframe(X, Y, Z, colors='r',linestyles='-')
		ax.plot_surface(X, Y, Z, linewidth=0, antialiased=False)
		ax.plot_wireframe(X, Y, Z_interp, colors='r',linestyles='-')

		ax.set_xlabel(r'$log(T_e) [eV]$')
		ax.set_ylabel(r'$log(n_e) [m^{-3}]$')
		ax.set_zlabel(r'$log(R) [m^3s^{-1}]$')

		plt.show()

