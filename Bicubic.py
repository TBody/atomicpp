import numpy as np
from atomicpp.Bicubic import PyBicubicSpline
from atomicpp.Bicubic import PyBilinearSpline
from atomicpp.Bicubic import PyBivariateBSpline
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator

import math
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, Normalize

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

	if  set(data_dict.keys()) != {'numpy_ndarrays', 'charge', 'help', 'log_density', 'number_of_charge_states', 'log_temperature', 'element', 'log_coeff', 'name', 'class'}\
	and set(data_dict.keys()) != {'numpy_ndarrays', 'help', 'Ne', 'Tnorm', 'P', 'Nnorm', 'Nn'}:
		warn('Imported JSON file {} does not have the expected set of keys - could result in an error'.format(file_name))

	# Convert jsonified numpy.ndarrays back from nested lists
	data_dict_dejsonified = deepcopy(data_dict)

	for key in data_dict['numpy_ndarrays']:
		data_dict_dejsonified[key] = np.array(data_dict_dejsonified[key])

	return data_dict_dejsonified

class wrapRGI():
	# Wraps RegularGridInterpolator so that it may be called with a standard method
	def __init__(self, x, y, z):
		self.interpolator = RegularGridInterpolator((x, y), z, method='linear', bounds_error=True)
	def __call__(self, x, y):
		return self.interpolator.__call__((x,y))

def gen_func(x, y):
	# z = np.cos(x**2 + y**2) #Wave ripple pattern
	# z = x**2 + (y-1)**2
	z = x**2 + (y-1)**2 + x*y
	# z = x
	# z = y
	# z = x+y
	# z = math.sin(x)
	# z = math.sin(y)
	# z = math.sin(x+y)
	# z = math.sin(x) + math.cos(y)

	return z

def generate_from_func_uniform(func, x_length, y_length, x_min, x_max, y_min, y_max):
	x_values = np.linspace(x_min, x_max, x_length)
	y_values = np.linspace(y_min, y_max, y_length)
	
	z_values = np.zeros((len(y_values), len(x_values)))

	for i in range(x_length):
		for j in range(y_length):
			x = x_values[i]
			y = y_values[j]

			try:
				z_value = func(x, y)
			except RuntimeError:
				# Edge values can throw RuntimeError in C++ interpolator
				z_value = 0
				if x == min(x_values):
					x_shift = (x_values[i+1]-x_values[i])*1e-6
				if y == min(y_values):
					y_shift = (y_values[j+1]-y_values[j])*1e-6

				if x == max(x_values):
					x_shift = (x_values[i-1]-x_values[i])*1e-6
				if y == max(y_values):
					y_shift = (y_values[j-1]-y_values[j])*1e-6

				z_value = func(x + x_shift, y + y_shift)

			# z_values[y_length -1 - j][i] = z_value
			z_values[j][i] = z_value

	return {'x': x_values, 'y': y_values, 'z': z_values}

def generate_from_func_vectors(func, x_values, y_values):
	z_values = np.zeros((len(y_values), len(x_values)))

	for i in range(len(x_values)):
		for j in range(len(y_values)):
			x = x_values[i]
			y = y_values[j]

			try:
				z_value = func(x, y)
			except RuntimeError:
				# Edge values can throw RuntimeError in C++ interpolator
				z_value = 0
				if x == min(x_values):
					x_shift = (x_values[i+1]-x_values[i])*1e-12
				if y == min(y_values):
					y_shift = (y_values[j+1]-y_values[j])*1e-12

				if x == max(x_values):
					x_shift = (x_values[i-1]-x_values[i])*1e-12
				if y == max(y_values):
					y_shift = (y_values[j-1]-y_values[j])*1e-12

				z_value = func(x + x_shift, y + y_shift)

			# z_values[y_length -1 - j][i] = z_value
			z_values[j][i] = z_value

	return {'x': x_values, 'y': y_values, 'z': z_values}

def upscale(func, x_length, y_length, _x_values, _y_values):
	x_values = np.linspace(min(_x_values), max(_x_values), x_length)
	y_values = np.linspace(min(_y_values), max(_y_values), y_length)
	
	z_values = np.zeros((len(y_values), len(x_values)))

	for i in range(x_length):
		for j in range(y_length):
			x = x_values[i]
			y = y_values[j]

			try:
				z_value = func(x, y)
			except RuntimeError:
				# Edge values can throw RuntimeError in C++ interpolator
				z_value = 0
				if x == min(x_values):
					x_shift = (x_values[i+1]-x_values[i])*1e-12
				if y == min(y_values):
					y_shift = (y_values[j+1]-y_values[j])*1e-12

				if x == max(x_values):
					x_shift = (x_values[i-1]-x_values[i])*1e-12
				if y == max(y_values):
					y_shift = (y_values[j-1]-y_values[j])*1e-12

				z_value = func(x + x_shift, y + y_shift)

			# z_values[y_length -1 - j][i] = z_value
			z_values[j][i] = z_value

	return {'x': x_values, 'y': y_values, 'z': z_values}

def plot_compare(low_res, hi_res, cpp_hi_res, pyx_hi_res):
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')

	ax1.pcolor(low_res['x'], low_res['y'], low_res['z'],cmap='RdBu')
	ax1.set_title('Low res')

	ax2.pcolor(hi_res['x'], hi_res['y'], hi_res['z'],cmap='RdBu')
	ax2.set_title('High res')

	ax3.pcolor(cpp_hi_res['x'], cpp_hi_res['y'], cpp_hi_res['z'],cmap='RdBu')
	ax3.set_title('C++')

	ax4.pcolor(pyx_hi_res['x'], pyx_hi_res['y'], pyx_hi_res['z'],cmap='RdBu')
	ax4.set_title('Python')

	plt.show()

def plot_difference(source, candidate):
	difference = source['z'] - candidate['z']

	f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

	ax1.pcolor(source['x'], source['y'], source['z'],cmap='RdBu')
	ax1.set_title('Source')

	im = ax2.pcolor(source['x'], source['y'], difference,cmap='RdBu')
	ax2.set_title('Difference')
	
	f.subplots_adjust(right=0.8)
	cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
	f.colorbar(im, cax=cbar_ax)

	plt.show()

def plot_ratio(source, candidate):
	ratio = candidate['z']/source['z']

	f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)

	ax1.pcolor(source['x'], source['y'], source['z'],cmap='RdBu')
	ax1.set_title('Source')

	im = ax2.pcolor(source['x'], source['y'], ratio,cmap='RdBu')
	ax2.set_title('Ratio')
	
	f.subplots_adjust(right=0.8)
	cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
	f.colorbar(im, cax=cbar_ax)

	plt.show()

def inspect_grid_coeff(source, candidate):

	plt.pcolor(source['x'], source['y'], source['z'],cmap='RdBu')

	x_values = np.array(candidate.get_x_values())
	y_values = np.array(candidate.get_y_values())

	grid_coeff = np.array(candidate.get_grid_coeff())
	grid_coeff = grid_coeff.transpose()
	data_shape = grid_coeff.shape

	print("Data shape is {} by {}".format(data_shape[0], data_shape[1]))
	
	f_grid     = np.zeros_like(grid_coeff)
	fdx_grid   = np.zeros_like(grid_coeff)
	fdy_grid   = np.zeros_like(grid_coeff)
	fdxdy_grid = np.zeros_like(grid_coeff)

	for i in range(data_shape[0]):
		for j in range(data_shape[1]):
			# coord = grid_coeff[i][j]['coord']
			# print("({}, {}) = ({}, {})".format(i, j, coord[0], coord[1]))
			datapoint = grid_coeff[i][j]['datapoint']
			f_grid[i][j]     = grid_coeff[i][j]['f']
			fdx_grid[i][j]   = grid_coeff[i][j]['fdx']
			fdy_grid[i][j]   = grid_coeff[i][j]['fdy']
			fdxdy_grid[i][j] = grid_coeff[i][j]['fdxdy']
			print("({:<+5.2f}, {:<+5.2f}) = ({:<+5.2f}, {:<+5.2f}) -> ({:<+5.2f}, {:<+5.2f}, {:<+5.2f}, {:<+5.2f})".format(x_values[j], y_values[i], datapoint[0], datapoint[1], f_grid[i][j],fdx_grid[i][j], fdy_grid[i][j], fdxdy_grid[i][j]))

	X, Y = np.meshgrid(x_values, y_values)

	if False:
		# Plot the function values
		# If the same colourmap is used between pcolor and scatter (which is default,cmap='RdBu')
		#   you shouldn't be able to see the points at all
		plt.scatter(X, Y, c=f_grid, marker = '.')
		plt.scatter(1, 1, c='r', marker='x')

	elif True:
		# Plot dx and dy as a quiver grid -- should look like a 'grad function'
		U = fdx_grid
		V = fdy_grid

		Q = plt.quiver(X.tolist(), Y.tolist(), U.tolist(), V.tolist(), pivot='mid',angles='xy')
		# plt.streamplot(X, Y, U, V)
	else:
		# What should this look like? http://www.math.harvard.edu/archive/21a_fall_08/exhibits/fxy/
		U = fdx_grid
		V = fdy_grid
		W = fdxdy_grid

		Q1 = plt.quiver(X, Y, U, 0, pivot='tail')
		# Q2 = plt.quiver(X, Y, 0, W, pivot='tail')



	plt.show()

def x_compare(hi_res, pyxspline, cppspline, y_const, trim_x):

	y_const_index = np.searchsorted(hi_res['y'], y_const)
	
	actual = hi_res['z'][y_const_index][:]
	pyx = pyxspline['z'][y_const_index][:]
	cpp = cppspline['z'][y_const_index][:]

	# plt.plot(hi_res['x'], actual,'g')
	# plt.plot(hi_res['x'], pyx,'b')
	# plt.plot(hi_res['x'], cpp,'r')

	start_in = np.searchsorted(hi_res['x'], trim_x[0])
	end_in = np.searchsorted(hi_res['x'], trim_x[1])

	plt.plot(hi_res['x'][start_in:end_in], (pyx - actual)[start_in:end_in],'b')
	plt.plot(hi_res['x'][start_in:end_in], (cpp - actual)[start_in:end_in],'r')

	plt.show()

def y_compare(hi_res, pyxspline, cppspline, x_const, trim_y):

	x_const_index = np.searchsorted(hi_res['y'], x_const)
	
	actual = hi_res['z'][:][x_const_index]
	pyx = pyxspline['z'][:][x_const_index]
	cpp = cppspline['z'][:][x_const_index]

	# plt.plot(hi_res['y'], actual,'g')
	# plt.plot(hi_res['y'], pyx,'b')
	# plt.plot(hi_res['y'], cpp,'r')

	start_in = np.searchsorted(hi_res['y'], trim_y[0])
	end_in = np.searchsorted(hi_res['y'], trim_y[1])

	plt.plot(hi_res['y'][start_in:end_in], (pyx - actual)[start_in:end_in],'b')
	plt.plot(hi_res['y'][start_in:end_in], (cpp - actual)[start_in:end_in],'r')

	plt.show()

if __name__ == '__main__':

	# data_dict = retrieveFromJSON('json_database/json_data/scd96_c.json')
	# x_values = data_dict['log_temperature']
	# y_values = data_dict['log_density']
	# x_min = min(x_values)
	# x_max = max(x_values)
	# y_min = min(y_values)
	# y_max = max(y_values)
	# z_values = data_dict['log_coeff']

	x_lowres = 10
	y_lowres = 15
	x_highres = 50
	y_highres = 50

	x_min = -2
	x_max = 3
	y_min = -4
	y_max = 3

	x_in_min = x_min + (x_max-x_min)/(x_lowres-1) #Edge-of-grid internal values
	x_in_max = x_max - (x_max-x_min)/(x_lowres-1)
	y_in_min = y_min + (y_max-y_min)/(y_lowres-1)
	y_in_max = y_max - (y_max-y_min)/(y_lowres-1)

	low_res = generate_from_func_uniform(gen_func, x_lowres, y_lowres, x_min, x_max, y_min, y_max)
	# low_res = generate_from_func_vectors(gen_func, x_values, y_values)
	hi_res = generate_from_func_uniform(gen_func, x_highres, y_highres, x_min, x_max, y_min, y_max)

	x_values = np.array(low_res['x'])
	y_values = np.array(low_res['y'])
	z_values = np.array(low_res['z'])

	# plt.pcolor(x_values, y_values, z_values,cmap='RdBu',
	# 	cmap='RdBu')
	# plt.colorbar()

	# plt.show()

	if True:
		cpp_interp = PyBivariateBSpline(x_values, y_values, z_values.transpose())
		pyx_interp = RectBivariateSpline(x_values, y_values, z_values.transpose(),kx=3,ky=3)
	elif True:
		cpp_interp = PyBicubicSpline(x_values, y_values, z_values.transpose())
		pyx_interp = RectBivariateSpline(x_values, y_values, z_values.transpose(),kx=3,ky=3)
	else:
		cpp_interp = PyBilinearSpline(x_values, y_values, z_values.transpose())
		pyx_interp = wrapRGI(x_values, y_values, z_values.transpose())

	# inspect_grid_coeff(hi_res, cpp_interp)

	# cpp_hi_res = generate_from_func(cpp_interp.call0D, x_highres, y_highres, x_min, x_max, y_min, y_max)
	cpp_hi_res = upscale(cpp_interp.call0D, x_highres, y_highres, low_res['x'], low_res['y'])
	# pyx_hi_res = generate_from_func(pyx_interp, x_highres, y_highres, x_min, x_max, y_min, y_max)
	pyx_hi_res = upscale(pyx_interp, x_highres, y_highres, low_res['x'], low_res['y'])

	# plot_compare(low_res, hi_res, cpp_hi_res, pyx_hi_res)

	# plot_difference(pyx_hi_res, cpp_hi_res)

	# plot_ratio(pyx_hi_res, cpp_hi_res)

	# cpp_hi_res = generate_from_func(cpp_interp.call0D, len(x_values), len(y_values), x_values.min(), x_values.max(), y_values.min(), y_values.max())

	# x_compare(hi_res, pyx_hi_res, cpp_hi_res, 0, (x_in_min, x_in_max))
	# y_compare(hi_res, pyx_hi_res, cpp_hi_res, 0, (y_in_min, y_in_max))


	# plt.pcolor(z_values[0] - cpp_hi_res['z'],cmap='RdBu',
	# 	cmap = 'jet',
	# 	extent=(y_values.min(), y_values.max(), x_values.min(), x_values.max())
	# 	)
	# plt.colorbar()

	# plt.show()






