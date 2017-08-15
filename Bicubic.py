import numpy as np
from atomicpp.Bicubic import PyBicubicSpline
from atomicpp.Bicubic import PyBilinearSpline
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

def gen_func(x, y):
	# z = 10 * np.cos(x**2 + y**2) #Wave ripple pattern
	# z = x**2 + y**2
	z = math.sin(x+y)

	return z

def generate_from_func(func, x_length, y_length, x_min, x_max, y_min, y_max):
	x_values = np.linspace(x_min, x_max, x_length)
	y_values = np.linspace(y_min, y_max, y_length)
	
	z_values = np.zeros((len(x_values), len(y_values)))

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
			z_values[i][j] = z_value

	return {'x': x_values, 'y': y_values, 'z': z_values}

def upscale(func, x_length, y_length, _x_values, _y_values):
	x_values = np.linspace(min(_x_values), max(_x_values), x_length)
	y_values = np.linspace(min(_y_values), max(_y_values), y_length)
	
	z_values = np.zeros((len(x_values), len(y_values)))

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
			z_values[i][j] = z_value

	return {'x': x_values, 'y': y_values, 'z': z_values}

def plot_compare(low_res, hi_res, cpp_hi_res, pyx_hi_res):
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

	ax1.imshow(low_res['z'],
		aspect='auto',
		origin='lower',
		extent=(low_res['x'].min(),low_res['x'].max(),low_res['y'].min(),low_res['y'].max())
		)
	ax1.set_title('Low res')

	ax2.imshow(hi_res['z'],
		aspect='auto',
		origin='lower',
		extent=(hi_res['x'].min(),hi_res['x'].max(),hi_res['y'].min(),hi_res['y'].max())
		)
	ax2.set_title('High res')

	ax3.imshow(cpp_hi_res['z'],
		aspect='auto',
		origin='lower',
		extent=(cpp_hi_res['x'].min(),cpp_hi_res['x'].max(),cpp_hi_res['y'].min(),cpp_hi_res['y'].max())
		)
	ax3.set_title('C++')

	ax4.imshow(pyx_hi_res['z'],
		aspect='auto',
		origin='lower',
		extent=(pyx_hi_res['x'].min(),pyx_hi_res['x'].max(),pyx_hi_res['y'].min(),pyx_hi_res['y'].max())
		)
	ax4.set_title('Python')

	plt.show()

def plot_difference(source, candidate):
	difference = source['z'] - candidate['z']

	f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

	ax1.imshow(source['z'],
		aspect='auto',
		# cmap = 'jet',
		origin='lower',
		extent=(source['x'].min(),source['x'].max(),source['y'].min(),source['y'].max())
		)
	ax1.set_title('Source')

	im = ax2.imshow(difference,
		aspect='auto',
		cmap = 'jet',
		origin='lower',
		extent=(source['x'].min(),source['x'].max(),source['y'].min(),source['y'].max())
		)
	ax2.set_title('Difference')
	
	f.subplots_adjust(right=0.8)
	cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
	f.colorbar(im, cax=cbar_ax)

	plt.show()

def plot_ratio(source, candidate):
	ratio = candidate['z']/source['z']

	f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)

	ax1.imshow(source['z'],
		aspect='auto',
		# cmap = 'jet',
		origin='lower',
		extent=(source['x'].min(),source['x'].max(),source['y'].min(),source['y'].max())
		)
	ax1.set_title('Source')

	im = ax2.imshow(ratio,
		aspect='auto',
		cmap = 'jet',
		origin='lower',
		extent=(source['x'].min(),source['x'].max(),source['y'].min(),source['y'].max())
		)
	ax2.set_title('Ratio')
	
	f.subplots_adjust(right=0.8)
	cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
	f.colorbar(im, cax=cbar_ax)

	plt.show()

class wrapRGI():
	def __init__(self, x, y, z):
		self.interpolator = RegularGridInterpolator((x, y), z, method='linear', bounds_error=True)
	def __call__(self, x, y):
		return self.interpolator.__call__((x,y))



if __name__ == '__main__':

	# data_dict = retrieveFromJSON('json_database/json_data/scd96_c.json')

	# x_values = data_dict['log_temperature']
	# y_values = data_dict['log_density']
	# z_values = data_dict['log_coeff']

	x_lowres = 10
	y_lowres = 10
	x_highres = 100
	y_highres = 100


	x_min = -4
	x_max = +2
	y_min = -4
	y_max = +2

	low_res = generate_from_func(gen_func, x_lowres, y_lowres, x_min, x_max, y_min, y_max)
	
	if False:
		cpp_interp = PyBicubicSpline(low_res['x'], low_res['y'], low_res['z'])
		pyx_interp = RectBivariateSpline(low_res['x'], low_res['y'], low_res['z'])
	else:
		cpp_interp = PyBilinearSpline(low_res['x'], low_res['y'], low_res['z'])
		# x_grid, y_grid = np.meshgrid(low_res['x'], low_res['y'])
		pyx_interp = wrapRGI(low_res['x'], low_res['y'], low_res['z'])

	# cpp_hi_res = generate_from_func(cpp_interp.call0D, x_highres, y_highres, x_min, x_max, y_min, y_max)
	cpp_hi_res = upscale(cpp_interp.call0D, x_highres, y_highres, low_res['x'], low_res['y'])
	# pyx_hi_res = generate_from_func(pyx_interp, x_highres, y_highres, x_min, x_max, y_min, y_max)

	pyx_hi_res = upscale(pyx_interp, x_highres, y_highres, low_res['x'], low_res['y'])

	hi_res = generate_from_func(gen_func, x_highres, y_highres, x_min, x_max, y_min, y_max)

	# plot_compare(low_res, hi_res, cpp_hi_res, pyx_hi_res)

	plot_difference(pyx_hi_res, cpp_hi_res)

	plot_ratio(pyx_hi_res, cpp_hi_res)

	# cpp_hi_res = generate_from_func(cpp_interp.call0D, len(x_values), len(y_values), x_values.min(), x_values.max(), y_values.min(), y_values.max())



	# plt.imshow(z_values[0] - cpp_hi_res['z'],
	# 	aspect='auto',
	# 	cmap = 'jet',
	# 	origin='lower',
	# 	extent=(y_values.min(), y_values.max(), x_values.min(), x_values.max())
	# 	)
	# plt.colorbar()

	# plt.show()






