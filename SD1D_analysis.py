
from boutdata import collect
# from boutdata.data import BoutOptionsFile
import matplotlib.pyplot as plt
# from scipy import array, sqrt
import numpy as np
from glob import glob
import os
from atomicpp import atomicpy

def replace_guards(var):
	"""
	This in-place replaces the points in the guard cells with the points on the boundary
	
	"""
	var[0] = 0.5*(var[0] + var[1])
	var[-1] = 0.5*(var[-1] + var[-2])

def extract_data(directory):

	fullpaths = glob(directory+"/*/")

	scan_data = {}

	# n_target = []
	max_density_pos = []
	max_z_rad_pos = []
	max_h_rad_pos = []
	n_upstream = []
	# p_target = []
	# p_upstream = []
	# te_target = []
	# te_upstream = []
	# nvi_target = []

	# Loop through the list of paths 
	for path in fullpaths:
		print("Reading from: "+path)
		try:
			# Extract the position from the data
			# Changed slice from [0,1:-1] to [-1,1:-1] - doubt that position grid changes with time, but just in case
			dy = collect("dy", path=path, yguards=True)[-1,1:-1]
			n = len(dy)
			pos = np.zeros(n)
			# position at the centre of the grid cell
			pos[0] = -0.5*dy[1]
			pos[1] = 0.5*dy[1]
			for i in range(2,n):
				pos[i] = pos[i-1] + 0.5*dy[i-1] + 0.5*dy[i]
			replace_guards(pos)

			# Read the final time point
			n = collect("Ne", path=path, yguards=True, tind=-1)
			# p = collect("P", path=path, yguards=True, tind=-1)
			# nvi = collect("NVi", path=path, yguards=True, tind=-1)
			
			# Read normalisations
			nnorm = collect("Nnorm", path=path)  # m^-3
			tnorm = collect("Tnorm", path=path)  # eV
			pnorm = nnorm*1.602e-19*tnorm # Pressure normalisation [Pa]
			cs0 = collect("Cs0", path=path) # m/s
			timenorm = collect("Omega_ci", path=path)
			
			# Get value at upper boundary
			# between cells 1 (guard cell) and 2 (in domain)
			n_up = 0.5*(n[-1,0,1,0] + n[-1,0,2,0])
			# p_up = 0.5*(p[-1,0,1,0] + p[-1,0,2,0])
			
			n_upstream.append( n_up * nnorm )
			# p_upstream.append( p_up * pnorm )
			# te_upstream.append( tnorm * 0.5*p_up / n_up )
			
			# target value between cell -2 (guard cell) and -3 (in domain)
			# n_targ = 0.5*(n[-1,0,-3,0] + n[-1,0,-2,0])
			# p_targ = 0.5*(p[-1,0,-3,0] + p[-1,0,-2,0])
			# nvi_targ = nvi[-1,0,-2,0] # Note: NVi in guard cells is Nout*Vout
			
			# n_target.append( n_targ * nnorm )
			# p_target.append( p_targ * pnorm )
			# te_target.append( tnorm * 0.5*p_targ / n_targ )
			# nvi_target.append( cs0 * nnorm * nvi_targ )

			Enorm = 1.602e-19*tnorm*nnorm*timenorm
			# Extract the radiation profiles
			# Rrec  = collect("Rrec",  path=path, tind=-1) #Recombination radiation
			# Riz   = collect("Riz",   path=path, tind=-1) #Ionisation radiation?
			Rzrad = collect("Rzrad", path=path, tind=-1) #Impurity radiation
			Rex   = collect("Rex",   path=path, tind=-1) #Hydrogen-excitation radiation

			# Erec = collect("Erec", path=path, tind=-1) # Energy transfer between neutrals and plasma
			# Eiz  = collect("Eiz", path=path, tind=-1)
			# Ecx  = collect("Ecx", path=path, tind=-1)
			# Eel  = collect("Eel", path=path, tind=-1)

			# How to use the recovered data - if you wanted to plot as a function of position
			# plt.plot(pos, (data["Rrec"]+data["Erec"])*Enorm, label="Recombination (Rrec+Erec)")
			# plt.plot(pos, (data["Riz"]+data["Eiz"])*Enorm, label="Ionisation (Riz+Eiz)")
			# plt.plot(pos, data["Ecx"]*Enorm, label="Charge exchange (Ecx)")
			# plt.plot(pos, data["Rzrad"]*Enorm, label="Impurity radiation (Rzrad)")
			# if data["Rex"] is not None:
			#     plt.plot(pos, data["Rex"]*Enorm, label="Hydrogen excitation (Rex)")
			# if data["Eel"] is not None:
			#     plt.plot(pos, data["Eel"]*Enorm, label="Elastic scattering (Eel)")

			# N.b. argmax only returns the first maximum value - might need to check that this is unique
			# Index from both ends?
			max_density_index = np.argmax(n)
			max_z_rad_index = np.argmax(Rzrad)
			max_h_rad_index = np.argmax(Rex)
			# Find what position that corresponds to
			max_density_pos.append(pos[max_density_index])
			max_z_rad_pos.append(pos[max_z_rad_index])
			max_h_rad_pos.append(pos[max_h_rad_index])

		except IOError:
			print("Couldn't read data in '{0}'. Skipping".format(path))

	scan_data['n_upstream'] = n_upstream
	scan_data['max_density_pos'] = max_density_pos
	scan_data['max_z_rad_pos'] = max_z_rad_pos
	scan_data['max_h_rad_pos'] = max_h_rad_pos

	return scan_data

def sort_by_upstream(scan_data):
	n_upstream = scan_data['n_upstream']
	max_density_pos = scan_data['max_density_pos']
	max_z_rad_pos = scan_data['max_z_rad_pos']
	max_h_rad_pos = scan_data['max_h_rad_pos']

	n_upstream = np.array(n_upstream)
	max_density_pos = np.array(max_density_pos)
	max_z_rad_pos = np.array(max_z_rad_pos)
	max_h_rad_pos = np.array(max_h_rad_pos)

	inds = np.argsort(n_upstream)

	n_upstream = n_upstream[inds]
	max_density_pos = max_density_pos[inds]
	max_z_rad_pos = max_z_rad_pos[inds]
	max_h_rad_pos = max_h_rad_pos[inds]
	
	processed_data = {}
	processed_data['n_upstream'] = n_upstream
	processed_data['max_density_pos'] = max_density_pos
	processed_data['max_z_rad_pos'] = max_z_rad_pos
	processed_data['max_h_rad_pos'] = max_h_rad_pos


	# for key, value in scan_data.items():
	# 	processed_data[key] = np.array(value)

	# # Sort by n_upstream

	# print(inds)

	# n_upstream = n_upstream[inds]
	# max_density_pos = max_density_pos[inds]
	# max_z_rad_pos = max_z_rad_pos[inds]
	# max_h_rad_pos = max_h_rad_pos[inds]

	# # processed_data = {"n_upstream" : n_upstream,"max_density_pos" : max_density_pos,"max_z_rad_pos" : max_z_rad_pos,"max_h_rad_pos" : max_h_rad_pos}

	# for key, value in scan_data.items():
	# 	processed_data[key] = value[inds]
	# 	# value = value[inds] #is this deep or shallow copy?

	return processed_data

def plot_compare(plot_data, show=False):

	fig, ax = plt.subplots()

	# plot adas as solid lines
	ADAS_x = plot_data['ADAS']['n_upstream']
	ax.plot(ADAS_x, 30-plot_data['ADAS']['max_density_pos'], 'r.-', label='Density')
	ax.plot(ADAS_x, 30-plot_data['ADAS']['max_z_rad_pos'], 'g.-', label='Impurity rad.')
	ax.plot(ADAS_x, 30-plot_data['ADAS']['max_h_rad_pos'], 'b.-', label='Hydrogen rad.')

	# Plot hutchinson as dashed
	Hutchinson_x = plot_data['Hutchinson']['n_upstream']
	ax.plot(Hutchinson_x, 30-plot_data['Hutchinson']['max_density_pos'], 'rx')
	ax.plot(Hutchinson_x, 30-plot_data['Hutchinson']['max_z_rad_pos'], 'gx')
	ax.plot(Hutchinson_x, 30-plot_data['Hutchinson']['max_h_rad_pos'], 'bx')

	ax.set_xlabel(r'Upstream density [$m^{-3}$]')
	ax.set_ylabel('Distance to strike point [m]')

	ax.legend(loc=4)
	ax.invert_yaxis()
	ax.set_ylim(1,0)
	ax.grid()
	
	if show:
		plt.show()
	return fig

if __name__ == '__main__':
	
	import pickle	
	# directories = {'ADAS': 'ex-C-1.0/area-2.0/nloss-0.0/tn-0.5', 'Hutchinson': 'ex-CH-1.0/area-2.0/nloss-0.0/tn-0.5'}
	directories = {'ADAS': 'ADAS', 'Hutchinson': 'Hutchinson'}
	path_to_output = 'Figures/'

	show = True
	reread = False
	if reread:
		plot_data = {}

		for key, directory in directories.items():

			scan_data = extract_data(directory)

			scan_data = sort_by_upstream(scan_data)

			plot_data[key] = scan_data

		with open('python_results/SD1D_data_density_scan.pickle', 'wb') as handle:
				pickle.dump(plot_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/SD1D_data_density_scan.pickle', 'rb') as handle:
			plot_data = pickle.load(handle)

	density_scan_plot = plot_compare(plot_data)
	density_scan_plot.set_size_inches(6.268, 3.52575, forward=True)
	plt.tight_layout()
	density_scan_plot.savefig(path_to_output+"density_scan_SD1D.pdf",bbox_inches = 'tight',pad_inches = 0.03)

	if show:
		plt.show()








	