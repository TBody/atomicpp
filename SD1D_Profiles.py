import numpy as np
from atomicpp import atomicpy
from scipy.integrate import odeint #ODEPACK, for numerical integration
# from scipy.integrate import simps #Simpson's rule, for definite integrals
import pickle
from atomicpp import atomicpy
import matplotlib.pyplot as plt
from numpy import sqrt

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

		Vi = data_dict["Vi"]
		Vn = data_dict["Vn"]
		sound_speed = data_dict["Cs0"]


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
		self.ion_velocity 	  = np.squeeze(np.array(Vi*sound_speed))
		self.neutral_velocity = np.squeeze(np.array(Vn*sound_speed))

def alpha_e(k):
	return 0.71 * (k**2)

def beta_i(k, mu):
	return 3*(mu + 5*sqrt(2)*(k**2) * (1.1 * (mu**(5/2)) - 0.34 * (mu*(3/2)) -1 ))/(2.6 - 2*mu * 5.4 * (mu**2))

if __name__ == '__main__':

	# input_file = "sd1d-MAST-U-nups_0.2.json"
	# input_file = "0.2.json"
	input_file = "0.5.json"

	fixed_fraction = 1e-2
	Z = 6

	test_data = SD1DData(input_file)

	position = -(test_data.position - max(test_data.position))
	# position = test_data.position

	impurity_symbol = b'c'
	impurity = atomicpy.PyImpuritySpecies(impurity_symbol)

	mz = 12.0107 #amu, Carbon
	mi = 1.00794 #amu, Hydrogen
	mu = mz/(mz+mi)
	amu = 1.6605e-27 #kg

	impurity_derivatives = atomicpy.PyRateEquations(impurity)
	impurity_derivatives.setThresholdDensity(-1.0) #Don't use a threshold density at first
	impurity_derivatives.setDominantIonMass(1.0)

	plot_Prad    = []
	plot_P_stage = []
	plot_P_line  = []
	plot_P_cont  = []
	plot_P_cx    = []
	plot_Nzk 	 = []
	plot_dNzk 	 = []
	# plot_dPe 	 = []
	calc_NVzk 	 = []

	time_step = 1e-12

	for index in range(len(position)):
		Te = test_data.temperature[index]
		Ne = test_data.density[index]
		Nn = test_data.neutral_density[index]
		Nz = Ne * fixed_fraction

		Vi = test_data.ion_velocity[index]
		Vn = test_data.neutral_velocity[index]

		Nzk = impurity.calculateNzk(Te, Ne, Nz, Nn)
		plot_Nzk.append(Nzk)
		Vzk = np.zeros((Z+1,))

		collision_frequency_ii_PF = impurity_derivatives.calculateIonIonDragFactor(50, 1e19)

		if not(index == 0 or index == len(position)-1):
			# Central difference
			dPe = (test_data.pressure[index+1]-test_data.pressure[index-1])/(test_data.position[index+1]-test_data.position[index-1])
			dTe = (test_data.temperature[index+1]-test_data.temperature[index-1])/(test_data.position[index+1]-test_data.position[index-1])
		elif (index == 0):
			# Forward difference
			dPe = (test_data.pressure[index+1]-test_data.pressure[index])/(test_data.position[index+1]-test_data.position[index])
			dTe = (test_data.temperature[index+1]-test_data.temperature[index])/(test_data.position[index+1]-test_data.position[index])
		elif (index == len(position)-1):
			# Backward difference
			dPe = (test_data.pressure[index]-test_data.pressure[index-1])/(test_data.position[index]-test_data.position[index-1])
			dTe = (test_data.temperature[index]-test_data.temperature[index-1])/(test_data.position[index]-test_data.position[index-1])

		# print(Vi)
		for k in range(1,Z+1):
			Vzk[k] = Vi + (1/(mz*collision_frequency_ii_PF*k*k)) * ((1+k)*(-1/Ne)*dPe + (alpha_e(k)+beta_i(k, mu))*dTe)
			# print(Vzk[k])

		derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk)

		Pcool = np.array(derivative_struct["Pcool"])
		Prad = np.array(derivative_struct["Prad"])
		dNzk = np.array(derivative_struct["dNzk"])
		F_zk = np.array(derivative_struct["F_zk"])
		dNe = np.array(derivative_struct["dNe"])
		F_i = np.array(derivative_struct["F_i"])
		dNn = np.array(derivative_struct["dNn"])
		F_n = np.array(derivative_struct["F_n"])
		P_stage = np.array(derivative_struct["P_stage"])
		P_line = np.array(derivative_struct["P_line"])
		P_cont = np.array(derivative_struct["P_cont"])
		P_cx = np.array(derivative_struct["P_cx"])

		plot_Prad.append(Prad)
		plot_P_stage.append(P_stage)
		plot_P_line.append(	P_line )
		plot_P_cont.append(	P_cont )
		plot_P_cx.append(	P_cx   )
		# plot_Nzk.append(	Nzk 	)
		plot_dNzk.append(	dNzk 	)

		# calc_NVzk.append(Nzk*F_zk/(mz*amu)*time_step)


	# for index in range(len(position)):

	# 	if not(index == 0 or index == len(position)-1):
	# 		# Central difference
	# 		dNVzk = (calc_NVzk[index+1]-calc_NVzk[index-1])/(test_data.position[index+1]-test_data.position[index-1])
	# 	elif (index == 0):
	# 		# Forward difference
	# 		dNVzk = (calc_NVzk[index+1]-calc_NVzk[index])/(test_data.position[index+1]-test_data.position[index])
	# 	elif (index == len(position)-1):
	# 		# Backward difference
	# 		dNVzk = (calc_NVzk[index]-calc_NVzk[index-1])/(test_data.position[index]-test_data.position[index-1])

	# 	# plot_dNzk[index] -= dNVzk


	plot_Prad    = np.array(plot_Prad)
	plot_P_stage = np.array(plot_P_stage)
	plot_P_line  = np.array(plot_P_line)
	plot_P_cont  = np.array(plot_P_cont)
	plot_P_cx    = np.array(plot_P_cx)
	plot_Nzk 	 = np.array(plot_Nzk)
	plot_dNzk 	 = np.array(plot_dNzk)


	f, (ax1, ax2) = plt.subplots(2, sharex = True)

	ax1.plot(position, plot_dNzk)
	# ax2.plot(position, calc_NVzk)
	plt.show()


	if False:
		f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)

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
		ax2.semilogy(position, plot_P_line, label=r'Line')
		ax2.semilogy(position, plot_P_cont, label=r'Cont.')
		ax2.semilogy(position, plot_P_cx, label=r'C.X.')
		ax2.grid()
		ax2.legend()
		ax2.set_ylabel(r"P$_{rad}$ [$W/m^{3}$]")
		# ax3.set_ylim(1e6, 1e21)

		# Bottom plot - instantaneous change

		for k in range(Z+1):
			if k == 0:
				ax3.semilogy(position, plot_Nzk[:,k], label='g.s.'.format(k))
			else:
				ax3.semilogy(position, plot_Nzk[:,k], label='C{}+'.format(k))

		ax3.legend()
		ax3.set_ylim(1e6, 1e21)

		ax3.set_xlabel("Position [$m$]")

		# ax1.set_xbound(lower = min(position), upper = max(position))
		ax1.set_xbound(lower = 0, upper = 5)
		ax1.invert_xaxis()

		plt.show()












