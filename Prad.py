# Program name: atomic1D/Prad.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 15 July 2017
#
# Program function: output the radiated power (Prad)
#                   by using OpenADAS rates on output JSON from SD1D run
#
# Under active development: <<TODO>> indicates development goal

from atomic1D import sharedFunctions
from atomic1D import ImpuritySpecies #Load the impurity species class
from atomic1D import SD1DData
import numpy as np

# Mapping between processes and their ADAS codes
# Note that the code will try to generate a file for each process listed
# Will return an error if file not found, except for charge exchange (ccd and prc)
# which is skipped if .has_charge_exchange = False
datatype_abbrevs = {
        'ionisation'           : 'scd',
        'recombination'        : 'acd',
        'cx_rec'               : 'ccd',
        'continuum_power'      : 'prb',
        'line_power'           : 'plt',
        'cx_power'             : 'prc',
        'ionisation_potential' : 'ecd', #N.b. ionisation_potential is not a rate-coefficient, but most of the
                                        #methods are transferable
}

# Invert the mapping of datatype_abbrevs
inv_datatype_abbrevs = {v: k for k, v in datatype_abbrevs.items()}

def calculateCollRadEquilibrium(impurity, experiment):
	# Calculates the fractional distribution across ionisation stages, assuming generalised collisional-radiative
	# equilibrium (i.e. effective ionisation and recombination rates are the only significant terms, charge-exchange
	# is ignored)
	#
	# Inputs: 	impurity 				= ImpuritySpecies object (the impurity for which we are calculating the radiation)
	# 			experiment.density    	= SD1DData object with electron/ion density (in m^-3)
	# 			experiment.temperature 	= SD1DData object with electron/ion density (in eV)
	#

	Z = impurity.atomic_number

	# Calculating for 1D linked data, output from SD1D
	data_length = len(experiment.temperature)
	# Assumed that temperature and density data are of the same length - check this before continuing
	assert data_length == len(experiment.density)

	# Preallocate an empty array for
	y = np.zeros((Z + 1, data_length))

	# Set the ground state density to zero (arbitrary - will normalise later)
	y[0,:] = 1

	for k in range(Z): #will return k = 0, 1, ..., Z-1 where Z = impurity.atomic_number
	                                        #i.e. iterate over all charge states except bare nucleus
		iz_coeffs = impurity.rate_coefficients['ionisation'].call1D(k, experiment.temperature, experiment.density)
		rec_coeffs = impurity.rate_coefficients['recombination'].call1D(k, experiment.temperature, experiment.density)

		# The ratio of ionisation from the (k)th stage and recombination from the (k+1)th
		# sets the equilibrium densities of the (k+1)th stage in terms of the (k)th (since
		# R = n_z * n_e * rate_coefficient)
		# N.b. Since there is no ionisation from the bare nucleus, and no recombination onto the
		# neutral (ignoring anion formation) the 'k' value of ionisation coeffs is shifted down
		# by one relative to the recombination coeffs - therefore this evaluation actually gives the
		# balance

		y[k+1,:] = y[k,:] * (iz_coeffs/rec_coeffs)


	# Normalise such that the sum over all ionisation stages is '1' at all points
	y = y / y.sum(axis=0)
	# Make sure that the fractional abundance returned has the same shape as the input arrays
	# assert density.shape == y.shape #For [t,s] case
	assert data_length == len(y[0,:]) #For [t=0, s] case
	assert np.allclose(y.sum(axis=0), 1.0)

	return y

def plot_iz_stage_distribution(experiment, iz_stage_distribution):
	# plot the plasma temperature, density and ionisation-stage distribution as a function of position
	# For a single time step

	import matplotlib.pyplot as plt

	# Create iterator for distance axis
	s = range(experiment.data_shape)

	plt.plot(s, experiment.temperature/max(experiment.temperature),\
		'--',label='T/{:.2f}[eV]'.format(max(experiment.temperature)))

	plt.plot(s, experiment.density/max(experiment.density),\
		'--',label=r'$n_e$/{:.2e}[$m^{{-3}}$]'.format(max(experiment.density)))

	for k in range(iz_stage_distribution.shape[0]): #iterate over all charge states
		plt.plot(iz_stage_distribution[k,:],label='{}'.format(k))

	plt.xlabel('Distance (downstream, a.u.)')
	plt.ylabel('Fraction')
	plt.legend()

	plt.show()

def computeRadiatedPower(impurity, experiment, iz_stage_distribution):
	# Calculate the power radiated from each location
	#
	radiated_power = {}
	stage_integrated_power = {}

	# Select the processes which provide a 'L' coefficient that contributes to radiation
	radiative_processes = ['continuum_power','line_power']
	if impurity.has_charge_exchange:
		# Include charge exchange if the impurity has this attribute
		radiative_processes.append('cx_power')

	for physics_process in radiative_processes:
		# Find the coefficient object (i.e. L s.t. Prad = L n_1 n_2)
		coeff = impurity.rate_coefficients[physics_process] #will return a RateCoefficient object

		# Coefficients have Z data entries, corresponding to a charge range of either 0 to (Z-1) for electron-
		# bound target or 1 to Z for charged target
		coeff_evaluated = np.zeros((impurity.atomic_number, experiment.data_shape))
		# Want to determine the power evaluated for all charge states, including fully ionised state
		# Therefore add +1 to range to include endpoint
		power_evaluated = np.zeros((impurity.atomic_number+1, experiment.data_shape))

		# Iterate over 0 to Z-1, i.e. not inclusive of endpoint since RateCoefficients are only 6 entries long
		# This corresponds to different target charge states depending on the process
		# For ionic target processes (recombination, charge exchange, cx power, continuum rec/brems power)
		# 	k = 1+ to Z+
		# For bound-electron target processes (ionisation and line excitation/relaxation power)
		# 	k = 0+ (neutral) to (Z-1)+

		for k in range(impurity.atomic_number):
			# Evaluate the coefficient
			# Will return a 1D array since call1D has grid=False on interpolation evaluation
			#
			# Make a list of evaluated rate coefficients [W m^3], where the first index gives the
			# ionisation stage and the second index gives the position
			coeff_evaluated[k,:] = coeff.call1D(k, experiment.temperature, experiment.density)

			if physics_process is 'line_power':
				# range of k is 0 to (Z-1)+ (needs bound electrons)
				target_charge_state = k #electron-bound target

				# Prad = L * n_e * n_z^k+
				#      = L * scale
				n_e = experiment.density
				n_z_charge_state = experiment.impurity_density * iz_stage_distribution[target_charge_state]
				scale = n_e * n_z_charge_state
			elif physics_process is 'continuum_power':
				# range of k is 1+ to Z+ (needs charged target)
				target_charge_state = k + 1 #charged target

				# Prad = L * n_e * n_z^(k+1)
				#      = L * scale
				n_e = experiment.density
				n_z_charge_state = experiment.impurity_density * iz_stage_distribution[target_charge_state]
				scale = n_e * n_z_charge_state
			elif physics_process is 'cx_power':
				# range of k is 1+ to Z+ (needs charged target)
				target_charge_state = k + 1 #charged target

				# Prad = L * n_0 * n_z^(k+1)+
				#      = L * scale
				n_0 = experiment.density * experiment.neutral_fraction
				n_z_charge_state = experiment.impurity_density * iz_stage_distribution[target_charge_state]
				scale = n_0 * n_z_charge_state

			# Computed the power radiated (W/m^3)
			# P_rad = L * scale
			power_evaluated[target_charge_state,:] = coeff_evaluated[k,:] * scale

		# Append the power radiated to the radiated_power dictionary
		# (i.e. add the full [Z * data_length] list to the dictionary)
		radiated_power[physics_process] = power_evaluated
		# Sum across all the ionisation stages to find total contribution from each
		# physics_process
		stage_integrated_power[physics_process] = np.sum(power_evaluated,axis=0)

	# Performs piece-wise add: has the same shape as the source
	# 'total' has sum over all physics_processes
	radiated_power['total'] = sum(radiated_power.values())
	stage_integrated_power['total'] = sum(stage_integrated_power.values())
	# total_power has total over all physics_processes and ionisation stages
	total_power = sum(radiated_power['total'])
	# Summing dimensions in different orders should give the same result
	assert np.allclose(total_power,stage_integrated_power['total'])

	# Check data shape
	if type(experiment.data_shape) == int:
		assert experiment.data_shape == len(total_power)
	elif type(experiment.data_shape) == np.ndarray:
		assert experiment.data_shape == total_power.shape
	else:
		raise NotImplementedError('Error checking data_shape match')

	# Copy the results to the experiment object
	experiment.radiated_power = radiated_power
	experiment.stage_integrated_power = stage_integrated_power
	experiment.total_power = total_power

def plotRadiatedPowerProcess(experiment):
	# plot the plasma temperature, density and radiated-power as a function of position
	# For a single time step

	import matplotlib.pyplot as plt

	labels = []

	# Create iterator for distance axis
	s = range(experiment.data_shape)

	fig, ax1 = plt.subplots()

	ax2 = ax1.twinx()

	for physics_process in experiment.stage_integrated_power.keys():
		ax1.plot(experiment.stage_integrated_power[physics_process],label='{}'.format(physics_process))


	ax1.set_xlabel('Distance (downstream, a.u.)')
	ax1.set_ylabel(r'Prad from process $W/m^3$', color='k')
	ax1.set_yscale('log')
	ax1.legend(loc=3)

	ax2.set_ylabel('Experiment parameter', color='b')

	ax2.plot(s, experiment.temperature/max(experiment.temperature),\
		'--',label='T/{:.2f}[eV]'.format(max(experiment.temperature)))

	ax2.plot(s, experiment.density/max(experiment.density),\
		'--',label=r'$n_e$/{:.2e}[$m^{{-3}}$]'.format(max(experiment.density)))

	ax2.legend(loc=4)

	plt.show()

def plotRadiatedPowerStage(impurity, experiment):
	# plot the plasma temperature, density and radiated-power as a function of position
	# For a single time step

	import matplotlib.pyplot as plt

	labels = []

	# Create iterator for distance axis
	s = range(experiment.data_shape)

	fig, ax1 = plt.subplots()

	ax2 = ax1.twinx()

	# Include the fully-ionised state (endpoint of range, requires +1)
	for k in range(impurity.atomic_number+1):
		ax1.plot(experiment.radiated_power['total'][k,:],label='{}'.format(k))
	ax1.plot(experiment.total_power,label='total')

	ax1.set_xlabel('Distance (downstream, a.u.)')
	ax1.set_ylabel(r'Prad from stage $W/m^3$', color='k')
	ax1.set_yscale('log')
	ax1.legend(loc=3)

	ax2.set_ylabel('Experiment parameter', color='b')

	ax2.plot(s, experiment.temperature/max(experiment.temperature),\
		'--',label='T/{:.2f}[eV]'.format(max(experiment.temperature)))

	ax2.plot(s, experiment.density/max(experiment.density),\
		'--',label=r'$n_e$/{:.2e}[$m^{{-3}}$]'.format(max(experiment.density)))

	ax2.legend(loc=4)

	plt.show()

def plotRadiatedChargeExchange(impurity, experiment):
	# plot the plasma temperature, density and radiated-power as a function of position
	# For a single time step

	import matplotlib.pyplot as plt

	labels = []

	# Create iterator for distance axis
	s = range(experiment.data_shape)

	fig, ax1 = plt.subplots()

	ax2 = ax1.twinx()

	# Include the fully-ionised state (endpoint of range, requires +1)
	# Note that cx_power requires charged target
	for k in range(1,impurity.atomic_number+1):
		ax1.plot(experiment.radiated_power['cx_power'][k,:],label='{}'.format(k))

	ax1.set_xlabel('Distance (downstream, a.u.)')
	ax1.set_ylabel(r'Prad from cx $W/m^3$', color='k')
	ax1.set_yscale('log')
	ax1.legend(loc=3)

	ax2.set_ylabel('Experiment parameter', color='b')

	ax2.plot(s, experiment.temperature/max(experiment.temperature),\
		'--',label='T/{:.2f}[eV]'.format(max(experiment.temperature)))

	ax2.plot(s, experiment.density/max(experiment.density),\
		'--',label=r'$n_e$/{:.2e}[$m^{{-3}}$]'.format(max(experiment.density)))

	ax2.plot(s, experiment.neutral_fraction/max(experiment.neutral_fraction),\
		'--',label=r'$n_o/n_e$/{:.2e}[$m^{{-3}}$]'.format(max(experiment.neutral_fraction)))

	ax2.legend(loc=4)

	plt.show()

def computeDerivs(impurity, experiment, iz_stage_distribution):

	constant_position_index = 0;
	Z = impurity.atomic_number;

	# Te, Ne, Nn, Nik - copy these out to C++ function form
	Te = experiment.temperature[constant_position_index] #in eV
	Ne = experiment.density[constant_position_index] #in m^-3
	Nn = experiment.density[constant_position_index] * experiment.neutral_fraction[constant_position_index]
	Ni = experiment.impurity_density[constant_position_index] #summed over all stages

	Nik = Ni * iz_stage_distribution[:,constant_position_index]

	# print("{:>8}:{:.6e}".format("Te",Te))
	# print("{:>8}:{:.6e}".format("Ne",Ne))
	# print("{:>8}:{:.6e}".format("Nn",Nn))
	# for k in range(Z+1):
	# 	print("{:>5}^{}: {:.2e}".format("Ni",k,Nik[k]))
	# for physics_process in ['ionisation','recombination','continuum_power','line_power','cx_rec','cx_power']:
	# 	coeff = impurity.rate_coefficients[physics_process]
	# 	for k in range(Z):
	# 		coeff_evaluated = coeff.call1D(k, Te, Ne) #return scalar
	# 		print("{}({}) = {:.5e}".format(physics_process, k, coeff_evaluated))

	Pcool = 0 # Pcool = dydt[0]
	Prad = 0 # Prad  = dydt[1]
	dNik = np.zeros_like(Nik) # dNik  = dydt[2:Z+3]
	dNe  = 0 # dNe   = dydt[Z+3]
	dNn  = 0 # dNn   = dydt[Z+3+1]

	# Select the processes which provide a 'k' coefficient that contributes to population equation
	population_processes = ['ionisation','recombination']
	# Select the processes which provide a 'L' coefficient that contributes to radiation
	radiative_processes = ['continuum_power','line_power']


	for physics_process in population_processes:
		# Find the coefficient object (i.e. k s.t. dNi = k n_1 n_2)
		coeff = impurity.rate_coefficients[physics_process] #will return a RateCoefficient object

		# Iterate over 0 to Z-1, i.e. not inclusive of endpoint since RateCoefficients are only 6 entries long
		# This corresponds to different target charge states depending on the process
		# For ionic target processes (recombination, charge exchange, cx power, continuum rec/brems power)
		# 	k = 1+ to Z+
		# For bound-electron target processes (ionisation and line excitation/relaxation power)
		# 	k = 0+ (neutral) to (Z-1)+

		for k in range(Z):
			# Evaluate the coefficient
			# Will return a 1D array since call1D has grid=False on interpolation evaluation
			#
			# Make a list of evaluated rate coefficients [W m^3], where the first index gives the
			# ionisation stage and the second index gives the position
			coeff_evaluated = coeff.call1D(k, Te, Ne) #return scalar

			if physics_process is 'ionisation':
				# range of k is 0 to (Z-1)+ (needs bound electrons)
				source_charge_state = k #electron-bound target
				sink_charge_state = k + 1 #resulting charge state

			elif physics_process is 'recombination':
				# range of k is 1+ to Z+ (needs charged target)
				source_charge_state = k + 1 #charged target
				sink_charge_state = k #resulting charge state

			rate = coeff_evaluated * Ne * Nik[source_charge_state]

			# Compute the resulting change in the ionisation distribution
			dNik[source_charge_state] -= rate
			dNik[sink_charge_state] += rate
			if physics_process is 'ionisation':
				dNe += rate
			elif physics_process is 'recombination':
				dNe -= rate

			print("{} ({})        R = {:+1.20E}".format(physics_process, k, rate));
			print("{} ({}) dNik ({}) = {:+1.20E}".format(physics_process, k, source_charge_state, dNik[source_charge_state]));
			print("{} ({}) dNik ({}) = {:+1.20E}".format(physics_process, k, sink_charge_state, dNik[sink_charge_state]));
			print("{} ({})      dNe = {:+1.20E}".format(physics_process, k, dNe));

			# print("{}({}),   R = {:.5e}".format(physics_process, k, rate))
			# print("{}({}), dNe = {:.5e}".format(physics_process, k, dNe))

	# Calculate the electron cooling power due to ionisation balance
	# (This could be added to the above loop for speed, but is kept here for clarity)
	# Extract the ionisation potential (in J per transition)
	ionisation_potential = impurity.rate_coefficients['ionisation_potential']

	eV_to_J = 1.60217662e-19 #J per eV
	for k in range(Z):
		Pcool += eV_to_J * ionisation_potential.call1D(k,Te,Ne) * dNik[k]
	# Find that this is a very small effect


	for physics_process in radiative_processes:
		# Find the coefficient object (i.e. L s.t. Prad = L n_1 n_2)
		coeff = impurity.rate_coefficients[physics_process] #will return a RateCoefficient object

		# Iterate over 0 to Z-1, i.e. not inclusive of endpoint since RateCoefficients are only Z entries long
		# This corresponds to different target charge states depending on the process
		# For ionic target processes (recombination, charge exchange, cx power, continuum rec/brems power)
		# 	k = 1+ to Z+
		# For bound-electron target processes (ionisation and line excitation/relaxation power)
		# 	k = 0+ (neutral) to (Z-1)+

		for k in range(Z):
			# Evaluate the coefficient
			# Will return a scalr since call1D is called with scalar arguments
			#
			# Evaluate rate coefficients [W m^3], where the first index gives the
			# ionisation stage and the second index gives the position
			coeff_evaluated = coeff.call1D(k, Te, Ne)

			if physics_process is 'line_power':
				# range of k is 0 to (Z-1)+ (needs bound electrons)
				source_charge_state = k #electron-bound target

				# Prad = L * Ne * n_z^k+
				#      = L * scale
				scale = Ne * Nik[source_charge_state]
			elif physics_process is 'continuum_power':
				# range of k is 1+ to Z+ (needs charged target)
				source_charge_state = k + 1 #charged target

				# Prad = L * Ne * n_z^(k+1)
				#      = L * scale
				scale = Ne * Nik[source_charge_state]

			# Computed the power radiated (W/m^3)
			# P_rad = L * scale
			Prad += coeff_evaluated * scale

	# Verified that radiated power gives the same result as experiment.total_power
	# Add the total (no cx) Prad to Pcool
	Pcool += Prad

	# if impurity.has_charge_exchange:
	# 	# Include charge exchange if the impurity has this attribute
	# 	physics_processes = ['cx_rec','cx_power']
	# 	coeff = impurity.rate_coefficients[physics_process] #will return a RateCoefficient object

	# 	for k in range(Z):
	# 		coeff_evaluated = coeff.call1D(k, Te, Ne) #return scalar

	# 		if physics_process is 'cx_rec':
	# 				# range of k is 1+ to Z+ (needs charged target)
	# 				source_charge_state = k + 1 #charged target
	# 				sink_charge_state = k #resulting charge state

	# 				# dNi = k * Nn * n_z^(k+1)+
	# 				#     = k * scale
	# 				scale = Nn * Nik[source_charge_state]

	# 				dNik[source_charge_state] -= coeff_evaluated * scale
	# 				dNik[sink_charge_state] += coeff_evaluated * scale
	# 				dNn -= coeff_evaluated * scale

	# 		elif physics_process is 'cx_power':
	# 				# range of k is 1+ to Z+ (needs charged target)
	# 				source_charge_state = k + 1 #charged target

	# 				# Prad = L * Nn * n_z^(k+1)+
	# 				#      = L * scale
	# 				scale = Nn * Nik[source_charge_state]
	# 				Prad += coeff_evaluated * scale

	dydt = np.zeros(Z+4+1)

	dydt[0] 	= Pcool
	dydt[1] 	= Prad
	dydt[2:Z+3] = dNik
	dydt[Z+3] 	= dNe
	dydt[Z+3+1] = dNn

	return dydt


if __name__ == '__main__':

	# Process command line arguments to set the path to the input file (from SD1D)
	# and the JSON database (from make json_update)
	[input_file, JSON_database_path, impurity] = sharedFunctions.processCommandLineArguments()

	# Add the JSON files associated with this impurity to its .adas_files_dict attribute
	# where the key is the (extended) process name, which maps to a filename (string)
	# Check that these files exist in the JSON_database_path/json_data/ directory
	for physics_process,filetype_code in datatype_abbrevs.items():
		if impurity.has_charge_exchange or not(filetype_code in {'ccd', 'prc'}):
			impurity.addJSONFiles(physics_process,filetype_code,JSON_database_path)

	# Use the .adas_file_dict files to generate RateCoefficient objects for each process
	# Uses the same keys as .adas_file_dict
	impurity.makeRateCoefficients(JSON_database_path)

	# Inspections on the RateCoefficient object
	# Can supply the process of interest (i.e. 'ionisation') and the ionisation stage (i.e. 1), the
	# function will return a 3D plot of the data used, and also the interpolation generated.
	# impurity.rate_coefficients['ionisation'].compare_interp_with_plot(1)

	# Process the input_file to extract
	# 	density(t,s) 					= electron density (in m^-3)
	# 	temperature(t,s)				= electron/ion temperature (in eV)
	# 	neutral_fraction(t,s)			= neutral density/electron density (no units)
	# where t is time index, s is 1D distance index
	experiment = SD1DData(input_file)

	t = experiment.data_shape[-1]
	# Extract data for a single time-step
	experiment.selectSingleTime(t)

	# Set the fixed fraction for impurity concentration
	experiment.setImpurityFraction(1e-2)

	# Calculate the distribution across ionisation stages, assuming collisional-radiative equilibrium
	iz_stage_distribution = calculateCollRadEquilibrium(impurity, experiment)

	# Plot the ionisation stage distribution as a function of distance
	# plot_iz_stage_distribution(experiment, iz_stage_distribution)

	# Compute radiated power
	# Returns total_power, stage_integrated_power (sum over all ionisation stages), and
	# radiated_power (resolved into different ionisation stages, with 'total' giving sum over
	# all physics_processes)
	# 	stage_integrated_power and radiated_power are dictionaries with physics_process keys and
	# 	data_length shape for stage_integrated_power and [Z, data_length] shape for radiated_power
	# 	total_power is an array of shape data_length
	computeRadiatedPower(impurity, experiment, iz_stage_distribution)

	# Export results/plot
	# plotRadiatedPowerProcess(experiment)
	# plotRadiatedPowerStage(impurity,experiment)
	# plotRadiatedChargeExchange(impurity,experiment)

	dydt = computeDerivs(impurity, experiment, iz_stage_distribution)

	Z = impurity.atomic_number

	Pcool = dydt[0]
	Prad  = dydt[1]
	dNik  = dydt[2:Z+3]
	dNe   = dydt[Z+3]
	dNn   = dydt[Z+3+1]

	print("{:>7}: {:+.2e} [{}]".format("Pcool",Pcool,"W m^-3"))
	print("{:>7}: {:+.2e} [{}]".format("Prad",Prad,"W m^-3"))
	for k in range(Z+1):
		print("{:>5}^{}: {:+.2e} [{}]".format("dNi",k,dNik[k],"m^-3 s^-1"))
	print("{:>7}: {:+.2e} [{}]".format("dNe",dNe,"m^-3 s^-1"))
	print("{:>7}: {:+.2e} [{}]".format("dNn",dNn,"m^-3 s^-1"))




























