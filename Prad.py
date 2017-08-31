# Program name: atomicpp/Prad.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 11 August 2017
#
# Use the atomic++ module to evaluate the rate-coefficients from OpenADAS

import numpy as np
from atomicpp import atomicpy
from scipy.integrate import odeint #ODEPACK, for numerical integration
# from scipy.integrate import simps #Simpson's rule, for definite integrals
import pickle
import matplotlib.pyplot as plt
import random
random.seed(1) #To ensure results are reproducible

# Code to hide OPEPACK/lsoda warnings (stdout)
# From https://stackoverflow.com/questions/31681946/disable-warnings-originating-from-scipy
import os
import sys
import contextlib

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextlib.contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    """
    https://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

# OOP method for solving the differential equations, to allow for additional output
class AtomicSolver(object):

	def __init__(self, impurity_symbol):

		# ImpuritySpecies
		self.impurity_symbol = impurity_symbol
		self.impurity = atomicpy.PyImpuritySpecies(impurity_symbol)

		# Evaluation parameters
		self.t_values = np.logspace(-6, 2, 200)
		self.Te_values = np.logspace(-0.69, 3.99, 100) #eV, span the entire array for which there is ADAS data
		self.Te_const = 50
		self.Ne_values = np.logspace(13.7, 21.3, 100) #m^-3
		self.Ne_const = 1e19
		self.Ne_tau_values = [1e17, 1e16, 1e15, 1e14] #m^-3 s, values to return Prad(tau) for

		# RateEquations
		self.impurity_derivatives = atomicpy.PyRateEquations(self.impurity)
		self.impurity_derivatives.setThresholdDensity(-1.0) #Don't use a threshold density at first
		self.impurity_derivatives.setDominantIonMass(1.0)

		# Initial values
		self.Z = self.impurity.get_atomic_number()
		self.Te  = 50 #eV
		self.Ne  = 1e19 #m^-3
		self.Vi  = 0 #m/s
		self.Nn  = 0 #m^-3
		self.Vn  = 0 #m/s
		self.Nzk = np.zeros((self.Z+1,)) #m^-3
		self.Nzk[0] = 1e17 #m^-3 - start in g.s.
		self.Vzk = np.zeros((self.Z+1,)) #m/s

		# Additional output initialisation
		self.additional_out = {'Prad':[], 'Pcool':[], 'dNzk':[], 'F_zk':[], 'dNe':[], 'F_i':[], 'dNn':[], 'F_n':[]} #Blank lists to append onto
		self.additional_out_keys = ['Prad', 'Pcool', 'dNzk'] #Keys to record data for
	
	@staticmethod
	def evolveDensity(Nzk, t, self, Te, Ne):
		# Te  = self.Te
		# Ne  = self.Ne
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		# Nzk = self.Nzk
		Vzk = self.Vzk

		# Prevent negative densities
		# (these are possible if the time-step is large)
		for k in range(len(Nzk)):
			if(Nzk[k] < 0):
				Nzk[k] = 0

		derivative_struct = self.impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

		dNzk = derivative_struct["dNzk"]
		for key in self.additional_out_keys:
			self.additional_out[key].append(derivative_struct[key])

		return dNzk

	@staticmethod
	def evolveDensity_withRefuelling(Nzk, t, self, Te, Ne, refuelling_rate):
		# Te  = self.Te
		# Ne  = self.Ne
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		# Nzk = self.Nzk
		Vzk = self.Vzk

		# Prevent negative densities
		# (these are possible if the time-step is large)
		for k in range(len(Nzk)):
			if(Nzk[k] < 0):
				Nzk[k] = 0

		derivative_struct = self.impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);
		dNzk = derivative_struct["dNzk"]

		fraction_in_stage = Nzk/sum(Nzk)
		# Add neutrals at a rate of tau^-1
		dNzk[0] += sum(Nzk)*refuelling_rate
		# Remove other stages based on their density
		dNzk -= sum(Nzk)*refuelling_rate*fraction_in_stage

		for key in self.additional_out_keys:
			self.additional_out[key].append(derivative_struct[key])

		return dNzk

	def reset_additional_out(self):
		self.additional_out = {'Prad':[], 'Pcool':[], 'dNzk':[], 'F_zk':[], 'dNe':[], 'F_i':[], 'dNn':[], 'F_n':[]}

	def timeIntegrate(self, Te, Ne, refuelling_rate = 0):
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		Vzk = self.Vzk

		print("Te = {:.2e}eV, Ne = {:.2e}/m3, tau_inv = {:.2e}".format(Te, Ne, refuelling_rate))
		
		if refuelling_rate == 0:
			# No refuelling - coronal equilibrium case
			with stdout_redirected():
				(result, output_dictionary) = odeint(self.evolveDensity, self.Nzk, self.t_values, args=(self, Te, Ne), printmessg=False, full_output=True, mxhnil=0)
		else:
			# Refuelling case
			with stdout_redirected():
				(result, output_dictionary) = odeint(self.evolveDensity_withRefuelling, self.Nzk, self.t_values, args=(self, Te, Ne, refuelling_rate), printmessg=False, full_output=True, mxhnil=0)
			# Will change the result, but may be treated the same as the CR case

		feval_at_step = output_dictionary['nfe'] #function evaluations at the time-step
		time_at_step = output_dictionary['tcur'] #time at the time-step

		time_indices = np.searchsorted(self.t_values, time_at_step, side='left') #find how the time-steps are distributed. Usually close but not 1 to 1 with self.t_values

		for key, value in self.additional_out.items():
			if value: #if list isn't empty - i.e. additional output has been recorded for this key
				output_feval = value #copy the output evaluated at each time-step

				try:
					# If the additional_out has a length (i.e. is an array)
					output_values = np.zeros((len(self.t_values), len(output_feval[0]))) #Output values corresponding to the self.t_values

					# Fill the first few values from the first function evaluation (corresponding to __init__)
					output_values[0:time_indices[0]] = output_feval[0]

					# Fill the rest of the values by matching the time of the time=step to a self.t_values time
					for step in range(len(feval_at_step)-1):
						# Might need one feval to span multiple self.t_values
						output_values[time_indices[step]:time_indices[step+1]] = output_feval[feval_at_step[step]-1]

					self.additional_out[key] = output_values #copy the adjusted array back onto the additional_out attribute

				except TypeError:
					output_values = np.zeros(len(self.t_values)) #Output values corresponding to the self.t_values

					# Fill the first few values from the first function evaluation (corresponding to __init__)
					output_values[0:time_indices[0]] = output_feval[0]

					# Fill the rest of the values by matching the time of the time=step to a self.t_values time
					for step in range(len(feval_at_step)-1):
						# Might need one feval to span multiple self.t_values
						output_values[time_indices[step]:time_indices[step+1]] = output_feval[feval_at_step[step]-1]

					self.additional_out[key] = output_values #copy the adjusted array back onto the additional_out attribute

			if refuelling_rate > 0: 
				self.Nzk = np.zeros((self.Z+1,)) #m^-3 Always need to start system is g.s. for Prad_tau calculation 
				self.Nzk[0] = 1e17 #m^-3 - start in g.s. 
			else: 
				# self.Nzk = result[-1,:] # Can use previous result to try speed up evaluation if not calculating Prad(tau)
				pass # However, this is found to result in odd numerical behaviour
				


		return result

	def scanTempCREquilibrium(self, output_index = -1):

		additional_out = {}
		for key in self.additional_out.keys():
			additional_out[key] = []

		results = np.zeros((len(self.Te_values),self.Z+1))

		for Te_iterator in range(len(self.Te_values)):
			self.reset_additional_out()
			Te = self.Te_values[Te_iterator]

			print("Evaluating test {} of {}".format(Te_iterator, len(self.Te_values)))

			result = self.timeIntegrate(Te, self.Ne_const)

			results[Te_iterator,:] = result[output_index,:]

			for key in self.additional_out_keys:
				additional_out[key].append(self.additional_out[key][output_index]) #Take the last time slice if using default output_index

		self.additional_out = additional_out #Replace the additional_out with the end-values

		return results

	def scanTempRefuelling(self):
		additional_out = {}
		refuelling_out = {}
		if len(self.Ne_tau_values) == 0:
			raise RuntimeError("Need non-zero list of Ne_tau values to calculate refuelling-dependant results")

		for key in self.additional_out.keys():
			additional_out[key] = []
			refuelling_out[key] = []

		results = np.zeros((len(self.Te_values),len(self.Ne_tau_values),self.Z+1))

		for Te_iterator in range(len(self.Te_values)):
			Te = self.Te_values[Te_iterator]

			print("Evaluating test {} of {}".format(Te_iterator, len(self.Te_values)))

			refuelling_rates = (self.Ne_const/np.array(self.Ne_tau_values))

			for key in self.additional_out_keys:
					refuelling_out[key] = [] #Reset for each time slice

			for refuelling_index in range(len(refuelling_rates)):
				self.reset_additional_out()
				refuelling_rate = refuelling_rates[refuelling_index]

				result = self.timeIntegrate(Te, self.Ne_const, refuelling_rate=refuelling_rate)

				results[Te_iterator,refuelling_index,:] = result[-1,:]

				for key in self.additional_out_keys:
					refuelling_out[key].append(self.additional_out[key][-1]) #Take the last time slice
			
			for key in self.additional_out_keys:
				additional_out[key].append(refuelling_out[key]) #Append the list of refuelling-specific values

		self.additional_out = additional_out #Replace the additional_out with the end-values

		return results

	def scanDensityCREquilibrium(self):
		self.reset_additional_out()

		results = np.zeros((self.Z+1, len(self.Ne_values)))

		for Ne_iterator in range(len(self.Ne_values)):
			Ne = self.Ne_values[Ne_iterator]

			print("Evaluating test {} of {}".format(Ne_iterator, len(self.Ne_values)))

			result = self.timeIntegrate(self.Te_const, Ne, self.t_values)

			results[:,Ne_iterator] = result[-1,:]

			self.Nzk = result[-1,:]

		return results

# End class methods

# Plotting methods
def plotResultFromDensityEvolution(solver, result, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none", align_yticks = True, show=False):
	fig, ax1 = plt.subplots()
	for k in range(solver.Z+1):
		if k == 0:
			ax1.semilogx(solver.t_values, result[:,k], label="{}".format("g.s."))
		else:
			ax1.semilogx(solver.t_values, result[:,k], label="{}+".format(k))

	# ax1.semilogx(solver.t_values, np.sum(result[:,:],1), label="Total")
	ax1.set_ylim(0, 1e17)
	ax1.set_xlabel(r'Time (s)')
	ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
	# plt.title('Time evolution of ionisation stages')
	ax1.tick_params('y', colors = 'b')
	# ax1.legend()


	ax1.grid(which=grid, axis='both')

	if plot_power:
		ax2 = ax1.twinx()
		scaled_Prad = np.array(solver.additional_out['Prad'])*1e-3
		print(scaled_Prad[-1])
		ax2.semilogx(solver.t_values, scaled_Prad,'k-.',label=r'$P_{rad}$',linewidth=1)
		scaled_Pcool = np.array(solver.additional_out['Pcool'])*1e-3
		ax2.semilogx(solver.t_values, scaled_Pcool,'r-.',label=r'$P_{cool}$',linewidth=1)
		ax2.set_ylim(min(min(scaled_Prad),min(scaled_Pcool)), max(max(scaled_Prad),max(scaled_Pcool)))
		ax2.set_ylabel(r'$Power$ (kW $m^{-3}$)')
		ax2.tick_params('y', colors='k')
		# ax2.legend(loc=0)

	h1, l1 = ax1.get_legend_handles_labels()
	h2, l2 = ax2.get_legend_handles_labels()
	ax1.legend(h1+h2, l1+l2, loc=0)
	
	ax1.set_xscale(x_axis_scale)
	ax1.set_yscale(y_axis_scale)
	ax1.set_xlim(min(solver.t_values), max(solver.t_values))
	if plot_power:
		ax2.set_yscale(y_axis_scale)
		if align_yticks:
			ax2.set_yticks(np.linspace(ax2.get_yticks()[0],ax2.get_yticks()[-1],len(ax1.get_yticks())))

	if show:
		plt.show()
	return fig

def plotScanTempCR_Dens(solver, reevaluate_scan=False, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none", align_yticks = True, show=False):
	from scipy.interpolate import interp1d

	fig, (ax1,ax2) = plt.subplots(2,sharex=True)

	if reevaluate_scan:

		scan_temp = solver.scanTempCREquilibrium()

		with open('python_results/scanTempCREquilibrium({}_at_{},res+{})-INTEG_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(scan_temp, handle, protocol=pickle.HIGHEST_PROTOCOL)
		with open('python_results/scanTempCREquilibrium({}_at_{},res+{})-ADDIT_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(solver.additional_out, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/scanTempCREquilibrium({}_at_{},res+{})-INTEG_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			scan_temp = pickle.load(handle)
		with open('python_results/scanTempCREquilibrium({}_at_{},res+{})-ADDIT_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			solver.additional_out = pickle.load(handle)

	# Extract the density curves from Florido 2009 paper. Plot the carbon result against these
	FLORIDO_x_eval = np.logspace(-1,3,100) #Points to return function values for
	FLORIDO_y_mean = []
	FLORIDO_y_diff  = []

	for k in range(solver.Z+1):

		# For each charge state, load the results from each model
		ABAKO_dens_data = np.loadtxt('python_results/ABAKO_C_{}.csv'.format(k),delimiter=', ')
		ABAKO_x = np.array(ABAKO_dens_data[:,0])
		ABAKO_y = np.array(ABAKO_dens_data[:,1])
		ABAKO_dens_interp = interp1d(ABAKO_x,ABAKO_y,bounds_error=False)

		ATOMIC_dens_data = np.loadtxt('python_results/ATOMIC_C_{}.csv'.format(k),delimiter=', ')
		ATOMIC_x = np.array(ATOMIC_dens_data[:,0])
		ATOMIC_y = np.array(ATOMIC_dens_data[:,1])
		ATOMIC_dens_interp = interp1d(ATOMIC_x,ATOMIC_y,bounds_error=False)

		# Linearly interpolate to the points specified
		ABAKO_y_interp = ABAKO_dens_interp(FLORIDO_x_eval)
		ATOMIC_y_interp = ATOMIC_dens_interp(FLORIDO_x_eval)

		# print(np.mean([ABAKO_y_interp,ATOMIC_y_interp],0))
		# print(np.std([ABAKO_y_interp,ATOMIC_y_interp],0))
		# Will give NaN if off grid
		k_state_mean = np.mean([ABAKO_y_interp,ATOMIC_y_interp],0)
		FLORIDO_y_mean.append(k_state_mean)
		# Find the std at the points specified
		k_state_diff = np.absolute(ABAKO_y_interp-ATOMIC_y_interp)/2
		FLORIDO_y_diff.append(k_state_diff)

		if k == 0:
			ax2.semilogx(FLORIDO_x_eval, ABAKO_y_interp, '{}-.'.format('C'+str(k)))
			ax2.semilogx(FLORIDO_x_eval, ATOMIC_y_interp, '{}-'.format('C'+str(k)), label="{}".format("g.s."))
		else:
			ax2.semilogx(FLORIDO_x_eval, ABAKO_y_interp, '{}-.'.format('C'+str(k)))
			ax2.semilogx(FLORIDO_x_eval, ATOMIC_y_interp, '{}-'.format('C'+str(k)), label=r"$C^{{{}}}$".format(str(k)+"+"))

	FLORIDO_y_mean = np.array(FLORIDO_y_mean)
	FLORIDO_y_diff = np.array(FLORIDO_y_diff)

	total_density = np.sum(scan_temp[-1,:],0)

	for k in range(solver.Z+1):
		if k == 0:
			ax1.semilogx(solver.Te_values, scan_temp[:,k]/total_density, label="{}".format("g.s."))
		else:
			ax1.semilogx(solver.Te_values, scan_temp[:,k]/total_density, label=r"$C^{{{}}}$".format(str(k)+"+"))
	
	# for k in range(solver.Z+1):
	# 	# ax1.errorbar(FLORIDO_x_eval, FLORIDO_y_mean[k,:]*1e17, yerr=FLORIDO_y_diff[k,:]*1e17, fmt='{}.'.format('C'+str(k)), ecolor='{}'.format('C'+str(k)), capthick=1, capsize=3)
	# 	ax1.errorbar(FLORIDO_x_eval, FLORIDO_y_mean[k,:]*1e17, yerr=FLORIDO_y_diff[k,:]*1e17, fmt='{}.'.format('C'+str(k)), ecolor='{}'.format('C'+str(k)), elinewidth=1 ,capthick=1, capsize=0)

	# plt.semilogx(solver.Te_values, np.sum(scan_temp[:,:],1), label="Total")

	# ax1.set_xlabel(r'Plasma temperature (eV)')
	ax2.set_xlabel(r'Plasma temperature (eV)')
	# ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
	# ax1.tick_params('y', colors = 'b')
	fig.text(0.0, 0.5, r'Relative density of stage', va='center', rotation='vertical')

	# ax1.set_xlim(min(solver.Te_values), max(solver.Te_values))
	ax1.set_xlim(10**-0.5, 10**3.5)
	ax1.set_ylim([0, 1])
	ax2.set_ylim([0, 1])

	ax1.grid(which=grid, axis='both')
	ax2.grid(which=grid, axis='both')

	if plot_power:
		ax1_twin = ax1.twinx()
		scaled_power = np.array(solver.additional_out['Prad'])*1e-3
		ax1_twin.semilogx(solver.Te_values, scaled_power,'k-.',label=r'$P_{rad}$',linewidth=1)
		ax1_twin.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
		ax1_twin.tick_params('y', colors='k')
		ax1_twin.set_ylim(0,)

		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax1_twin.get_legend_handles_labels()
		ax1.legend(h1+h2, l1+l2, loc=0)
	else:
		ax2.legend(loc=4)
	
	ax1.set_xscale(x_axis_scale)
	ax1.set_yscale(y_axis_scale)
	ax2.set_yscale(y_axis_scale)

	if plot_power:
		ax1_twin.set_yscale(y_axis_scale)
		if align_yticks:
			ax1_twin.set_yticks(np.linspace(ax1_twin.get_yticks()[0],ax1_twin.get_yticks()[-1],len(ax1.get_yticks())))

	fig.set_size_inches(6.268, 3.52575, forward=True)
	plt.tight_layout()
	vals1 = ax1.get_yticks()
	ax1.set_yticklabels(['{:3.0f}%'.format(y*100) for y in vals1])
	vals2 = ax2.get_yticks()
	ax2.set_yticklabels(['{:3.0f}%'.format(y*100) for y in vals2])

	if show:
		plt.show()
	return fig

def plotScanTempCR_Prad_tau(solver, x_axis_scale = "log", y_axis_scale = "log", grid = "none", show=False, ylim = [1e-37, 1e-30]):
	from scipy.interpolate import interp1d
	
	if reevaluate_scan:
		scan_temp_refuelling = solver.scanTempRefuelling()

		with open('python_results/scanTempRefuelling({}_at_{},res+{})-INTEG_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(scan_temp_refuelling, handle, protocol=pickle.HIGHEST_PROTOCOL)
		with open('python_results/scanTempRefuelling({}_at_{},res+{})-ADDIT_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(solver.additional_out, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/scanTempRefuelling({}_at_{},res+{})-INTEG_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			scan_temp_refuelling = pickle.load(handle)
		with open('python_results/scanTempRefuelling({}_at_{},res+{})-ADDIT_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			solver.additional_out = pickle.load(handle)

	# Post PSI radiation curves (4 datasets) are saved as Prad#.csv in python_results

	# Extract the radiation curves from Post PSI paper. Plot the carbon result against these
	POST_x_eval = np.logspace(0.1,3,20) #Points to return function values for
	POST_y_eval = []
	for dataset in range(0,4):
		POST_Prad_data = np.loadtxt('python_results/POST_Prad_data_{}.csv'.format(dataset+1),delimiter=', ')
		
		x = np.array(POST_Prad_data[:,0])
		y = np.array(POST_Prad_data[:,1])/1e12

		Prad_interp = interp1d(x,y,bounds_error=False)
		POST_y_eval.append(Prad_interp(POST_x_eval))

	POST_y_eval      = np.array(POST_y_eval)
	POST_y_mean      = np.nanmean(POST_y_eval,0)
	POST_y_std       = np.nanstd(POST_y_eval,0)
	POST_y_lower     = np.maximum(ylim[0], POST_y_mean - POST_y_std)
	POST_y_std_lower = POST_y_mean - POST_y_lower

	# Extract the radiation Ne_tau curves from Mavrin 2017 paper. Plot the carbon result against these
	MAVRIN_x_eval = np.logspace(0,4,25) #Points to return function values for
	MAVRIN_y_eval = {}
	for dataset in [16, 17, 18, 19]:
		MAVRIN_Ne_tau_data = np.loadtxt('python_results/MAVRIN_{}.csv'.format(dataset),delimiter=', ')
		
		x = np.array(MAVRIN_Ne_tau_data[:,0])
		y = np.array(MAVRIN_Ne_tau_data[:,1])

		Ne_tau_interp = interp1d(x,y,bounds_error=False)
		# print(dataset)
		# print(Ne_tau_interp(MAVRIN_x_eval))
		MAVRIN_y_eval[dataset] = np.array(Ne_tau_interp(MAVRIN_x_eval))

	# Construct the carbon cooling curve from I.Hutchinson
	HUTCHINSON_y = 2e-31 * (solver.Te_values/10.)**3 / ( 1.0 + (solver.Te_values/10.)**4.5 )

	# Plot the results for the specified ne_tau values
	fig, ax = plt.subplots()

	Prad = np.array(solver.additional_out['Prad'])

	for ne_tau_index in range(len(solver.Ne_tau_values)-1,-1,-1): #Plot in reverse
		ne_tau = solver.Ne_tau_values[ne_tau_index]
		if ne_tau == 1e19:
			# This line makes the resulting plot look cluttered, since it's basically CR
			continue
		colour_index = ne_tau_index -1
		if ne_tau_index == 4:
			colour_index = 0

		total_density = np.sum(scan_temp_refuelling[-1,ne_tau_index,:],0)

		# ax.loglog(MAVRIN_x_eval, MAVRIN_y_eval['low_T'], "oC{}".format(1),markersize=3,markerfacecolor='none')
		if ne_tau > 1e20:
			ax.loglog(solver.Te_values, Prad[:,ne_tau_index]/(total_density*solver.Ne_const),'k-.',label="ADAS-CR", linewidth=1)
		else:
			power = int(np.floor(np.log10(ne_tau)))
			factor = ne_tau/10**power
			if(factor == 1.0): #if ne_tau is a perfect power of 10
				ax.loglog(solver.Te_values, Prad[:,ne_tau_index]/(total_density*solver.Ne_const),"C{}".format(colour_index),label=r"$N_e\tau$=$10^{{{:d}}}$".format(power))
				ax.loglog(MAVRIN_x_eval, MAVRIN_y_eval[power], "oC{}".format(colour_index),markersize=3,markerfacecolor='none')
			else:
				ax.loglog(solver.Te_values, Prad[:,ne_tau_index]/(total_density*solver.Ne_const),label=r"$N_e\tau$={:.1f}$\times10^{{{:d}}}$".format(factor,power))
		# ax.loglog(solver.Te_values, Prad[:,ne_tau_index]/(total_density*solver.Ne_const),label=r"{:.2e}".format(ne_tau))

	ax.set_ylabel(r'$P_{rad}(\tau)/(N_e N_z)$ ($W m^3$)')
	ax.set_xlim(min(solver.Te_values), max(solver.Te_values))
	ax.set_ylim(ylim[0], ylim[1])

	# Plot the POST radiation curves
	ax.errorbar(POST_x_eval, POST_y_mean, yerr=[POST_y_std_lower, 2*POST_y_std], xerr=0, fmt='k.', ecolor='k', capthick=1, capsize=3, label='Post')
	# Plot the HUTCHINSON radiation curves
	ax.loglog(solver.Te_values, HUTCHINSON_y, 'r-.',label="Hutchinson",linewidth=1)

	ax.set_xlabel(r'Plasma temperature (eV)')
	plt.legend(loc=4)

	ax.grid(which=grid, axis='both')
	
	ax.set_xscale(x_axis_scale)
	ax.set_yscale(y_axis_scale)

	if show:
		plt.show()
	return fig

def plotScanTempCR_Prad_tau_noDiff(solver, x_axis_scale = "log", y_axis_scale = "log", grid = "none", show=False, ylim = [1e-37, 1e-30]):
	from scipy.interpolate import interp1d
	
	if reevaluate_scan:
		cropped_results = {}
		cropped_additional = {}
		for Ne_tau in solver.Ne_tau_values:
			print("For Ne_tau = {:.2e}".format(Ne_tau))
			tau = Ne_tau/solver.Ne_const

			time_index = np.searchsorted(solver.t_values, tau) #Find the output slice to return for
			if time_index == len(solver.t_values):
				time_index -= 1
				# Shift the CR time index back into the time values, since otherwise will get an out-of-bounds error

			results = solver.scanTempCREquilibrium(output_index=time_index)
			cropped_results[Ne_tau] = results
			cropped_additional[Ne_tau] = solver.additional_out

		with open('python_results/scanTemp_tau_noDiff({}_at_{},res+{})-INTEG_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(cropped_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
		with open('python_results/scanTemp_tau_noDiff({}_at_{},res+{})-ADDIT_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(cropped_additional, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/scanTemp_tau_noDiff({}_at_{},res+{})-INTEG_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			cropped_results = pickle.load(handle)
		with open('python_results/scanTemp_tau_noDiff({}_at_{},res+{})-ADDIT_results.pickle'.format(len(solver.Te_values),solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			cropped_additional = pickle.load(handle)

	# Post PSI radiation curves (4 datasets) are saved as Prad#.csv in python_results

	# Extract the radiation curves from Post PSI paper. Plot the carbon result against these
	POST_x_eval = np.logspace(0.1,3,20) #Points to return function values for
	POST_y_eval = []
	for dataset in range(0,4):
		POST_Prad_data = np.loadtxt('python_results/POST_Prad_data_{}.csv'.format(dataset+1),delimiter=', ')
		
		x = np.array(POST_Prad_data[:,0])
		y = np.array(POST_Prad_data[:,1])/1e12

		Prad_interp = interp1d(x,y,bounds_error=False)
		POST_y_eval.append(Prad_interp(POST_x_eval))

	POST_y_eval      = np.array(POST_y_eval)
	POST_y_mean      = np.nanmean(POST_y_eval,0)
	POST_y_std       = np.nanstd(POST_y_eval,0)
	POST_y_lower     = np.maximum(ylim[0], POST_y_mean - POST_y_std)
	POST_y_std_lower = POST_y_mean - POST_y_lower

	# Construct the carbon cooling curve from I.Hutchinson
	HUTCHINSON_y = 2e-31 * (solver.Te_values/10.)**3 / ( 1.0 + (solver.Te_values/10.)**4.5 )

	# Plot the results for the specified ne_tau values
	fig, ax = plt.subplots()

	for Ne_tau_index in range(len(solver.Ne_tau_values)-1,-1,-1):
		
		Ne_tau = solver.Ne_tau_values[Ne_tau_index]
		print("For Ne_tau = {:.2e}".format(Ne_tau))
		tau = Ne_tau/solver.Ne_const

		Prad = np.array(cropped_additional[Ne_tau]['Prad'])
		total_density = np.sum(cropped_results[Ne_tau][-1,:],0)

		if Ne_tau > 1e18:
			ax.semilogx(solver.Te_values, Prad[:]/(total_density*solver.Ne_const),'k-.',label="ADAS-CR", linewidth=1)
		else:
			power = int(np.floor(np.log10(Ne_tau)))
			factor = Ne_tau/10**power
			if(factor == 1.0): #if Ne_tau is a perfect power of 10
				ax.semilogx(solver.Te_values, Prad[:]/(total_density*solver.Ne_const),label=r"$N_e\tau$=$10^{{{:d}}}$".format(power))
			else:
				ax.semilogx(solver.Te_values, Prad[:]/(total_density*solver.Ne_const),label=r"$N_e\tau$={:.1f}$\times10^{{{:d}}}$".format(factor,power))
		# ax.semilogx(solver.Te_values, Prad[:,Ne_tau_index]/(total_density*solver.Ne_const),label=r"{:.2e}".format(Ne_tau))

	ax.set_ylabel(r'$P_{rad}(\tau)/(N_e N_z)$ ($W m^3$)')
	ax.set_xlim(min(solver.Te_values), max(solver.Te_values))
	ax.set_ylim(ylim[0], ylim[1])

	# Plot the POST radiation curves
	ax.errorbar(POST_x_eval, POST_y_mean, yerr=[POST_y_std_lower, 2*POST_y_std], xerr=0, fmt='k.', ecolor='k', capthick=1, capsize=3, label='Post')
	# Plot the HUTCHINSON radiation curves
	ax.semilogx(solver.Te_values, HUTCHINSON_y, 'r-.',label="Hutchinson",linewidth=1)

	ax.set_xlabel(r'Plasma temperature (eV)')
	plt.legend(loc=4)

	ax.grid(which=grid, axis='both')
	
	ax.set_xscale(x_axis_scale)
	ax.set_yscale(y_axis_scale)

	if show:
		plt.show()
	return fig

def plotTestTimeIntegrator(solver, reevaluate_scan=False, show=False):

	fig, (ax1, ax2) = plt.subplots(2, sharex=False)

	# Determine high-resolution comparison results to compare against
	t_values_hi_res = np.logspace(-10, 2, 10000)

	if reevaluate_scan:
		prev_t_values     = solver.t_values

		solver.t_values   = t_values_hi_res #Use high resolution t values for this evaluation
		error_analysis    = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
		comparison_values = error_analysis[-1,:]
		comparison_values = np.append(comparison_values,solver.additional_out['Prad'][-1])

		solver.t_values   = prev_t_values #Reset to original t values
		with open('python_results/error_analysis({},{})(res+{})-comparison_values.pickle'.format(solver.Te_const,solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(comparison_values, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/error_analysis({},{})(res+{})-comparison_values.pickle'.format(solver.Te_const,solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			comparison_values = pickle.load(handle)

	# Test whether specified time resolution affects the result

	time_test_values = np.round(np.logspace(1, 3, 20))
	if reevaluate_scan:
		time_test_results = [];
		prev_t_values = solver.t_values
		for time_iterator in range(len(time_test_values)):
			solver.reset_additional_out()

			time_test                             = time_test_values[time_iterator]
			print("Evaluating for time-resolution = {}pts".format(time_test))
			solver.t_values                       = np.logspace(-6, 2, time_test)
			error_analysis                        = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
			test_values                           = error_analysis[-1,:]
			test_values                           = np.append(test_values,solver.additional_out['Prad'][-1])
			time_test_results.append(comparison_values-test_values)

		solver.t_values   = prev_t_values #Reset to original t values	
		with open('python_results/error_analysis({},{})(res+{})-time_test_results.pickle'.format(solver.Te_const,solver.Ne_const,len(time_test_values)), 'wb') as handle:
			pickle.dump(time_test_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/error_analysis({},{})(res+{})-time_test_results.pickle'.format(solver.Te_const,solver.Ne_const,len(time_test_values)), 'rb') as handle:
			time_test_results = pickle.load(handle)

	ax1.plot(time_test_values, time_test_results)
	ax1.set_xlabel("Specified time steps")
	
	# Test whether the specified start time affects the result

	shift_test_values = np.linspace(-10,4, num=100)
	original_test_length = len(shift_test_values)
	if reevaluate_scan:	
		prev_t_values      = solver.t_values
		shift_test_values  = shift_test_values.tolist()
		shift_test_results = [];
		failed_test_values = [];
		for shift_iterator in range(len(shift_test_values)):
			solver.reset_additional_out()

			shift_test = shift_test_values[shift_iterator]
			print("Evaluating for shift = {}".format(shift_test))
			try:
				solver.t_values = np.logspace(shift_test, 5, 200)
				error_analysis  = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
				test_values     = error_analysis[-1,:]
				test_values     = np.append(test_values,solver.additional_out['Prad'][-1])
				shift_test_results.append((comparison_values - test_values)/comparison_values)
			except:
				print("Evaluation failed for shift = {}".format(shift_test))
				failed_test_values.append(shift_test)

		for shift_test in failed_test_values:
			shift_test_values.remove(shift_test)

		shift_test_results         = np.absolute(shift_test_results)
		shift_test_values          = np.array(shift_test_values)
		start_times                = np.power(10*np.ones_like(shift_test_values), shift_test_values)

		shift_test_data            = {}
		shift_test_data['results'] = shift_test_results
		shift_test_data['times']   = start_times
		solver.t_values            = prev_t_values #Reset to original t values
		with open('python_results/error_analysis({},{})(res+{})-shift_test_data.pickle'.format(solver.Te_const,solver.Ne_const,original_test_length), 'wb') as handle:
			pickle.dump(shift_test_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/error_analysis({},{})(res+{})-shift_test_data.pickle'.format(solver.Te_const,solver.Ne_const,original_test_length), 'rb') as handle:
			shift_test_data = pickle.load(handle)

	for k in range(solver.Z+2):
		if k == 0:
			ax2.loglog(shift_test_data['times'], shift_test_data['results'][:,k], label="{}".format("g.s."))
		elif k == solver.Z+1:
			ax2.loglog(shift_test_data['times'], shift_test_data['results'][:,k], label="{}".format(r"$P_{rad}$"))
		else:
			ax2.loglog(shift_test_data['times'], shift_test_data['results'][:,k], label="{}+".format(k))
	
	ax2.set_xlabel('Start time for evaluation (s)')
	fig.text(0.0, 0.5, r'Relative deviation from expected answer ($\Delta x/x$)', va='center', rotation='vertical')
	
	# ax2.legend(loc=0)
	plt.subplots_adjust(hspace=0.3, left=0.15)

	if show:
		plt.show()
	return fig

def findStddev(solver, reevaluate_scan=False):
	solver.reset_additional_out()
	# Determine high-resolution comparison results to compare against
	t_values_hi_res = np.logspace(-10, 2, 10000)

	if reevaluate_scan:
		prev_t_values     = solver.t_values

		solver.t_values   = t_values_hi_res #Use high resolution t values for this evaluation
		error_analysis    = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
		comparison_values = error_analysis[-1,:]
		comparison_values = np.append(comparison_values,solver.additional_out['Prad'][-1])

		solver.t_values   = prev_t_values #Reset to original t values
		with open('python_results/error_analysis({},{})(res+{})-comparison_values.pickle'.format(solver.Te_const,solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(comparison_values, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/error_analysis({},{})(res+{})-comparison_values.pickle'.format(solver.Te_const,solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			comparison_values = pickle.load(handle)

	samples_per_point = 50
	if reevaluate_scan:
		random_results = []
		store_Nzk = solver.Nzk
		for iterator in range(samples_per_point):
			solver.reset_additional_out()
			Nzk = np.zeros((solver.Z+1,))
			for k in range(solver.Z+1):
				Nzk[k] = random.random()*(10**random.uniform(1, 17))
			solver.Nzk = 1e17*Nzk/sum(Nzk)

			solver.t_values = np.logspace(-6, 2, 200)
			random_init = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
			random_values = random_init[-1,:]
			random_values = np.append(random_values,solver.additional_out['Prad'][-1])
			random_results.append(random_values)
		
		random_results = np.array(random_results)
		solver.Nzk = store_Nzk
		with open('python_results/error_analysis({},{})(res+{})-random_results.pickle'.format(solver.Te_const,solver.Ne_const,len(solver.t_values)), 'wb') as handle:
			pickle.dump(random_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open('python_results/error_analysis({},{})(res+{})-random_results.pickle'.format(solver.Te_const,solver.Ne_const,len(solver.t_values)), 'rb') as handle:
			random_results = pickle.load(handle)

	for k in range(solver.Z+2):
		mean = np.mean(random_results[:,k])
		stdev = np.std(random_results[:,k])
		stdev_norm = stdev/mean
		diff = (mean - comparison_values[k])/comparison_values[k]
		if k == 0:
			print("{:5} -> mean = {:.2e}, stdev = {:.2e}, stdev_norm = {:.2e}, mean_diff = {:.2e}".format("g.s.", mean, stdev, stdev_norm, diff))
		elif k == solver.Z+1:
			print("{:5} -> mean = {:.2e}, stdev = {:.2e}, stdev_norm = {:.2e}, mean_diff = {:.2e}".format("P_rad", mean, stdev, stdev_norm, diff))
		else:
			print("{:5} -> mean = {:.2e}, stdev = {:.2e}, stdev_norm = {:.2e}, mean_diff = {:.2e}".format(k, mean, stdev, stdev_norm, diff))

def plotErrorPropagation(solver, reevaluate_scan=False, show=False, plot='both', show_species=[], find_regression=False):

	stdev_Te = np.linspace(0,solver.Te_const/2,num=20)
	stdev_Ne = np.linspace(0,solver.Ne_const/2,num=20)
	samples_per_point = 50
	if reevaluate_scan:
		solver.t_values = np.logspace(-6, 2, 200)

		if plot in ['Te','both']:
			stdev_norm_Te = []
			for sigma in stdev_Te:
				random_results = []
				for iterator in range(samples_per_point):
					try:
						solver.reset_additional_out()
						Te = random.normalvariate(solver.Te_const, sigma)
						Nzk = np.zeros((solver.Z+1,))
						for k in range(solver.Z+1):
							Nzk[k] = random.random()*(10**random.uniform(1, 17))
						solver.Nzk = 1e17*Nzk/sum(Nzk)

						random_te = solver.timeIntegrate(Te, solver.Ne_const, 0)
						random_values = random_te[-1,:]
						random_values = np.append(random_values,solver.additional_out['Prad'][-1])
						if not(np.isnan(random_values).any()):
							random_results.append(random_values)
						else:
							print("NaN for Te = {}".format(Te))
					except:
						print("Error for Te = {}".format(Te))

				mean = np.mean(random_results,0)
				stdev = np.std(random_results,0)
				stdev_norm_Te.append(np.absolute(stdev/mean))

			stdev_norm_Te = np.array(stdev_norm_Te)
			with open('python_results/error_analysis({},{})(res+{})-stdev_norm_Te.pickle'.format(solver.Te_const,solver.Ne_const,len(stdev_Te)), 'wb') as handle:
				pickle.dump(stdev_norm_Te, handle, protocol=pickle.HIGHEST_PROTOCOL)

		if plot in ['Ne','both']:
			stdev_norm_Ne = []
			for sigma in stdev_Ne:
				random_results = []
				for iterator in range(samples_per_point):
					try:
						solver.reset_additional_out()
						Ne = random.normalvariate(solver.Ne_const, sigma)
						Nzk = np.zeros((solver.Z+1,))
						for k in range(solver.Z+1):
							Nzk[k] = random.random()*(10**random.uniform(1, 17))
						solver.Nzk = 1e17*Nzk/sum(Nzk)

						random_te = solver.timeIntegrate(solver.Te_const, Ne, 0)
						random_values = random_te[-1,:]
						random_values = np.append(random_values,solver.additional_out['Prad'][-1])
						if not(np.isnan(random_values).any()):
							random_results.append(random_values)
						else:
							print("NaN for Ne = {}".format(Ne))
					except:
						print("Error for Ne = {}".format(Ne))

				mean = np.mean(random_results,0)
				stdev = np.std(random_results,0)
				stdev_norm_Ne.append(np.absolute(stdev/mean))

			stdev_norm_Ne = np.array(stdev_norm_Ne)
			with open('python_results/error_analysis({},{})(res+{})-stdev_norm_Ne.pickle'.format(solver.Te_const,solver.Ne_const,len(stdev_Ne)), 'wb') as handle:
				pickle.dump(stdev_norm_Ne, handle, protocol=pickle.HIGHEST_PROTOCOL)

	else:
		with open('python_results/error_analysis({},{})(res+{})-stdev_norm_Te.pickle'.format(solver.Te_const,solver.Ne_const,len(stdev_Te)), 'rb') as handle:
				stdev_norm_Te = pickle.load(handle)
		with open('python_results/error_analysis({},{})(res+{})-stdev_norm_Ne.pickle'.format(solver.Te_const,solver.Ne_const,len(stdev_Ne)), 'rb') as handle:
				stdev_norm_Ne = pickle.load(handle)

	if plot is 'Te':
		fig, ax = plt.subplots()	
		for k in range(solver.Z+2):
			if k == 0 and 0 in show_species:
				ax.plot(stdev_Te/solver.Te_const, stdev_norm_Te[:,k], label="{}".format("g.s."))
			elif k == solver.Z+1:
				ax.plot(stdev_Te/solver.Te_const, stdev_norm_Te[:,k], label="{}".format(r"$P_{rad}$"))
			else:
				if(k in show_species):
					ax.plot(stdev_Te/solver.Te_const, stdev_norm_Te[:,k], label="{}+".format(k))
		ax.legend()
		ax.set_xlabel(r'Relative error in $T_e$')
		ax.set_ylabel(r'Relative error in parameter ($\sigma/\mu$)')
		vals = ax.get_yticks()
		ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in vals])
		ax.grid()
		ax.set_xlim(min(stdev_Te/solver.Te_const), max(stdev_Te/solver.Te_const))

		vals = ax.get_xticks()
		ax.set_xticklabels(['{:3.0f}%'.format(x*100) for x in vals])
		vals = ax.get_yticks()
		ax.set_yticklabels(['{:3.0f}%'.format(y*100) for y in vals])

	if plot is 'Ne':
		fig, ax = plt.subplots()	
		for k in range(solver.Z+2):
			if k == 0 and 0 in show_species:
				ax.plot(stdev_Ne/solver.Ne_const, stdev_norm_Ne[:,k], label="{}".format("g.s."))
			elif k == solver.Z+1:
				ax.plot(stdev_Ne/solver.Ne_const, stdev_norm_Ne[:,k], label="{}".format(r"$P_{rad}$"))
			else:
				if(k in show_species):
					ax.plot(stdev_Ne/solver.Ne_const, stdev_norm_Ne[:,k], label="{}+".format(k))
		ax.legend()
		ax.set_xlabel(r'Relative error in $N_e$')
		ax.set_ylabel(r'Relative error in parameter ($\sigma/\mu$)')
		vals = ax.get_yticks()
		ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in vals])
		ax.grid()
		ax.set_xlim(min(stdev_Ne/solver.Ne_const), max(stdev_Ne/solver.Ne_const))

		vals = ax.get_xticks()
		ax.set_xticklabels(['{:3.0f}%'.format(x*100) for x in vals])
		vals = ax.get_yticks()
		ax.set_yticklabels(['{:3.0f}%'.format(y*100) for y in vals])


	if plot is 'both':
		fig, (ax1, ax2) = plt.subplots(2, sharex = False)

		for k in range(solver.Z+2):
			if k == 0 and 0 in show_species:
				ax1.plot(stdev_Te/solver.Te_const, stdev_norm_Te[:,k], label="{}".format("g.s."))
			elif k == solver.Z+1:
				ax1.plot(stdev_Te/solver.Te_const, stdev_norm_Te[:,k], 'k-.', label="{}".format(r"$P_{rad}$"), linewidth=1)
			else:
				if k in show_species:
					ax1.plot(stdev_Te/solver.Te_const, stdev_norm_Te[:,k], label="{}+".format(k))
		ax1.legend()
		ax1.set_xlabel(r'Relative error in $T_e$')
		ax1.grid()

		for k in range(solver.Z+2):
			if k == 0 and 0 in show_species:
				ax2.plot(stdev_Ne/solver.Ne_const, stdev_norm_Ne[:,k], label="{}".format("g.s."))
			elif k == solver.Z+1:
				ax2.plot(stdev_Ne/solver.Ne_const, stdev_norm_Ne[:,k], 'k-.', label="{}".format(r"$P_{rad}$"), linewidth=1)
			else:
				if k in show_species:
					ax2.plot(stdev_Ne/solver.Ne_const, stdev_norm_Ne[:,k], label="{}+".format(k))
		ax2.legend()
		ax2.set_xlabel(r'Relative error in $N_e$')
		ax2.grid()

		ax1.set_xlim(min(stdev_Te/solver.Te_const), max(stdev_Te/solver.Te_const))
		ax2.set_xlim(min(stdev_Ne/solver.Ne_const), max(stdev_Ne/solver.Ne_const))

		if find_regression:
			crop_ind_Te = np.searchsorted(stdev_Te/solver.Te_const, 0.25)
			p_Te = np.polyfit((stdev_Te/solver.Te_const)[:crop_ind_Te], stdev_norm_Te[:crop_ind_Te,solver.Z+1],1)
			fit_regress_Te = np.polyval(p_Te, (stdev_Te/solver.Te_const)[:crop_ind_Te])
			ax1.plot((stdev_Te/solver.Te_const)[:crop_ind_Te], fit_regress_Te)
			# ax1.set_ylim(0,0.5)
			print("dPrad = {:.2f}dTe + {:.2f}".format(p_Te[0],p_Te[1]))

			crop_ind_Ne = np.searchsorted(stdev_Ne/solver.Ne_const, 0.5)
			p_Ne = np.polyfit((stdev_Ne/solver.Ne_const)[:crop_ind_Ne], stdev_norm_Ne[:crop_ind_Ne,solver.Z+1],1)
			fit_regress_Ne = np.polyval(p_Ne, (stdev_Ne/solver.Ne_const)[:crop_ind_Ne])
			ax2.plot((stdev_Ne/solver.Ne_const)[:crop_ind_Ne], fit_regress_Ne)
			print("dPrad = {:.2f}dNe + {:.2f}".format(p_Ne[0],p_Ne[1]))

		vals = ax1.get_xticks()
		ax1.set_xticklabels(['{:3.0f}%'.format(x*100) for x in vals])
		vals = ax2.get_xticks()
		ax2.set_xticklabels(['{:3.0f}%'.format(x*100) for x in vals])

		vals = ax1.get_yticks()
		ax1.set_yticklabels(['{:3.0f}%'.format(y*100) for y in vals])
		vals = ax2.get_yticks()
		ax2.set_yticklabels(['{:3.0f}%'.format(y*100) for y in vals])

		
		fig.text(0.0, 0.5, r'Relative error in parameter ($\sigma/\mu$)', va='center', rotation='vertical')

		plt.subplots_adjust(hspace=0.3, left=0.15)

	if show:
		plt.show()
	return fig

if __name__ == "__main__":

	# Control booleans
	reevaluate_scan           = False
	plot_solver_evolution     = False
	find_stddev               = False
	plot_test_time_integrator = False
	plot_error_propagation    = False
	plot_scan_temp_dens       = True
	plot_scan_temp_prad_tau   = False
	plot_scan_temp_prad_tau_noDiff = False #no diffusion

	SMALL_SIZE = 10
	MEDIUM_SIZE = 12
	BIGGER_SIZE = 14

	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

	impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

	solver = AtomicSolver(impurity_symbol)

	path_to_output = 'Figures/'
	# tight_layout() can take keyword arguments of pad, w_pad and h_pad. These control the extra padding around the figure border and between subplots. The pads are specified in fraction of fontsize.

	if plot_solver_evolution:
		solver_evolution = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
		plot_solver_evolution = plotResultFromDensityEvolution(solver, solver_evolution, plot_power = True, grid="major", show=False, y_axis_scale="linear")

		plot_solver_evolution.set_size_inches(6.268, 3.52575, forward=True)
		plt.tight_layout()
		plot_solver_evolution.savefig(path_to_output+"solver_evolution.pdf",bbox_inches = 'tight',pad_inches = 0.03,transparent=True)

	if find_stddev:
		findStddev(solver, reevaluate_scan = reevaluate_scan)

	if plot_test_time_integrator:
		plot_test_time_integrator = plotTestTimeIntegrator(solver, reevaluate_scan = reevaluate_scan)

		plot_test_time_integrator.set_size_inches(6.268, 3.52575, forward=True)
		plt.tight_layout()
		plot_test_time_integrator.savefig(path_to_output+"test_time_integrator.pdf",bbox_inches = 'tight',pad_inches = 0.03,transparent=True)

	if plot_error_propagation:
		plot_error_propagation = plotErrorPropagation(solver, show_species=[4, 5], reevaluate_scan = reevaluate_scan, find_regression=False)

		plot_error_propagation.set_size_inches(6.268, 3.52575, forward=True)
		plt.tight_layout()
		plot_error_propagation.savefig(path_to_output+"error_propagation.pdf",bbox_inches = 'tight',pad_inches = 0.03,transparent=True)

	if plot_scan_temp_dens:
		if not(impurity_symbol is b'c'):
			raise NotImplementedError('Prad_tau plot comparison data is for Carbon. Will need to add data for species {}'.format(str(impurity_symbol,'utf-8')))
		plot_scan_temp_dens = plotScanTempCR_Dens(solver, grid="major", plot_power=False, reevaluate_scan = reevaluate_scan)

		plot_scan_temp_dens.set_size_inches(6.268, 3.52575, forward=True)
		plt.tight_layout()
		plot_scan_temp_dens.savefig(path_to_output+"plot_scan_temp_dens.pdf",bbox_inches = 'tight',pad_inches = 0.03,transparent=True)

	if plot_scan_temp_prad_tau:
		solver.Ne_tau_values = [1e30, 1e19, 1e18, 1e17, 1e16] #m^-3 s, values to return Prad(tau) for, for comparison against Mavrin 2017
		if not(impurity_symbol is b'c'):
			raise NotImplementedError('Prad_tau plot comparison data is for Carbon. Will need to add data for species {}'.format(str(impurity_symbol,'utf-8')))
		plot_scan_temp_prad_tau = plotScanTempCR_Prad_tau(solver, grid="major")

		plot_scan_temp_prad_tau.set_size_inches(6.268, 3.52575, forward=True)
		plt.tight_layout()
		plot_scan_temp_prad_tau.savefig(path_to_output+"plot_scan_temp_prad_tau.pdf",bbox_inches = 'tight',pad_inches = 0.03,transparent=True)

	if plot_scan_temp_prad_tau_noDiff:
		solver.Ne_tau_values = [1e30, 1e17, 1e16, 1e15, 1e14] #m^-3 s, values to return Prad(tau) for, for comparison against Stangeby 2000
		if not(impurity_symbol is b'c'):
			raise NotImplementedError('Prad_tau plot comparison data is for Carbon. Will need to add data for species {}'.format(str(impurity_symbol,'utf-8')))
		plot_scan_temp_prad_tau_noDiff = plotScanTempCR_Prad_tau_noDiff(solver, grid="major")

		plot_scan_temp_prad_tau_noDiff.set_size_inches(6.268, 3.52575, forward=True)
		plt.tight_layout()
		plot_scan_temp_prad_tau_noDiff.savefig(path_to_output+"plot_scan_temp_prad_tau_noDiff.pdf",bbox_inches = 'tight',pad_inches = 0.03,transparent=True)

	plt.show()
	

























