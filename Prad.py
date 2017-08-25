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

		# Control booleans
		self.plot_solver_evolution   = False
		self.reevaluate_scan_temp    = False
		self.plot_scan_temp_dens     = False
		self.plot_scan_temp_prad     = False
		self.plot_scan_temp_prad_tau = False

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
		self.additional_out_keys = ['Prad', 'dNzk'] #Keys to record data for
	
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

		# print("Refuel: ", sum(Nzk)*refuelling_rate)
		# print("No refuelling", dNzk)
		# dNzk_store = dNzk
		fraction_in_stage = Nzk/sum(Nzk)
		# Add neutrals at a rate of tau^-1
		dNzk[0] += sum(Nzk)*refuelling_rate
		# Remove other stages based on their density
		dNzk -= sum(Nzk)*refuelling_rate*fraction_in_stage
		# print("Refuelling", dNzk)
		# print("Ratio", dNzk/dNzk_store)

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
			(result, output_dictionary) = odeint(self.evolveDensity, self.Nzk, self.t_values, args=(self, Te, Ne), printmessg=False, full_output=True)
		else:
			# Refuelling case
			(result, output_dictionary) = odeint(self.evolveDensity_withRefuelling, self.Nzk, self.t_values, args=(self, Te, Ne, refuelling_rate), printmessg=False, full_output=True)
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
				# Can use previous result to try speed up evaluation if not calculating Prad(tau) 
				self.Nzk = result[-1,:] 

		return result

	# def calculateTimeIntegratedPower(self, Ne, integration_endpoints = []):
		# # To be called after timeIntegrate (or another function which returns 'Prad' additional_out
		# # corresponding to self.t_values)
		# time_integrated_power = []

		# for ne_tau in integration_endpoints:
		# 	# Reset sliced arrays to original
		# 	Prad = self.additional_out['Prad']
		# 	Prad_times = self.t_values
		# 	# Extract next value of tau
		# 	tau = ne_tau/Ne
		# 	# Find which index is closest to the tau value
		# 	time_index = np.searchsorted(self.t_values, tau, side='left')
		# 	time_left = self.t_values[time_index-1]
		# 	time_right = self.t_values[time_index]
		# 	time_diff = time_right - time_left
		# 	# Linearly interpolate between two time values
		# 	t_next_weighting = (tau - time_left)/(time_right - time_left)
		# 	Prad_t_exact = Prad[time_index-1] * (1 - t_next_weighting) + Prad[time_index] * t_next_weighting
		# 	# Slice the arrays, and then append the interpolated values
		# 	Prad_times = Prad_times[:time_index-1]
		# 	np.append(Prad_times, tau)
		# 	Prad = Prad[:time_index-1]
		# 	np.append(Prad, Prad_t_exact)
		# 	#Calculate the definite integral via Simpson's rule
		# 	time_integrated_power.append(simps(Prad, Prad_times))

		# return time_integrated_power

# Scan methods
	def scanTempCREquilibrium(self):

		additional_out = {}
		for key in self.additional_out.keys():
			additional_out[key] = []

		results = np.zeros((len(self.Te_values),self.Z+1))

		for Te_iterator in range(len(self.Te_values)):
			self.reset_additional_out()
			Te = self.Te_values[Te_iterator]

			print("Evaluating test {} of {}".format(Te_iterator, len(self.Te_values)))

			result = self.timeIntegrate(Te, self.Ne_const)

			results[Te_iterator,:] = result[-1,:]

			for key in self.additional_out_keys:
				additional_out[key].append(self.additional_out[key][-1]) #Take the last time slice

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

# Plotting methods
	def plotResultFromDensityEvolution(self, result, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none"):
		fig, ax1 = plt.subplots()
		for k in range(solver.Z+1):
			if k == 0:
				ax1.semilogx(self.t_values, result[:,k], label="{}".format("g.s."))
			else:
				ax1.semilogx(self.t_values, result[:,k], label="{}+".format(k))

		ax1.semilogx(self.t_values, np.sum(result[:,:],1), label="Total")

		ax1.set_xlabel(r'Time (s)')
		ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
		# plt.title('Time evolution of ionisation stages')
		ax1.tick_params('y', colors = 'b')
		ax1.legend()

		ax1.set_xlim(min(self.t_values), max(self.t_values))

		ax1.grid(which=grid, axis='both')

		if plot_power:
			ax2 = ax1.twinx()
			scaled_power = np.array(self.additional_out['Prad'])*1e-3
			ax2.semilogx(self.t_values, scaled_power,'k--',label=r'$P_{rad}$')
			ax2.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
			ax2.tick_params('y', colors='k')
			ax2.legend(loc=0)
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)
		if plot_power:
			ax2.set_yscale(y_axis_scale)

		plt.show()

	def plotScanTempCR_Dens(self, results, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none"):
		
		fig, ax1 = plt.subplots()

		for k in range(self.Z+1):
			if k == 0:
				ax1.semilogx(self.Te_values, results[:,k], label="{}".format("g.s."))
			else:
				ax1.semilogx(self.Te_values, results[:,k], label="{}+".format(k))

		plt.semilogx(self.Te_values, np.sum(results[:,:],1), label="Total")

		total_density = np.sum(results[-1,:],0)
		ax1.set_ylim([1e-3*total_density, total_density])
		ax1.set_xlabel(r'Plasma temperature (eV)')
		ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
		# plt.title(r'Ionisation stage at C.R. as $f(T_e)$')
		ax1.tick_params('y', colors = 'b')
		plt.legend()

		ax1.set_xlim(min(self.Te_values), max(self.Te_values))

		ax1.grid(which=grid, axis='both')

		if plot_power:
			ax2 = ax1.twinx()
			scaled_power = np.array(self.additional_out['Prad'])*1e-3
			ax2.semilogx(self.Te_values, scaled_power,'k--',label=r'$P_{rad}$')
			ax2.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
			ax2.tick_params('y', colors='k')
			ax2.legend(loc=0)
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)
		if plot_power:
			ax2.set_yscale(y_axis_scale)

		plt.show()

	def plotScanTempCR_Prad(self, results, x_axis_scale = "log", y_axis_scale = "log", grid = "none"):
		
		fig, ax1 = plt.subplots()

		total_density = np.sum(results[-1,:],0)
		scaled_power = np.array(self.additional_out['Prad'])

		compare_to_Post_PSI = True
		if compare_to_Post_PSI:
			convertCM = 1e-6 #convert from cm^3 to m^3

			ax1.semilogx(self.Te_values, scaled_power/(total_density*self.Ne_const*convertCM),'k--',label=r'$L$')
			
			ax1.set_ylabel(r'$L$ (W $cm^3$)')

			ax1.set_xlim(1, 1e3)
			ax1.set_ylim(1e-29, 1e-24)

		else:
			ax1.semilogx(self.Te_values, scaled_power/(total_density*self.Ne_const),'k--',label=r'$L$')

			ax1.set_ylabel(r'$L$ (W $m^3$)')
			
			ax1.set_xlim(min(self.Te_values), max(self.Te_values))

		ax1.set_xlabel(r'Plasma temperature (eV)')
		# plt.legend()

		ax1.grid(which=grid, axis='both')
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)

		plt.show()

	def plotScanTempCR_Prad_tau(self, refuelling_results, x_axis_scale = "log", y_axis_scale = "log", grid = "none"):
		
		fig, ax1 = plt.subplots()

		
		Prad = np.array(self.additional_out['Prad'])

		for ne_tau_index in range(len(self.Ne_tau_values)):
			ne_tau = self.Ne_tau_values[ne_tau_index]

			total_density = np.sum(refuelling_results[-1,ne_tau_index,:],0)

			ax1.semilogx(self.Te_values, Prad[:,ne_tau_index]/(total_density*self.Ne_const),label="{:.1e}".format(ne_tau))

			ax1.set_ylabel(r'$P_{rad}(\tau)/(n_e n_z)$ ($W m^3$)')
			
			ax1.set_xlim(min(self.Te_values), max(self.Te_values))

		ax1.set_xlabel(r'Plasma temperature (eV)')
		plt.legend()

		ax1.grid(which=grid, axis='both')
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)

		plt.show()

if __name__ == "__main__":

	import matplotlib.pyplot as plt

	impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

	solver = AtomicSolver(impurity_symbol)

	solver.plot_solver_evolution   = False
	solver.reevaluate_scan_temp    = True
	solver.plot_scan_temp_dens     = False
	solver.plot_scan_temp_prad     = False
	solver.plot_scan_temp_prad_tau = True

	solver.Ne_tau_values = [1e30, 1e17, 1e16, 1e15, 1e14] #m^-3 s, values to return Prad(tau) for
	# solver.Ne_tau_values = np.logspace(12, 20, 7) #m^-3 s, values to return Prad(tau) for
	solver.Ne_tau_values = np.append(solver.Ne_tau_values, 1e30)

	# solver.Te_values = np.logspace(0, 3, 20) #eV

	if solver.plot_solver_evolution:
		solver_evolution = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 1e14)
		# solver.plotResultFromDensityEvolution(solver_evolution, plot_power = True, grid="major")
		print(solver.additional_out['Prad'][-1])
	
	if solver.plot_scan_temp_dens or solver.plot_scan_temp_prad:
		if solver.reevaluate_scan_temp:

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

		if solver.plot_scan_temp_dens:
			solver.plotScanTempCR_Dens(scan_temp, grid="major")
		if solver.plot_scan_temp_prad:
			solver.plotScanTempCR_Prad(scan_temp, grid="major")

	if solver.plot_scan_temp_prad_tau:
		if solver.reevaluate_scan_temp:

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

		if solver.plot_scan_temp_prad_tau:
			solver.plotScanTempCR_Prad_tau(scan_temp_refuelling, grid="major")



	

























