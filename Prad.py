# Program name: atomicpp/Prad.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 11 August 2017
#
# Use the atomic++ module to evaluate the rate-coefficients from OpenADAS

import numpy as np
from atomicpp import atomicpy
from scipy.integrate import odeint #ODEPACK, for numerical integration
from scipy.integrate import simps #Simpson's rule, for definite integrals
import pickle


class AtomicSolver(object):

	def __init__(self, impurity_symbol):

		self.impurity = atomicpy.PyImpuritySpecies(impurity_symbol)

		self.Z = self.impurity.get_atomic_number()

		self.impurity_derivatives = atomicpy.PyRateEquations(self.impurity)
		self.impurity_derivatives.setThresholdDensity(-1.0) #Don't use a threshold density at first
		self.impurity_derivatives.setDominantIonMass(1.0)

		self.Te  = 50 #eV
		self.Ne  = 1e19 #m^-3
		self.Vi  = 0 #m/s
		self.Nn  = 0 #m^-3
		self.Vn  = 0 #m/s
		self.Nzk = np.zeros((self.Z+1,)) #m^-3
		self.Nzk[0] = 1e17 #m^-3 - start in g.s.
		# self.Nzk = np.array([1e17, 0, 0, 0, 0, 0, 0]) 
		self.Vzk = np.zeros((self.Z+1,)) #m/s

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

		# self.additional_out['Prad'].append(derivative_struct['Prad'])
		# self.additional_out['dNzk'].append(derivative_struct['dNzk'])

		return dNzk

	def reset_additional_out(self):
		include_Prad_tau = ('Prad_tau' in self.additional_out.keys())
		self.additional_out = {'Prad':[], 'Pcool':[], 'dNzk':[], 'F_zk':[], 'dNe':[], 'F_i':[], 'dNn':[], 'F_n':[]}
		if include_Prad_tau:
			self.additional_out['Prad_tau'] = []

	def solveCollisionalRadiativeEquilibrium(self, Te, Ne, t_values = np.logspace(-6, 2, 200), ne_tau_values = []):
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		Vzk = self.Vzk

		print("Te = {}eV, Ne = {}/m3".format(Te, Ne))
		
		(result, output_dictionary) = odeint(self.evolveDensity, self.Nzk, t_values, args=(self,Te, Ne), printmessg=False, full_output=True)

		feval_at_step = output_dictionary['nfe'] #function evaluations at the time-step
		time_at_step = output_dictionary['tcur'] #time at the time-step

		time_indices = np.searchsorted(t_values, time_at_step, side='left') #find how the time-steps are distributed. Usually close but not 1 to 1 with t_values

		for key, value in self.additional_out.items():
			if value and not(key is 'Prad_tau'): #if list isn't empty - i.e. additional output has been recorded for this key
				output_feval = value #copy the output evaluated at each time-step

				try:
					# If the additional_out has a length (i.e. is an array)
					output_values = np.zeros((len(t_values), len(output_feval[0]))) #Output values corresponding to the t_values

					# Fill the first few values from the first function evaluation (corresponding to __init__)
					output_values[0:time_indices[0]] = output_feval[0]

					# Fill the rest of the values by matching the time of the time=step to a t_values time
					for step in range(len(feval_at_step)-1):
						# Might need one feval to span multiple t_values
						output_values[time_indices[step]:time_indices[step+1]] = output_feval[feval_at_step[step]-1]

					self.additional_out[key] = output_values #copy the adjusted array back onto the additional_out attribute

				except TypeError:
					output_values = np.zeros(len(t_values)) #Output values corresponding to the t_values

					# Fill the first few values from the first function evaluation (corresponding to __init__)
					output_values[0:time_indices[0]] = output_feval[0]

					# Fill the rest of the values by matching the time of the time=step to a t_values time
					for step in range(len(feval_at_step)-1):
						# Might need one feval to span multiple t_values
						output_values[time_indices[step]:time_indices[step+1]] = output_feval[feval_at_step[step]-1]

					self.additional_out[key] = output_values #copy the adjusted array back onto the additional_out attribute

			if len(ne_tau_values) > 0:
				Prad_tau = solver.calculatePradTau(Ne, t_values=t_values, ne_tau_values=ne_tau_values)
				# Must be initialised to [] in the calling function
				self.additional_out['Prad_tau'] = Prad_tau
				# N.b. will need to take all indices from Prad_tau, not just -1 (already have sliced the array)
				self.Nzk = np.zeros((self.Z+1,)) #m^-3 Always need to start system is g.s. for Prad_tau calculation
				self.Nzk[0] = 1e17 #m^-3 - start in g.s.
			else:
				# Can use previous result to try speed up evaluation if not calculating Prad(tau)
				solver.Nzk = result[-1,:]

		return result

	def calculatePradTau(self, Ne, t_values = np.logspace(-6, 2, 200), ne_tau_values = []):
		# To be called after solveCollisionalRadiativeEquilibrium (or another function which returns 'Prad' additional_out
		# corresponding to t_values)
		Prad_tau = []

		for ne_tau in ne_tau_values:
			# Reset sliced arrays to original
			Prad = self.additional_out['Prad']
			Prad_times = t_values
			# Extract next value of tau
			tau = ne_tau/Ne
			# Find which index is closest to the tau value
			time_index = np.searchsorted(t_values, tau, side='left')
			time_left = t_values[time_index-1]
			time_right = t_values[time_index]
			time_diff = time_right - time_left
			# Linearly interpolate between two time values
			t_next_weighting = (tau - time_left)/(time_right - time_left)
			Prad_t_exact = Prad[time_index-1] * (1 - t_next_weighting) + Prad[time_index] * t_next_weighting
			# Slice the arrays, and then append the interpolated values
			Prad_times = Prad_times[:time_index-1]
			np.append(Prad_times, tau)
			Prad = Prad[:time_index-1]
			np.append(Prad, Prad_t_exact)
			#Calculate the definite integral via Simpson's rule
			Prad_tau.append(simps(Prad, Prad_times))

		return Prad_tau

	def plotResultFromDensityEvolution(self, result, t_values, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none"):
		fig, ax1 = plt.subplots()
		for k in range(solver.Z+1):
			if k == 0:
				ax1.semilogx(t_values, result[:,k], label="{}".format("g.s."))
			else:
				ax1.semilogx(t_values, result[:,k], label="{}+".format(k))

		ax1.semilogx(t_values, np.sum(result[:,:],1), label="Total")

		ax1.set_xlabel(r'Time (s)')
		ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
		# plt.title('Time evolution of ionisation stages')
		ax1.tick_params('y', colors = 'b')
		ax1.legend()

		ax1.set_xlim(min(t_values), max(t_values))

		ax1.grid(which=grid, axis='both')

		if plot_power:
			ax2 = ax1.twinx()
			scaled_power = np.array(self.additional_out['Prad'])*1e-3
			ax2.semilogx(t_values, scaled_power,'k--',label=r'$P_{rad}$')
			ax2.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
			ax2.tick_params('y', colors='k')
			ax2.legend(loc=0)
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)
		if plot_power:
			ax2.set_yscale(y_axis_scale)

		plt.show()

	def scanTempCREquilibrium(self, Te_values, Ne_const, t_values = np.logspace(-6, 2, 200), ne_tau_values = []):

		additional_out = {}
		for key in self.additional_out.keys():
			additional_out[key] = []
		if len(ne_tau_values) > 0:
			self.additional_out['Prad_tau'] = []
			additional_out['Prad_tau'] = []

		results = np.zeros((len(Te_values),self.Z+1))

		for Te_iterator in range(len(Te_values)):
			self.reset_additional_out()
			Te = Te_values[Te_iterator]

			print("Evaluating test {} of {}".format(Te_iterator, len(Te_values)))

			result = self.solveCollisionalRadiativeEquilibrium(Te, Ne_const, t_values=t_values, ne_tau_values=ne_tau_values)

			results[Te_iterator,:] = result[-1,:]

			for key in self.additional_out_keys: #Never put Prad_tau in additional_out_keys!
				additional_out[key].append(self.additional_out[key][-1]) #Take the last time slice
			if len(ne_tau_values) > 0:
				additional_out['Prad_tau'].append(self.additional_out['Prad_tau']) #Have already sliced the array
				# Indices of Prad_tau correspond to ne_tau values

			self.Nzk = result[-1,:]

		self.additional_out = additional_out #Replace the additional_out with the end-values
			

		return results

	def plotScanTempCR_Dens(self, results, Te_values, Ne_const, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none"):
		
		fig, ax1 = plt.subplots()

		for k in range(self.Z+1):
			if k == 0:
				ax1.semilogx(Te_values, results[:,k], label="{}".format("g.s."))
			else:
				ax1.semilogx(Te_values, results[:,k], label="{}+".format(k))

		plt.semilogx(Te_values, np.sum(results[:,:],1), label="Total")

		total_density = np.sum(results[-1,:],0)
		ax1.set_ylim([1e-3*total_density, total_density])
		ax1.set_xlabel(r'Plasma temperature (eV)')
		ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
		# plt.title(r'Ionisation stage at C.R. as $f(T_e)$')
		ax1.tick_params('y', colors = 'b')
		plt.legend()

		ax1.set_xlim(min(Te_values), max(Te_values))

		ax1.grid(which=grid, axis='both')

		if plot_power:
			ax2 = ax1.twinx()
			scaled_power = np.array(self.additional_out['Prad'])*1e-3
			ax2.semilogx(Te_values, scaled_power,'k--',label=r'$P_{rad}$')
			ax2.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
			ax2.tick_params('y', colors='k')
			ax2.legend(loc=0)
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)
		if plot_power:
			ax2.set_yscale(y_axis_scale)

		plt.show()

	def plotScanTempCR_Prad(self, results, Te_values, Ne_const, x_axis_scale = "log", y_axis_scale = "log", grid = "none"):
		
		fig, ax1 = plt.subplots()

		total_density = np.sum(results[-1,:],0)
		scaled_power = np.array(self.additional_out['Prad'])

		compare_to_Post_PSI = True
		if compare_to_Post_PSI:
			convertCM = 1e-6 #convert from cm^3 to m^3

			ax1.semilogx(Te_values, scaled_power/(total_density*Ne_const*convertCM),'k--',label=r'$L$')
			
			ax1.set_ylabel(r'$L$ (W $cm^3$)')

			ax1.set_xlim(1, 1e3)
			ax1.set_ylim(1e-29, 1e-24)

		else:
			ax1.semilogx(Te_values, scaled_power/(total_density*Ne_const),'k--',label=r'$L$')

			ax1.set_ylabel(r'$L$ (W $m^3$)')
			
			ax1.set_xlim(min(Te_values), max(Te_values))

		ax1.set_xlabel(r'Plasma temperature (eV)')
		# plt.legend()

		ax1.grid(which=grid, axis='both')
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)

		plt.show()

	def plotScanTempCR_Prad_tau(self, results, Te_values, Ne_const, ne_tau_values=[], x_axis_scale = "log", y_axis_scale = "log", grid = "none"):
		
		fig, ax1 = plt.subplots()

		total_density = np.sum(results[-1,:],0)
		
		Prad_tau = np.array(self.additional_out['Prad_tau'])

		for ne_tau_index in range(len(ne_tau_values)):
			ne_tau = ne_tau_values[ne_tau_index]

			ax1.semilogx(Te_values, Prad_tau[:,ne_tau_index]/(total_density*Ne_const),label="{:.1e}".format(ne_tau))

			ax1.set_ylabel(r'$P_{rad}(\tau)/(n_e n_z)$ ($W m^3$)')
			
			ax1.set_xlim(min(Te_values), max(Te_values))

		ax1.set_xlabel(r'Plasma temperature (eV)')
		plt.legend()

		ax1.grid(which=grid, axis='both')
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)

		plt.show()

	def scanDensityCREquilibrium(self, Ne_values, Te_const, t_values = np.logspace(-6, 2, 200)):
		self.reset_additional_out()

		results = np.zeros((self.Z+1, len(Ne_values)))

		for Ne_iterator in range(len(Ne_values)):
			Ne = Ne_values[Ne_iterator]

			print("Evaluating test {} of {}".format(Ne_iterator, len(Ne_values)))

			result = self.solveCollisionalRadiativeEquilibrium(Te_const, Ne, t_values)

			results[:,Ne_iterator] = result[-1,:]

			self.Nzk = result[-1,:]

		return results

if __name__ == "__main__":

	import matplotlib.pyplot as plt

	impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

	solver = AtomicSolver(impurity_symbol)

	solver.Te = 50
	solver.Nzk = np.array([1e17, 0, 0, 0, 0, 0, 0])

	t_values = np.logspace(-6, 2, 200)
	Te_values = np.logspace(-0.69, 3.99, 100) #Span the entire array for which there is ADAS data
	Te_const = 50
	Ne_values = np.logspace(13.7, 21.3, 100)
	Ne_const = 1e19
	Ne_tau_values = [1e17, 1e16, 1e15, 1e14] #Values to return Prad(tau) for

	plot_solver_evolution = False
	reevaluate_scan_temp = False
	plot_scan_temp_dens = False
	plot_scan_temp_prad = False
	plot_scan_temp_prad_tau = True

	if plot_solver_evolution:
		solver_evolution = solver.solveCollisionalRadiativeEquilibrium(Te_const, Ne_const, t_values)
		solver.plotResultFromDensityEvolution(solver_evolution, t_values, plot_power = False, grid="major")
	
	if plot_scan_temp_dens or plot_scan_temp_prad or plot_scan_temp_prad_tau or reevaluate_scan_temp:
		if reevaluate_scan_temp:

			scan_temp = solver.scanTempCREquilibrium(Te_values, Ne_const, t_values=t_values, ne_tau_values=Ne_tau_values)

			with open('python_results/results_Te(scan{})_Ne({})_res({}).pickle'.format(len(Te_values),Ne_const,len(t_values)), 'wb') as handle:
				pickle.dump(scan_temp, handle, protocol=pickle.HIGHEST_PROTOCOL)
			with open('python_results/additional_out_Te(scan{})_Ne({})_res({}).pickle'.format(len(Te_values),Ne_const,len(t_values)), 'wb') as handle:
				pickle.dump(solver.additional_out, handle, protocol=pickle.HIGHEST_PROTOCOL)
		else:
			with open('python_results/results_Te(scan{})_Ne({})_res({}).pickle'.format(len(Te_values),Ne_const,len(t_values)), 'rb') as handle:
				scan_temp = pickle.load(handle)
			with open('python_results/additional_out_Te(scan{})_Ne({})_res({}).pickle'.format(len(Te_values),Ne_const,len(t_values)), 'rb') as handle:
				solver.additional_out = pickle.load(handle)

		if plot_scan_temp_dens:
			solver.plotScanTempCR_Dens(scan_temp, Te_values, Ne_const, grid="major")
		if plot_scan_temp_prad:
			solver.plotScanTempCR_Prad(scan_temp, Te_values, Ne_const, grid="major")
		if plot_scan_temp_prad_tau:	
			solver.plotScanTempCR_Prad_tau(scan_temp, Te_values, Ne_const, ne_tau_values=Ne_tau_values, grid="major")


	

























