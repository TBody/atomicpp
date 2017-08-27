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
		self.reevaluate_scan         = False
		self.error_estimation        = False
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
				# self.Nzk = result[-1,:] # Can use previous result to try speed up evaluation if not calculating Prad(tau)
				pass # However, this is found to result in odd numerical behaviour
				


		return result

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
	def plotResultFromDensityEvolution(self, result, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none", show=False):
		fig, ax1 = plt.subplots()
		for k in range(self.Z+1):
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
			ax2.set_ylim(min(scaled_power), max(scaled_power))
			ax2.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
			ax2.tick_params('y', colors='k')
			ax2.legend(loc=0)
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)
		if plot_power:
			ax2.set_yscale(y_axis_scale)

		if show:
			plt.show()
		return fig

	def plotScanTempCR_Dens(self, results, plot_power = False, x_axis_scale = "log", y_axis_scale = "linear", grid = "none", show=False):
		
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

		if show:
			plt.show()
		return fig

	def plotScanTempCR_Prad(self, results, x_axis_scale = "log", y_axis_scale = "log", grid = "none", show=False):
		
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

		if show:
			plt.show()
		return fig

	def plotScanTempCR_Prad_tau(self, refuelling_results, x_axis_scale = "log", y_axis_scale = "log", grid = "none", show=False):
		from scipy.interpolate import interp1d
		

		# Extract the radiation curves from Post PSI paper. Plot the carbon result against these
		x_eval_points = np.logspace(0.1,3,15)
		x_eval_points_hires = np.logspace(0.1,3,100)
		y_data_points = []
		y_data_points_hires = []
		for dataset in range(0,4):
			Prad = np.loadtxt('python_results/Prad{}.csv'.format(dataset+1),delimiter=', ')
			x = np.array(Prad[:,0])
			y = np.array(Prad[:,1])/1e12
			Prad_interp = interp1d(x,y,bounds_error=False)
			Prad_interp_hires = interp1d(x,y,bounds_error=False)
			# plt.loglog(x,y,label=(dataset+1))
			y_data_points.append(Prad_interp(x_eval_points))
			y_data_points_hires.append(Prad_interp(x_eval_points_hires))
		y_data_points = np.array(y_data_points)
		y_mean = np.nanmean(y_data_points,0)
		y_mean_hires = np.nanmean(y_data_points_hires,0)
		y_std = np.nanstd(y_data_points,0)
		ylower = np.maximum(1e-50, y_mean - y_std)
		y_std_lower = y_mean - ylower

		fig, ax1 = plt.subplots()

		Prad = np.array(self.additional_out['Prad'])

		for ne_tau_index in range(len(self.Ne_tau_values)-1,-1,-1):
			ne_tau = self.Ne_tau_values[ne_tau_index]

			total_density = np.sum(refuelling_results[-1,ne_tau_index,:],0)

			if ne_tau > 1e18:
				ax1.semilogx(self.Te_values, Prad[:,ne_tau_index]/(total_density*self.Ne_const),'r-.',label="CR")
			else:
				ax1.semilogx(self.Te_values, Prad[:,ne_tau_index]/(total_density*self.Ne_const),label="{:.1e}".format(ne_tau))

		ax1.set_ylabel(r'$P_{rad}(\tau)/(N_e N_z)$ ($W m^3$)')
		ax1.set_xlim(min(self.Te_values), max(self.Te_values))
		ax1.set_ylim(1e-37, 1e-30)

		ax1.semilogx(x_eval_points_hires, y_mean_hires, 'k--', label='CR expected')
		ax1.errorbar(x_eval_points, y_mean, yerr=[y_std_lower, 2*y_std], xerr=0, fmt='k.', ecolor='k', capthick=1, capsize=3)

		# Plot the carbon cooling curve from I.Hutchinson
		Te = self.Te_values
		P = 2e-31 * (Te/10.)**3 / ( 1.0 + (Te/10.)**4.5 )
		ax1.semilogx(self.Te_values, P, 'b-.',label="Hutchinson")


		ax1.set_xlabel(r'Plasma temperature (eV)')
		plt.legend(loc=0)

		ax1.grid(which=grid, axis='both')
		
		ax1.set_xscale(x_axis_scale)
		ax1.set_yscale(y_axis_scale)

		if show:
			plt.show()
		return fig

	def error_estimations(self, show=False, plot_time_test = False, plot_start_test = False, find_stddev = False, plot_stddev_Ne = False, plot_stddev_Te = False, plot_together = False):
		self.t_values = np.logspace(-10, 2, 10000)
		if self.reevaluate_scan:
			error_analysis = self.timeIntegrate(self.Te_const, self.Ne_const, 0)
			comparison_values = error_analysis[-1,:]
			comparison_values = np.append(comparison_values,self.additional_out['Prad'][-1])

			with open('python_results/error_analysis({},{})(res+{})-comparison_values.pickle'.format(self.Te_const,self.Ne_const,len(self.t_values)), 'wb') as handle:
				pickle.dump(comparison_values, handle, protocol=pickle.HIGHEST_PROTOCOL)
		else:
			with open('python_results/error_analysis({},{})(res+{})-comparison_values.pickle'.format(self.Te_const,self.Ne_const,len(self.t_values)), 'rb') as handle:
				comparison_values = pickle.load(handle)
		
		if plot_time_test:
			time_test_values = np.round(np.logspace(1, 3, 20))
			time_test_results = [];
			for time_it in range(len(time_test_values)):
				self.reset_additional_out()

				time_test = time_test_values[time_it]
				print("Evaluating for time-resolution = {}pts".format(time_test))
				self.t_values = np.logspace(-6, 2, time_test)
				error_analysis = self.timeIntegrate(self.Te_const, self.Ne_const, 0)
				# plot_self_evolution = self.plotResultFromDensityEvolution(error_analysis, plot_power = True, grid="major", show=True)
				test_values = error_analysis[-1,:]
				test_values = np.append(test_values,self.additional_out['Prad'][-1])
				time_test_results.append(comparison_values-test_values)
			plt.plot(time_test_values, time_test_results)

		if plot_start_test:
			start_time_log = -6
			shift_test_values = np.linspace(-3,10, num=100)
			original_test_length = len(shift_test_values)
			if self.reevaluate_scan:	
				shift_test_values = shift_test_values.tolist()
				shift_test_results = [];
				failed_test_values = [];
				for shift_it in range(len(shift_test_values)):
					self.reset_additional_out()

					shift_test = shift_test_values[shift_it]
					print("Evaluating for shift = {}".format(shift_test))
					try:
						self.t_values = np.logspace(start_time_log+shift_test, 5, 200)
						error_analysis = self.timeIntegrate(self.Te_const, self.Ne_const, 0)
						# print("COMPARISON: ", comparison_values)
						# plot_self_evolution = self.plotResultFromDensityEvolution(error_analysis, plot_power = True, grid="major", show=True)
						test_values = error_analysis[-1,:]
						test_values = np.append(test_values,self.additional_out['Prad'][-1])
						shift_test_results.append((comparison_values - test_values)/comparison_values)
					except:
						print("Evaluation failed for shift = {}".format(shift_test))
						failed_test_values.append(shift_test)

				for shift_test in failed_test_values:
					shift_test_values.remove(shift_test)

				shift_test_results = np.absolute(shift_test_results)
				shift_test_values = np.array(shift_test_values)+start_time_log
				start_times = np.power(10*np.ones_like(shift_test_values), shift_test_values)

				shift_test_data = {}
				shift_test_data['results'] = shift_test_results
				shift_test_data['times'] = start_times

				with open('python_results/error_analysis({},{})(res+{})-shift_test_data.pickle'.format(self.Te_const,self.Ne_const,original_test_length), 'wb') as handle:
					pickle.dump(shift_test_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
			else:
				with open('python_results/error_analysis({},{})(res+{})-shift_test_data.pickle'.format(self.Te_const,self.Ne_const,original_test_length), 'rb') as handle:
					shift_test_data = pickle.load(handle)

			for k in range(self.Z+2):
				if k == 0:
					plt.loglog(shift_test_data['times'], shift_test_data['results'][:,k], label="{}".format("g.s."))
				elif k == self.Z+1:
					plt.loglog(shift_test_data['times'], shift_test_data['results'][:,k], label="{}".format(r"$P_{rad}$"))
				else:
					plt.loglog(shift_test_data['times'], shift_test_data['results'][:,k], label="{}+".format(k))
			
			plt.xlabel('Start time for evaluation (s)')
			plt.ylabel('Deviation from expected answer (relative error)')

			plt.legend()

		if find_stddev:
			if self.reevaluate_scan:
				import random
				random.seed(1)
				random_results = []
				for iterator in range(100):
					self.reset_additional_out()
					Nzk = np.zeros((self.Z+1,))
					for k in range(self.Z+1):
						Nzk[k] = random.random()*(10**random.uniform(1, 17))
					self.Nzk = 1e17*Nzk/sum(Nzk)

					self.t_values = np.logspace(-6, 2, 200)
					random_init = self.timeIntegrate(self.Te_const, self.Ne_const, 0)
					random_values = random_init[-1,:]
					random_values = np.append(random_values,self.additional_out['Prad'][-1])
					random_results.append(random_values)
				
				random_results = np.array(random_results)
				with open('python_results/error_analysis({},{})(res+{})-random_results.pickle'.format(self.Te_const,self.Ne_const,len(self.t_values)), 'wb') as handle:
					pickle.dump(random_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
			else:
				with open('python_results/error_analysis({},{})(res+{})-random_results.pickle'.format(self.Te_const,self.Ne_const,len(self.t_values)), 'rb') as handle:
					random_results = pickle.load(handle)


			for k in range(self.Z+2):
				mean = np.mean(random_results[:,k])
				stdev = np.std(random_results[:,k])
				stdev_norm = stdev/mean
				diff = (mean - comparison_values[k])/comparison_values[k]
				if k == 0:
					print("{:5} -> mean = {:.2e}, stdev = {:.2e}, stdev_norm = {:.2e}, mean_diff = {:.2e}".format("g.s.", mean, stdev, stdev_norm, diff))
				elif k == self.Z+1:
					print("{:5} -> mean = {:.2e}, stdev = {:.2e}, stdev_norm = {:.2e}, mean_diff = {:.2e}".format("P_rad", mean, stdev, stdev_norm, diff))
				else:
					print("{:5} -> mean = {:.2e}, stdev = {:.2e}, stdev_norm = {:.2e}, mean_diff = {:.2e}".format(k, mean, stdev, stdev_norm, diff))

			self.Nzk = np.zeros((self.Z+1,)) #m^-3
			self.Nzk[0] = 1e17 #m^-3 - start in g.s.

		if plot_stddev_Te:
			self.t_values = np.logspace(-6, 2, 200)
			stdev_Te = np.linspace(0,self.Te_const/2,num=20)
			stdev_norm = []
			if self.reevaluate_scan:
				import random
				random.seed(1)
				for sigma in stdev_Te:
					random_results = []
					for iterator in range(50):
						try:
							self.reset_additional_out()
							Te = random.normalvariate(self.Te_const, sigma)
							Nzk = np.zeros((self.Z+1,))
							for k in range(self.Z+1):
								Nzk[k] = random.random()*(10**random.uniform(1, 17))
							self.Nzk = 1e17*Nzk/sum(Nzk)

							random_te = self.timeIntegrate(Te, self.Ne_const, 0)
							random_values = random_te[-1,:]
							random_values = np.append(random_values,self.additional_out['Prad'][-1])
							if not(np.isnan(random_values).any()):
								random_results.append(random_values)
							else:
								print("NaN for Te = {}".format(Te))
						except:
							print("Error for Te = {}".format(Te))

					mean = np.mean(random_results,0)
					stdev = np.std(random_results,0)
					stdev_norm.append(np.absolute(stdev/mean))

					
				stdev_norm = np.array(stdev_norm)
				with open('python_results/error_analysis(rand{},{})(res+{})-stdev_norm.pickle'.format(self.Te_const,self.Ne_const,len(stdev_Te)), 'wb') as handle:
					pickle.dump(stdev_norm, handle, protocol=pickle.HIGHEST_PROTOCOL)
			else:
				with open('python_results/error_analysis(rand{},{})(res+{})-stdev_norm.pickle'.format(self.Te_const,self.Ne_const,len(stdev_Te)), 'rb') as handle:
					stdev_norm = pickle.load(handle)

			for k in range(self.Z+2):
				if k == 0:
					plt.plot(stdev_Te/self.Te_const, stdev_norm[:,k], label="{}".format("g.s."))
				elif k == self.Z+1:
					plt.plot(stdev_Te/self.Te_const, stdev_norm[:,k], label="{}".format(r"$P_{rad}$"))
				else:
					plt.plot(stdev_Te/self.Te_const, stdev_norm[:,k], label="{}+".format(k))
			plt.legend()
			plt.xlabel(r'Relative error in $T_e$')
			plt.ylabel(r'Relative error in parameter ($\sigma/\mu$)')

		if plot_stddev_Ne:
			self.t_values = np.logspace(-6, 2, 200)
			stdev_Ne = np.linspace(0,self.Ne_const/2,num=20)
			stdev_norm = []
			if self.reevaluate_scan:
				import random
				random.seed(1)
				for sigma in stdev_Ne:
					random_results = []
					for iterator in range(50):
						try:
							self.reset_additional_out()
							Ne = random.normalvariate(self.Ne_const, sigma)
							Nzk = np.zeros((self.Z+1,))
							for k in range(self.Z+1):
								Nzk[k] = random.random()*(10**random.uniform(1, 17))
							self.Nzk = 1e17*Nzk/sum(Nzk)

							random_te = self.timeIntegrate(self.Te_const, Ne, 0)
							random_values = random_te[-1,:]
							random_values = np.append(random_values,self.additional_out['Prad'][-1])
							if not(np.isnan(random_values).any()):
								random_results.append(random_values)
							else:
								print("NaN for Ne = {}".format(Ne))
						except:
							print("Error for Ne = {}".format(Ne))

					mean = np.mean(random_results,0)
					stdev = np.std(random_results,0)
					stdev_norm.append(np.absolute(stdev/mean))

					
				stdev_norm = np.array(stdev_norm)
				with open('python_results/error_analysis({},rand{})(res+{})-stdev_norm.pickle'.format(self.Te_const,self.Ne_const,len(stdev_Ne)), 'wb') as handle:
					pickle.dump(stdev_norm, handle, protocol=pickle.HIGHEST_PROTOCOL)
			else:
				with open('python_results/error_analysis({},rand{})(res+{})-stdev_norm.pickle'.format(self.Te_const,self.Ne_const,len(stdev_Ne)), 'rb') as handle:
					stdev_norm = pickle.load(handle)

			for k in range(self.Z+2):
				if k == 0:
					plt.plot(stdev_Ne/self.Ne_const, stdev_norm[:,k], label="{}".format("g.s."))
				elif k == self.Z+1:
					plt.plot(stdev_Ne/self.Ne_const, stdev_norm[:,k], label="{}".format(r"$P_{rad}$"))
				else:
					plt.plot(stdev_Ne/self.Ne_const, stdev_norm[:,k], label="{}+".format(k))
					# if(k==4 or k==5):
					# 	plt.plot(stdev_Ne/self.Ne_const, stdev_norm[:,k], label="{}+".format(k))
			plt.legend()
			plt.xlabel(r'Relative error in $N_e$')
			plt.ylabel(r'Relative error in parameter ($\sigma/\mu$)')

		if plot_together:
			fig, (ax1, ax2) = plt.subplots(2, sharex = True)

			stdev_Te = np.linspace(0,self.Te_const/2,num=20)

			with open('python_results/error_analysis(rand{},{})(res+{})-stdev_norm.pickle'.format(self.Te_const,self.Ne_const,len(stdev_Te)), 'rb') as handle:
				stdev_norm_Te = pickle.load(handle)

			for k in range(self.Z+2):
				if k == 0:
					pass
					# ax1.plot(stdev_Te/self.Te_const, stdev_norm_Te[:,k], label="{}".format("g.s."))
				elif k == self.Z+1:
					ax1.plot(stdev_Te/self.Te_const, stdev_norm_Te[:,k], label="{}".format(r"$P_{rad}$"))
				else:
					# ax1.plot(stdev_Te/self.Te_const, stdev_norm_Te[:,k], label="{}+".format(k))
					if(k in [4,5]):
						ax1.plot(stdev_Te/self.Te_const, stdev_norm_Te[:,k], label="{}+".format(k))
			ax1.legend()
			ax1.set_xlabel(r'Relative error in $T_e$')
			# ax1.set_ylabel(r'Relative error in parameter ($\sigma/\mu$)')

			stdev_Ne = np.linspace(0,self.Ne_const/2,num=20)
			with open('python_results/error_analysis({},rand{})(res+{})-stdev_norm.pickle'.format(self.Te_const,self.Ne_const,len(stdev_Ne)), 'rb') as handle:
					stdev_norm_Ne = pickle.load(handle)
			for k in range(self.Z+2):
				if k == 0:
					pass
					# ax2.plot(stdev_Ne/self.Ne_const, stdev_norm_Ne[:,k], label="{}".format("g.s."))
				elif k == self.Z+1:
					ax2.plot(stdev_Ne/self.Ne_const, stdev_norm_Ne[:,k], label="{}".format(r"$P_{rad}$"))
				else:
					# ax2.plot(stdev_Ne/self.Ne_const, stdev_norm_Ne[:,k], label="{}+".format(k))
					if(k in [4,5]):
						ax2.plot(stdev_Ne/self.Ne_const, stdev_norm_Ne[:,k], label="{}+".format(k))
			ax2.legend()
			ax2.set_xlabel(r'Relative error in $N_e$')
			# ax2.set_ylabel(r'Relative error in parameter ($\sigma/\mu$)')
			fig.text(0.04, 0.5, r'Relative error in parameter ($\sigma/\mu$)', va='center', rotation='vertical')

		self.Nzk = np.zeros((self.Z+1,)) #m^-3
		self.Nzk[0] = 1e17 #m^-3 - start in g.s.
		self.t_values = np.logspace(-6, 2, 200)

		if show:
			plt.show()
		return fig


if __name__ == "__main__":

	import matplotlib.pyplot as plt

	impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

	solver = AtomicSolver(impurity_symbol)

	path_to_output = 'Figures/'

	solver.plot_solver_evolution   = False
	solver.reevaluate_scan         = False
	solver.error_estimation        = False
	solver.plot_scan_temp_dens     = False
	solver.plot_scan_temp_prad     = False
	solver.plot_scan_temp_prad_tau = True

	# solver.Ne_tau_values = [1e30, 1e17, 1e16, 1e15, 1e14] #m^-3 s, values to return Prad(tau) for
	solver.Ne_tau_values = [1e30, 1e17, 1e16, 1e15] #m^-3 s, values to return Prad(tau) for

	if solver.plot_solver_evolution:
		solver_evolution = solver.timeIntegrate(solver.Te_const, solver.Ne_const, 0)
		plot_solver_evolution = solver.plotResultFromDensityEvolution(solver_evolution, plot_power = True, grid="major", show=False, y_axis_scale="linear")
		plot_solver_evolution.savefig(path_to_output+"solver_evolution.pdf")

	if solver.error_estimation:
		error_estimation = solver.error_estimations(show=True, plot_together=True)
		error_estimation.savefig(path_to_output+"error_estimation.pdf")

	if solver.plot_scan_temp_dens or solver.plot_scan_temp_prad:
		if solver.reevaluate_scan:

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
			plot_scan_temp_dens = solver.plotScanTempCR_Dens(scan_temp, grid="major", show=False)
			plot_scan_temp_dens.savefig(path_to_output+"plot_scan_temp_dens.pdf")
		if solver.plot_scan_temp_prad:
			plot_scan_temp_prad = solver.plotScanTempCR_Prad(scan_temp, grid="major", show=False)
			plot_scan_temp_prad.savefig(path_to_output+"plot_scan_temp_prad.pdf")

	if solver.plot_scan_temp_prad_tau:
		if solver.reevaluate_scan:

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
			plot_scan_temp_prad_tau = solver.plotScanTempCR_Prad_tau(scan_temp_refuelling, grid="major", show=True)
			plot_scan_temp_prad_tau.savefig(path_to_output+"plot_scan_temp_prad_tau.pdf")


	# plt.show()



	

























