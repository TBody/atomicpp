# Program name: atomicpp/Prad.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 11 August 2017
#
# Use the atomic++ module to evaluate the rate-coefficients from OpenADAS

import numpy as np
from atomicpp import atomicpy
from scipy.integrate import odeint


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
		self.Nzk = np.array([1e17, 0, 0, 0, 0, 0, 0]) #m^-3
		self.Vzk = np.zeros((self.Z+1,)) #m/s

		self.additional_out = {'Prad':[], 'Pcool':[], 'dNzk':[], 'F_zk':[], 'dNe':[], 'F_i':[], 'dNn':[], 'F_n':[]}
	
	@staticmethod
	def evolveDensity(Nzk, t, self, Te, Ne, additional_out_keys = []):
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
		for key in additional_out_keys:
			self.additional_out[key].append(derivative_struct[key])
		# self.additional_out['Prad'].append(derivative_struct['Prad'])
		# self.additional_out['dNzk'].append(derivative_struct['dNzk'])

		return dNzk

	def reset_additional_out(self):
		self.additional_out = {'Prad':[], 'Pcool':[], 'dNzk':[], 'F_zk':[], 'dNe':[], 'F_i':[], 'dNn':[], 'F_n':[]}

	def solveCollisionalRadiativeEquilibrium(self, Te, Ne, t_values = np.logspace(-6, 2, 200), additional_out_keys = []):
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		Vzk = self.Vzk

		print("Te = {}eV, Ne = {}/m3".format(Te, Ne))
		
		(result, output_dictionary) = odeint(self.evolveDensity, self.Nzk, t_values, args=(self,Te, Ne, additional_out_keys), printmessg=False, full_output=True)

		feval_at_step = output_dictionary['nfe'] #function evaluations at the time-step
		time_at_step = output_dictionary['tcur'] #time at the time-step

		time_indices = np.searchsorted(t_values, time_at_step, side='left') #find how the time-steps are distributed. Usually close but not 1 to 1 with t_values

		for key, value in self.additional_out.items():
			if value: #if list isn't empty - i.e. additional output has been recorded for this key
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
			
				# plt.semilogx(t_values, output_values)
				# plt.show()


		# for k in range(self.Z+1):
		# 	print("{:20} Nz^{} = {}".format("Equilibrium found for",k,result[-1,k]))

		solver.Nzk = result[-1,:]

		return result

	def scanTempCREquilibrium(self, Te_values, Ne_const, t_values = np.logspace(-6, 2, 200)):
		self.reset_additional_out()

		results = np.zeros((self.Z+1, len(Te_values)))

		for Te_iterator in range(len(Te_values)):
			Te = Te_values[Te_iterator]

			print("Evaluating test {} of {}".format(Te_iterator, len(Te_values)))

			result = self.solveCollisionalRadiativeEquilibrium(Te, Ne_const, t_values)

			results[:,Te_iterator] = result[-1,:]

			self.Nzk = result[-1,:]

		return results

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

	def plotResultFromDensityEvolution(self, result, t_values, plot_power = False):
		fig, ax1 = plt.subplots()
		for k in range(solver.Z+1):
			ax1.semilogx(t_values, result[:,k], label="{}".format(k))

		ax1.semilogx(t_values, np.sum(result[:,:],1), label="Total")

		ax1.set_xlabel(r'Time (s)')
		ax1.set_ylabel(r'Density of stage ($m^{-3}$)')
		plt.title('Time evolution of ionisation stages')
		ax1.tick_params('y', colors = 'b')
		ax1.legend()

		if plot_power:
			ax2 = ax1.twinx()
			scaled_power = np.array(self.additional_out['Prad'])*1e-3
			ax2.semilogx(t_values, scaled_power,'k--',label=r'$P_{rad}$')
			ax2.set_ylabel(r'$P_{rad}$ (KW $m^{-3}$)')
			ax2.tick_params('y', colors='k')
			ax2.legend(loc=0)

		plt.show()

if __name__ == "__main__":

	import matplotlib.pyplot as plt

	impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

	solver = AtomicSolver(impurity_symbol)

	# solver.Te = 50
	# solver.Nzk = np.array([1e17, 0, 0, 0, 0, 0, 0])

	t = np.logspace(-6, 2, 200)

	# result = solver.solveCollisionalRadiativeEquilibrium(1, 1e19, t, ['Prad'])
	# solver.plotResultFromDensityEvolution(result, t)

	Te_values = np.logspace(-0.69, 3.99, 100) #Span the entire array for which there is ADAS data
	Ne_values = np.logspace(13.7, 21.3, 100)
	results = solver.scanTempCREquilibrium(Te_values, 1e19, t)

	for k in range(solver.Z+1):
		plt.loglog(Te_values, results[k,:], label="{}".format(k))

	plt.loglog(Te_values, np.sum(results[:,:],0), label="Total")
	total_density = np.sum(results[:,-1],0)
	plt.ylim([1e-3*total_density, total_density])
	plt.xlabel(r'Plasma temperature (eV)')
	plt.ylabel(r'Density of stage ($m^{-3}$)')
	plt.title(r'Ionisation stage at C.R. as $f(T_e)$')
	plt.legend()
	plt.show()

























