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

	@staticmethod
	def evolveDensity(Nzk, t, self):
		Te  = self.Te
		Ne  = self.Ne
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		# Nzk = self.Nzk
		Vzk = self.Vzk

		# Prevent negative densities
		# (these are possible if the time-step is reasonably large)
		# for k in range(len(Nzk)):
		# 	if(Nzk[k] < 0):
		# 		Nzk[k] = 0

		derivative_struct = self.impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

		dNzk = derivative_struct["dNzk"]

		return dNzk

	def solveCollisionalRadiativeEquilibrium(self, Te, Ne, t_values = np.logspace(-6, 2, 200)):
		Vi  = self.Vi
		Nn  = self.Nn
		Vn  = self.Vn
		Vzk = self.Vzk

		print("Te = {}eV, Ne = {}/m3".format(Te, Ne))
		
		result = odeint(self.evolveDensity, self.Nzk, t_values, args=(self,), printmessg=False)

		# for k in range(self.Z+1):
		# 	print("{:20} Nz^{} = {}".format("Equilibrium found for",k,result[-1,k]))

		return result

	def scanCollisionalRadiativeEquilibrium(self, Te_values, Ne_values, t_values = np.logspace(-6, 2, 200)):

		result_grid = np.zeros((self.Z+1, len(Te_values), len(Ne_values)))

		for Te_iterator in range(len(Te_values)):
			Te = Te_values[Te_iterator]
			for Ne_iterator in range(len(Ne_values)):
				Ne = Ne_values[Ne_iterator]

				print("Evaluating test {} of {}".format(Te_iterator*len(Te_values)+Ne_iterator, len(Te_values)*len(Ne_values)))

				result = self.solveCollisionalRadiativeEquilibrium(Te, Ne, t_values)

				result_grid[:,Te_iterator,Ne_iterator] = result[-1,:]

		return result_grid

				



if __name__ == "__main__":

	import matplotlib.pyplot as plt

	impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

	solver = AtomicSolver(impurity_symbol)

	# solver.Te = 50
	# solver.Nzk = np.array([1e17, 0, 0, 0, 0, 0, 0])

	t = np.logspace(-6, 2, 200)
	result = solver.solveCollisionalRadiativeEquilibrium(50, 1e19, t)

	for k in range(solver.Z+1):
		plt.semilogx(t, result[:,k], label="{}".format(k))
		
	plt.semilogx(t, np.sum(result[:,:],1), label="Total")
	plt.xlabel(r'Time (s)')
	plt.ylabel(r'Density of stage ($m^{-3}$)')
	plt.title('Time evolution of ionisation stages')
	plt.legend()
	plt.show()




	# Te_values = np.logspace(-0.69, 3.99, 100) #Span the entire array for which there is ADAS data
	# Ne_values = np.logspace(13.7, 21.3, 100)

	# t = np.logspace(-6, 2, 2000)

	# Ne = 1e19; #m^-3

	# Nzk_init = np.ones((Z+1,)) * 1e17/Z
	# # Nzk_init[0] = 1e17

	# Ionisation_stage_distibution = np.zeros((Z+1,len(Te_values)))

	# for Te_iterator in range(len(Te_values)):
	# 	Te = Te_values[Te_iterator]

	# 	result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))

	# 	# for k in range(Z+1):
	# 		# print("Nz^{} = {}".format(k,result[-1,k]))

	# 	Ionisation_stage_distibution[:,Te_iterator] = result[-1,:]
	# 	# Nzk_init = result[-1,:]

	# for k in range(Z+1):
	# 	plt.loglog(Te_values, Ionisation_stage_distibution[k,:], label="{}".format(k))

	# plt.loglog(Te_values, np.sum(Ionisation_stage_distibution[:,:],0), label="Total")
	# plt.ylim([1e-3*np.sum(Ionisation_stage_distibution[:,-1],0), 1*np.sum(Ionisation_stage_distibution[:,-1],0)])
	# plt.xlabel(r'Plasma temperature (eV)')
	# plt.ylabel(r'Density of stage ($m^{-3}$)')
	# plt.title(r'Ionisation stage at C.R. as $f(T_e)$')
	# plt.legend()
	# plt.show()

























