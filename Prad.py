# Program name: atomicpp/Prad.py
# Author: Thomas Body
# Author email: tajb500@york.ac.uk
# Date of creation: 11 August 2017
#
# Use the atomic++ module to evaluate the rate-coefficients from OpenADAS

import numpy as np
from atomicpp import atomicpy

impurity_symbol = b'c' #need to include b (bytes) before the string for it to be sent as a std::string to C++

impurity = atomicpy.PyImpuritySpecies(impurity_symbol)

impurity_derivatives = atomicpy.PyRateEquations(impurity);
impurity_derivatives.setThresholdDensity(0.0);
impurity_derivatives.setDominantIonMass(1.0);

Te = 50; #eV
Ne = 1e19; #m^-3
Vi = 0; #m/s
Nn = 0; #m^-3
Vn = 0; #m/s

Z=impurity.get_atomic_number()

# From collisional radiative equilibrium, start values for the Carbon impurity densities
Nzk_init = np.array([1.747803e-01, 1.366167e+05, 8.865589e+09, 6.294431e+13, 9.049412e+16, 9.440710e+15, 3.206463e+13])
# Don't treat momentum yet
Vzk = np.zeros((Z+1,))

from scipy.integrate import odeint
from scipy.optimize import fsolve
# Two different methods for finding stable values

def evolve_density_TD(Nzk, t, Te, Ne, Vi, Nn, Vn, Vzk):
	# Time-dependant, for odeint

	# Prevent negative densities
	# (these are possible if the time-step is reasonably large)
	# for k in range(len(Nzk)):
	# 	if(Nzk[k] < 0):
	# 		Nzk[k] = 0

	derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

	dNzk = derivative_struct["dNzk"]

	return derivative_struct["dNzk"]

if __name__ == "__main__":

	import matplotlib.pyplot as plt

	# Te = 1

	# t = np.logspace(-6, 2, 200)
	# result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk), printmessg=True)

	# for k in range(Z+1):
	# 	plt.semilogx(t, result[:,k], label="{}".format(k))
	# 	# print("Nz^{} = {}".format(k,result[-1,k]))
	# plt.semilogx(t, np.sum(result[:,:],1), label="Total")
	# plt.xlabel(r'Time (s)')
	# plt.ylabel(r'Density of stage ($m^{-3}$)')
	# plt.title('Time evolution of ionisation stages')
	# plt.legend()
	# plt.show()




	Te_values = np.logspace(-0.69, 3.99, 100) #Span the entire array for which there is ADAS data
	Ne_values = np.logspace(13.7, 21.3, 100)

	t = np.logspace(-6, 2, 2000)

	Ne = 1e19; #m^-3

	Nzk_init = np.ones((Z+1,)) * 1e17/Z
	# Nzk_init[0] = 1e17

	Ionisation_stage_distibution = np.zeros((Z+1,len(Te_values)))

	for Te_iterator in range(len(Te_values)):
		Te = Te_values[Te_iterator]

		result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))

		# for k in range(Z+1):
			# print("Nz^{} = {}".format(k,result[-1,k]))

		Ionisation_stage_distibution[:,Te_iterator] = result[-1,:]
		# Nzk_init = result[-1,:]

	for k in range(Z+1):
		plt.loglog(Te_values, Ionisation_stage_distibution[k,:], label="{}".format(k))
	
	plt.loglog(Te_values, np.sum(Ionisation_stage_distibution[:,:],0), label="Total")
	plt.ylim([1e-3*np.sum(Ionisation_stage_distibution[:,-1],0), 1*np.sum(Ionisation_stage_distibution[:,-1],0)])
	plt.xlabel(r'Plasma temperature (eV)')
	plt.ylabel(r'Density of stage ($m^{-3}$)')
	plt.title(r'Ionisation stage at C.R. as $f(T_e)$')
	plt.legend()
	plt.show()

























