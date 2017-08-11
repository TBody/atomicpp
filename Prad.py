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
impurity_derivatives.setThresholdDensity(1e9);
impurity_derivatives.setDominantIonMass(1.0);

Te = 50; #eV
Ne = 1e19; #m^-3
Vi = 0; #m/s
Nn = 0; #m^-3
Vn = 0; #m/s

Z=impurity.get_atomic_number()

print(Z)

# From collisional radiative equilibrium, start values for the Carbon impurity densities
Nzk_init = np.array([1.747803e-01, 1.366167e+05, 8.865589e+09, 6.294431e+13, 9.049412e+16, 9.440710e+15, 3.206463e+13])
Vzk = np.zeros((Z+1,))

from scipy.integrate import odeint
from scipy.optimize import fsolve
# Two different methods for finding stable values

def evolve_density_TD(Nzk, t, Te, Ne, Vi, Nn, Vn, Vzk):
	# Time-dependant, for odeint
	derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

	return derivative_struct["dNzk"]


def evolve_density_TI(Nzk, Te, Ne, Vi, Nn, Vn, Vzk):
	# Time-independent, for optimise
	derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

	# print(derivative_struct["dNzk"]) #Can see that this is still changing reasonably quickly ~50 at final step
	return derivative_struct["dNzk"]

if __name__ == "__main__":

	# import matplotlib.pyplot as plt

	# t = np.logspace(0.0, 10, 30)
	# result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# for k in range(Z+1):
	# 	plt.loglog(t, result[:,k], label="{}".format(k))	
	# plt.legend()
	# plt.show()

	# eq = fsolve(evolve_density_TI,Nzk_init, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# for k in range(Z+1):
	# 	print("{:>5}^{}: {:.2e}".format("Ni",k,result[-1,k]))
	# 	print("{:>5}^{}: {:.2e}".format("Ni",k,eq[k]))

	Te_tests = np.logspace(-0.6, 3.8, 500);
	
	Nz_results = np.zeros((len(Te_tests),Z+1))

	# Nz_results = np.zeroslike()

	for Te_index in range(len(Te_tests)):
		Te = Te_tests[Te_index]
		print("{:.2e}".format(Te))
		eq = fsolve(evolve_density_TI,Nzk_init, args=(Te, Ne, Vi, Nn, Vn, Vzk))
		Nz_results[Te_index,:] = eq/sum(eq)

	import matplotlib.pyplot as plt
	for k in range(Z+1):
		plt.semilogx(Te_tests, Nz_results[:,k], label="{}".format(k))	
	plt.legend()
	plt.show()
























