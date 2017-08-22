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

# print(Z)

# From collisional radiative equilibrium, start values for the Carbon impurity densities
Nzk_init = np.array([1.747803e-01, 1.366167e+05, 8.865589e+09, 6.294431e+13, 9.049412e+16, 9.440710e+15, 3.206463e+13])
# Nzk_init = np.array([100, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2])
Vzk = np.zeros((Z+1,))

from scipy.integrate import odeint
from scipy.optimize import fsolve
# Two different methods for finding stable values

def evolve_density_TD(Nzk, t, Te, Ne, Vi, Nn, Vn, Vzk):
	# Time-dependant, for odeint
	# for k in range(len(Nzk)):
	# 	if(Nzk[k] < 0):
	# 		Nzk[k] = 0

	derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

	dNzk = derivative_struct["dNzk"]
	# for k in range(len(Nzk)):
	# 	if(Nzk[k] < 0):
	# 		dNzk[k] = 0
	
	return derivative_struct["dNzk"]


def evolve_density_TI(Nzk, Te, Ne, Vi, Nn, Vn, Vzk):
	# Time-independent, for optimise
	derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

	

	# print(derivative_struct["dNzk"]) #Can see that this is still changing reasonably quickly ~50 at final step
	# print(sum(dNzk))

	return derivative_struct["dNzk"]

if __name__ == "__main__":

	import matplotlib.pyplot as plt


	# Te_tests = np.logspace(-0.6, 3.8, 100)
	Te = 1.6990
	# Te = 50
	# Nz_results = np.zeros((len(Te_tests),Z+1))

	# for Te_it in range(len(Te_tests)):
	# Te = Te_tests[Te_it]

	t = np.logspace(-6, 0, 2000)
	result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	for k in range(Z+1):
		plt.semilogx(t, result[:,k], label="{}".format(k))
		# plt.semilogx(t, result[:,k]-result[-1,k], label="{}".format(k))
		print("Nz^{} = {}".format(k,result[-1,k]))
	plt.semilogx(t, np.sum(result[:,:],1), label="Total")

	plt.legend()
	plt.show()




























	# # t = np.logspace(0.0, 10, 30)
	# # result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# # for k in range(Z+1):
	# # 	plt.loglog(t, result[:,k]-result[-1,k], label="{}".format(k))	
	# # 	print(result[-1,k])
	# # plt.legend()
	# # plt.show()	

	# # eq = fsolve(evolve_density_TI,Nzk_init, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# # for k in range(Z+1):
	# # 	print("{:>5}^{}: {:.2e}".format("Ni",k,result[-1,k]))
	# # 	print("{:>5}^{}: {:.2e}".format("Ni",k,eq[k]))


	# # quit()

	# if True:
	# 	Te_tests = np.logspace(-0.6, 3.8, 500);
		
	# 	Nz_results = np.zeros((len(Te_tests),Z+1))

	# 	# Nz_results = np.zeroslike()

	# 	for Te_index in range(len(Te_tests)):
	# 		Te = Te_tests[Te_index]



	# 		t = np.logspace(0.0, 10, 30)
	# 		result = odeint(evolve_density_TD, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# 		# for k in range(Z+1):
	# 			# plt.loglog(t, result[:,k]-result[-1,k], label="{}".format(k))
	# 			# print(result[-1,k])
	# 		# plt.legend()
	# 		# plt.show()


	# 		# print("{:.2e}".format(Te))
	# 		eq = fsolve(evolve_density_TI,Nzk_init, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# 		# Nz_results[Te_index,:] = eq/sum(eq)
	# 		# Nz_results[Te_index,:] = eq
	# 		Nz_results[Te_index,:] = result[-1,:]

	# 	# import matplotlib.pyplot as plt
	# 	for k in range(Z+1):
	# 		plt.semilogx(Te_tests, Nz_results[:,k], label="{}".format(k))

	# 	# plt.ylim([0, 1])
	# 	# plt.xlim([1e-1, 1e2])
	# 	plt.legend()
	# 	plt.show()

	# if False:
	# 	# Te_tests = np.logspace(-0.6, 3.8, 500);
	# 	Ne_tests = np.logspace(14,21,500)
		
	# 	Nz_results = np.zeros((len(Ne_tests),Z+1))

	# 	# Nz_results = np.zeroslike()

	# 	for Ne_index in range(len(Ne_tests)):
	# 		Ne = Ne_tests[Ne_index]
	# 		print("{:.2e}".format(Ne))
	# 		eq = fsolve(evolve_density_TI,Nzk_init, args=(Te, Ne, Vi, Nn, Vn, Vzk))
	# 		Nz_results[Ne_index,:] = eq/sum(eq)

	# 	import matplotlib.pyplot as plt
	# 	for k in range(Z+1):
	# 		plt.loglog(Ne_tests, Nz_results[:,k], label="{}".format(k))

	# 	plt.ylim([1e-3, 1])
	# 	# plt.xlim([1e-1, 1e2])
	# 	plt.legend()
	# 	plt.show()
























