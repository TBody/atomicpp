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

Z=6

# From collisional radiative equilibrium, start values for the Carbon impurity densities
Nzk_init = np.array([1.747803e-01, 1.366167e+05, 8.865589e+09, 6.294431e+13, 9.049412e+16, 9.440710e+15, 3.206463e+13])
Vzk = np.zeros((7,))

from scipy.integrate import odeint

def evolve_density(Nzk, t, Te, Ne, Vi, Nn, Vn, Vzk):
	derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);

	return derivative_struct["dNzk"]

if __name__ == "__main__":
	t = np.logspace(0.0, 1, 20)
	result = odeint(evolve_density, Nzk_init, t, args=(Te, Ne, Vi, Nn, Vn, Vzk))

	import matplotlib.pyplot as plt

	for i in range(Z+1):
		plt.loglog(t, result[:,i], label="{}".format(i))
		
	plt.legend()
	plt.show()