# distutils: language = c++
# distutils: sources = [ImpuritySpecies.cpp, RateCoefficient.cpp, sharedFunctions.cpp, RateEquations.cpp]

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
from libcpp.string cimport string #to pass strings to/from C++
from libcpp.map cimport map
from libcpp.memory cimport shared_ptr

# Note - only need to declare the functions that you will use. Internal functions and attributes do not 

cdef extern from "ImpuritySpecies.hpp":
	cdef cppclass ImpuritySpecies:
		ImpuritySpecies(string& impurity_symbol_supplied) except +

cdef extern from "RateEquations.hpp":
	cdef cppclass RateEquations:
		RateEquations(ImpuritySpecies& impurity, double Nthres_set = 1e9, double mi_in_amu = 1) except +
		void setThresholdDensity(double density_threshold)
		void setDominantIonMass(double mi_in_amu)
		std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > computeDerivs(double Te, double Ne, double Vi, double Nn, double Vn, std::vector<double>& Nzk, std::vector<double>& Vzk)

