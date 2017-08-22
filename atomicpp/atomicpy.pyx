# distutils: language = c++
# distutils: sources = [ImpuritySpecies.cpp, RateCoefficient.cpp, sharedFunctions.cpp, RateEquations.cpp, Spline/BilinearSpline.cpp, Spline/BicubicSpline.cpp, Spline/BivariateBSpline.cpp]

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
from libcpp.string cimport string #to pass strings to/from C++
from libcpp.map cimport map
from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

# Note - only need to declare the functions that you will use. Internal functions and attributes do not need to be declared.
# See https://dmtn-013.lsst.io/ for the use of unique_ptr in classes

cdef extern from "ImpuritySpecies.hpp" namespace "atomicpp":
	cdef cppclass ImpuritySpecies:
		ImpuritySpecies(string& impurity_symbol_supplied) except +
		string get_name()
		int get_atomic_number()

cdef extern from "RateEquations.hpp" namespace "atomicpp":
	cdef struct DerivStruct:
		double Pcool
		double Prad
		vector[double] dNzk
		vector[double] F_zk
		double dNe
		double F_i
		double dNn
		double F_n

cdef extern from "RateEquations.hpp" namespace "atomicpp":
	cdef cppclass RateEquations:
		RateEquations(ImpuritySpecies& impurity) except + #different case to Python class
		void setThresholdDensity(double density_threshold)
		void setDominantIonMass(double mi_in_amu)
		DerivStruct computeDerivs(double Te, double Ne, double Vi, double Nn, double Vn, vector[double]& Nzk, vector[double]& Vzk)

cdef class PyImpuritySpecies:
	cdef unique_ptr[ImpuritySpecies] ImpuritySpeciesPtr
	def __init__(self, string impurity_symbol_supplied):
		self.ImpuritySpeciesPtr.reset(new ImpuritySpecies(impurity_symbol_supplied))

	def get_name(self):
		return deref(self.ImpuritySpeciesPtr).get_name()
	def get_atomic_number(self):
		return deref(self.ImpuritySpeciesPtr).get_atomic_number()

cdef class PyRateEquations:
	cdef unique_ptr[RateEquations] RateEquationsPtr
	def __init__(self, PyImpuritySpecies impurity):
		self.RateEquationsPtr.reset(new RateEquations(deref(impurity.ImpuritySpeciesPtr)))

	def setThresholdDensity(self, double density_threshold):
		deref(self.RateEquationsPtr).setThresholdDensity(density_threshold)
	def setDominantIonMass(self, double mi_in_amu):
		deref(self.RateEquationsPtr).setDominantIonMass(mi_in_amu)
	def computeDerivs(self, double Te, double Ne, double Vi, double Nn, double Vn, vector[double]& Nzk, vector[double]& Vzk):
		return deref(self.RateEquationsPtr).computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk)







