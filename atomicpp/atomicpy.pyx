# distutils: language = c++
# distutils: sources = [ImpuritySpecies.cpp, RateCoefficient.cpp, sharedFunctions.cpp, RateEquations.cpp]

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
from libcpp.string cimport string #to pass strings to/from C++
from libcpp.map cimport map
from libcpp.memory cimport shared_ptr

# Note - only need to declare the functions that you will use. Internal functions and attributes do not need to be declared.

cdef extern from "ImpuritySpecies.hpp" namespace "atomicpp":
	cdef cppclass ImpuritySpecies:
		ImpuritySpecies(string& impurity_symbol_supplied) except +

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
	cdef ImpuritySpecies *ImpuritySpeciesPtr
	def __cinit__(self, string impurity_symbol_supplied):
		self.ImpuritySpeciesPtr = new ImpuritySpecies(impurity_symbol_supplied)
	def __dealloc__(self):
		del self.ImpuritySpeciesPtr

cdef class PyRateEquations:
	cdef RateEquations *RateEquationsPtr
	def __cinit__(self, PyImpuritySpecies impurity):
		self.RateEquationsPtr = new RateEquations(impurity.ImpuritySpeciesPtr[0]) #will have to play about to get this working with default values
	def __dealloc__(self):
		del self.RateEquationsPtr
	def setThresholdDensity(self, double density_threshold):
		self.RateEquationsPtr.setThresholdDensity(density_threshold)
	def setDominantIonMass(self, double mi_in_amu):
		self.RateEquationsPtr.setDominantIonMass(mi_in_amu)
	def computeDerivs(self, double Te, double Ne, double Vi, double Nn, double Vn, vector[double]& Nzk, vector[double]& Vzk):
		return self.RateEquationsPtr.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk)







