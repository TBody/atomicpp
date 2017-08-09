# distutils: language = c++
# distutils: sources = Adder.cpp

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
from libcpp.string cimport string #to pass strings to/from C++

# c++ interface to cython
cdef extern from "Adder.hpp":
	cdef cppclass Adder:
		Adder()
		Adder(vector[double] Input)
		vector[double] ReturnVector()
		void PlusOne()
		void PlusVector(vector[double] vector_to_add)
		string Print()
		vector[double] internal

# creating a cython wrapper class
cdef class PyAdder:
	cdef Adder *AdderPtr
	def __cinit__(self, vector[double] Input):
		self.AdderPtr = new Adder(Input)
	def __dealloc__(self):
		del self.AdderPtr
	def ReturnVector(self):
		return self.AdderPtr.ReturnVector()
	def PlusOne(self):
		self.AdderPtr.PlusOne()
	def PlusVector(self, vector_to_add):
		self.AdderPtr.PlusVector(vector_to_add)
	def Print(self):
		return self.AdderPtr.Print()
	