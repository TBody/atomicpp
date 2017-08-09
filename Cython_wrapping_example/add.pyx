# distutils: language = c++
# distutils: sources = Adder.cpp

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
from libcpp.string cimport string #to pass strings to/from C++

cdef extern from "Adder.hpp":
	cdef struct twoInts:
		int a
		int b

# c++ interface to cython
cdef extern from "Adder.hpp":
	cdef cppclass Adder:
		Adder()
		Adder(vector[double] Input)
		vector[double] ReturnVector() #passing vectors out
		# void PlusOne()
		void PlusTwo() #can call this method, even though the internal method PlusOne() is not declared to the Python interface
		void PlusVector(vector[double] vector_to_add) #passing vectors in
		string Print()
		string sayHello() #can call this method, even though the internal attribute privatestring is not declared to the Python interface
		vector[double] internal
		twoInts returntwoInts()

# creating a cython wrapper class
cdef class PyAdder:
	cdef Adder *AdderPtr
	def __cinit__(self, vector[double] Input):
		self.AdderPtr = new Adder(Input)
	def __dealloc__(self):
		del self.AdderPtr
	def ReturnVector(self):
		return self.AdderPtr.ReturnVector()
	# def PlusOne(self):
	# 	self.AdderPtr.PlusOne()
	def PlusTwo(self):
		self.AdderPtr.PlusTwo()
	def PlusVector(self, vector_to_add):
		self.AdderPtr.PlusVector(vector_to_add)
	def Print(self):
		return self.AdderPtr.Print()
	def sayHello(self):
		return self.AdderPtr.sayHello()
	def returntwoInts(self):
		return self.AdderPtr.returntwoInts()
	