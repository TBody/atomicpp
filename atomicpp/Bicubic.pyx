# distutils: language = c++
# distutils: sources = BicubicSpline.cpp

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
# from libcpp.string cimport string #to pass strings to/from C++
# from libcpp.map cimport map
from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

# Note - only need to declare the functions that you will use. Internal functions and attributes do not need to be declared.
# See https://dmtn-013.lsst.io/ for the use of unique_ptr in classes

cdef extern from "BicubicSpline.hpp" namespace "atomicpp":
	cdef cppclass BicubicSpline:
		# BicubicSpline()
		BicubicSpline( vector[double]& _x_values, vector[double]& _y_values, vector[ vector[double] ]& _z_values )
		double call0D(double eval_x, double eval_y)
		vector[ vector[double] ] get_z_values()
		vector[double] get_x_values()
		vector[double] get_y_values()

cdef class PyBicubicSpline:
	cdef unique_ptr[BicubicSpline] BicubicSplinePtr
	def __init__(self, vector[double]& _x_values, vector[double]& _y_values, vector[ vector[double] ]& _z_values ):
		self.BicubicSplinePtr.reset(new BicubicSpline(_x_values, _y_values, _z_values))
	def call0D(self, double eval_x, double eval_y):
		return deref(self.BicubicSplinePtr).call0D(eval_x, eval_y)