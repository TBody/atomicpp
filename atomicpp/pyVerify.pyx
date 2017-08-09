# distutils: language = c++
# distutils: sources = [ImpuritySpecies.cpp, RateCoefficient.cpp, sharedFunctions.cpp, RateEquations.cpp]

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector #to pass vectors to/from C++
from libcpp.string cimport string #to pass strings to/from C++
from libcpp.map cimport map
from libcpp.memory cimport shared_ptr

cdef extern from "ImpuritySpecies.hpp":
	cdef cppclass ImpuritySpecies:
		ImpuritySpecies(string& impurity_symbol_supplied)
		void addJSONFiles(string& physics_process, string& filetype_code, string& json_database_path, int year_fallback = 1996) except +
		void makeRateCoefficients()
		string get_symbol()
		string get_name()
		int get_year()
		bool get_has_charge_exchange()
		int get_atomic_number()
		double get_mass()
		map[string,string] get_adas_files_dict()
		map[string,std::shared_ptr<RateCoefficient> ] get_rate_coefficients()
		bool get_has_shared_interpolation()
		void add_to_rate_coefficients(string key, std::shared_ptr<RateCoefficient> value)
		std::shared_ptr<RateCoefficient> get_rate_coefficient(string& key)
		void initialiseSharedInterpolation()
	
		string symbol
		string name
		int year
		bool has_charge_exchange
		int atomic_number
		double mass
		std::map<string,string> adas_files_dict
		std::map<string,std::shared_ptr<RateCoefficient> > rate_coefficients
		bool has_shared_interpolation


		
string get_json_database_path()
string get_impurity_user_input()

double eV_to_J = 1.60217662e-19
double amu_to_kg = 1.66054e-27

