// Program name: atomic++/Prad.cpp
// Author: Thomas Body
// Author email: tajb500@york.ac.uk
// Date of creation: 17 July 2017
//
// Program function: output the radiated power (Prad)
//                   by using OpenADAS rates on output JSON from SD1D run
//
// Based on the TBody/atomic1D code, which is in turn based on the cfe316/atomic code
//
// Under active development: <<TODO>> indicates development goal

// Include declarations
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "json.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

// N.b. perform all plotting routines as Python post-processing

map<string,string>datatype_abbrevs={
	{"ionisation",           "scd"},
	{"recombination",        "acd"},
	{"cx_recc",              "ccd"},
	{"continuum_power",      "prb"},
	{"line_power",           "plt"},
	{"cx_power",             "prc"},
	{"ionisation_potential", "ecd"} //N.b. ionisation_potential is not a rate-coefficient, but most of the methods are transferable
};
const string input_file="sd1d-case-05.json";
const string json_database_path="json_database/json_data";
const string impurity_symbol="C";

void retrieveFromJSON(
	);

class ImpuritySpecies{
		// # For storing OpenADAS data related to a particular impurity species
		// # Loosely based on cfe316/atomic/atomic_data.py/AtomicData class (although with much less code since
		// # all of the F77 importing is done in the seperate <<make json_update>> code since BOUT++ protocol
		// # requires fortran code be isolated from main operation)
	public:
		// For constructor and member function declarations
		ImpuritySpecies(string symbol);
		void addJSONFiles();
		void makeRateCoefficients();

	private:
		// Data fields
		string symbol;
		string name;
		int year;
		bool has_charge_exchange;
		int atomic_number;
		map<string,string> adas_files_dict;
		// map<string,RateCoefficient> rate_coefficients;
	};
	ImpuritySpecies::ImpuritySpecies(string symbol){
		cout << "Constructing ImpuritySpecies object for " << symbol << "\n";
		symbol="";
		name="";
		year=0;
		has_charge_exchange=true;
		atomic_number=0;
		// adas_files_dict={};
	};
	void ImpuritySpecies::addJSONFiles(){
		// # 1. Make the filename string expected for the json adas file
		// # 2. Check that this file exists in the JSON_database_path/json_data directory
		// # 3. Add this file to the atomic data .adas_files_dict attribute
	};
	void ImpuritySpecies::makeRateCoefficients(){
		// # Calls the RateCoefficient constructor method for each entry in the .adas_files_dict
		// # Generates a dictionary of RateCoefficient objects as .rate_coefficients
	};

class RateCoefficient{
	// # For storing the RateCoefficients encoded in an OpenADAS data file
	// # Intended to be called from the .makeRateCoefficients method of an ImpuritySpecies object
	// #
	// # Closely based on the cfe316/atomic/atomic_data.py/RateCoefficient class
	// #
	// # Interpolation tables for the rate of some physical process.
	// # Contains one 2D spline interpolator for each charge state of an element,
	// # per  process like 'ionisation', 'recombination'.

	// # Attributes:
	// #     atomic_number (int) : The element's Z.
	// #     element (str)       : Short element name like 'c'
	// #     adf11_file (str)    : The /full/filename it came from (link to .json, not .dat)
	// #     log_temperature     : np.array of log10 of temperature values
	// #     log_density         : np.array of log10 of density values
	// #     log_coeff           : a 3D np.array with shape (Z, temp, dens)
	// #     splines             : list of scipy.interpolate.fitpack2.RectBivariateSpline
	// #         The list has length Z and is interpolations of log_coeff.
	public:
		RateCoefficient(ImpuritySpecies& impurity, string filename);
		void compute_interpolating_splines(); //Could consider fixing length, since it will always be the same shape
		vector<double> call1D(int k, double Te, double ne);
	private:
		int atomic_number;
		string element;
		string adf11_file;
		// double log_temperature [x] [y]; //Need to set x, y. Dynamic sizing of arrays? Size is set from JSON
		// double log_density [x] [y]; //
		// double log_coeff [x] [y]; //
		// splines (interpolation functions on 2D grid which can be called to return value)
		// see https://en.wikipedia.org/wiki/List_of_numerical_libraries#C.2B.2B for C++ math libraries
	};
	RateCoefficient::RateCoefficient(ImpuritySpecies& impurity, string filename){
	};
	void RateCoefficient::compute_interpolating_splines(){
	}; //Could consider fixing length, since it will always be the same shape
	vector<double> RateCoefficient::call1D(int k, double Te, double ne){
	};

class SD1DData{
	// # For storing the data output from SD1D. To create the required JSON run the function
	// # data_dict_export.py in an I/O (case) folder in SD1D.
	public:
		SD1DData(string input_file);
		void setImpurityFraction(float impurity_fraction);
		void setImpurityDensity(float impurity_density);
		void selectSingleTime(float t);
	private:
		vector<double> temperature;
		vector<double> density;
		vector<double> neutral_fraction;
		vector<double> impurity_density;
		double data_shape[2];
		double impurity_fraction;
	};
	SD1DData::SD1DData(string input_file){
	};
	void SD1DData::setImpurityFraction(float impurity_fraction){
	};
	void SD1DData::setImpurityDensity(float impurity_density){
	};
	void SD1DData::selectSingleTime(float t){
	};


int main(){

	cout<<"Hello world\n";

	ImpuritySpecies("D");


	cout<<datatype_abbrevs["cx_power"]<<"\n";
}





























