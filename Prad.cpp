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
#include <fstream>

#include "atomicpp/ImpuritySpecies.hpp"
#include "atomicpp/RateCoefficient.hpp"
#include "atomicpp/SD1DData.hpp"

// Look at
// http://kluge.in-chemnitz.de/opensource/spline/
// for spline interpolation

#include "atomicpp/json.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

// How to use auto keyword: http://www.acodersjourney.com/2016/02/c-11-auto/
// (with some good code examples)

// N.b. perform all plotting routines as Python post-processing

// Shared variables
map<string,string>datatype_abbrevs={
	{"ionisation",           "scd"},
	{"recombination",        "acd"},
	{"cx_recc",              "ccd"},
	{"continuum_power",      "prb"},
	{"line_power",           "plt"},
	{"cx_power",             "prc"},
	{"ionisation_potential", "ecd"} //N.b. ionisation_potential is not a rate-coefficient, but most of the methods are transferable
};

// iz_stage_distribution = calculateCollRadEquilibrium(impurity, experiment)
// computeRadiatedPower(impurity, experiment, iz_stage_distribution)



int main(){
	const string user_file="user_input.json";
	const string input_file="sd1d-case-05.json";
	const string json_database_path="json_database/json_data";
	const string impurity_symbol="c"; //Must be lower-case! Haven't yet written .tolower into the code

	// make an ImpuritySpecies object 'impurity' from the user_input.json file and the impurity_symbol variable
	ImpuritySpecies impurity(impurity_symbol, user_file);

	// cout << "The atomic_number of c is " << j_object["c"]["atomic_number"] << "\n";
	cout << "The atomic_number of c is " << impurity.get_year()+ 1<< "\n";
	// std::cout << "value at key \"c\": " << *it_carbon << '\n';

	cout<<datatype_abbrevs["cx_power"]<<"\n";

	// # Calculate the distribution across ionisation stages, assuming collisional-radiative equilibrium
	// iz_stage_distribution = calculateCollRadEquilibrium(impurity, experiment)

	// # Compute radiated power
	// # Returns total_power, stage_integrated_power (sum over all ionisation stages), and
	// # radiated_power (resolved into different ionisation stages, with 'total' giving sum over
	// # all physics_processes)
	// # 	stage_integrated_power and radiated_power are dictionaries with physics_process keys and
	// # 	data_length shape for stage_integrated_power and [Z, data_length] shape for radiated_power
	// # 	total_power is an array of shape data_length
	// computeRadiatedPower(impurity, experiment, iz_stage_distribution)
}





























