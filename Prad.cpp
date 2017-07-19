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

#include <memory>

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

vector<vector<double> > calculateCollRadEquilibrium(ImpuritySpecies& impurity, SD1DData& experiment){
	vector<vector<double> > iz_stage_distribution;

	int Z = impurity.get_atomic_number();

	int data_length = experiment.get_temperature().size();
	assert(data_length == experiment.get_density().size());

	// Perform this calculation for each charge state
	// vector<double> k_state_density(data_length);
	// k_state_density[0] = 1;

	// for(int k=0; k<Z; ++k){
	// }
	
	// Load from experiment for contiguous access?
	vector<double> temperature = experiment.get_temperature();
	vector<double> density = experiment.get_density();
	// Take log10, pipelining?
	vector<double> log10_temperature(data_length);
	vector<double> log10_density(data_length);
	vector<double> log10_iz_coeffs(data_length);
	vector<double> log10_rec_coeffs(data_length);
	vector<double> log10_coeffs(data_length);
	for(int i=0; i<data_length; ++i){
		// Could split loops for contiguous access -- loop overheard versus unit march
		log10_temperature[i] = log10(temperature[i]);
		log10_density[i] = log10(density[i]);
	}
	// Should this be i<data_length or <= data_length? temperature[data_length] is off the end of the array (uncomment below to see)
	// cout << temperature[data_length-1] << endl;
	// cout << log10_temperature[data_length-1] << endl;
	
	int k = 1;
	shared_ptr<RateCoefficient> iz_rate_coefficient = impurity.get_rate_coefficient("ionisation");
	// Call1D. k & data_length are const int, temperature and density are const vector<double>&, coeff is vector<double>
	iz_rate_coefficient->call1D(k, data_length, log10_temperature, log10_density, log10_iz_coeffs);
	shared_ptr<RateCoefficient> rec_rate_coefficient = impurity.get_rate_coefficient("recombination");
	rec_rate_coefficient->call1D(k, data_length, log10_temperature, log10_density, log10_rec_coeffs);
	for(int s=0; s<data_length; ++s){

	}




	// // Define and initialize a vector with 2D array
	// vector<vector<double>> iz_stage_distribution = {
	// 	vector<double>(begin(arr[0]), end(arr[0])),
	// 	vector<double>(begin(arr[1]), end(arr[1]))
	// };

	return iz_stage_distribution;
}

int main(){
	const string user_file="user_input.json";
	const string input_file="sd1d-case-05.json";
	const string json_database_path="json_database/json_data";
	string impurity_symbol="c";

	impurity_symbol[0] = tolower(impurity_symbol[0]);

	// make an ImpuritySpecies object 'impurity' from the user_input.json file and the impurity_symbol variable
	ImpuritySpecies impurity(impurity_symbol, user_file);

	// # Add the JSON files associated with this impurity to its .adas_files_dict attribute
	// # where the key is the (extended) process name, which maps to a filename (string)
	// # <<TODO>> Check that these files exist in the JSON_database_path/json_data/ directory
	vector<string> cx_processes{"ccd","prc"};
	string physics_process;
	string filetype_code;
	for (auto& kv : datatype_abbrevs) {
		physics_process = kv.first;
		filetype_code = kv.second;
		bool code_is_cx = find(begin(cx_processes), end(cx_processes), filetype_code) != end(cx_processes);
		// See if the datatype code is in the list of cx_processes
		if (impurity.get_has_charge_exchange() or not(code_is_cx)){
		    impurity.addJSONFiles(physics_process, filetype_code, json_database_path);
		};
	}

	// Print get_adas_files_dict to check copy
	// for (auto& kv : impurity.get_adas_files_dict()) {
	// 	cout << kv.first << ": " << kv.second << "\n";
	// }

	// # Use the .adas_file_dict files to generate RateCoefficient objects for each process
	// # Uses the same keys as .adas_file_dict
	impurity.makeRateCoefficients();

	// # Process the input_file to extract
	// # 	density(s) 					= electron density (in m^-3)
	// # 	temperature(s)				= electron/ion temperature (in eV)
	// # 	neutral_fraction(s)			= neutral density/electron density (no units)
	// # where s is 1D distance index. Time is already contracted (using final time-step)
	// Second argument is impurity fraction
	SD1DData experiment(input_file, 1e-2);

	// # Calculate the distribution across ionisation stages, assuming collisional-radiative equilibrium
	vector<vector<double> > iz_stage_distribution = calculateCollRadEquilibrium(impurity, experiment);

	// # Compute radiated power
	// # Returns total_power, stage_integrated_power (sum over all ionisation stages), and
	// # radiated_power (resolved into different ionisation stages, with 'total' giving sum over
	// # all physics_processes)
	// # 	stage_integrated_power and radiated_power are
	//  dictionaries with physics_process keys and
	// # 	data_length shape for stage_integrated_power and [Z, data_length] shape for radiated_power
	// # 	total_power is an array of shape data_length
	// computeRadiatedPower(impurity, experiment, iz_stage_distribution)
}





























