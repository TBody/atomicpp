#include <map>
#include <string>
#include <ostream>
#include <stdexcept> //For error-throwing
#include <memory> //For smart pointers
#include <set>

#include "ExternalModules/json.hpp"
using json = nlohmann::json;
#include "ImpuritySpecies.hpp"
#include "sharedFunctions.hpp"
#include "RateCoefficient.hpp"

using namespace atomicpp;
ImpuritySpecies::ImpuritySpecies(std::string& impurity_symbol_supplied){
	// Constructor for the ImpuritySpecies class
	// Input is a std::string, which typically should be of length 1 (i.e. "c" for carbon, "n" for nitrogen)
	// Will return an ImpuritySpecies object
	
	// Convert impurity_symbol_supplied to lowercase (this is what impurity_user_input uses as a key)
	impurity_symbol_supplied[0] = tolower(impurity_symbol_supplied[0]);

	// Declare the path to the JSON file which contains hard-coded parameters for the ADAS data corresponding to the impurities
	// std::string user_file = get_impurity_user_input();
	std::string user_file = "impurity_user_input.json";

	// Declare the path the directory which contains the JSON files created from OpenADAS .dat files
	// To convert .dat to .json, see https://github.com/TBody/OpenADAS_to_JSON (should be included as a subdirectory of SD1D)
	// std::string json_database_path = get_json_database_path();
	std::string json_database_path = "json_database/json_data";

	// Define a std::map between a recognisable key for various physics process, and the 3-letter code used by OpenADAS to represent this data
	std::map<std::string,std::string>datatype_abbrevs={
		{"ionisation",           "scd"},
		{"recombination",        "acd"},
		{"cx_rec",              "ccd"},
		{"continuum_power",      "prb"},
		{"line_power",           "plt"},
		{"cx_power",             "prc"},
		{"ionisation_potential", "ecd"} //N.b. ionisation_potential is not a rate-coefficient, but most of the methods are transferable
	};

	// Save from impurity_symbol_supplied and the hard-coded data within user_file
	symbol              = impurity_symbol_supplied;
	// retrieve from JSON and std::set attributes
	
	json j_object = retrieveFromJSON(user_file);
	auto check_symbol_in_file = j_object.find(impurity_symbol_supplied);
	if ((check_symbol_in_file != j_object.end())){
		name                = j_object[impurity_symbol_supplied]["name"];
		year                = j_object[impurity_symbol_supplied]["year"];
		has_charge_exchange = j_object[impurity_symbol_supplied]["has_charge_exchange"];
		atomic_number       = j_object[impurity_symbol_supplied]["atomic_number"];
		double mass_in_amu  = j_object[impurity_symbol_supplied]["mass"]; //in amu
		mass                = mass_in_amu;
	} else {
		std::cout << "Error: "<< impurity_symbol_supplied << " not found in " << user_file << "\n";
		throw std::invalid_argument( "Key not found in user input file" );
	};

	// Add the JSON files associated with this impurity to its .adas_files_dict attribute where the key is the (extended) process name, which std::maps to a filename (std::string)
	// Mark the physics processes which rely on charge exchange
	std::set<std::string> cx_processes{datatype_abbrevs["cx_rec"],datatype_abbrevs["cx_power"]};
	// Iterate over the key -> value pairs in datatype_abbrevs
	for (auto& kv : datatype_abbrevs) {
		std::string physics_process = kv.first;
		std::string filetype_code = kv.second;
		// Check to see if physics_process is a charge exchange process
		bool code_is_cx = find(begin(cx_processes), end(cx_processes), filetype_code) != end(cx_processes);
		// If it is, and the element doesn't have charge exchange data, then skip this process
		if (get_has_charge_exchange() or not(code_is_cx)){
			// Otherwise, add the filenames of the json files corresponding to the rate-coefficients for the physics process to .adas_files_dict
		    addJSONFiles(physics_process, filetype_code, json_database_path);
		};
	}

	// Use the .adas_file_dict files to generate RateCoefficient objects for each process
	// Uses the same keys as .adas_file_dict
	makeRateCoefficients();

};
void ImpuritySpecies::addJSONFiles(const std::string& physics_process, const std::string& filetype_code, const std::string& json_database_path, const int year_fallback /* = 1996*/){
	// # 1. Make the filename std::string expected for the json adas file
	// # 2. Check that this file exists in the JSON_database_path/json_data directory
	// # 3. Add this file to the atomic data .adas_files_dict attribute
	std::string filename;
	std::string year_to_string = std::to_string(year);

	filename = json_database_path + "/" + filetype_code + year_to_string.substr(2,4) + "_" + symbol + ".json";

	if (test_file_exists(filename)) {
		adas_files_dict[physics_process] = filename;
	} else {
		std::printf("%s not found in the specified JSON database\n", filename.c_str());

		std::string filename_fallback;
		std::string year_to_string_fallback = std::to_string(year_fallback);
		filename_fallback = json_database_path + "/" + filetype_code + year_to_string_fallback.substr(2,4) + "_" + symbol + ".json";
		if (test_file_exists(filename_fallback)) {
			std::printf("Warning: Fallback (%s) will be used\n", filename_fallback.c_str());
			adas_files_dict[physics_process] = filename_fallback;
		} else {
			std::printf("Error: Fallback to %s failed\n", year_to_string_fallback.c_str());
			throw std::runtime_error("File not found error (in ImpuritySpecies::addJSONFiles)");
		}
	}

};
void ImpuritySpecies::makeRateCoefficients(){
	// # Calls the RateCoefficient constructor method for each entry in the .adas_files_dict
	// # Generates a dictionary smart pointer to RateCoefficient objects as .rate_coefficients
	// See http://umich.edu/~eecs381/handouts/C++11_smart_ptrs.pdf for information on smart pointers (memory-managed)

	for (auto& kv : get_adas_files_dict()) {
		std::string physics_process = kv.first;
		std::string filename = kv.second;
		// Make a new RateCoefficient object by calling the RateCoefficient constructor on 'filename'
		// Create a smart pointer 'RC' that points to this object
		std::shared_ptr<RateCoefficient> RC(new RateCoefficient(filename));
		// Add 'RC' to the rate_coefficients attribute of ImpuritySpecies
		// (n.b. this is a std::map from a std::string 'physics_process' to a smart pointer which points to a RateCoefficient object)
		rate_coefficients[physics_process] = RC;
	}
};
// Accessor functions
	std::string ImpuritySpecies::get_symbol(){
		return symbol;
	};
	std::string ImpuritySpecies::get_name(){
		return name;
	};
	int ImpuritySpecies::get_year(){
		return year;
	};
	bool ImpuritySpecies::get_has_charge_exchange(){
		return has_charge_exchange;
	};
	int ImpuritySpecies::get_atomic_number(){
		return atomic_number;
	};
	double ImpuritySpecies::get_mass(){
		return mass;
	};
	std::map<std::string,std::string> ImpuritySpecies::get_adas_files_dict(){
		return adas_files_dict;
	};
	std::map<std::string,std::shared_ptr<RateCoefficient> > ImpuritySpecies::get_rate_coefficients(){
		return rate_coefficients;
	};
	void ImpuritySpecies::add_to_rate_coefficients(std::string key, std::shared_ptr<RateCoefficient> value){
		rate_coefficients[key] = value;
	};
	std::shared_ptr<RateCoefficient> ImpuritySpecies::get_rate_coefficient(const std::string& key){
		return rate_coefficients[key];
	};

	std::vector<double> ImpuritySpecies::calculateNzk(double Te, double Ne, double Nz, double Nn){
		// Calculates the relative distribution across ionisation stages of the impurity by assuming collisional-radiative
		// equilbrium. This is then used to calculate the density within each state, allowing the total power at a point
		// to be evaluated
		// std::cout << "Called for Te = " << Te << ", Ne = " << Ne << ", Nz = " << Nz << ", Nn = " << Nn << std::endl;

		int Z = atomic_number;
		std::vector<double> iz_stage_distribution(Z+1);

		// std::set GS density equal to 1 (arbitrary)
		iz_stage_distribution[0] = 1;
		double sum_iz = 1;

		// Loop over 0, 1, ..., Z-1
		// Each charge state is std::set in terms of the density of the previous
		for(int k=0; k<Z; ++k){
			// Ionisation
			// Get the RateCoefficient from the rate_coefficient std::map (atrribute of impurity)
			std::shared_ptr<RateCoefficient> iz_rate_coefficient = rate_coefficients["ionisation"];
			// Evaluate the RateCoefficient at the point
			double k_iz_evaluated = iz_rate_coefficient->call0D(k, Te, Ne);

			// Recombination
			// Get the RateCoefficient from the rate_coefficient std::map (atrribute of impurity)
			std::shared_ptr<RateCoefficient> rec_rate_coefficient = rate_coefficients["recombination"];
			// Evaluate the RateCoefficient at the point
			double k_rec_evaluated = rec_rate_coefficient->call0D(k, Te, Ne);

			std::shared_ptr<RateCoefficient> cx_rate_coefficient = rate_coefficients["cx_rec"];
			// Evaluate the RateCoefficient at the point
			double k_cx_evaluated = cx_rate_coefficient->call0D(k, Te, Ne);

			// The ratio of ionisation from the (k)th stage and recombination from the (k+1)th std::sets the equilibrium densities
			// of the (k+1)th stage in terms of the (k)th (since R = Nz * Ne * rate_coefficient) N.b. Since there is no
			// ionisation from the bare nucleus, and no recombination onto the neutral (ignoring anion formation) the 'k'
			// value of ionisation coeffs is shifted down  by one relative to the recombination coeffs - therefore this
			// evaluation actually gives the balance

			iz_stage_distribution[k+1] = iz_stage_distribution[k] * (k_iz_evaluated*Ne/(k_rec_evaluated*Ne + k_cx_evaluated*Nn));
			sum_iz += iz_stage_distribution[k+1];
		}

		// # Normalise such that the sum over all ionisation stages is '1' at all points
		for(int k=0; k<=Z; ++k){
			iz_stage_distribution[k] = iz_stage_distribution[k] / sum_iz;
		}

		std::vector<double> Nzk(Z+1);
		for(int k=0; k<=Z; ++k){
			Nzk[k] = iz_stage_distribution[k]*Nz;
		}
		return Nzk;
	}