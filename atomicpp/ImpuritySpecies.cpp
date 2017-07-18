// 1. Classes
// 2. Constants
// 3. Non-member functions
// 4. Global variables

#include <map>
#include <string>
#include <iostream>
#include <stdexcept> //For error-throwing

#include "json.hpp"
using json = nlohmann::json;
#include "ImpuritySpecies.hpp"
#include "sharedFunctions.hpp"
#include "RateCoefficient.hpp"

using namespace std;

ImpuritySpecies::ImpuritySpecies(const string& symbol_key_supplied, const string& user_file){
	// cout << "Constructing ImpuritySpecies object for " << symbol_key_supplied << "\n";
	symbol              = symbol_key_supplied;

	json j_object = retrieveFromJSON(user_file);

	auto check_symbol_in_file = j_object.find(symbol_key_supplied);
	if ((check_symbol_in_file != j_object.end())){
		name                = j_object[symbol_key_supplied]["name"];
		year                = j_object[symbol_key_supplied]["year"];
		has_charge_exchange = j_object[symbol_key_supplied]["has_charge_exchange"];
		atomic_number       = j_object[symbol_key_supplied]["atomic_number"];
	} else {
		cout << "Error: "<< symbol_key_supplied << " not found in " << user_file << "\n";
		throw invalid_argument( "Key not found in user input file" );
	};

	// adas_files_dict={};
};
void ImpuritySpecies::addJSONFiles(const string& physics_process, const string& filetype_code, const string& json_database_path){
	// # 1. Make the filename string expected for the json adas file
	// # 2. Check that this file exists in the JSON_database_path/json_data directory
	// # 3. Add this file to the atomic data .adas_files_dict attribute
	string filename;
	string year_to_string = to_string(year);

	filename = json_database_path + "/" + filetype_code + year_to_string.substr(2,4) + "_" + symbol + ".json";

	// cout << filename + "\n";
	
	// cout << boolalpha;
	// cout << "File found: " << test_file_exists(filename) << "\n";

	adas_files_dict[physics_process] = filename;
};
void ImpuritySpecies::makeRateCoefficients(){
	// # Calls the RateCoefficient constructor method for each entry in the .adas_files_dict
	// # Generates a dictionary of RateCoefficient objects as .rate_coefficients
	
	cout << "makeRateCoefficients called \n";

	string filename;
	filename = adas_files_dict["ionisation"];
	cout << filename << "\n";

	RateCoefficient r(filename);


	// for physics_process, filename in self.adas_files_dict.items():
	// 	full_path = '{}/json_data/{}'.format(JSON_database_path,filename)
	// 	self.rate_coefficients[physics_process] = RateCoefficient(full_path)
};
// Accessor functions
	string ImpuritySpecies::get_symbol(){
		return symbol;
	};
	string ImpuritySpecies::get_name(){
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
	map<string,string> ImpuritySpecies::get_adas_files_dict(){
		return adas_files_dict;
	};