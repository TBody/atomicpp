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

using namespace std;

ImpuritySpecies::ImpuritySpecies(string symbol, string user_file){
	cout << "Constructing ImpuritySpecies object for " << symbol << "\n";
	symbol              = symbol;

	json j_object = retrieveFromJSON(user_file);

	auto check_symbol_in_file = j_object.find(symbol);
	if ((check_symbol_in_file != j_object.end())){
		name                = j_object[symbol]["name"];
		year                = j_object[symbol]["year"];
		has_charge_exchange = j_object[symbol]["has_charge_exchange"];
		atomic_number       = j_object[symbol]["atomic_number"];
	} else {
		cout << "Error: "<< symbol << " not found in " << user_file << "\n";
		throw std::invalid_argument( "Key not found in user input file" );
	};

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