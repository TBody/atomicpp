#ifndef SHAREDFUNCTIONS_H //Preprocessor directives to prevent multiple definitions
#define SHAREDFUNCTIONS_H

	#include <string>
	#include <fstream>

	#include "json.hpp"
	using json = nlohmann::json;

	using namespace std; //saves having to prepend std:: onto common functions
	// for convenience

	json retrieveFromJSON(string path_to_file);
	// Reads a .json file given at path_to_file
	// Uses the json module at https://github.com/nlohmann/json/
	// This relies upon the "json.hpp" header which must be included in the same folder as the source
#endif