#include <string>
#include <fstream>
#include <stdexcept> //For error-throwing
#include <sstream>

#include "sharedFunctions.hpp"

#include "ExternalModules/json.hpp"
using json = nlohmann::json;

using namespace atomicpp;

json atomicpp::retrieveFromJSON(std::string path_to_file){
	// Do not pass path_to_file by reference - results in error!
	// Reads a .json file given at path_to_file
	// Uses the json module at https://github.com/nlohmann/json/
	// This relies upon the "json.hpp" header which must be included in the same folder as the source
	
	// Open a file-stream at path_to_file
	if (atomicpp::test_file_exists(path_to_file)){
		std::ifstream json_file(path_to_file);
		// Initialise a json file object at j_object
		json j_object;
		json_file >> j_object;
		return j_object;
	} else {
		std::stringstream errMsg;
		errMsg << "atomicpp::sharedFunctions::retrieveFromJSON error: file (" << path_to_file << ") not found";
		throw std::runtime_error(errMsg.str());
	}
};
bool atomicpp::test_file_exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
};