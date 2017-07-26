#ifndef SHAREDFUNCTIONS_H //Preprocessor directives to prevent multiple definitions
#define SHAREDFUNCTIONS_H

	#include <string>
	#include <fstream>

	#include "json.hpp"
	using json = nlohmann::json;

	/**
	 * @brief Returns a JSON object from a specified file location
	 * @details Reads a .json file given at path_to_file
	 * Uses the json module at https://github.com/nlohmann/json/
	 * This is supplied as a header file "json.hpp" which must be included in the same folder as the source
	 * N.b. do NOT pass path_to_file by reference - results in error (Undefined symbols)!
	 * @param path_to_file std::string, giving relative or absolute path to file
	 * @return JSON object - basically a std::map/dictionary which returns the stored value associated with a given (std::string) key
	 * The JSON object tries to convert data from a custom class to whatever you ask it to assign to - most of the time
	 * it works as intuition would suggest
	 */
	json retrieveFromJSON(std::string path_to_file);
	/**
	 * @brief A simple program to check if a file exists
	 * 
	 * @param name std::string, giving relative or absolute path to file
	 * @return boolean - true if file found, false if not
	 */
	bool test_file_exists (const std::string& name);
#endif