#ifndef IMPURITYSPECIES_H //Preprocessor directives to prevent multiple definitions
#define IMPURITYSPECIES_H

	#include <map>
	#include <string>
	#include "RateCoefficient.hpp"

	#include <memory>
	
	namespace atomicpp{
	class ImpuritySpecies{
		// # For storing OpenADAS data related to a particular impurity species
		// # Loosely based on cfe316/atomic/atomic_data.py/AtomicData class (although with much less code since
		// # all of the F77 importing is done in the seperate <<make json_update>> code since BOUT++ protocol
		// # requires fortran code be isolated from main operation)
	public:
		/**
		 * @brief ImpuritySpecies Constructor
		 * @details Uses the "impurity_user_input" file to std::set hard-coded impurity
		 * parameters corresponding to the impurity represented by the symbol
		 * Determines the OpenADAS file-dependencies and adds them to .adas_data_dict
		 * Makes RateCoefficient objects corresponding to the files and adds smart
		 * pointers for them to .rate_coefficients
		 * 
		 * @param impurity_symbol_supplied typically should be of length 1 (i.e. "c" for carbon, "n" for nitrogen)
		 */
		ImpuritySpecies(std::string& impurity_symbol_supplied);
		/**
		 * @brief Determines the OpenADAS json files which the impurity data is given 
		 * in and adds them to .adas_data_dict 
		 * Will throw a runtime error if the file isn't found in the database
		 * 
		 * @param physics_process a std::string corresponding to a physics process
		 * @param filetype_code the code used by OpenADAS to represent this process
		 * @param json_database_path relative or absolute path to where the json data files from OpenADAS are located
		 * @param year_fallback if data isn't found for the specified year, try another year (default 1996) rather than immediately throw an error
		 */
		void addJSONFiles(const std::string& physics_process, const std::string& filetype_code, const std::string& json_database_path, const int year_fallback = 1996);
		/**
		 * @brief Uses the OpenADAS files determined in addJSONFiles to construct a
		 * std::map between a physics_process std::string and a smart-pointer to a corresponding
		 * RateCoefficient object
		 * 	
		 */
		void makeRateCoefficients();
		std::string get_symbol();
		std::string get_name();
		int get_year();
		bool get_has_charge_exchange();
		int get_atomic_number();
		double get_mass();
		std::map<std::string,std::string> get_adas_files_dict();
		std::map<std::string,std::shared_ptr<RateCoefficient> > get_rate_coefficients();
		/**
		 * @brief Adds the key, value pair to the rate_coefficients map attribute
		 * 
		 * @param key, an std::string to access the value from the map
		 * @param value, a smart pointer to a RateCoefficient object
		 */
		void add_to_rate_coefficients(std::string key, std::shared_ptr<RateCoefficient> value);
		/**
		 * @brief Accesses the value of the rate_coefficient std::map corresponding
		 * to the supplied std::string key
		 * 
		 * @param key std::string corresponding to a physics_process
		 * @return shared (smart) pointer to a RateCoefficient object
		 */
		std::shared_ptr<RateCoefficient> get_rate_coefficient(const std::string& key);
		std::vector<double> calculateNzk(double Te, double Ne, double Nz, double Nn);
	private:
		// Data fields
		std::string symbol;
		std::string name;
		int year;
		bool has_charge_exchange;
		int atomic_number; //in elementary charges
		double mass; //in amu
		std::map<std::string,std::string> adas_files_dict;
		std::map<std::string,std::shared_ptr<RateCoefficient> > rate_coefficients;
	};

	const double eV_to_J = 1.60217662e-19; //Conversion factor between electron-volts and joules (effective units J/eV)
	const double amu_to_kg = 1.66054e-27; ////Conversion factor between atomic-mass-units and kilograms (effective units kg/amu)
	}
#endif