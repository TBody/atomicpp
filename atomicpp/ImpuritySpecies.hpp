#ifndef IMPURITYSPECIES_H //Preprocessor directives to prevent multiple definitions
#define IMPURITYSPECIES_H
  // 1. Classes
  // 2. Constants
  // 3. Non-member functions
  // 4. Global variables


	#include <map>
	#include <string>
	#include "RateCoefficient.hpp"
	
	using namespace std;

	class ImpuritySpecies{
		// # For storing OpenADAS data related to a particular impurity species
		// # Loosely based on cfe316/atomic/atomic_data.py/AtomicData class (although with much less code since
		// # all of the F77 importing is done in the seperate <<make json_update>> code since BOUT++ protocol
		// # requires fortran code be isolated from main operation)
	public:
		// For constructor and member function declarations
		ImpuritySpecies(const string& symbol, const string& user_file);
		void addJSONFiles(const string& physics_process, const string& filetype_code, const string& json_database_path);
		void makeRateCoefficients();
		string get_symbol();
		string get_name();
		int get_year();
		bool get_has_charge_exchange();
		int get_atomic_number();
		map<string,string> get_adas_files_dict();

	private:
		// Data fields
		string symbol;
		string name;
		int year;
		bool has_charge_exchange;
		int atomic_number;
		map<string,string> adas_files_dict;
		map<string,RateCoefficient> rate_coefficients;
	};

#endif