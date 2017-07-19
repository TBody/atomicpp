#ifndef SD1DDATA_H
#define SD1DDATA_H

	#include <string>
	#include <vector>
	#include <fstream>
	#include "json.hpp"
	using namespace std; //saves having to prepend std:: onto common functions

	// for convenience
	using json = nlohmann::json;

	class SD1DData{
		// # For storing the data output from SD1D. To create the required JSON run the function
		// # data_dict_export.py in an I/O (case) folder in SD1D.
		public:
			SD1DData(const string& input_file, double impurity_fraction);
		private:
			vector<double> temperature;
			vector<double> density;
			vector<double> neutral_fraction;
			double impurity_fraction;
			vector<double> impurity_density;
			// double data_shape[2];
		};
#endif