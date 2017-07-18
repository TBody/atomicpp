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
			SD1DData(string input_file);
			void setImpurityFraction(float impurity_fraction);
			void setImpurityDensity(float impurity_density);
			void selectSingleTime(float t);
		private:
			vector<double> temperature;
			vector<double> density;
			vector<double> neutral_fraction;
			vector<double> impurity_density;
			double data_shape[2];
			double impurity_fraction;
		};
#endif