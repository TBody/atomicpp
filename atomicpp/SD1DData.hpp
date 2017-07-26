#ifndef SD1DDATA_H
#define SD1DDATA_H

	#include <string>
	#include <vector>
	#include <fstream>
	#include "json.hpp"

	using json = nlohmann::json;

	class SD1DData{
		// # For storing the data output from SD1D. To create the required JSON run the function
		// # data_dict_export.py in an I/O (case) folder in SD1D.
		public:
			SD1DData(const std::string& expt_results_json, double impurity_fraction);
			std::vector<double> get_temperature();
			std::vector<double> get_density();
			std::vector<double> get_neutral_fraction();
			double get_impurity_fraction();
			std::vector<double> get_impurity_density();
		private:
			std::vector<double> temperature;
			std::vector<double> density;
			std::vector<double> neutral_fraction;
			double impurity_fraction;
			std::vector<double> impurity_density;
		};
#endif