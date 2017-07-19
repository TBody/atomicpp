#ifndef RATECOEFFICIENT_H
#define RATECOEFFICIENT_H
	#include <string>
	#include <vector>
	#include <fstream>
	#include "json.hpp"

	using namespace std; //saves having to prepend std:: onto common functions

	// for convenience
	using json = nlohmann::json;
	class RateCoefficient{
		// # For storing the RateCoefficients encoded in an OpenADAS data file
		// # Intended to be called from the .makeRateCoefficients method of an ImpuritySpecies object
		// #
		// # Closely based on the cfe316/atomic/atomic_data.py/RateCoefficient class
		// #
		// # Interpolation tables for the rate of some physical process.
		// # Contains one 2D spline interpolator for each charge state of an element,
		// # per  process like 'ionisation', 'recombination'.

		// # Attributes:
		// #     atomic_number (int) : The element's Z.
		// #     element (str)       : Short element name like 'c'
		// #     adf11_file (str)    : The /full/filename it came from (link to .json, not .dat)
		// #     log_temperature     : np.array of log10 of temperature values
		// #     log_density         : np.array of log10 of density values
		// #     log_coeff           : a 3D np.array with shape (Z, temp, dens)
		// #         The list has length Z and is interpolations of log_coeff.
		public:
			RateCoefficient(const string& filename);
			double call0D(const int k, const double log10_Te, const double log10_ne);
			void call1D(const int k, const int data_length, const vector<double>& log10_Te, const vector<double>& log10_ne, vector<double>& log10_coeffs);
			friend ostream& operator<<(ostream& os, const RateCoefficient& RC); //Define the __str__ return to cout
			int get_atomic_number();
			string get_element();
			string get_adf11_file();
			vector<vector< vector<double> > > get_log_coeff();
			vector<double> get_log_temperature();
			vector<double> get_log_density();
		private:
			int atomic_number;
			string element;
			string adf11_file;
			vector<vector< vector<double> > > log_coeff;
			vector<double> log_temperature;
			vector<double> log_density;
		};
#endif