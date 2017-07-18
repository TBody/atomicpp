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
		// #     splines             : list of scipy.interpolate.fitpack2.RectBivariateSpline
		// #         The list has length Z and is interpolations of log_coeff.
		public:
			RateCoefficient(string filename);
			void compute_interpolating_splines(); //Could consider fixing length, since it will always be the same shape
			vector<double> call1D(int k, double Te, double ne);
		private:
			int atomic_number;
			string element;
			string adf11_file;
			// double log_temperature [x] [y]; //Need to set x, y. Dynamic sizing of arrays? Size is set from JSON
			// double log_density [x] [y]; //
			// double log_coeff [x] [y]; //
			// splines (interpolation functions on 2D grid which can be called to return value)
			// see https://en.wikipedia.org/wiki/List_of_numerical_libraries#C.2B.2B for C++ math libraries
		};
#endif