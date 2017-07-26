#ifndef RATECOEFFICIENT_H
#define RATECOEFFICIENT_H
	#include <string>
	#include <vector>
	#include <fstream>
	#include "json.hpp"
	using json = nlohmann::json;
	
	class RateCoefficient{
		// # For storing the RateCoefficients encoded in an OpenADAS data file
		// # Intended to be called from the .makeRateCoefficients method of an ImpuritySpecies object
		// #
		// # Closely based on the cfe316/atomic/atomic_data.py/RateCoefficient class

		// # Attributes:
		// #     atomic_number (int) : The element's Z.
		// #     element (str)       : Short element name like 'c'
		// #     adf11_file (str)    : The /full/filename it came from (link to .json, not .dat)
		// #     log_temperature     : std::vector<double> of log10 of temperature values for building interpolation grid
		// #     log_density         : std::vector<double> of log10 of density values for building interpolation grid
		// #     log_coeff           : nested 3D std::vector<double> with shape (Z, temp, dens)
		// #         The list has length Z and is interpolations of log_coeff.
		public:
			/**
			 * @brief RateCoefficient constructor
			 * 
			 * @param filename JSON file from OpenADAS which supplies the rate coefficient data
			 */
			RateCoefficient(const std::string& filename);
			/**
			 * @brief blank RateCoefficient constructor
			 * @details Copies all characteristics except log_coeff
			 * 
			 * @param source_rc a smart pointer to a sample RateCoefficent object to source from
			 */
			RateCoefficient(const std::shared_ptr<RateCoefficient> source_rc);
			/**
			 * @brief Returns the rate coefficient for a (scalar) Te and Ne supplied
			 * @details Performs a simple bivariate (multilinear) interpolation to return the rate coefficient
			 * at the supplied Te and Ne values. N.b. will throw a std::runtime_error if the supplied Te or Ne
			 * value are not on the interpolating grid (otherwise you'll get a seg fault)
			 * 
			 * @param k The charge state index process (actually k+=1 for charged-target processes, but we don't implement this here)
			 * @param eval_Te electron temperature (Te) at a point (in eV)
			 * @param eval_Ne electron density (Ne) at a point (in m^-3)
			 * @return eval_coeff evaluated rate coefficient in m^3/s
			 */
			double call0D(const int k, const double eval_Te, const double eval_Ne);
			double call0DSharedInterpolation(const int k, const std::pair<int, double> Te_interp, const std::pair<int, double> Ne_interp);
			friend std::ostream& operator<<(std::ostream& os, const RateCoefficient& RC); //Define the __str__ return to std::cout
			int get_atomic_number();
			std::string get_element();
			std::string get_adf11_file();
			std::vector<std::vector< std::vector<double> > > get_log_coeff();
			std::vector<double> get_log_temperature();
			std::vector<double> get_log_density();
		private:
			int atomic_number;
			std::string element;
			std::string adf11_file;
			std::vector<std::vector< std::vector<double> > > log_coeff;
			std::vector<double> log_temperature;
			std::vector<double> log_density;
		};
#endif