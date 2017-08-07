#ifndef RATEEEQUATIONS_H //Preprocessor directives to prevent multiple definitions
#define RATEEEQUATIONS_H
	
	// Include declarations
	#include <ostream>
	#include <map>
	#include <string>
	#include <vector>
	#include <tuple>
	#include <fstream>
	#include <set>
	#include <stdexcept> //For error-throwing
	#include <memory> //For smart pointers
	#include <cstdio> //For print formatting (printf, fprintf, sprintf, snprintf)
	#include <cmath> //For abs

	#include <limits>
	#include <cstdint>
	#include <cinttypes>

	#include "ImpuritySpecies.hpp"
	#include "RateCoefficient.hpp"

	class RateEquations{
	public:
	RateEquations(ImpuritySpecies& impurity, const double Nthres_set = 1e9);
	std::pair<int, double> findSharedInterpolation(const std::vector<double>& log_grid, const double eval);
	void calculate_ElectronImpact_PopulationEquation(
		const double Ne,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk,
		const std::pair<int, double>& Te_interp,
		const std::pair<int, double>& Ne_interp
		);
	std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > computeDerivs(
		const double Te,
		const double Ne,
		const double Vi,
		const double Nn,
		const double Vn,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk);
	private:
		//Map of RateCoefficient objects, copied from an ImpuritySpecies object
		std::map<std::string,std::shared_ptr<RateCoefficient> > rate_coefficients2;
		//If all the rate coefficients have the same log_temperature and log_density then can use the same scaling...
		//...values from a single bilinear interpolation, to save shared computation. Copied from an ImpuritySpecies object.
		bool use_shared_interpolation;
		//Whether charge exchange should be considered
		bool use_charge_exchange;
		//Nuclear charge of the impurity, in elementary charge units
		int Z; 
		//Mass of the impurity, in amu
		double mz; 
		// Threshold density for impurity stages, below which the time evolution of this stage is ignored. Default is 1e9 (constant),
		// although it is recommended that a time-step dependence be added in the calling code (overloaded call to computeDerivs).
		double Nthres;

		// Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
		double Pcool; 												// = std::get<0>(derivative_tuple);  
		// Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
		double Prad; 												// = std::get<1>(derivative_tuple);  
		
		//Change in each ionisation stage of the impurity population, in particles m^-3 s^-1
		//The index corresponds to the charge of the ionisation stage
		//	i.e. the elements are N_z^0+, N_z^1+, ... N_z^Z+ where Z is the nuclear charge
		std::vector<double> dNzk; 									// = std::get<2>(derivative_tuple); 
		//Force on each particle of ionisation stage k of the impurity population, in N
		//The index corresponds to the charge of the ionisation stage
		//	i.e. the elements are F on 0+ stage, F on 1+ stage, ..., F on Z+ stage where Z is the nuclear charge
		std::vector<double> F_zk; 									// = std::get<3>(derivative_tuple); 
		// The underscore in the name doesn't really mean anything - it's just for spacing since easier to read aligned text
		// Corrections for Kahan-Neumaier summation
		std::vector<double> dNzk_correction; 						// = std::get<2>(derivative_tuple); 
		std::vector<double> F_zk_correction; 						// = std::get<3>(derivative_tuple); 
		
		// Perturbation change in the electron density (in particles m^-3 s^-1)...
		double dNe; 												// = std::get<4>(derivative_tuple); 
		// ... and perturbation force (in N) on the electron population due to atomic processes
		double F_i; 												// = std::get<5>(derivative_tuple); 
		// Perturbation change in the neutral density (in particles m^-3 s^-1)...
		double dNn; 												// = std::get<6>(derivative_tuple); 
		// ... and perturbation force (in N) on the neutral population due to atomic processes
		double F_n; 												// = std::get<7>(derivative_tuple); 
	};

	//Public function
	std::pair<double, double> neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction = 0.0);
#endif