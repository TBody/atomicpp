// Program name: atomic++/Prad.cpp
// Author: Thomas Body
// Author email: tajb500@york.ac.uk
// Date of creation: 17 July 2017
//
// Program function: output the radiated power (Prad)
//                   by using OpenADAS rates on output JSON from SD1D run
//
// Based on the TBody/atomic1D code, which is in turn based on the cfe316/atomic code
//

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

#include "atomicpp/ImpuritySpecies.hpp"
#include "atomicpp/RateCoefficient.hpp"
#include "atomicpp/SD1DData.hpp"
#include "atomicpp/RateEquations.hpp"

#include "atomicpp/json.hpp"
using json = nlohmann::json;

// extern const double eV_to_J; //Conversion factor between electron-volts and joules (effective units J/eV)
// extern const double amu_to_kg; ////Conversion factor between atomic-mass-units and kilograms (effective units kg/amu)

//Only used for getting guess values for the impurity ionisation-stage densities -- not required for final code
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// NOT FOR FINAL SIMULATION CODE //////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> computeIonisationDistribution(ImpuritySpecies& impurity, double Te, double Ne, double Nz, double Nn){

	int Z = impurity.get_atomic_number();
	std::vector<double> iz_stage_distribution(Z+1);

	// std::set GS density equal to 1 (arbitrary)
	iz_stage_distribution[0] = 1;
	double sum_iz = 1;

	// Loop over 0, 1, ..., Z-1
	// Each charge state is std::set in terms of the density of the previous
	for(int k=0; k<Z; ++k){
		// Ionisation
		// Get the RateCoefficient from the rate_coefficient std::map (atrribute of impurity)
		std::shared_ptr<RateCoefficient> iz_rate_coefficient = impurity.get_rate_coefficient("ionisation");
		// Evaluate the RateCoefficient at the point
		double k_iz_evaluated = iz_rate_coefficient->call0D(k, Te, Ne);

		// Recombination
		// Get the RateCoefficient from the rate_coefficient std::map (atrribute of impurity)
		std::shared_ptr<RateCoefficient> rec_rate_coefficient = impurity.get_rate_coefficient("recombination");
		// Evaluate the RateCoefficient at the point
		double k_rec_evaluated = rec_rate_coefficient->call0D(k, Te, Ne);

		// The ratio of ionisation from the (k)th stage and recombination from the (k+1)th std::sets the equilibrium densities
		// of the (k+1)th stage in terms of the (k)th (since R = Nz * Ne * rate_coefficient) N.b. Since there is no
		// ionisation from the bare nucleus, and no recombination onto the neutral (ignoring anion formation) the 'k'
		// value of ionisation coeffs is shifted down  by one relative to the recombination coeffs - therefore this
		// evaluation actually gives the balance

		iz_stage_distribution[k+1] = iz_stage_distribution[k] * (k_iz_evaluated/k_rec_evaluated);
		sum_iz += iz_stage_distribution[k+1];
	}

	// # Normalise such that the sum over all ionisation stages is '1' at all points
	for(int k=0; k<=Z; ++k){
		iz_stage_distribution[k] = iz_stage_distribution[k] / sum_iz;
	}
	return iz_stage_distribution;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// NOT FOR FINAL SIMULATION CODE //////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup the ImpuritySpecies object
	const std::string expt_results_json="sd1d-case-05.json";
	std::string impurity_symbol="c";
	std::string hydrogen_symbol="h";

	ImpuritySpecies impurity(impurity_symbol);
	ImpuritySpecies hydrogen(hydrogen_symbol);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process the expt_results_json to extract
	// 	density(s) 					= electron density (in m^-3)
	// 	temperature(s)				= electron/ion temperature (in eV)
	// 	neutral_fraction(s)			= neutral density/electron density (no units)
	// where s is 1D distance index. Time is already contracted (using final time-step)
	// Second argument is impurity fraction
	SD1DData experiment(expt_results_json, 1e-2);
	// N.b. This is only for training the data

	//Cast the SD1D data into a form which is like how the function will be called by SD1D
	//This is all not required for SD1D -- it's just to train the code with reasonable numbers
		const double mz = impurity.get_mass(); //amu

		int constant_position_index = 0;
		double Te = experiment.get_temperature()[constant_position_index];
		double Ne = experiment.get_density()[constant_position_index];
		double Vi = 3.10080599e-01 * 69205.6142373; //Picked random values from SD1D output. Using rho_s0 * Omega_ci to normalise (probably wrong!!)
		double Vn = 1.32440666e-01 * 69205.6142373; //Picked random values from SD1D output. Using rho_s0 * Omega_ci to normalise (probably wrong!!)
		double neutral_fraction = experiment.get_neutral_fraction()[constant_position_index];
		double Nn = Ne * neutral_fraction;
		double Nz = experiment.get_impurity_density()[constant_position_index];
		// Compute the iz-stage-distribution to create the Nzk (charged-resolved impurity) density std::vector
		std::vector<double> iz_stage_distribution = computeIonisationDistribution(impurity, Te, Ne, Nz, Nn);
		std::vector<double> Nzk(impurity.get_atomic_number()+1);
		for(int k=0; k<=impurity.get_atomic_number(); ++k){
			// Use the plasma temperature, and then add a scaling based on the charge so that there's different velocities for each charge
			// (so that momentum transfer results in an observable effect)
			Nzk[k] = Nz * iz_stage_distribution[k];
		}
		std::vector<double> Vzk(impurity.get_atomic_number()+1);
		for(int k=0; k<=impurity.get_atomic_number(); ++k){
			Vzk[k] = Vi;
			// sqrt(2 * Te * eV_to_J / (mz * amu_to_kg));
			// std::printf("Vz_i^(%i):  %+.2e [m/s]\n",k ,Vzk[k]);
		}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Time dependant solver code
	
	RateEquations atomic_derivatives(impurity); //Organised as a RateEquations object for cleanliness
	atomic_derivatives.setThresholdDensity(1e9); //Density threshold - ignore ionisation stages which don't have at least this density
	atomic_derivatives.setDominantIonMass(1.0); //Dominant ion mass in amu, for the stopping time calculation

	auto derivative_tuple = atomic_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);
	
	std::printf("\nDerivatives for %s\n",impurity.get_name().c_str());
	atomic_derivatives.printDerivativeTuple(derivative_tuple);
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Comparison to Post PSI
	// Ne = 1e18;
	// Te = 6;
	// Nn = 0;
	// total_power = computeRadiatedPower(impurity, Te, Ne, Nz, Nn);
	// std::cout << "Comparison to PSI paper:" << total_power << "W/m3" << std::endl;
	// std::cout << "Comparison to PSI paper:" << total_power/(Nz*Ne) << "Wm3" << std::endl;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// std::pair<int, double> Te_interp, Ne_interp;
	// Te_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_temperature(), Te);
	// Ne_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_density(), Ne);

	// int Z = impurity.get_atomic_number();
	// const double eV_to_J = 1.60217662e-19; //J per eV

	// if (impurity.get_has_charge_exchange()){
	// 	std::vector<double> dNn_from_stage(impurity.get_atomic_number(), 0.0);
	// 	std::vector<double>   cx_rec_to_below(Z+1, 0.0);
	// 	std::vector<double> cx_rec_from_above(Z+1, 0.0);

	// 	std::shared_ptr<RateCoefficient> cx_recombination_coefficient = impurity.get_rate_coefficient("cx_rec");
	// 	std::shared_ptr<RateCoefficient> cx_power_coefficient = impurity.get_rate_coefficient("cx_power");

	// 	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
	// 		// m^-3 s^-1
	// 		double cx_recombination_coefficient_evaluated = cx_recombination_coefficient->call0D(k, Te_interp, Ne_interp);
	// 		double cx_recombination_rate = cx_recombination_coefficient_evaluated * Nn * Nzk[k+1];

	// 		// W m^-3
	// 		double cx_power_coefficient_evaluated = cx_power_coefficient->call0D(k, Te_interp, Ne_interp);
	// 		double cx_power_rate = cx_power_coefficient_evaluated * Nn * Nzk[k+1] / eV_to_J;

	// 		std::printf("cx: %e [m^-3 s^-1], cx_power: %e [eV m^-3 s^-1], per transition: %e [eV] \n", cx_recombination_rate, cx_power_rate, cx_power_rate/cx_recombination_rate);
	// 	}
	// }

}





























