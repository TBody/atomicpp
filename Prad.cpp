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

#include "atomicpp/json.hpp"
using json = nlohmann::json;


/**
 * @brief Calculate the total radiated power assuming collisional-radiative equilibrium
 * @details  Assumes collisional-radiative equilibrium (i.e. infinite impurity retention time,
 * no charge-exchange recombination) to provide a simple method for calculating the ionisation stage distribution
 * for the impurity. This is then used to calculate the power due to each considered physics process (line power,
 * continuum power and charge-exchange power) for each ionisation stage. The total power (in W/m^3) is returned.
 * 
 * @param impurity ImpuritySpecies object, which contains OpenADAS data on relevant atomic-physics rate-coefficients
 * @param Te electron temperature in eV
 * @param Ne electron density in m^-3
 * @param Ni impurity density in m^-3, summed over all ionisation stages
 * @param Nn neutral density in m^-3
 * @return Total power in W/m^3
 */
double computeRadiatedPower(ImpuritySpecies& impurity, double Te, double Ne, double Ni, double Nn){
	// Calculates the relative distribution across ionisation stages of the impurity by assuming collisional-radiative
	// equilbrium. This is then used to calculate the density within each state, allowing the total power at a point
	// to be evaluated
	// std::cout << "Called for Te = " << Te << ", Ne = " << Ne << ", Ni = " << Ni << ", Nn = " << Nn << std::endl;

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

	std::set<std::string> radiative_processes = {"line_power","continuum_power"};
	if (impurity.get_has_charge_exchange()){
		radiative_processes.insert("cx_power");
	}

	double total_power = 0;

	for(int k=0; k< Z; ++k){
		double k_power = 0;
		for(std::set<std::string>::iterator iter = radiative_processes.begin();iter != radiative_processes.end();++iter){
				
			std::shared_ptr<RateCoefficient> rate_coefficient = impurity.get_rate_coefficient(*iter);
			double k_evaluated = rate_coefficient->call0D(k, Te, Ne);

			double scale;
			int target_charge_state;

			if (*iter == "line_power"){
				//# range of k is 0 to (Z-1)+ (needs bound electrons)
				target_charge_state = k; //#electron-bound target
				//# Prad = L * Ne * Nz^k+
				//#      = L * scale
				// N.b. Ne is function input
				double Nz_charge_state = Ni * iz_stage_distribution[target_charge_state];
				scale = Ne * Nz_charge_state;
			} else if (*iter == "continuum_power"){
				//# range of k is 1+ to Z+ (needs charged target)
				target_charge_state = k + 1; //#charged target
				//# Prad = L * Ne * Nz^(k+1)
				//#      = L * scale
				// N.b. Ne is function input
				double Nz_charge_state = Ni * iz_stage_distribution[target_charge_state];
				scale = Ne * Nz_charge_state;
			} else if (*iter == "cx_power"){
				//# range of k is 1+ to Z+ (needs charged target)
				target_charge_state = k + 1; //#charged target
				//# Prad = L * n_0 * Nz^(k+1)+
				//#      = L * scale
				// N.b. Nn is function input
				double Nz_charge_state = Ni * iz_stage_distribution[target_charge_state];
				scale = Nn * Nz_charge_state;
			} else {
				throw std::invalid_argument( "radiative_process not recognised (in computeRadiatedPower)" );
			}
		double power = scale * k_evaluated;
		// N.b. These won't quite give the power from the kth charge state. Instead they give the
		// power from the kth element on the rate coefficient, which may be kth or (k+1)th charge state
		// std::cout << "Power due to "<< *iter << " from k="<<k<<" is "<<power<<" [W/m3]"<<std::endl;
		k_power += power;
		}
		// N.b. These won't quite give the power from the kth charge state. Instead they give the
		// power from the kth element on the rate coefficient, which may be kth or (k+1)th charge state
		// std::cout << "Total power from all procs. from k="<<k<<" is "<<k_power<<" [W/m3]\n"<<std::endl;
		total_power += k_power;
	}

	return total_power;
}

std::vector<double> computeIonisationDistribution(ImpuritySpecies& impurity, double Te, double Ne, double Ni, double Nn){
	//Identical to main code -- used to train time-dependent solver

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
/**
 * @brief Check that shared interpolation (for speed) can be used
 * @details Checks that the log_temperature and log_density attributes of the 
 * RateCoefficients in the impurity.rate_coefficient map are identical. Also
 * adds a "blank" RateCoefficient object that doesn't have any coefficients - 
 * hopefully throws an error if you try to do something incorrectly.
 * 
 * @param impurity An ImpuritySpecies object
 */
void initialiseSharedInterpolation(ImpuritySpecies& impurity){
	// Make a blank RateCoefficient object by calling the RateCoefficient constructor on another RateCoefficient object
	//   (Choose ionisation as source since this is one of the processes always included in the model)
	//   (Might consider pushing this into seperate method and constructor, but this works for now)
	// Create a smart pointer 'RC' that points to this object
	std::shared_ptr<RateCoefficient> blank_RC(new RateCoefficient(impurity.get_rate_coefficients()["ionisation"]));
	// Add 'blank_RC' to the rate_coefficients attribute of ImpuritySpecies
	// (n.b. this is a std::map from a std::string 'physics_process' to a smart pointer which points to a RateCoefficient object)
	impurity.add_to_rate_coefficients("blank", blank_RC);

	for (auto& kv : impurity.get_rate_coefficients()) {
		std::string physics_process = kv.first;
		std::shared_ptr<RateCoefficient> RC_to_compare = kv.second;
		// Seems to implicitly compare based on a small tolerance -- works for now
		if (not(blank_RC->get_log_temperature() == RC_to_compare->get_log_temperature())){
			std::cout << "\n log_temperature doesn't match between ionisation and " << physics_process << ". Can't use shared interpolation." << std::endl;
			throw std::runtime_error("non-shared interpolation method not written - contact developer.");
		}
		if (not(blank_RC->get_log_density() == RC_to_compare->get_log_density())){
			std::cout << "\n log_density doesn't match between ionisation and " << physics_process << ". Can't use shared interpolation." << std::endl;
			throw std::runtime_error("non-shared interpolation method not written - contact developer.");
		}
	}
}
std::pair<int, double> findSharedInterpolation(const std::vector<double>& log_grid, const double eval){
	// Perform a basic interpolation based on linear distance
	// values to search for
	double eval_log10 = log10(eval);
	// Look through the grid to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int interp_gridpoint = lower_bound(log_grid.begin(), log_grid.end(), eval_log10) - log_grid.begin() - 1;

	// Bounds checking -- make sure you haven't dropped off the end of the array
	if ((interp_gridpoint == log_grid.size()-1) or (interp_gridpoint == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Te called to point off the grid for which it was defined (will give seg fault)");
	};

	int next_gridpoint = interp_gridpoint + 1;

	double grid_norm = 1/(log_grid[next_gridpoint] - log_grid[interp_gridpoint]); //Spacing between grid points

	double interp_fraction = (eval_log10 - log_grid[interp_gridpoint])*grid_norm;

	std::pair<int, double> interp_pair(interp_gridpoint, interp_fraction);
	return interp_pair;

}
/**
 * @brief Uses Neumaier algorithm to add the elements of a list
 * @details Extension on Kahan summation algorithm for an unsorted list
 * Uses a compensated sum to improve precision when summing numbers of 
 * significantly different magnitude
 * 
 * @param list_to_sum The list of numbers to sum
 * @return neumaier_pair the uncompensated sum and the compensation
 * Compensated sum is sum + correction - however, this is left external
 * in case summation is restarted
 */
std::pair<double, double> neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction = 0.0){
    double sum = 0.0;

    double correction = previous_correction;                 // A running compensation for lost low-order bits. Use previous result to restart summation

    for(int i=0; i < list_to_sum.size(); ++i){
        double temporary_sum = sum + list_to_sum[i];
        if (abs(sum) >= abs(list_to_sum[i])){
            correction += (sum - temporary_sum) + list_to_sum[i]; // If sum is bigger, low-order digits of list_to_sum[i] are lost.
        } else {
            correction += (list_to_sum[i] - temporary_sum) + sum; // Else low-order digits of sum are lost
        }
        sum = temporary_sum;
	}

	std::pair<double, double> neumaier_pair(sum, correction);

	return neumaier_pair;
}
/**
 * @brief Calculates the rate of change (input units per second) for plasma parameters due to OpenADAS atomic physics processes
 * @details Still under development
 * 
 * * @param impurity ImpuritySpecies object, which contains OpenADAS data on relevant atomic-physics rate-coefficients
 * @param Te electron temperature in eV
 * @param Ne electron density in m^-3
 * @param Nn neutral density in m^-3
 * @param Nik impurity density in m^-3, std::vector of densities of the form [Ni^0, Ni^1+, Ni^2+, ..., Ni^Z+]
 * @param Nthres threshold density for impurity stages, below which the time evolution of this stage is ignored. Default is 1e9,
 * although it is recommended that a time-step dependance be added in the calling code.
 * return dydt;
 * //where the derivative std::vector may be unpacked as
 *   double Pcool = dydt[0]; //Electron-cooling power - rate at which energy is lost from the electron population - in W/m^3
 *   double Prad  = dydt[1]; //Radiated power - rate at which energy is dissipated as radiation (for diagnostics) - in W/m^3
 *   std::vector<double> dNik(impurity.get_atomic_number()+1); //Density change for each ionisation stage of the impurity - in 1/(m^3 s)
 *   for(int k=0; k<=impurity.get_atomic_number(); ++k){
 *   	int dydt_index = k + 2;
 *   	dNik[k] = dydt[dydt_index];
 *   }
 *   double dNe   = dydt[(impurity.get_atomic_number()+2) + 1]; //Density change for electrons due to impurity-atomic processes (perturbation) - in 1/(m^3 s)
 *   double dNn   = dydt[(impurity.get_atomic_number()+2) + 2]; //Density change for neutrals due to impurity-atomic processes (perturbation) - in 1/(m^3 s)
 */
std::vector<double> computeDerivs(ImpuritySpecies& impurity, const double Te, const double Ne, const double Nn,
	const std::vector<double>& Nik, const double Nthres = 1e9){
	std::vector<double> dydt(Nik.size()+4);
	
	int Z = impurity.get_atomic_number();
	const double eV_to_J = 1.60217662e-19; //J per eV

	double Pcool = 0.0; // Pcool = dydt[0]
	double Prad = 0.0; // Prad  = dydt[1]
	std::vector<double> dNik(impurity.get_atomic_number()+1, 0.0); // dNik  = dydt[2:Z+3]
	std::vector<double> dNik_c(impurity.get_atomic_number()+1, 0.0); // corrections for neumaierSum
	double dNe  = 0.0; // dNe   = dydt[Z+3]
	double dNn  = 0.0; // dNn   = dydt[Z+3+1]

	// Find the points on the grid which correspond to Te and Ne. Since we have determined in initialiseSharedInterpolation that the grids are identical
	// we can use the same interpolated points for each
	std::pair<int, double> Te_interp, Ne_interp;
	Te_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_temperature(), Te);
	Ne_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_density(), Ne);

	// Initialise vectors as all zeros (this is default, but it doesn't hurt to be explicit)
	// These will be summed with Kahan summation
	std::vector<double>    iz_to_above(Z+1, 0.0);
	std::vector<double>  iz_from_below(Z+1, 0.0);
	std::vector<double>   rec_to_below(Z+1, 0.0);
	std::vector<double> rec_from_above(Z+1, 0.0);


	std::shared_ptr<RateCoefficient> ionisation_coefficient = impurity.get_rate_coefficient("ionisation");
	std::shared_ptr<RateCoefficient> recombination_coefficient = impurity.get_rate_coefficient("recombination");
	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
		double ionisation_coefficient_evaluated = ionisation_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
		double recombination_coefficient_evaluated = recombination_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
		
		// Note that recombination coefficients are indexed from 1 to Z, while ionisation is indexed from 0 to Z-1 (consider what 
		// charge the target must have in each case)

		double ionisation_rate = ionisation_coefficient_evaluated * Ne * Nik[k];
		// std::printf("ionisation(%i)    K = %e, Ne = %e, Nik = %e, R= %e\n",k,ionisation_coefficient_evaluated, Ne, Nik[k],ionisation_rate);
		double recombination_rate = recombination_coefficient_evaluated * Ne * Nik[k+1];
		// std::printf("recombination(%i) K = %e, Ne = %e, Nik = %e, R= %e\n",k+1,recombination_coefficient_evaluated, Ne, Nik[k+1],recombination_rate);

		// Want both the target and source densities to be above the Nthres density threshold
		// If we allow target to be below the source density, then won't get particle balance if ignoring stage or alternatively
		// will artificially pump a low density stage if we consider only source terms
		if((Nik[k+1] > Nthres) and (Nik[k] > Nthres)){
			iz_to_above[k] = -ionisation_rate;
			// std::printf("iz_to_above[%i]    : %+e\n",k,iz_to_above[k]);
			iz_from_below[k+1] = ionisation_rate;
			// std::printf("iz_from_below[%i]  : %+e\n",k+1,iz_from_below[k+1]);
			rec_from_above[k] = recombination_rate;
			// std::printf("rec_from_above[%i] : %+e\n",k,rec_from_above[k]);
			rec_to_below[k+1] = -recombination_rate;
			// std::printf("rec_to_below[%i]   : %+e\n",k+1,rec_to_below[k+1]);
		}
		// std::printf("\n");
		// Otherwise, the rates stay as default = 0
	}

	// std::printf("\n-------------------------------------------------------------\n");
	std::vector<double> dNe_from_stage(impurity.get_atomic_number(), 0.0);
	std::shared_ptr<RateCoefficient> ionisation_potential = impurity.get_rate_coefficient("ionisation_potential");

	for(int k=0; k < Z; ++k){//N.b. bare nucleus will not contribute to electron density nor have an ionisation potential -- treat it seperately
		// N.b. rates will be zero for edge cases (from initialisation)
		// std::printf("iz_to_above[%i]    : %+e\n",k,iz_to_above[k]);
		// std::printf("iz_from_below[%i]  : %+e\n",k,iz_from_below[k]);
		// std::printf("rec_to_below[%i]   : %+e\n",k,rec_to_below[k]);
		// std::printf("rec_from_above[%i] : %+e\n",k,rec_from_above[k]);
		// std::printf("\n");
		std::vector<double> rates_for_stage = {iz_to_above[k],iz_from_below[k],rec_to_below[k],rec_from_above[k]};
		std::pair<double, double> neumaier_pair_rates = neumaierSum(rates_for_stage);
		dNik[k] = neumaier_pair_rates.first;
		dNik_c[k] = neumaier_pair_rates.second; //Store compensation seperately for later evaluation

		double ionisation_potential_evaluated = ionisation_potential->call0DSharedInterpolation(k, Te_interp, Ne_interp);
		Pcool += eV_to_J * ionisation_potential_evaluated * (iz_to_above[k] - rec_from_above[k]);
		dNe_from_stage[k] = (iz_to_above[k] - rec_from_above[k]);
	}
	std::pair<double, double> neumaier_pair_dNe = neumaierSum(dNe_from_stage);
	dNe = neumaier_pair_dNe.first + neumaier_pair_dNe.second;
	// Bare nucleus case
	// std::printf("iz_to_above[%i]    : %+e\n",Z,iz_to_above[Z]);
	// std::printf("iz_from_below[%i]  : %+e\n",Z,iz_from_below[Z]);
	// std::printf("rec_to_below[%i]   : %+e\n",Z,rec_to_below[Z]);
	// std::printf("rec_from_above[%i] : %+e\n",Z,rec_from_above[Z]);
	// std::printf("\n");
	std::vector<double> rates_for_bare_nucleus = {iz_to_above[Z],iz_from_below[Z],rec_to_below[Z],rec_from_above[Z]};
	std::pair<double, double> neumaier_pair_dNik = neumaierSum(rates_for_bare_nucleus);
	dNik[Z] = neumaier_pair_dNik.first;
	dNik_c[Z] = neumaier_pair_dNik.second; //Store compensation seperately for later evaluation

	// Consider charge exchange after calculating Pcool
	if (impurity.get_has_charge_exchange()){
		std::vector<double> dNn_from_stage(impurity.get_atomic_number(), 0.0);
		std::vector<double>   cx_rec_to_below(Z+1, 0.0);
		std::vector<double> cx_rec_from_above(Z+1, 0.0);

		std::shared_ptr<RateCoefficient> cx_recombination_coefficient = impurity.get_rate_coefficient("cx_rec");
		for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
			double cx_recombination_coefficient_evaluated = cx_recombination_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
			
			// Note that cx_recombination coefficients are indexed from 1 to Z

			double cx_recombination_rate = cx_recombination_coefficient_evaluated * Nn * Nik[k+1];
			// std::printf("cx_recombination(%i) K = %e, Nn = %e, Nik = %e, R= %e\n",k+1,cx_recombination_coefficient_evaluated, Nn, Nik[k+1],cx_recombination_rate);

			// Want both the target and source densities to be above the Nthres density threshold
			if((Nik[k+1] > Nthres) and (Nik[k] > Nthres)){
				cx_rec_from_above[k] = cx_recombination_rate;
				// std::printf("cx_rec_from_above[%i] : %+e\n",k,cx_rec_from_above[k]);
				cx_rec_to_below[k+1] = -cx_recombination_rate;
				// std::printf("cx_rec_to_below[%i]   : %+e\n",k+1,cx_rec_to_below[k+1]);
				dNn_from_stage[k] = -cx_recombination_rate;
			}
			// std::printf("\n");
			// Otherwise, the rates stay as default = 0
		}

		for(int k=0; k <= Z; ++k){//Consider all states at once
			std::vector<double> rates_for_stage = {cx_rec_to_below[k],cx_rec_from_above[k]};
			// std::printf("cx_rec_to_below[%i]   : %+e\n",k,cx_rec_to_below[k]);
			// std::printf("cx_rec_from_above[%i] : %+e\n",k,cx_rec_from_above[k]);
			// std::printf("\n");
			std::pair<double, double> neumaier_pair_dNik = neumaierSum(rates_for_stage,dNik_c[k]); //Extend on previous compensation
			dNik[k] += neumaier_pair_dNik.first; //Add cx to the sum
			dNik_c[k] = neumaier_pair_dNik.second; //Overwrite compensation with updated value
		}
		std::pair<double, double> neumaier_pair_dNn = neumaierSum(dNn_from_stage);
		dNn = neumaier_pair_dNn.first + neumaier_pair_dNn.second;
	}

	// Verify that the sum over all elements equals zero (or very close to)
	std::vector<double> dNik_corrected(Z+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		 dNik_corrected[k] = dNik[k] + dNik_c[k];
		 // std::printf("dNik[%i]: %e\n", k, dNik_corrected[k]);
	}
	std::pair<double, double> neumaier_pair_total_dNik = neumaierSum(dNik_corrected);
	// std::printf("Total sum: %e\n", neumaier_pair_total_dNik.first + neumaier_pair_total_dNik.second);
	if(abs(neumaier_pair_total_dNik.first + neumaier_pair_total_dNik.second) > 1){
		std::printf("Warning: total sum of dNik elements is non-zero (=%e) - may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Kahan-Neumaier summation)\n", neumaier_pair_total_dNik.first + neumaier_pair_total_dNik.second);
	}

	//Calculate the power as well - doesn't need as high precision since everything is positive
	std::shared_ptr<RateCoefficient> line_power_coefficient = impurity.get_rate_coefficient("line_power");
	std::shared_ptr<RateCoefficient> continuum_power_coefficient = impurity.get_rate_coefficient("continuum_power");
	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
		double line_power_coefficient_evaluated = line_power_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
		double continuum_power_coefficient_evaluated = continuum_power_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
		
		// Note that continuum_power coefficients are indexed from 1 to Z, while line_power is indexed from 0 to Z-1 (consider what 
		// charge the target must have in each case)

		double line_power_rate = line_power_coefficient_evaluated * Ne * Nik[k];
		// std::printf("line_power(%i)    K = %e, Ne = %e, Nik = %e, R= %e\n",k,line_power_coefficient_evaluated, Ne, Nik[k],line_power_rate);
		double continuum_power_rate = continuum_power_coefficient_evaluated * Ne * Nik[k+1];
		// std::printf("continuum_power(%i) K = %e, Ne = %e, Nik = %e, R= %e\n",k+1,continuum_power_coefficient_evaluated, Ne, Nik[k+1],continuum_power_rate);

		Prad  += line_power_rate + continuum_power_rate;
		Pcool += line_power_rate + continuum_power_rate;
	}
	if (impurity.get_has_charge_exchange()){
		std::shared_ptr<RateCoefficient> cx_power_coefficient = impurity.get_rate_coefficient("cx_power");
		for(int k=0; k < Z; ++k){
			double cx_power_coefficient_evaluated = cx_power_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
			double cx_power_rate = cx_power_coefficient_evaluated * Nn * Nik[k+1];
			// std::printf("cx_power(%i) K = %e, Nn = %e, Nik = %e, R= %e\n",k+1,cx_power_coefficient_evaluated, Nn, Nik[k+1],cx_power_rate);
			Prad  += cx_power_rate;
		}
	}

	dydt[0] 	= Pcool;
	dydt[1] 	= Prad;
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		int dydt_index = k + 2;
		dydt[dydt_index] = dNik[k] + dNik_c[k];
	}
	dydt[impurity.get_atomic_number()+3] 	= dNe;
	dydt[impurity.get_atomic_number()+3+1] 	= dNn;

	return dydt;
}
std::vector<double> computeDerivs_old(ImpuritySpecies& impurity, const double Te, const double Ne, const double Nn, const std::vector<double>& Nik){

	int Z = impurity.get_atomic_number();
	const double eV_to_J = 1.60217662e-19; //J per eV

	double Pcool = 0; // Pcool = dydt[0]
	double Prad = 0; // Prad  = dydt[1]
	std::vector<double> dNik(impurity.get_atomic_number()+1); // dNik  = dydt[2:Z+3]
	for(int k=0; k<dNik.size(); ++k){
		dNik[k] = 0;
	}
	double dNe  = 0; // dNe   = dydt[Z+3]
	double dNn  = 0; // dNn   = dydt[Z+3+1]
	
	// Find the points on the grid which correspond to Te and Ne. Since we have determined in initialiseSharedInterpolation that the grids are identical
	// we can use the same interpolated points for each
	std::pair<int, double> Te_interp, Ne_interp;
	Te_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_temperature(), Te);
	Ne_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_density(), Ne);

	std::cout << Te << std::endl;
	std::cout << Ne << std::endl;
	std::cout << Nn << std::endl;
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		std::printf("Nz^%i: %+1.5E\n", k, Nik[k]);
		// std::cout << "Nz^" << k << ": " << Nik[k] << std::endl;
	}
	std::set<std::string> all_processes = {"ionisation","recombination","continuum_power","line_power","cx_rec","cx_power"};
	for(std::set<std::string>::iterator iter = all_processes.begin();iter != all_processes.end();++iter){
		std::shared_ptr<RateCoefficient> rate_coefficient = impurity.get_rate_coefficient(*iter);
		for(int k=0; k < Z; ++k){
			double coeff_evaluated = rate_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
			std::cout << *iter << "(" << k << ")" << " = " << coeff_evaluated << std::endl;
		}
	}

	// Select all the non-cx processes - always need these, and also they are the only ones which affect the electron population and energy
	std::set<std::string> non_cx_processes = {"ionisation","recombination","line_power","continuum_power"};
	std::shared_ptr<RateCoefficient> ionisation_potential = impurity.get_rate_coefficient("ionisation_potential");

	for(std::set<std::string>::iterator iter = non_cx_processes.begin();iter != non_cx_processes.end();++iter){
		std::shared_ptr<RateCoefficient> rate_coefficient = impurity.get_rate_coefficient(*iter);
		for(int k=0; k < Z; ++k){
			double coeff_evaluated = rate_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
			double iz_potential_evaluated = eV_to_J * ionisation_potential->call0DSharedInterpolation(k, Te_interp, Ne_interp);
			if (*iter == "ionisation"){
				int source_charge_state = k; 
				int sink_charge_state = k + 1;
				double rate = coeff_evaluated * Ne * Nik[source_charge_state]; //Returns the rate in reactions m^-3 s^-1

				
				// Compute the resulting change in the ionisation distribution
				dNik[source_charge_state] -= rate; //Loss from the source state (i.e. reactant)
				dNik[sink_charge_state] += rate; //Gain by the sink state (i.e. product)
				dNe += rate; //Electron is liberated
				Pcool -= iz_potential_evaluated * rate; //Energy loss from the electrons (iz_potential_evaluation is in J per reaction)
				// N.b. Pcool due to iz/rec is found to be very small
				// std::printf("%s (%i)        R = %+1.20E\n", iter->c_str(), k, rate);
				// std::printf("%s (%i) dNik (%i) = %+1.20E\n", iter->c_str(), k, source_charge_state, dNik[source_charge_state]);
				// std::printf("%s (%i) dNik (%i) = %+1.20E\n", iter->c_str(), k, sink_charge_state, dNik[sink_charge_state]);
				// std::printf("%s (%i)      dNe = %+1.20E\n", iter->c_str(), k, dNe);

			} else if (*iter == "recombination"){
				int source_charge_state = k + 1;
				int sink_charge_state = k; //i.e. source_charge_state - 1 
				double rate = coeff_evaluated * Ne * Nik[source_charge_state]; //Returns the rate in reactions m^-3 s^-1
				
				// Compute the resulting change in the ionisation distribution
				dNik[source_charge_state] -= rate;
				dNik[sink_charge_state] += rate;
				dNe -= rate;
				Pcool += iz_potential_evaluated * rate;

				// std::printf("%s (%i)        R = %+1.20E\n", iter->c_str(), k, rate);
				// std::printf("%s (%i) dNik (%i) = %+1.20E\n", iter->c_str(), k, source_charge_state, dNik[source_charge_state]);
				// std::printf("%s (%i) dNik (%i) = %+1.20E\n", iter->c_str(), k, sink_charge_state, dNik[sink_charge_state]);
				// std::printf("%s (%i)      dNe = %+1.20E\n", iter->c_str(), k, dNe);

			} else if (*iter == "line_power"){
				int source_charge_state = k;
				double rate = coeff_evaluated * Ne * Nik[source_charge_state]; //Returns the rate in reactions m^-3 s^-1

				// Compute the power radiated
				Prad += rate;
				Pcool += rate;

			} else if (*iter == "continuum_power"){
				int source_charge_state = k + 1;
				double rate = coeff_evaluated * Ne * Nik[source_charge_state]; //Returns the rate in reactions m^-3 s^-1

				// Compute the power radiated
				Prad += rate;
				Pcool += rate;

			} else {
				throw std::invalid_argument( "not_cx_process not recognised (in computeDerivs)" );
			}
		}
	}

	// if (impurity.get_has_charge_exchange()){
	// 	std::set<std::string> cx_processes = {"cx_rec","cx_power"};
		
	// 	for(std::set<std::string>::iterator iter = cx_processes.begin();iter != cx_processes.end();++iter){
	// 		std::shared_ptr<RateCoefficient> rate_coefficient = impurity.get_rate_coefficient(*iter);
	// 		for(int k=0; k < Z; ++k){
	// 			double coeff_evaluated = rate_coefficient->call0DSharedInterpolation(k, Te_interp, Ne_interp);
	// 			if (*iter == "cx_rec"){
	// 				int source_charge_state = k + 1;
	// 				int sink_charge_state = k; //i.e. source_charge_state - 1 
	// 				double rate = coeff_evaluated * Ne * Nik[source_charge_state]; //Returns the rate in reactions m^-3 s^-1
					
	// 				// Compute the resulting change in the ionisation distribution
	// 				dNik[source_charge_state] -= rate;
	// 				dNik[sink_charge_state] += rate;
	// 				dNn -= rate;

	// 			} else if (*iter == "cx_power"){
	// 				int source_charge_state = k;
	// 				double rate = coeff_evaluated * Ne * Nik[source_charge_state]; //Returns the rate in reactions m^-3 s^-1

	// 				// Compute the power radiated
	// 				Prad += rate;

	// 			} else {
	// 				throw std::invalid_argument( "cx_process not recognised (in computeDerivs)" );
	// 			}
	// 		}
	// 	}
	// }

	std::vector<double> dydt(Nik.size()+4);
	dydt[0] 	= Pcool;
	dydt[1] 	= Prad;
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		int dydt_index = k + 2;
		dydt[dydt_index] = dNik[k];
	}
	dydt[impurity.get_atomic_number()+3] 	= dNe;
	dydt[impurity.get_atomic_number()+3+1] 	= dNn;

	return dydt;
}


int main(){
	const std::string expt_results_json="sd1d-case-05.json";
	std::string impurity_symbol="c";

	ImpuritySpecies impurity(impurity_symbol);

	// # Process the expt_results_json to extract
	// # 	density(s) 					= electron density (in m^-3)
	// # 	temperature(s)				= electron/ion temperature (in eV)
	// # 	neutral_fraction(s)			= neutral density/electron density (no units)
	// # where s is 1D distance index. Time is already contracted (using final time-step)
	// Second argument is impurity fraction
	// 
	// N.b. This is only for training the data
	SD1DData experiment(expt_results_json, 1e-2);

	int constant_position_index = 0;

	double Te = experiment.get_temperature()[constant_position_index];
	double Ne = experiment.get_density()[constant_position_index];
	double neutral_fraction = experiment.get_neutral_fraction()[constant_position_index];
	double Nn = Ne * neutral_fraction;
	double Ni = experiment.get_impurity_density()[constant_position_index];

	// BOUT++/SD1D form for compute radiated power
	// double total_power = computeRadiatedPower(impurity, Te, Ne, Ni, Nn);

	// std::cout << "Total power from all stages is "<<total_power<<" [W/m3]\n"<<std::endl;

	// Time dependant solver code
	// Repeat the iz-stage-distribution calculation to create the Nik (charged-resolved impurity) density std::vector
	std::vector<double> iz_stage_distribution = computeIonisationDistribution(impurity, Te, Ne, Ni, Nn);
	std::vector<double> Nik(impurity.get_atomic_number()+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		Nik[k] = Ni * iz_stage_distribution[k];
		// std::cout << "k = " << k << ", fraction = " << iz_stage_distribution[k] << ", Nz^" << k << ": " << Nik[k] << std::endl;
	}

	initialiseSharedInterpolation(impurity);

	Nik[0] = 2.18178912e-01;
	Nik[1] = 1.62785542e+05;
	Nik[2] = 1.00794644e+10;
	Nik[3] = 6.77750113e+13;
	Nik[4] = 9.53000511e+16;
	Nik[5] = 7.69525043e+15;
	Nik[6] = 1.85892574e+13;
	std::vector<double> dydt = computeDerivs(impurity, 4.987422e+01, 1.030817e+19, 4.226823e+12, Nik);
	
	double Pcool = dydt[0];
	double Prad  = dydt[1];
	std::vector<double> dNik(impurity.get_atomic_number()+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		int dydt_index = k + 2;
		dNik[k] = dydt[dydt_index];
	}
	double dNe   = dydt[(impurity.get_atomic_number()+2) + 1];
	double dNn   = dydt[(impurity.get_atomic_number()+2) + 2];

	std::printf("Pcool:       %+.2e [J m^-3 s^-1]\n", Pcool);
	std::printf("Prad:        %+.2e [J m^-3 s^-1]\n" , Prad);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
	std::printf("dNz^(%i)/dt:  %+.2e [p m^-3 s^-1]\n",k ,dNik[k]);
	}
	std::printf("dNe/dt:      %+.2e [p m^-3 s^-1]\n",dNe);
	std::printf("dNn/dt:      %+.2e [p m^-3 s^-1]\n",dNn);


	// Comparison to Post PSI
	// Ne = 1e18;
	// Te = 6;
	// Nn = 0;
	// total_power = computeRadiatedPower(impurity, Te, Ne, Ni, Nn);
	// std::cout << "Comparison to PSI paper:" << total_power << "W/m3" << std::endl;
	// std::cout << "Comparison to PSI paper:" << total_power/(Ni*Ne) << "Wm3" << std::endl;

}





























