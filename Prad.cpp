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

#include "atomicpp/json.hpp"
using json = nlohmann::json;

//Only used for getting guess values for the impurity ionisation-stage densities -- not required for final code
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






/**
 * @brief find the lower-bound gridpoint and fraction within the grid for the given point at which to interpolate
 * @details Using bilinear interpolation, the scaling factors for interpolating the rate coefficients are the same
 * regardless of which process is called (since the underlying log_temperature and log_density grids are the same).
 * Therefore, the grid-point and fraction pair may be shared by any rate-coefficient. Have overloaded call0D such that,
 * if a <int, double> pair is supplied as an argument then the shared interpolation method will be called
 * 
 * @param log_grid grid-points for which data is given
 * @param eval point at which the interpolation should be performed
 * 
 * @return <int, double> pair where int is the lower-bound grid-point and fraction is the scaling factor (fractional distance
 * between the lower and upper-bound gridpoints)
 */
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
 * @param impurity ImpuritySpecies object, which contains OpenADAS data on relevant atomic-physics rate-coefficients
 * @param Te electron temperature in eV
 * @param Ne electron density in m^-3
 * @param Vi ion velocity in m/s
 * @param Nn neutral density in m^-3
 * @param Vn neutral velocity in m/s
 * @param Nzk impurity density in m^-3, std::vector of densities of the form [Nz^0, Nz^1+, Nz^2+, ..., Nz^Z+]
 * @param Vzk impurity velocity in m/s, std::vector of densities of the form [Vz^0, Vz^1+, Vz^2+, ..., Vz^Z+]
 * @param Nthres threshold density for impurity stages, below which the time evolution of this stage is ignored. Default is 1e9,
 * although it is recommended that a time-step dependance be added in the calling code.
 * return derivative_tuple;
 * where the return values are unpacked as
 * 	double Pcool = std::get<0>(derivative_tuple);	//Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
 * 	double Prad  = std::get<1>(derivative_tuple);	//Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
 * 	std::vector<double> dNzk = std::get<2>(derivative_tuple);	//Change in each ionisation stage of the impurity population, in particles m^-3 s^-1
 * 	std::vector<double> F_zk = std::get<3>(derivative_tuple);	//Force on each particle of ionisation stage k of the impurity population, in N
 * (returned perturbation values, not essential for modelling)
 * 	double dNe = std::get<4>(derivative_tuple); 	//Perturbation change in the electron density (in particles m^-3 s^-1) and
 * 	double F_i  = std::get<5>(derivative_tuple);	//  perturbation force (in N) on the electron population due to atomic processes
 * 	double dNn = std::get<6>(derivative_tuple); 	//Perturbation change in the neutral density (in particles m^-3 s^-1) and 
 * 	double F_n  = std::get<7>(derivative_tuple);	// 	perturbation force (in N) on the neutral population due to atomic processes
 */
std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double >
	computeDerivs(	ImpuritySpecies&impurity,
					const double Te,
					const double Ne,
					const double Vi,
					const double Nn,
					const double Vn,
					const std::vector<double>& Nzk,
					const std::vector<double>& Vzk,
					const double Nthres = 1e9){
	
	// Set parameters that are useful for multiple functions
		//Nuclear charge of the impurity, in elementary charge units
		int Z = impurity.get_atomic_number();
		//Mass of the impurity, in kilograms
		const double mz = impurity.get_mass();
		//Conversion factor between electron-volts and joules (effective units J/eV)
		const double eV_to_J = 1.60217662e-19;

	// Initialise all of the output variables
		//Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
		double Pcool = 0.0; // = std::get<0>(derivative_tuple); 

		//Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
		double Prad  = 0.0; // = std::get<1>(derivative_tuple); 

		//Change in each ionisation stage of the impurity population, in particles m^-3 s^-1
		//The index corresponds to the charge of the ionisation stage
		//	i.e. the elements are N_z^0+, N_z^1+, ... N_z^Z+ where Z is the nuclear charge
		std::vector<double> dNzk(Z+1, 0.0); // = std::get<2>(derivative_tuple);

		//Force on each particle of ionisation stage k of the impurity population, in N
		//The index corresponds to the charge of the ionisation stage
		//	i.e. the elements are F on 0+ stage, F on 1+ stage, ..., F on Z+ stage where Z is the nuclear charge
		std::vector<double> F_zk(Z+1, 0.0); // = std::get<3>(derivative_tuple);
		// The underscore in the name doesn't really mean anything - it's just for spacing since easier to read aligned text
		
		//Correction vectors for neumairSum
		std::vector<double> dNzk_correction(Z+1, 0.0); // = std::get<2>(derivative_tuple);
		std::vector<double> F_zk_correction(Z+1, 0.0); // = std::get<3>(derivative_tuple);

		//Perturbation change in the electron density (in particles m^-3 s^-1) and perturbation force (in N) on the electron population due to atomic processes
		double dNe = 0.0; // = std::get<4>(derivative_tuple);
		double F_i = 0.0; // = std::get<5>(derivative_tuple);

		//Perturbation change in the neutral density (in particles m^-3 s^-1) and perturbation force (in N) on the neutral population due to atomic processes
		double dNn = 0.0; // = std::get<6>(derivative_tuple);
		double F_n = 0.0; // = std::get<7>(derivative_tuple);

	// Perform 'sharedInterpolation'
	std::pair<int, double> Te_interp, Ne_interp;
	if (impurity.get_shared_interpolation()){
		Te_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_temperature(), Te);
		Ne_interp = findSharedInterpolation(impurity.get_rate_coefficient("blank")->get_log_density(), Ne);
	} else {
		throw std::runtime_error("Non-shared interpolation method requires switching of method. Declare Te_interp and Ne_interp as doubles instead of <int, double> pairs.");
		// //Pass Te_interp and Ne_interp as doubles instead of pairs and the program will auto-switch to the full interpolation method.
		// double Te_interp = Te;
		// double Ne_interp = Ne;
	}



	// Initialise vectors as all zeros (this is default, but it doesn't hurt to be explicit)
	// These will be summed with Kahan summation
	std::vector<double>    iz_to_above(Z+1, 0.0);
	std::vector<double>  iz_from_below(Z+1, 0.0);
	std::vector<double>   rec_to_below(Z+1, 0.0);
	std::vector<double> rec_from_above(Z+1, 0.0);

	std::vector<double>    iz_p_to_above(Z+1, 0.0); //Momentum (kg m/s s^-1 = N)
	std::vector<double>  iz_p_from_below(Z+1, 0.0); //Momentum (kg m/s s^-1 = N)
	std::vector<double>   rec_p_to_below(Z+1, 0.0); //Momentum (kg m/s s^-1 = N)
	std::vector<double> rec_p_from_above(Z+1, 0.0); //Momentum (kg m/s s^-1 = N)


	std::shared_ptr<RateCoefficient> ionisation_coefficient = impurity.get_rate_coefficient("ionisation");
	std::shared_ptr<RateCoefficient> recombination_coefficient = impurity.get_rate_coefficient("recombination");
	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
		double ionisation_coefficient_evaluated = ionisation_coefficient->call0D(k, Te_interp, Ne_interp);
		double recombination_coefficient_evaluated = recombination_coefficient->call0D(k, Te_interp, Ne_interp);
		
		// Note that recombination coefficients are indexed from 1 to Z, while ionisation is indexed from 0 to Z-1 (consider what 
		// charge the target must have in each case)
		int k_rec = k + 1; //The target charge state for recombination -- makes it a bit easier to understand

		double ionisation_rate = ionisation_coefficient_evaluated * Ne;
		// std::printf("ionisation(%i)    K = %e, Ne = %e, Nzk = %e, R= %e\n",k,ionisation_coefficient_evaluated, Ne, Nzk[k],ionisation_rate);
		double recombination_rate = recombination_coefficient_evaluated * Ne;
		// std::printf("recombination(%i) K = %e, Ne = %e, Nzk = %e, R= %e\n",k+1,recombination_coefficient_evaluated, Ne, Nzk[k+1],recombination_rate);

		// Want both the target and source densities to be above the Nthres density threshold
		// If we allow target to be below the source density, then won't get particle balance if ignoring stage or alternatively
		// will artificially pump a low density stage if we consider only source terms
		if((Nzk[k+1] > Nthres) and (Nzk[k] > Nthres)){
			iz_to_above[k] = -ionisation_rate * Nzk[k];
			iz_from_below[k+1] = ionisation_rate * Nzk[k];
			rec_from_above[k_rec-1] = recombination_rate * Nzk[k_rec];
			rec_to_below[k_rec] = -recombination_rate * Nzk[k_rec];

			iz_p_to_above[k] = -ionisation_rate * mz * Vzk[k];
			iz_p_from_below[k+1] = ionisation_rate * mz * Vzk[k+1];
			rec_p_from_above[k_rec-1] = recombination_rate * mz * Vzk[k_rec];
			rec_p_to_below[k_rec] = -recombination_rate * mz * Vzk[k_rec];
		}
		// Otherwise, the rates stay as default = 0
	}

	std::vector<double> dNe_from_stage(impurity.get_atomic_number(), 0.0);
	std::shared_ptr<RateCoefficient> ionisation_potential = impurity.get_rate_coefficient("ionisation_potential");

	for(int k=0; k < Z; ++k){//N.b. bare nucleus will not contribute to electron density nor have an ionisation potential -- treat it seperately
		// N.b. rates will be zero for edge cases (from initialisation)
		std::vector<double> rates_for_stage = {iz_to_above[k],iz_from_below[k],rec_to_below[k],rec_from_above[k]};
		std::pair<double, double> neumaier_pair_rates = neumaierSum(rates_for_stage);
		dNzk[k] = neumaier_pair_rates.first;
		dNzk_correction[k] = neumaier_pair_rates.second; //Store compensation separately for later evaluation

		std::vector<double> momentum_for_stage = {iz_p_to_above[k],iz_p_from_below[k],rec_p_to_below[k],rec_p_from_above[k]};
		std::pair<double, double> neumaier_pair_momentum = neumaierSum(momentum_for_stage);
		F_zk[k] = neumaier_pair_momentum.first;
		F_zk_correction[k] = neumaier_pair_momentum.second; //Store compensation separately for later evaluation

		double ionisation_potential_evaluated = ionisation_potential->call0D(k, Te_interp, Ne_interp);
		Pcool += eV_to_J * ionisation_potential_evaluated * (iz_to_above[k] - rec_from_above[k]);
		dNe_from_stage[k] = (iz_to_above[k] - rec_from_above[k]);
	}
	// Perturbation on electron density
	std::pair<double, double> neumaier_pair_dNe = neumaierSum(dNe_from_stage);
	dNe = neumaier_pair_dNe.first + neumaier_pair_dNe.second;
	// Bare nucleus case
	std::vector<double> rates_for_bare_nucleus = {iz_to_above[Z],iz_from_below[Z],rec_to_below[Z],rec_from_above[Z]};
	std::pair<double, double> neumaier_pair_dNzk = neumaierSum(rates_for_bare_nucleus);
	dNzk[Z] = neumaier_pair_dNzk.first;
	dNzk_correction[Z] = neumaier_pair_dNzk.second; //Store compensation separately for later evaluation

	std::vector<double> momentum_for_stage = {iz_p_to_above[Z],iz_p_from_below[Z],rec_p_to_below[Z],rec_p_from_above[Z]};
	std::pair<double, double> neumaier_pair_momentum = neumaierSum(momentum_for_stage);
	F_zk[Z] = neumaier_pair_momentum.first;
	F_zk_correction[Z] = neumaier_pair_momentum.second; //Store compensation separately for later evaluation

	// Consider charge exchange after calculating Pcool
	if (impurity.get_has_charge_exchange()){
		std::vector<double> dNn_from_stage(impurity.get_atomic_number(), 0.0);
		std::vector<double>   cx_rec_to_below(Z+1, 0.0);
		std::vector<double> cx_rec_from_above(Z+1, 0.0);

		std::shared_ptr<RateCoefficient> cx_recombination_coefficient = impurity.get_rate_coefficient("cx_rec");
		for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
			double cx_recombination_coefficient_evaluated = cx_recombination_coefficient->call0D(k, Te_interp, Ne_interp);
			
			// Note that cx_recombination coefficients are indexed from 1 to Z

			double cx_recombination_rate = cx_recombination_coefficient_evaluated * Nn * Nzk[k+1];

			// Want both the target and source densities to be above the Nthres density threshold
			if((Nzk[k+1] > Nthres) and (Nzk[k] > Nthres)){
				cx_rec_from_above[k] = cx_recombination_rate;
				cx_rec_to_below[k+1] = -cx_recombination_rate;
				dNn_from_stage[k] = -cx_recombination_rate;
			}
			// Otherwise, the rates stay as default = 0
		}

		for(int k=0; k <= Z; ++k){//Consider all states at once
			std::vector<double> rates_for_stage = {dNzk[k], cx_rec_to_below[k], cx_rec_from_above[k]};
			std::pair<double, double> neumaier_pair_dNzk = neumaierSum(rates_for_stage,dNzk_correction[k]); //Extend on previous compensation
			dNzk[k] = neumaier_pair_dNzk.first; //Add cx to the sum
			dNzk_correction[k] = neumaier_pair_dNzk.second; //Overwrite compensation with updated value
		}
		std::pair<double, double> neumaier_pair_dNn = neumaierSum(dNn_from_stage);
		dNn = neumaier_pair_dNn.first + neumaier_pair_dNn.second;
	}

	// Verify that the sum over all elements equals zero (or very close to)
	std::vector<double> dNzk_correctionorrected(Z+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		 dNzk_correctionorrected[k] = dNzk[k] + dNzk_correction[k];
	}
	std::pair<double, double> neumaier_pair_total_dNzk = neumaierSum(dNzk_correctionorrected);
	// std::printf("Total sum: %e\n", neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second);
	if(abs(neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second) > 1){
		std::printf("Warning: total sum of dNzk elements is non-zero (=%e) - may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Kahan-Neumaier summation)\n", neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second);
	}

	// Verify that the sum over all elements equals zero (or very close to)
	std::vector<double> F_zk_correctionorrected(Z+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		 F_zk_correctionorrected[k] = F_zk[k] + F_zk_correction[k];
	}
	std::pair<double, double> neumaier_pair_total_F_zk = neumaierSum(F_zk_correctionorrected);
	// std::printf("Total sum: %e\n", neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second);
	if(abs(neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second) > 1){
		std::printf("Warning: total sum of F_zk elements is non-zero (=%e) - may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Kahan-Neumaier summation)\n", neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second);
	}

	//Calculate the power as well - doesn't need as high precision since everything is positive
	std::shared_ptr<RateCoefficient> line_power_coefficient = impurity.get_rate_coefficient("line_power");
	std::shared_ptr<RateCoefficient> continuum_power_coefficient = impurity.get_rate_coefficient("continuum_power");
	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
		double line_power_coefficient_evaluated = line_power_coefficient->call0D(k, Te_interp, Ne_interp);
		double continuum_power_coefficient_evaluated = continuum_power_coefficient->call0D(k, Te_interp, Ne_interp);
		
		// Note that continuum_power coefficients are indexed from 1 to Z, while line_power is indexed from 0 to Z-1 (consider what 
		// charge the target must have in each case)

		double line_power_rate = line_power_coefficient_evaluated * Ne * Nzk[k];
		double continuum_power_rate = continuum_power_coefficient_evaluated * Ne * Nzk[k+1];

		Prad  += line_power_rate + continuum_power_rate;
		Pcool += line_power_rate + continuum_power_rate;
	}
	if (impurity.get_has_charge_exchange()){
		std::shared_ptr<RateCoefficient> cx_power_coefficient = impurity.get_rate_coefficient("cx_power");
		for(int k=0; k < Z; ++k){
			double cx_power_coefficient_evaluated = cx_power_coefficient->call0D(k, Te_interp, Ne_interp);
			double cx_power_rate = cx_power_coefficient_evaluated * Nn * Nzk[k+1];
			Prad  += cx_power_rate;
		}
	}

	// dydt[0] 	= Pcool;
	// dydt[1] 	= Prad;
	// dydt[impurity.get_atomic_number()+3] 	= dNe;
	// dydt[impurity.get_atomic_number()+3+1] 	= dNn;
	
	//Apply neumairSum corrections
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		dNzk[k] += dNzk_correction[k];
		F_zk[k] += F_zk_correction[k];
	}


	std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double >derivative_tuple = 
	std::make_tuple(
		Pcool,
		Prad,
		dNzk,
		F_zk,
		dNe,
		F_i,
		dNn,
		F_n
	);


	// double Pcool = std::get<0>(derivative_tuple);
	// double Prad  = std::get<1>(derivative_tuple);
	// std::vector<double> dNzk = std::get<2>(derivative_tuple);
	// std::vector<double> F_zk = std::get<3>(derivative_tuple);

	// double dNe = std::get<4>(derivative_tuple);
	// double F_i  = std::get<5>(derivative_tuple);
	// double dNn = std::get<6>(derivative_tuple);
	// double F_n  = std::get<7>(derivative_tuple);

	return derivative_tuple;
}

int main(){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup the ImpuritySpecies object
	const std::string expt_results_json="sd1d-case-05.json";
	std::string impurity_symbol="c";

	ImpuritySpecies impurity(impurity_symbol);

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
	const int Z = impurity.get_atomic_number();
	const double eV_to_J = 1.60217662e-19; //J per eV
	const double mz = impurity.get_mass(); //Kilograms

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
		Nzk[k] = Nz * iz_stage_distribution[k] * k/Z;
	}
	std::vector<double> Vzk(impurity.get_atomic_number()+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		Vzk[k] = sqrt(2*Te*eV_to_J/mz) * k/Z;
		// std::printf("Vz_i^(%i):  %+.2e [m/s]\n",k ,Vzk[k]);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Time dependant solver code
	double Nthres = 1e9; //Density threshold - ignore ionisation stages which don't have at least this density
	auto derivative_tuple = computeDerivs(impurity, Te, Ne, Vi, Nn, Vn, Nzk, Vzk, Nthres);
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Unpacking the return from computeDerivs
	double Pcool = std::get<0>(derivative_tuple);
	double Prad  = std::get<1>(derivative_tuple);
	std::vector<double> dNzk = std::get<2>(derivative_tuple);
	std::vector<double> F_zk = std::get<3>(derivative_tuple);

	double dNe = std::get<4>(derivative_tuple);
	double F_i  = std::get<5>(derivative_tuple);
	double dNn = std::get<6>(derivative_tuple);
	double F_n  = std::get<7>(derivative_tuple);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Print-verifying the return from computeDerivs
	std::printf("Pcool:       %+.2e [J m^-3 s^-1]\n", Pcool);
	std::printf("Prad:        %+.2e [J m^-3 s^-1]\n" , Prad);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
	std::printf("dNz^(%i)/dt:  %+.2e [p m^-3 s^-1]\n",k ,dNzk[k]);
	}
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
	std::printf("Fz^(%i):      %+.2e [N]\n",k ,F_zk[k]);
	}
	std::printf("dNe/dt:      %+.2e [p m^-3 s^-1]\n",dNe);
	std::printf("F_i/dt:      %+.2e [N]\n",F_i);
	std::printf("dNn/dt:      %+.2e [p m^-3 s^-1]\n",dNn);
	std::printf("F_n/dt:      %+.2e [N]\n",F_n);

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





























