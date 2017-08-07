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
#include "RateEquations.hpp"

//Load global variables from ImpuritySpecies.hpp
extern const double eV_to_J; //Conversion factor between electron-volts and joules (effective units J/eV)
extern const double amu_to_kg; ////Conversion factor between atomic-mass-units and kilograms (effective units kg/amu)

RateEquations::RateEquations(ImpuritySpecies& impurity, const double Nthres_set /*= 1e9*/){
	// Set parameters that are useful for multiple functions
	//Nuclear charge of the impurity, in elementary charge units
	Z = impurity.get_atomic_number();
	//Mass of the impurity, in kilograms
	mz = impurity.get_mass();

	Nthres = Nthres_set;

	// Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
	Pcool = 0.0;
	// Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
	Prad  = 0.0;
	// Perturbation change in the electron density (in particles m^-3 s^-1) and perturbation force (in N) on the electron population due to atomic processes
	dNe = 0.0;
	F_i = 0.0;
	// Perturbation change in the neutral density (in particles m^-3 s^-1) and perturbation force (in N) on the neutral population due to atomic processes
	dNn = 0.0;
	F_n = 0.0;
};
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
std::pair<int, double> RateEquations::findSharedInterpolation(const std::vector<double>& log_grid, const double eval){
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
};
/**
 * @brief calculates the effects of electron-impact collisions on the impurity-species populations
 * @details Uses Neumaier summation to prevent floating-point rounding error when taking difference of
 * values with significantly varied magnitudes. See computeDerivs for description of parameters.
 * @param[in] Z 
 * @param[in] mz 
 * @param[in] eV_to_J 
 * @param[in] Nthres 
 * @param[in] Ne 
 * @param[in] Nzk 
 * @param[in] Vzk 
 * @param[in] Te_interp 
 * @param[in] Ne_interp 
 * @param[out] dNzk
 * @param[out] F_zk
 * @param[out] dNzk_correction
 * @param[out] F_zk_correction
 * @param[out] dNe
 * @param[out] Pcool
 */
void RateEquations::calculate_ElectronImpact_PopulationEquation(
	ImpuritySpecies& impurity,
	const int Z,
	const double mz,
	const double eV_to_J,
	const double Ne,
	const std::vector<double>& Nzk,
	const std::vector<double>& Vzk,
	const std::pair<int, double>& Te_interp,
	const std::pair<int, double>& Ne_interp,
	std::vector<double>& dNzk,
	std::vector<double>& F_zk,
	std::vector<double>& dNzk_correction,
	std::vector<double>& F_zk_correction,
	double dNe,
	double Pcool
	){
	
	// Initialise vectors as all zeros (this is default, but it doesn't hurt to be explicit)
	// These will be summed with Neumaier summation
	//Rates for population equation
		std::vector<double>    iz_to_above(Z+1, 0.0); //Rate (m^-3 s^-1)
		std::vector<double>  iz_from_below(Z+1, 0.0);
		std::vector<double>   rec_to_below(Z+1, 0.0);
		std::vector<double> rec_from_above(Z+1, 0.0);

	//Rates for momentum equation
		std::vector<double>    iz_p_to_above(Z+1, 0.0); //Momentum (kg m/s s^-1 = N)
		std::vector<double>  iz_p_from_below(Z+1, 0.0);
		std::vector<double>   rec_p_to_below(Z+1, 0.0);
		std::vector<double> rec_p_from_above(Z+1, 0.0);

	std::shared_ptr<RateCoefficient> ionisation_coefficient    = impurity.get_rate_coefficient("ionisation");
	std::shared_ptr<RateCoefficient> recombination_coefficient = impurity.get_rate_coefficient("recombination");
	
	//N.b. iterating over all data indices of the rate coefficient, hence the <
	for(int k=0; k < Z; ++k){
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
		// 	will artificially pump a low density stage if we consider only source terms
		// 	N.b. Nzk[k+1] = Nzk[k_rec], Nzk[k] = Nzk[k_rec - 1], so can test just one condition and it will work for ionisation and recombination
		if((Nzk[k+1] > Nthres) and (Nzk[k] > Nthres)){
			//Rates 'to' a state are loss terms from that state. They are paired with Rates 'from'
			//which are source terms for other states
			iz_to_above[k]            = -ionisation_rate * Nzk[k];
			iz_from_below[k+1]        = +ionisation_rate * Nzk[k];
			rec_to_below[k_rec]       = -recombination_rate * Nzk[k_rec];
			rec_from_above[k_rec-1]   = +recombination_rate * Nzk[k_rec];

			iz_p_to_above[k]          = -ionisation_rate * mz * Vzk[k];
			iz_p_from_below[k+1]      = +ionisation_rate * mz * Vzk[k+1];
			rec_p_to_below[k_rec]     = -recombination_rate * mz * Vzk[k_rec];
			rec_p_from_above[k_rec-1] = +recombination_rate * mz * Vzk[k_rec];
		}
		// Otherwise, the rates stay as default = 0
	}

	//For Neumaier summation of the change in Ne from the different rates 
	std::vector<double> dNe_from_stage(Z, 0.0);
	
	// Ionisation potential (in eV per transition) is not a rate coefficient, but is treated as one since the interpolation methods and
	// data shape are identical. This is used to calculate the effect of iz and rec on the electron energy (in Pcool).
	std::shared_ptr<RateCoefficient> ionisation_potential = impurity.get_rate_coefficient("ionisation_potential");

	for(int k=0; k < Z; ++k){
		// N.b. bare nucleus will not contribute to electron density nor have an ionisation potential -- treat it separately
		// Rates will be zero for edge cases (from initialisation)
		std::vector<double> rates_for_stage              = {iz_to_above[k],iz_from_below[k],rec_to_below[k],rec_from_above[k]};
		std::pair<double, double> neumaier_pair_rates    = neumaierSum(rates_for_stage);
		dNzk[k]                                          = neumaier_pair_rates.first;
		dNzk_correction[k]                               = neumaier_pair_rates.second; //Store compensation separately for later evaluation

		std::vector<double> momentum_for_stage           = {iz_p_to_above[k],iz_p_from_below[k],rec_p_to_below[k],rec_p_from_above[k]};
		std::pair<double, double> neumaier_pair_momentum = neumaierSum(momentum_for_stage);
		F_zk[k]                                          = neumaier_pair_momentum.first;
		F_zk_correction[k]                               = neumaier_pair_momentum.second; //Store compensation separately for later evaluation

		double ionisation_potential_evaluated            = ionisation_potential->call0D(k, Te_interp, Ne_interp);
		Pcool                                            += eV_to_J * ionisation_potential_evaluated * (iz_to_above[k] - rec_from_above[k]);
		dNe_from_stage[k]                                = (iz_to_above[k] - rec_from_above[k]);
	}

	// Perturbation on electron density
	std::pair<double, double> neumaier_pair_dNe      = neumaierSum(dNe_from_stage);
	dNe                                              = neumaier_pair_dNe.first + neumaier_pair_dNe.second;

	// Bare nucleus case
	std::vector<double> rates_for_bare_nucleus       = {iz_to_above[Z],iz_from_below[Z],rec_to_below[Z],rec_from_above[Z]};
	std::pair<double, double> neumaier_pair_dNzk     = neumaierSum(rates_for_bare_nucleus);
	dNzk[Z]                                          = neumaier_pair_dNzk.first;
	dNzk_correction[Z]                               = neumaier_pair_dNzk.second;

	std::vector<double> momentum_for_stage           = {iz_p_to_above[Z],iz_p_from_below[Z],rec_p_to_below[Z],rec_p_from_above[Z]};
	std::pair<double, double> neumaier_pair_momentum = neumaierSum(momentum_for_stage);
	F_zk[Z]                                          = neumaier_pair_momentum.first;
	F_zk_correction[Z]                               = neumaier_pair_momentum.second;
};
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
std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > RateEquations::computeDerivs(
	ImpuritySpecies& impurity,
	const double Te,
	const double Ne,
	const double Vi,
	const double Nn,
	const double Vn,
	const std::vector<double>& Nzk,
	const std::vector<double>& Vzk){

	// // Set parameters that are useful for multiple functions
	// 	//Nuclear charge of the impurity, in elementary charge units
	// 	const int Z = impurity.get_atomic_number();
	// 	//Mass of the impurity, in kilograms
	// 	const double mz = impurity.get_mass();

	// Initialise all of the output variables
		//Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
		// double Pcool = 0.0; // = std::get<0>(derivative_tuple); 

		//Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
		// double Prad  = 0.0; // = std::get<1>(derivative_tuple); 

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
		// double dNe = 0.0; // = std::get<4>(derivative_tuple);
		// double F_i = 0.0; // = std::get<5>(derivative_tuple);

		//Perturbation change in the neutral density (in particles m^-3 s^-1) and perturbation force (in N) on the neutral population due to atomic processes
		// double dNn = 0.0; // = std::get<6>(derivative_tuple);
		// double F_n = 0.0; // = std::get<7>(derivative_tuple);

	// Perform 'sharedInterpolation' - find the lower-bound gridpoint and fraction into the grid for both Te and Ne
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


	// calculate_ElectronImpact_PopulationEquation(
	// 	impurity, Z, mz, eV_to_J,
	// 	Ne, Nzk, Vzk, Te_interp, Ne_interp,
	// 	dNzk, F_zk, dNzk_correction, F_zk_correction, dNe, Pcool
	// 	);

	// Consider charge exchange after calculating Pcool
	// if (impurity.get_has_charge_exchange()){
	// 	std::vector<double> dNn_from_stage(Z, 0.0);
	// 	std::vector<double>   cx_rec_to_below(Z+1, 0.0);
	// 	std::vector<double> cx_rec_from_above(Z+1, 0.0);

	// 	std::shared_ptr<RateCoefficient> cx_recombination_coefficient = impurity.get_rate_coefficient("cx_rec");
	// 	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
	// 		double cx_recombination_coefficient_evaluated = cx_recombination_coefficient->call0D(k, Te_interp, Ne_interp);
			
	// 		// Note that cx_recombination coefficients are indexed from 1 to Z

	// 		double cx_recombination_rate = cx_recombination_coefficient_evaluated * Nn * Nzk[k+1];

	// 		// Want both the target and source densities to be above the Nthres density threshold
	// 		if((Nzk[k+1] > Nthres) and (Nzk[k] > Nthres)){
	// 			cx_rec_from_above[k] = cx_recombination_rate;
	// 			cx_rec_to_below[k+1] = -cx_recombination_rate;
	// 			dNn_from_stage[k] = -cx_recombination_rate;
	// 		}
	// 		// Otherwise, the rates stay as default = 0
	// 	}

	// 	for(int k=0; k <= Z; ++k){//Consider all states at once
	// 		std::vector<double> rates_for_stage = {dNzk[k], cx_rec_to_below[k], cx_rec_from_above[k]};
	// 		std::pair<double, double> neumaier_pair_dNzk = neumaierSum(rates_for_stage,dNzk_correction[k]); //Extend on previous compensation
	// 		dNzk[k] = neumaier_pair_dNzk.first; //Add cx to the sum
	// 		dNzk_correction[k] = neumaier_pair_dNzk.second; //Overwrite compensation with updated value
	// 	}
	// 	std::pair<double, double> neumaier_pair_dNn = neumaierSum(dNn_from_stage);
	// 	dNn = neumaier_pair_dNn.first + neumaier_pair_dNn.second;
	// }

	// Verify that the sum over all elements equals zero (or very close to)
	std::vector<double> dNzk_correctionorrected(Z+1);
	for(int k=0; k<=Z; ++k){
		 dNzk_correctionorrected[k] = dNzk[k] + dNzk_correction[k];
	}
	std::pair<double, double> neumaier_pair_total_dNzk = neumaierSum(dNzk_correctionorrected);
	// std::printf("Total sum: %e\n", neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second);
	if(abs(neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second) > 1){
		std::printf("Warning: total sum of dNzk elements is non-zero (=%e) - may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Neumaier-Neumaier summation)\n", neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second);
	}

	// Verify that the sum over all elements equals zero (or very close to)
	std::vector<double> F_zk_correctionorrected(Z+1);
	for(int k=0; k<=Z; ++k){
		 F_zk_correctionorrected[k] = F_zk[k] + F_zk_correction[k];
	}
	std::pair<double, double> neumaier_pair_total_F_zk = neumaierSum(F_zk_correctionorrected);
	// std::printf("Total sum: %e\n", neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second);
	if(abs(neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second) > 1){
		std::printf("Warning: total sum of F_zk elements is non-zero (=%e) - may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Neumaier-Neumaier summation)\n", neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second);
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
	
	//Apply neumairSum corrections
	for(int k=0; k<=Z; ++k){
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
};
/**
 * @brief Uses Neumaier algorithm to add the elements of a list
 * @details Extension on Neumaier summation algorithm for an unsorted list
 * Uses a compensated sum to improve precision when summing numbers of 
 * significantly different magnitude
 * 
 * @param list_to_sum The list of numbers to sum
 * @return neumaier_pair the uncompensated sum and the compensation
 * Compensated sum is sum + correction - however, this is left external
 * in case summation is restarted
 */
std::pair<double, double> neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction /* = 0.0*/){
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
};