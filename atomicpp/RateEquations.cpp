// Include declarations
#include <ostream>
#include <vector>
#include <tuple>
#include <stdexcept> //For error-throwing
#include <memory> //For smart pointers
#include <cstdio> //For print formatting (printf, fprintf, sprintf, snprintf)
#include <cmath> //For std::fabs

#include "ImpuritySpecies.hpp"
#include "RateCoefficient.hpp"
#include "RateEquations.hpp"

#include "ExternalModules/prettyprint.hpp"

using namespace atomicpp;

//Load global variables from ImpuritySpecies.hpp
extern const double eV_to_J; //Conversion factor between electron-volts and joules (effective units J/eV)
extern const double amu_to_kg; ////Conversion factor between atomic-mass-units and kilograms (effective units kg/amu)

/**
 * @brief Initialises a RateEquation object
 * @details Reads the supplied ImpuritySpecies and copies useful attributes to its own methods. Sets the Nthres and mi attributes
 * to default values unless another value is given. Calculates the tau_s_ii stopping time for the ion-ion friction term. Sets all
 * evaluation values to zero
 * 
 * @param impurity an ImpuritySpecies object for which the radiation is calculated 
 * @param density_threshold sets Nthres, a value below which the evolution of a ionisation stage is not monitored (for numerical stability)
 * @param mi_in_amu the dominant ion mass in amu (i.e. 1 for protium, 2 for deuterium)
 * @return [description]
 */
RateEquations::RateEquations(ImpuritySpecies& impurity, const double density_threshold /*= 1e9*/, const double mi_in_amu /*= 1*/) : 
	dNzk(impurity.get_atomic_number()+1, 0.0),
	F_zk(impurity.get_atomic_number()+1, 0.0),
	dNzk_correction(impurity.get_atomic_number()+1, 0.0),
	F_zk_correction(impurity.get_atomic_number()+1, 0.0),
	Pstage(impurity.get_atomic_number()+1, 0.0)
	{

	// Set parameters that are useful for multiple functions
	rate_coefficients        = impurity.get_rate_coefficients();
	use_charge_exchange      = impurity.get_has_charge_exchange();
	Z                        = impurity.get_atomic_number();
	mz                       = impurity.get_mass();

	Nthres                   = density_threshold;
	mi                       = mi_in_amu;

	//Initialise the 
	Pcool                    = 0.0;
	Prad                     = 0.0;
	dNe                      = 0.0;
	F_i                      = 0.0;
	dNn                      = 0.0;
	F_n                      = 0.0;

	Pline = 0.0;
	Pcont = 0.0;
	Pcx = 0.0;

	//Calculate the factor that doesn't change throughout evaluation
	calculateStoppingTimeConstantFactor();
};
void RateEquations::setThresholdDensity(const double density_threshold){
	Nthres = density_threshold;
};
void RateEquations::setDominantIonMass(const double mi_in_amu){
	//External set for dominant ion mass in amu
	mi = mi_in_amu;
	//Recalculate stopping time
	calculateStoppingTimeConstantFactor();
	// double tau_s_ii_pointfactor = calculateStoppingTimePointFactor(50.0, 1e19);
};
void RateEquations::calculateStoppingTimeConstantFactor(const double coulomb_logarithm){
	//Stangeby "The Plasma Boundary of Magnetic Fusion Devices" (2000), eqn 6.35
	tau_s_ii_CF = ( 1.47e13 * mz * sqrt(1/mi) ) / ( (1 + mi/mz) * coulomb_logarithm );
	//s is 'stopping'
	//ii is 'ion ion'
	//CF stands for 'Constant Factor'
	collision_frequency_s_ii_CF = 1/tau_s_ii_CF;
	//Use the collision frequency instead of the characteristic time, since multiplication is faster than division
};
double RateEquations::calculateIonIonDragFactor(const double Ti, const double Ni){
	//Stangeby "The Plasma Boundary of Magnetic Fusion Devices" (2000), eqn 6.35
	double collision_frequency_ii_PF = collision_frequency_s_ii_CF * Ni/sqrt(Ti);
	//ii is 'ion ion'
	//PF stands for 'Point Factor' -- i.e. this factor is constant for a single location point at a single time
	return collision_frequency_ii_PF;
	//collision_frequency_s_ii on Impurity^k+ = collision_frequency_ii_PF * k*k where k is the impurity change
};
DerivStruct RateEquations::computeDerivs(
	const double Te,
	const double Ne,
	const double Vi,
	const double Nn,
	const double Vn,
	const std::vector<double>& Nzk,
	const std::vector<double>& Vzk){

	resetDerivatives(); //Reset all the derivatives to zero, since the object has a memory of the previous step

	calculateElectronImpactPopulationEquation(Ne, Nzk, Vzk, Te);

	if (use_charge_exchange){
		calculateChargeExchangePopulationEquation(Nn, Nzk, Vzk, Te, Ne);
	}

	verifyNeumaierSummation(Te, Ne);

	calculateIonIonDrag(Ne, Te, Vi, Nzk, Vzk);
	
	calculateElectronImpactPowerEquation(Ne, Nzk, Te);
	
	if (use_charge_exchange){
		calculateChargeExchangePowerEquation(Nn, Nzk, Te, Ne);
	}

	// Apply neumairSum corrections
	for(int k=0; k<=Z; ++k){
		dNzk[k] += dNzk_correction[k];
		F_zk[k] += F_zk_correction[k];
	}
	
	DerivStruct derivative_struct = makeDerivativeStruct();
	return derivative_struct;
};
DerivStruct RateEquations::computeDerivsHydrogen(
	const double Te,
	const double Ne,
	const std::vector<double>& Nhk,
	const std::vector<double>& Vhk){

	resetDerivatives(); //Reset all the derivatives to zero, since the object has a memory of the previous step

	calculateElectronImpactPopulationEquation(Ne, Nhk, Vhk, Te);

	//Apply neumairSum corrections
	for(int k=0; k<=Z; ++k){
		dNzk[k] += dNzk_correction[k]; //Uses 'z' notation since these are the attributes stored in RateEquations object
		F_zk[k] += F_zk_correction[k]; //
	}

	verifyNeumaierSummation(Te, Ne);
	
	calculateElectronImpactPowerEquation(Ne, Nhk, Te);

	// auto derivative_tuple = makeDerivativeTuple();
	// return derivative_tuple;
	DerivStruct derivative_struct = makeDerivativeStruct();
	return derivative_struct;
};
void RateEquations::resetDerivatives(){
	Pcool                    = 0.0;
	Prad                     = 0.0;
	dNe                      = 0.0;
	F_i                      = 0.0;
	dNn                      = 0.0;
	F_n                      = 0.0;

	Pline = 0.0;
	Pcont = 0.0;
	Pcx = 0.0;
	
	for(int k = 0; k<= Z; ++k){
		dNzk[k] = 0.0;
		F_zk[k] = 0.0;
		dNzk_correction[k] = 0.0;
		F_zk_correction[k] = 0.0;
		Pstage[k] = 0.0;
	};
};
void RateEquations::calculateElectronImpactPopulationEquation(
	const double Ne,
	const std::vector<double>& Nzk,
	const std::vector<double>& Vzk,
	const double Te
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

	std::shared_ptr<RateCoefficient> ionisation_coefficient    = rate_coefficients["ionisation"];
	std::shared_ptr<RateCoefficient> recombination_coefficient = rate_coefficients["recombination"];
	
	//N.b. iterating over all data indices of the rate coefficient, hence the < since these are 0,...,Z-1  or 1,...,Z indexed
	for(int k=0; k < Z; ++k){
		double ionisation_coefficient_evaluated = ionisation_coefficient->call0D(k, Te, Ne);
		double recombination_coefficient_evaluated = recombination_coefficient->call0D(k, Te, Ne);
		
		// Note that recombination coefficients are indexed from 1 to Z, while ionisation is indexed from 0 to Z-1 (consider what 
		// charge the target must have in each case)
		int k_rec = k + 1; //The target charge state for recombination -- makes it a bit easier to understand

		double ionisation_rate_factor = ionisation_coefficient_evaluated * Ne;
		// std::printf("ionisation(%i)    K = %e, Ne = %e, Nzk = %e, R= %e\n",k,ionisation_coefficient_evaluated, Ne, Nzk[k],ionisation_rate_factor);
		double recombination_rate_factor = recombination_coefficient_evaluated * Ne;
		// std::printf("recombination(%i) K = %e, Ne = %e, Nzk = %e, R= %e\n",k+1,recombination_coefficient_evaluated, Ne, Nzk[k+1],recombination_rate_factor);

		// Want both the target and source densities to be above the Nthres density threshold
		// If we allow target to be below the source density, then won't get particle balance if ignoring stage or alternatively
		// 	will artificially pump a low density stage if we consider only source terms
		// 	N.b. Nzk[k+1] = Nzk[k_rec], Nzk[k] = Nzk[k_rec - 1], so can test just one condition and it will work for ionisation and recombination
		if((Nzk[k+1] > Nthres) and (Nzk[k] > Nthres)){
			//Rates 'to' a state are loss terms from that state. They are paired with Rates 'from'
			//which are source terms for other states
			iz_to_above.at(k)               = -ionisation_rate_factor * Nzk.at(k);
			iz_from_below.at(k+1)        = +ionisation_rate_factor * Nzk.at(k);
			rec_to_below.at(k_rec)       = -recombination_rate_factor * Nzk.at(k_rec);
			rec_from_above.at(k_rec-1)   = +recombination_rate_factor * Nzk.at(k_rec);

			iz_p_to_above.at(k)          = -ionisation_rate_factor * mz * Vzk.at(k);
			iz_p_from_below.at(k+1)      = +ionisation_rate_factor * mz * Vzk.at(k);
			rec_p_to_below.at(k_rec)     = -recombination_rate_factor * mz * Vzk.at(k_rec);
			rec_p_from_above.at(k_rec-1) = +recombination_rate_factor * mz * Vzk.at(k_rec);
		}
		// Otherwise, the rates stay as default = 0
	}

	for(int k=0; k <= Z; ++k){
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
	}

	//For Neumaier summation of the change in Ne from the different rates 
	std::vector<double> dNe_from_stage(Z, 0.0);
	// Ionisation potential (in eV per transition) is not a rate coefficient, but is treated as one since the interpolation methods and
	// data shape are identical. This is used to calculate the effect of iz and rec on the electron energy (in Pcool).
	std::shared_ptr<RateCoefficient> ionisation_potential = rate_coefficients["ionisation_potential"];
	for(int k=0; k < Z; ++k){
		// N.b. bare nucleus will not contribute to electron density nor have an ionisation potential
		// Rates will be zero for edge cases (from initialisation)
		double ionisation_potential_evaluated            = ionisation_potential->call0D(k, Te, Ne);
		Pcool                                            += eV_to_J * ionisation_potential_evaluated * (iz_to_above[k] - rec_from_above[k]);
		dNe_from_stage[k]                                = -(iz_to_above[k] + rec_from_above[k]); //N.b. rates already have signs
	}
	// Perturbation on electron density
	std::pair<double, double> neumaier_pair_dNe      = neumaierSum(dNe_from_stage);
	dNe                                              = neumaier_pair_dNe.first + neumaier_pair_dNe.second;
};
void RateEquations::calculateChargeExchangePopulationEquation(
	const double Nn,
	const std::vector<double>& Nzk,
	const std::vector<double>& Vzk,
	const double Te,
	const double Ne
	){
	// Consider charge exchange after calculating Pcool
	
	std::vector<double>     cx_rec_to_below(Z+1, 0.0); //Rate (m^-3 s^-1)
	std::vector<double>   cx_rec_from_above(Z+1, 0.0);
	std::vector<double>   cx_rec_p_to_below(Z+1, 0.0); //Momentum (kg m/s s^-1 = N)
	std::vector<double> cx_rec_p_from_above(Z+1, 0.0);

	std::shared_ptr<RateCoefficient> cx_recombination_coefficient = rate_coefficients["cx_rec"];
	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
		double cx_recombination_coefficient_evaluated = cx_recombination_coefficient->call0D(k, Te, Ne);
		
		// Note that cx_recombination coefficients are indexed from 1 to Z
		int k_rec = k + 1; //The target charge state for recombination -- makes it a bit easier to understand

		double cx_recombination_rate_factor = cx_recombination_coefficient_evaluated * Nn;

		// Want both the target and source densities to be above the Nthres density threshold
		if((Nzk[k_rec] > Nthres) and (Nzk[k_rec-1] > Nthres)){
			cx_rec_to_below[k_rec]     = -cx_recombination_rate_factor * Nzk[k_rec];
			cx_rec_from_above[k_rec-1] = +cx_recombination_rate_factor * Nzk[k_rec];

			cx_rec_p_to_below[k_rec]     = -cx_recombination_rate_factor * mz * Vzk[k_rec];
			cx_rec_p_from_above[k_rec-1] = +cx_recombination_rate_factor * mz * Vzk[k_rec];
		}
		// Otherwise, the rates stay as default = 0
	}

	//For Neumaier summation of the change in Nn from the different rates 
	std::vector<double> dNn_from_stage(Z+1, 0.0); //Only need Z indices, but adding an extra to the end means that the calculation
	// can be performed in a single loop
	for(int k=0; k <= Z; ++k){//Consider all states at once
		std::vector<double> rates_for_stage          = {dNzk[k], cx_rec_to_below[k], cx_rec_from_above[k]}; //Add cx to the sum
		std::pair<double, double> neumaier_pair_dNzk = neumaierSum(rates_for_stage,dNzk_correction[k]); //Extend on previous compensation
		dNzk[k]                                      = neumaier_pair_dNzk.first;
		dNzk_correction[k]                           = neumaier_pair_dNzk.second; //Overwrite compensation with updated value

		std::vector<double> momentum_for_stage           = {F_zk[k], cx_rec_p_to_below[k], cx_rec_p_from_above[k]}; //Add cx to the sum
		std::pair<double, double> neumaier_pair_momentum = neumaierSum(momentum_for_stage,F_zk_correction[k]);
		F_zk[k]                                          = neumaier_pair_momentum.first; 
		F_zk_correction[k]                               = neumaier_pair_momentum.second; //Overwrite compensation with updated value
		dNn_from_stage[k] = -cx_rec_from_above[k]; //N.b. cx_rec_from_above[Z] = 0.0 always, included so we can use a single loop
	}
	
	std::pair<double, double> neumaier_pair_dNn = neumaierSum(dNn_from_stage);
	dNn = neumaier_pair_dNn.first + neumaier_pair_dNn.second;
};
void RateEquations::verifyNeumaierSummation(const double Te, const double Ne, const double NeumaierTolerance){
	// Verify that the sum over all elements equals zero (or very close to) 
	
	std::vector<double> check_dNzk(2*(Z+1), 0.0);
	std::vector<double> check_F_zk(2*(Z+1), 0.0);
	for(int k=0; k<=Z; ++k){
		check_dNzk[2*k] = dNzk[k];
		check_dNzk[2*k+1] = dNzk_correction[k];
		check_F_zk[2*k] = F_zk[k];
		check_F_zk[2*k+1] = F_zk_correction[k];
	}
	std::pair<double, double> neumaier_pair_total_dNzk = neumaierSum(check_dNzk);
	std::pair<double, double> neumaier_pair_total_F_zk = neumaierSum(check_F_zk); 

	// std::printf("Total sum: %e\n", neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second); 
	if(std::fabs(neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second) > NeumaierTolerance){ 
		std::printf("Warning: total sum of dNzk elements is non-zero (=%e) for Te = %f, Ne = %e- may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Kahan-Neumaier summation)\n", neumaier_pair_total_dNzk.first + neumaier_pair_total_dNzk.second, Te, Ne);
	}

	// std::printf("Total sum: %e\n", neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second); 
	if(std::fabs(neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second) > NeumaierTolerance){ 
		std::printf("Warning: total sum of F_zk elements is non-zero (=%e) for Te = %f, Ne = %e- may result in error\n>>>in Prad.cpp/computeDerivs (May be an error with Kahan-Neumaier summation)\n", neumaier_pair_total_F_zk.first + neumaier_pair_total_F_zk.second, Te, Ne);
	}
};
void RateEquations::calculateIonIonDrag(
	const double Ne,
	const double Te,
	const double Vi,
	const std::vector<double>& Nzk,
	const std::vector<double>& Vzk
	){
	
	//Assume that Te = Ti (isothermal), Ne = Ni (quasi-neutrality for Z=1 dominant ion)
	double collision_frequency_ii_PF = calculateIonIonDragFactor(Te, Ne);
	//ii is 'ion ion'
	//PF is 'Point Factor' -- i.e. a constant value shared for this point

	for(int k=1; k <= Z; ++k){//Don't consider the ground state -- not an ion so must be treated separately
		if(Nzk[k] > Nthres){
			double collision_frequency_s_ii = collision_frequency_ii_PF * k*k;
			double IonIonDrag_FF = mz * (Vi - Vzk[k]) * collision_frequency_s_ii;
			F_zk[k] += IonIonDrag_FF; //Calculate the force on the impurity ion due to ion-ion drag
			F_i     -= IonIonDrag_FF; //Calculate the force on the dominant ion due to ion-ion drag
		};
	}
};
void RateEquations::calculateElectronImpactPowerEquation(
	const double Ne,
	const std::vector<double>& Nzk,
	const double Te
	){
	//Calculate the power - doesn't need as high precision since everything is positive
	std::shared_ptr<RateCoefficient> line_power_coefficient = rate_coefficients["line_power"];
	std::shared_ptr<RateCoefficient> continuum_power_coefficient = rate_coefficients["continuum_power"];
	for(int k=0; k < Z; ++k){//N.b. iterating over all data indicies of the rate coefficient, hence the <
		double line_power_coefficient_evaluated = line_power_coefficient->call0D(k, Te, Ne);
		double continuum_power_coefficient_evaluated = continuum_power_coefficient->call0D(k, Te, Ne);
		
		// Note that continuum_power coefficients are indexed from 1 to Z, while line_power is indexed from 0 to Z-1 (consider what 
		// charge the target must have in each case)

		double line_power_rate = line_power_coefficient_evaluated * Ne * Nzk[k];
		double continuum_power_rate = continuum_power_coefficient_evaluated * Ne * Nzk[k+1];

		Pstage[k] += line_power_rate + continuum_power_rate;
		Pline += line_power_rate;
		Pcont += continuum_power_rate;
		Prad  += line_power_rate + continuum_power_rate;
		Pcool += line_power_rate + continuum_power_rate;
	}
};
void RateEquations::calculateChargeExchangePowerEquation(
	const double Nn,
	const std::vector<double>& Nzk,
	const double Te,
	const double Ne
	){
	//Calculate the power - doesn't need as high precision since everything is positive
	
	std::shared_ptr<RateCoefficient> cx_power_coefficient = rate_coefficients["cx_power"];
	for(int k=0; k < Z; ++k){
		double cx_power_coefficient_evaluated = cx_power_coefficient->call0D(k, Te, Ne);
		double cx_power_rate = cx_power_coefficient_evaluated * Nn * Nzk[k+1];
		Prad  += cx_power_rate;
		Pstage[k] += cx_power_rate;
		Pcx += cx_power_rate;
	}
};
DerivStruct RateEquations::makeDerivativeStruct(){
	DerivStruct derivative_struct;

	derivative_struct.Pcool = Pcool;
	derivative_struct.Prad  = Prad;
	derivative_struct.dNzk  = dNzk;
	derivative_struct.F_zk  = F_zk;
	derivative_struct.dNe   = dNe;
	derivative_struct.F_i   = F_i;
	derivative_struct.dNn   = dNn;
	derivative_struct.F_n   = F_n;

	derivative_struct.Pstage = Pstage;
	derivative_struct.Pline = Pline;
	derivative_struct.Pcont = Pcont;
	derivative_struct.Pcx = Pcx;

	return derivative_struct;
};
void RateEquations::printDerivativeStruct(DerivStruct& derivative_struct){
	double Pcool             = derivative_struct.Pcool;
	double Prad              = derivative_struct.Prad;
	std::vector<double> dNzk = derivative_struct.dNzk;
	std::vector<double> F_zk = derivative_struct.F_zk;
	double dNe               = derivative_struct.dNe;
	double F_i               = derivative_struct.F_i;
	double dNn               = derivative_struct.dNn;
	double F_n               = derivative_struct.F_n;

	//Print-verifying the return from computeDerivs
	std::printf("Pcool:       %+.2e [J m^-3 s^-1]\n", Pcool);
	std::printf("Prad:        %+.2e [J m^-3 s^-1]\n" , Prad);
	for(int k=0; k<=Z; ++k){
	std::printf("dNz^(%i)/dt:  %+.2e [p m^-3 s^-1]\n",k ,dNzk[k]);
	}
	for(int k=0; k<=Z; ++k){
	std::printf("Fz^(%i):      %+.2e [N]\n",k ,F_zk[k]);
	}
	std::printf("dNe/dt:      %+.2e [p m^-3 s^-1]\n",dNe);
	std::printf("F_i:         %+.2e [N]\n",F_i);
	std::printf("dNn/dt:      %+.2e [p m^-3 s^-1]\n",dNn);
	std::printf("F_n:         %+.2e [N]\n",F_n);
};
std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > RateEquations::makeDerivativeTuple(){
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
	return derivative_tuple;
	//Code to unpack
		// double Pcool = std::get<0>(derivative_tuple);
		// double Prad  = std::get<1>(derivative_tuple);
		// std::vector<double> dNzk = std::get<2>(derivative_tuple);
		// std::vector<double> F_zk = std::get<3>(derivative_tuple);

		// double dNe = std::get<4>(derivative_tuple);
		// double F_i  = std::get<5>(derivative_tuple);
		// double dNn = std::get<6>(derivative_tuple);
		// double F_n  = std::get<7>(derivative_tuple);
};
void RateEquations::printDerivativeTuple(std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > derivative_tuple){
	// //Unpacking the return from computeDerivs
	double Pcool             = std::get<0>(derivative_tuple);
	double Prad              = std::get<1>(derivative_tuple);
	std::vector<double> dNzk = std::get<2>(derivative_tuple);
	std::vector<double> F_zk = std::get<3>(derivative_tuple);

	double dNe               = std::get<4>(derivative_tuple);
	double F_i               = std::get<5>(derivative_tuple);
	double dNn               = std::get<6>(derivative_tuple);
	double F_n               = std::get<7>(derivative_tuple);

	//Print-verifying the return from computeDerivs
	std::printf("Pcool:       %+.2e [J m^-3 s^-1]\n", Pcool);
	std::printf("Prad:        %+.2e [J m^-3 s^-1]\n" , Prad);
	for(int k=0; k<=Z; ++k){
	std::printf("dNz^(%i)/dt:  %+.2e [p m^-3 s^-1]\n",k ,dNzk[k]);
	}
	for(int k=0; k<=Z; ++k){
	std::printf("Fz^(%i):      %+.2e [N]\n",k ,F_zk[k]);
	}
	std::printf("dNe/dt:      %+.2e [p m^-3 s^-1]\n",dNe);
	std::printf("F_i:         %+.2e [N]\n",F_i);
	std::printf("dNn/dt:      %+.2e [p m^-3 s^-1]\n",dNn);
	std::printf("F_n:         %+.2e [N]\n",F_n);
};
std::pair<double, double> atomicpp::neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction /* = 0.0*/){
    double sum = 0.0;

    double correction = previous_correction;                 // A running compensation for lost low-order bits. Use previous result to restart summation

    for(int i=0; i < (int)(list_to_sum.size()); ++i){
        double temporary_sum = sum + list_to_sum[i];
        if (std::fabs(sum) >= std::fabs(list_to_sum[i])){
            correction += (sum - temporary_sum) + list_to_sum[i]; // If sum is bigger, low-order digits of list_to_sum[i] are lost.
        } else {
            correction += (list_to_sum[i] - temporary_sum) + sum; // Else low-order digits of sum are lost
        }
        sum = temporary_sum;
	}

	std::pair<double, double> neumaier_pair(sum, correction);

	return neumaier_pair;
};











