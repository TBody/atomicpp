#ifndef RATEEEQUATIONS_H //Preprocessor directives to prevent multiple definitions
#define RATEEEQUATIONS_H

	// Include declarations
	#include <vector>
	#include <tuple>
	#include <memory> //For smart pointers

	#include "ImpuritySpecies.hpp"
	#include "RateCoefficient.hpp"

	namespace atomicpp{

	struct DerivStruct{
		double Pcool;
		double Prad;
		std::vector<double> dNzk;
		std::vector<double> F_zk;
		double dNe;
		double F_i;
		double dNn;
		double F_n;
		std::vector<double> Pstage;
		double Pline;
		double Pcont;
		double Pcx;
	};

	class RateEquations{
	public:
	RateEquations(ImpuritySpecies& impurity, const double Nthres_set = 1e9, const double mi_in_amu = 1);

	// RateEquations(ImpuritySpecies impurity);

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
	 * 	double Pcool             = std::get<0>(derivative_tuple);	//Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
	 * 	double Prad              = std::get<1>(derivative_tuple);	//Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
	 * 	std::vector<double> dNzk = std::get<2>(derivative_tuple);	//Change in each ionisation stage of the impurity population, in particles m^-3 s^-1
	 * 	std::vector<double> F_zk = std::get<3>(derivative_tuple);	//Force on each particle of ionisation stage k of the impurity population, in N
	 * (returned perturbation values, not essential for modelling)
	 * 	double dNe               = std::get<4>(derivative_tuple); 	//Perturbation change in the electron density (in particles m^-3 s^-1) and
	 * 	double F_i               = std::get<5>(derivative_tuple);	//  perturbation force (in N) on the electron population due to atomic processes
	 * 	double dNn               = std::get<6>(derivative_tuple); 	//Perturbation change in the neutral density (in particles m^-3 s^-1) and
	 * 	double F_n               = std::get<7>(derivative_tuple);	// 	perturbation force (in N) on the neutral population due to atomic processes
	 */
	DerivStruct computeDerivs(
		const double Te,
		const double Ne,
		const double Vi,
		const double Nn,
		const double Vn,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk);

	/**
	 * @brief Similar to computeDerivs, except simplified for the dominant ion/neutral species (i.e. Hydrogen)
	 * @details No cx, cx_power or ion-ion drag
	 *
	 * @param Te electron temperature in eV
	 * @param Ne electron density in m^-3
	 * @param Nhk Hydrogen density in m^-3, std::vector of densities of the form [Nh^0, Nh^1+, Nh^2+, ..., Nh^Z+]
	 * @param Vhk Hydrogen velocity in m/s, std::vector of densities of the form [Vh^0, Vh^1+, Vh^2+, ..., Vh^Z+]
	 * @return Same as computeDerivs
	 */
	DerivStruct computeDerivsHydrogen(
		const double Te,
		const double Ne,
		const std::vector<double>& Nhk,
		const std::vector<double>& Vhk);

	/**
	 * @brief Set the Nthres value
	 * @param density_threshold in m^-3
	 */
	void setThresholdDensity(const double density_threshold);

	/**
	 * @brief Set the mi value
	 *
	 * @param mi_in_amu in atomic mass units
	 */
	void setDominantIonMass(const double mi_in_amu);

	/**
	 * @brief Calculates the invariant term of the stopping time
	 *
	 * @param coulomb_logarithm optional argument to set the Coulomb algorithm - default is 15
	 */
	void calculateStoppingTimeConstantFactor(const double coulomb_logarithm = 15);

	/**
	 * @brief Calculate the stopping time, except for the factor of Z^-2 since this depends on the impurity charge
	 *
	 * @param Ti in eV
	 * @param Ni in m^-3
	 */
	double calculateIonIonDragFactor(const double Ti, const double Ni);
	void resetDerivatives();

	/**
	 * @brief calculates the effects of electron-impact collisions on the impurity-species populations
	 * @details Uses Neumaier summation to prevent floating-point rounding error when taking difference of
	 * values with significantly varied magnitudes. See computeDerivs for description of parameters.
	 * Also computes the effect of ionisation and recombination on the electron cooling function (via the
	 * ionisation potential) and computes the transfer of momentum between states due to iz/rec transfer.
	 * @param[in] Ne
	 * @param[in] Nzk
	 * @param[in] Vzk
	 * @param[in] Te_interp
	 * @param[in] Ne_interp
	 */
	void calculateElectronImpactPopulationEquation(
		const double Ne,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk,
		const double Te
		);

	/**
	 * @brief calculates the effects of charge-exchange on the impurity-species populations
	 * @details Uses Neumaier summation to prevent floating-point rounding error when taking difference of
	 * values with significantly varied magnitudes. See computeDerivs for description of parameters.
	 * Also computes the transfer of momentum between states due to cx_rec transfer. N.b. cx will not affect
	 * the electron population
	 * @param[in] Nn
	 * @param[in] Nzk
	 * @param[in] Vzk
	 * @param[in] Te_interp
	 * @param[in] Ne_interp
	 */
	void calculateChargeExchangePopulationEquation(
		const double Nn,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk,
		const double Te,
		const double Ne
		);

	/**
	 * @brief Checks that the Neumaier sum for dNzk and F_zk is close to zero (no net particle source or force from transfer equations)
	 */
	void verifyNeumaierSummation(const double Te, const double Ne, const double NeumaierTolerance = 1);

	void calculateIonIonDrag(
		const double Ne,
		const double Te,
		const double Vi,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk
	);

	/**
	 * @brief Calculates the radiation rates for line power and continuum power
	 */
	void calculateElectronImpactPowerEquation(
		const double Ne,
		const std::vector<double>& Nzk,
		const double Te
	);

	/**
	 * @brief Calculates the radiation rates for cx power
	 */
	void calculateChargeExchangePowerEquation(
		const double Nn,
		const std::vector<double>& Nzk,
		const double Te,
		const double Ne
	);

	DerivStruct makeDerivativeStruct();

	void printDerivativeStruct(DerivStruct& derivative_struct);

	/**
	 * @brief makeDerivativeTuple
	 * @details packs the calculated derivatives into a tuple to return to the main code
	 * @return See commented out source code for a method to unpack return
	 */
	std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > makeDerivativeTuple();

	/**
	 * @brief print check for the returned derivative tuple
	 */
	void printDerivativeTuple(std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > derivative_tuple);

	private:
	//Map of RateCoefficient objects, copied from an ImpuritySpecies object
	std::map<std::string,std::shared_ptr<RateCoefficient> > rate_coefficients;
	//Whether charge exchange should be considered
	bool use_charge_exchange;
	//Nuclear charge of the impurity, in elementary charge units
	int Z;
	//Mass of the impurity, in amu
	double mz;
	// Threshold density for impurity stages, below which the time evolution of this stage is ignored. Default is 1e9 (constant),
	// although it is recommended that a time-step dependence be added in the calling code (overloaded call to computeDerivs).
	double Nthres;
	//Mass of the dominant ion, in amu
	double mi;
	//FF ion-ion friction-force characteristic time which is constant for the whole evaluation
	double tau_s_ii_CF;
	// FF ion-ion friction-force collision frequency which is constant for the whole evaluation = 1/tau_s_ii_CF
	double collision_frequency_s_ii_CF;

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
	std::vector<double> Pstage;
	double Pline;
	double Pcont;
	double Pcx;
	};

	//Public functions

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
	std::pair<double, double> neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction = 0.0);
	}
#endif