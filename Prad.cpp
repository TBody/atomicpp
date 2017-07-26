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

#include <memory>

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
 * @brief Calculates the rate of change (input units per second) for plasma parameters due to OpenADAS atomic physics processes
 * @details Still under development
 * 
 * * @param impurity ImpuritySpecies object, which contains OpenADAS data on relevant atomic-physics rate-coefficients
 * @param Te electron temperature in eV
 * @param Ne electron density in m^-3
 * @param Nn neutral density in m^-3
 * @param Nik impurity density in m^-3, std::vector of densities of the form [Ni^0, Ni^1+, Ni^2+, ..., Ni^Z+]
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
std::vector<double> computeDerivs(ImpuritySpecies& impurity, const double Te, const double Ne, const double Nn, const std::vector<double>& Nik){
	// std::cout << "State std::vector in" << std::endl;
	// std::cout << "Te: " << state_std::vector[0] << std::endl;
	// std::cout << "Ne: " << state_std::vector[1] << std::endl;
	// std::cout << "Nn: " << state_std::vector[2] << std::endl;
	// for(int k=0; k<=impurity.get_atomic_number(); ++k){
	// 	int state_std::vector_index = k + 3;
	// 	std::cout << "Nz^" << k << ": " << state_std::vector[state_std::vector_index] << std::endl;
	// }

	std::vector<double> dydt(Nik.size()+4);
	
	// Start with perfectly steady-state (dydt = 0 for all elements)
	for(int i=0; i<dydt.size(); ++i){
		dydt[i] = 1;
	}
	// Can replace this with derivative calculation
	int Z = impurity.get_atomic_number();
	// Derivatives of the charge states
	for(int k=0; k<=Z; ++k){
		
		if (k==0){
			// std::cout << "Nz^" << k << ": " << Nik[k] << std::endl;
		} else if (k == Z){
			// std::cout << "Nz^" << k << ": " << Nik[k] << std::endl;
		} else {
			// std::cout << "Nz^" << k << ": " << Nik[k] << std::endl;
		};


		// // Ionisation
		// // Get the RateCoefficient from the rate_coefficient std::map (atrribute of impurity)
		// std::shared_ptr<RateCoefficient> iz_rate_coefficient = impurity.get_rate_coefficient("ionisation");
		// // Evaluate the RateCoefficient at the point
		// double k_iz_evaluated = iz_rate_coefficient->call0D(k, Te, Ne);

		// // Recombination
		// // Get the RateCoefficient from the rate_coefficient std::map (atrribute of impurity)
		// std::shared_ptr<RateCoefficient> rec_rate_coefficient = impurity.get_rate_coefficient("recombination");
		// // Evaluate the RateCoefficient at the point
		// double k_rec_evaluated = rec_rate_coefficient->call0D(k, Te, Ne);
		// dydt[i] = 

		// std::cout << "Nz^" << k << ": " << state_std::vector[state_std::vector_index] << std::endl;
	}

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

	int constant_position_index = 10;

	double Te = experiment.get_temperature()[constant_position_index];
	double Ne = experiment.get_density()[constant_position_index];
	double neutral_fraction = experiment.get_neutral_fraction()[constant_position_index];
	double Nn = Ne * neutral_fraction;
	double Ni = experiment.get_impurity_density()[constant_position_index];

	// BOUT++/SD1D form for compute radiated power
	double total_power = computeRadiatedPower(impurity, Te, Ne, Ni, Nn);

	// std::cout << "Total power from all stages is "<<total_power<<" [W/m3]\n"<<std::endl;

	// Time dependant solver code
	// Repeat the iz-stage-distribution calculation to create the Nik (charged-resolved impurity) density std::vector
	std::vector<double> iz_stage_distribution = computeIonisationDistribution(impurity, Te, Ne, Ni, Nn);
	std::vector<double> Nik(impurity.get_atomic_number()+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		Nik[k] = Ni * iz_stage_distribution[k];
		// std::cout << "k = " << k << ", fraction = " << iz_stage_distribution[k] << ", Nz^" << k << ": " << Nik[k] << std::endl;
	}

	std::vector<double> dydt = computeDerivs(impurity, Te, Ne, Nn, Nik);
	double Pcool = dydt[0];
	double Prad  = dydt[1];
	std::vector<double> dNik(impurity.get_atomic_number()+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		int dydt_index = k + 2;
		dNik[k] = dydt[dydt_index];
	}
	double dNe   = dydt[(impurity.get_atomic_number()+2) + 1];
	double dNn   = dydt[(impurity.get_atomic_number()+2) + 2];

	// std::cout << "\nRate-of-change" << std::endl;
	// std::cout << "Pcool: " << Pcool << std::endl;
	// std::cout << "Prad: "  << Prad << std::endl;
	// for(int k=0; k<=impurity.get_atomic_number(); ++k){
	// 	std::cout << "dNz^" << k << "/dt: " << dNik[k] << std::endl;
	// }
	// std::cout << "dNe/dt: " << dNe << std::endl;
	// std::cout << "dNn/dt: " << dNn << std::endl;


	// Comparison to Post PSI
	// Ne = 1e18;
	// Te = 6;
	// Nn = 0;
	// total_power = computeRadiatedPower(impurity, Te, Ne, Ni, Nn);
	// std::cout << "Comparison to PSI paper:" << total_power << "W/m3" << std::endl;
	// std::cout << "Comparison to PSI paper:" << total_power/(Ni*Ne) << "Wm3" << std::endl;

	std::cout << "Running" << std::endl;
}





























