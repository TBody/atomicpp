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
#include <iostream>
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
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

double computeRadiatedPower(ImpuritySpecies& impurity, double Te, double Ne, double Ni, double Nn){
	// Calculates the relative distribution across ionisation stages of the impurity by assuming collisional-radiative
	// equilbrium. This is then used to calculate the density within each state, allowing the total power at a point
	// to be evaluated
	cout << "Called for Te = " << Te << ", Ne = " << Ne << ", Ni = " << Ni << ", Nn = " << Nn << endl;

	int Z = impurity.get_atomic_number();
	vector<double> iz_stage_distribution(Z+1);

	// Set GS density equal to 1 (arbitrary)
	iz_stage_distribution[0] = 1;
	double sum_iz = 1;

	// Loop over 0, 1, ..., Z-1
	// Each charge state is set in terms of the density of the previous
	for(int k=0; k<Z; ++k){
		// Ionisation
		// Get the RateCoefficient from the rate_coefficient map (atrribute of impurity)
		shared_ptr<RateCoefficient> iz_rate_coefficient = impurity.get_rate_coefficient("ionisation");
		// Evaluate the RateCoefficient at the point
		double k_iz_evaluated = iz_rate_coefficient->call0D(k, Te, Ne);

		// Recombination
		// Get the RateCoefficient from the rate_coefficient map (atrribute of impurity)
		shared_ptr<RateCoefficient> rec_rate_coefficient = impurity.get_rate_coefficient("recombination");
		// Evaluate the RateCoefficient at the point
		double k_rec_evaluated = rec_rate_coefficient->call0D(k, Te, Ne);

		// The ratio of ionisation from the (k)th stage and recombination from the (k+1)th sets the equilibrium densities
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

	set<string> radiative_processes = {"line_power","continuum_power"};
	if (impurity.get_has_charge_exchange()){
		radiative_processes.insert("cx_power");
	}

	double total_power = 0;

	for(int k=0; k< Z; ++k){
		double k_power = 0;
		for(set<string>::iterator iter = radiative_processes.begin();iter != radiative_processes.end();++iter){
				
			shared_ptr<RateCoefficient> rate_coefficient = impurity.get_rate_coefficient(*iter);
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
				throw invalid_argument( "radiative_process not recognised (in computeRadiatedPower)" );
			}
		double power = scale * k_evaluated;
		// N.b. These won't quite give the power from the kth charge state. Instead they give the
		// power from the kth element on the rate coefficient, which may be kth or (k+1)th charge state
		// cout << "Power due to "<< *iter << " from k="<<k<<" is "<<power<<" [W/m3]"<<endl;
		k_power += power;
		}
		// N.b. These won't quite give the power from the kth charge state. Instead they give the
		// power from the kth element on the rate coefficient, which may be kth or (k+1)th charge state
		// cout << "Total power from all procs. from k="<<k<<" is "<<k_power<<" [W/m3]\n"<<endl;
		total_power += k_power;
	}

	return total_power;
}

vector<double> computeIonisationDistribution(ImpuritySpecies& impurity, double Te, double Ne, double Ni, double Nn){
	//Identical to main code -- used to train time-dependent solver

	int Z = impurity.get_atomic_number();
	vector<double> iz_stage_distribution(Z+1);

	// Set GS density equal to 1 (arbitrary)
	iz_stage_distribution[0] = 1;
	double sum_iz = 1;

	// Loop over 0, 1, ..., Z-1
	// Each charge state is set in terms of the density of the previous
	for(int k=0; k<Z; ++k){
		// Ionisation
		// Get the RateCoefficient from the rate_coefficient map (atrribute of impurity)
		shared_ptr<RateCoefficient> iz_rate_coefficient = impurity.get_rate_coefficient("ionisation");
		// Evaluate the RateCoefficient at the point
		double k_iz_evaluated = iz_rate_coefficient->call0D(k, Te, Ne);

		// Recombination
		// Get the RateCoefficient from the rate_coefficient map (atrribute of impurity)
		shared_ptr<RateCoefficient> rec_rate_coefficient = impurity.get_rate_coefficient("recombination");
		// Evaluate the RateCoefficient at the point
		double k_rec_evaluated = rec_rate_coefficient->call0D(k, Te, Ne);

		// The ratio of ionisation from the (k)th stage and recombination from the (k+1)th sets the equilibrium densities
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

vector<double> computeStateVectorDerivs(ImpuritySpecies& impurity, vector<double>& state_vector){
	cout << "State vector in" << endl;
	cout << "Te: " << state_vector[0] << endl;
	cout << "Ne: " << state_vector[1] << endl;
	cout << "Nn: " << state_vector[2] << endl;
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		int state_vector_index = k + 3;
		cout << "Nz^" << k << ": " << state_vector[state_vector_index] << endl;
	}

	// Start with perfectly steady-state
	vector<double> dydt(state_vector.size());
	for(int i=0; i<state_vector.size(); ++i){
		dydt[i] = 0;
	}
	// Can replace this with derivative calculation


	cout << "\nRate-of-change of state vector out" << endl;
	cout << "dTe/dt: " << dydt[0] << endl;
	cout << "dNe/dt: " << dydt[1] << endl;
	cout << "dNn/dt: " << dydt[2] << endl;
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		int dydt_index = k + 3;
		cout << "dNz^" << k << "/dt: " << dydt[dydt_index] << endl;
	}

	return dydt;
}

int main(){
	const string expt_results_json="sd1d-case-05.json";
	string impurity_symbol="c";

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

	int constant_position_index = 199;

	double Te = experiment.get_temperature()[constant_position_index];
	double Ne = experiment.get_density()[constant_position_index];
	double neutral_fraction = experiment.get_neutral_fraction()[constant_position_index];
	double Nn = Ne * neutral_fraction;
	double Ni = experiment.get_impurity_density()[constant_position_index];

	// BOUT++/SD1D form for compute radiated power
	double total_power = computeRadiatedPower(impurity, Te, Ne, Ni, Nn);

	cout << "Total power from all stages is "<<total_power<<" [W/m3]\n"<<endl;

	// Time dependant solver code
	// Repeat the iz-stage-distribution calculation to create the Nik (charged-resolved impurity) density vector
	vector<double> iz_stage_distribution = computeIonisationDistribution(impurity, Te, Ne, Ni, Nn);
	vector<double> Nik(impurity.get_atomic_number()+1);
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		Nik[k] = Ni * iz_stage_distribution[k];
		// cout << "k = " << k << ", fraction = " << iz_stage_distribution[k] << ", Nz^" << k << ": " << Nik[k] << endl;
	}

	// Create the 'state vector' for the differential equation
	vector<double> state_vector = {Te, Ne, Nn};
	for(int k=0; k<=impurity.get_atomic_number(); ++k){
		state_vector.push_back(Nik[k]);
	}
	cout << "state_vector with elements (Te, Ne, Nn, Nz^0, Nz^1, ..., Nz^k) initialised" << endl;
	cout << "state_vector length: " << state_vector.size() << ", monitoring " << impurity.get_atomic_number() + 1 << " charge states of " << impurity.get_name() << " (including g.s.)" << endl;

	vector<double> dydt = computeStateVectorDerivs(impurity, state_vector);
}





























