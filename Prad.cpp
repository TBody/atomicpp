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
// Under active development: <<TODO>> indicates development goal

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

// Look at
// http://kluge.in-chemnitz.de/opensource/spline/
// for spline interpolation

#include "atomicpp/json.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

// How to use auto keyword: http://www.acodersjourney.com/2016/02/c-11-auto/
// (with some good code examples)

// N.b. perform all plotting routines as Python post-processing

vector<double> calculateCollRadEquilibriumOD(ImpuritySpecies& impurity, SD1DData& experiment){
	// Returns the relative abundance across ionisation stages (assumes collisional radiation equilibrium)
	// 
	int constant_position_index = 0;

	int Z = impurity.get_atomic_number();

	vector<double> iz_stage_distribution(Z+1);
	
	// Do for a single point since this is what BOUT wants
	double expt_log_temperature = log10(experiment.get_temperature()[constant_position_index]);
	double expt_log_density = log10(experiment.get_density()[constant_position_index]);

	// Set GS density equal to 1 (arbitrary)
	iz_stage_distribution[0] = 1;
	double sum_iz = 1;

	for(int k=0; k<Z; ++k){

		shared_ptr<RateCoefficient> iz_rate_coefficient = impurity.get_rate_coefficient("ionisation");
		double log_k_iz_evaluated = iz_rate_coefficient->call0D(k, expt_log_temperature, expt_log_density);
		double k_iz_evaluated  = pow(10, log_k_iz_evaluated );
		shared_ptr<RateCoefficient> rec_rate_coefficient = impurity.get_rate_coefficient("recombination");
		double log_k_rec_evaluated = rec_rate_coefficient->call0D(k, expt_log_temperature, expt_log_density);
		double k_rec_evaluated = pow(10, log_k_rec_evaluated);

		// # The ratio of ionisation from the (k)th stage and recombination from the (k+1)th
		// # sets the equilibrium densities of the (k+1)th stage in terms of the (k)th (since
		// # R = n_z * n_e * rate_coefficient)
		// # N.b. Since there is no ionisation from the bare nucleus, and no recombination onto the
		// # neutral (ignoring anion formation) the 'k' value of ionisation coeffs is shifted down 
		// # by one relative to the recombination coeffs - therefore this evaluation actually gives the
		// # balance

		iz_stage_distribution[k+1] = iz_stage_distribution[k] * (k_iz_evaluated/k_rec_evaluated);
		sum_iz += iz_stage_distribution[k+1];
	}

	// # Normalise such that the sum over all ionisation stages is '1' at all points
	for(int k=0; k<=Z; ++k){
		iz_stage_distribution[k] = iz_stage_distribution[k] / sum_iz;
	}

	return iz_stage_distribution;
}
double computeRadiatedPower(ImpuritySpecies& impurity, SD1DData& experiment, vector<double>& iz_stage_distribution){
	int constant_position_index = 0;

	int Z = impurity.get_atomic_number();
	double neutral_fraction = 1e-2;
	set<string> radiative_processes = {"line_power","continuum_power"};
	if (impurity.get_has_charge_exchange()){
		radiative_processes.insert("cx_power");
	}


	// Do for a single point since this is what BOUT wants
	double expt_log_temperature = log10(experiment.get_temperature()[0]);
	double expt_log_density = log10(experiment.get_density()[0]);

	double total_power = 0;

	for(int k=0; k< Z; ++k){
		double k_power = 0;
		for(set<string>::iterator iter = radiative_processes.begin();iter != radiative_processes.end();++iter){
				
			shared_ptr<RateCoefficient> rate_coefficient = impurity.get_rate_coefficient(*iter);
			double log_k_evaluated = rate_coefficient->call0D(k, expt_log_temperature, expt_log_density);
			double k_evaluated = pow(10, log_k_evaluated);
			// cout << k_evaluated << endl;

			double scale;
			int target_charge_state;

			if (*iter == "line_power"){
				//# range of k is 0 to (Z-1)+ (needs bound electrons)
				target_charge_state = k; //#electron-bound target

				//# Prad = L * n_e * n_z^k+
				//#      = L * scale
				double n_e = experiment.get_density()[constant_position_index];
				double n_z_charge_state = experiment.get_impurity_density()[constant_position_index] * iz_stage_distribution[target_charge_state];
				scale = n_e * n_z_charge_state;
			} else if (*iter == "continuum_power"){
				//# range of k is 1+ to Z+ (needs charged target)
				target_charge_state = k + 1; //#charged target

				//# Prad = L * n_e * n_z^(k+1)
				//#      = L * scale
				double n_e = experiment.get_density()[constant_position_index];
				double n_z_charge_state = experiment.get_impurity_density()[constant_position_index] * iz_stage_distribution[target_charge_state];
				scale = n_e * n_z_charge_state;
			} else if (*iter == "cx_power"){
				//# range of k is 1+ to Z+ (needs charged target)
				target_charge_state = k + 1; //#charged target

				//# Prad = L * n_0 * n_z^(k+1)+
				//#      = L * scale
				double n_0 = experiment.get_density()[constant_position_index] * neutral_fraction;
				double n_z_charge_state = experiment.get_impurity_density()[constant_position_index] * iz_stage_distribution[target_charge_state];
				scale = n_0 * n_z_charge_state;
			} else {
				throw invalid_argument( "radiative_processe not recognised (in TBody/atomicpp/Prad.cpp/computeRadiatedPower)" );
			}
		double power = scale * k_evaluated;
		cout << "Power due to "<< *iter << " from k="<<k<<" is "<<power<<" [W/m3]"<<endl;
		k_power += power;
		}
		cout << "Total power from all procs. from k="<<k<<" is "<<k_power<<" [W/m3]\n"<<endl;
		total_power += k_power;
	}

	return total_power;
}

int main(){
	const string expt_results_json="sd1d-case-05.json";
	const string json_database_path="json_database/json_data";
	string impurity_symbol="c";

	ImpuritySpecies impurity(impurity_symbol);

	// # Process the expt_results_json to extract
	// # 	density(s) 					= electron density (in m^-3)
	// # 	temperature(s)				= electron/ion temperature (in eV)
	// # 	neutral_fraction(s)			= neutral density/electron density (no units)
	// # where s is 1D distance index. Time is already contracted (using final time-step)
	// Second argument is impurity fraction
	SD1DData experiment(expt_results_json, 1e-2);

	// # Calculate the distribution across ionisation stages, assuming collisional-radiative equilibrium
	vector<double> iz_stage_distribution = calculateCollRadEquilibriumOD(impurity, experiment);

	// for(int k=0; k<=impurity.get_atomic_number(); ++k){
	// 	cout <<"k: "<< k <<" val: " << iz_stage_distribution[k] << endl;
	// }

	// # Compute radiated power
	// # Returns total_power, stage_integrated_power (sum over all ionisation stages), and
	// # radiated_power (resolved into different ionisation stages, with 'total' giving sum over
	// # all physics_processes)
	// # 	stage_integrated_power and radiated_power are
	//  dictionaries with physics_process keys and
	// # 	data_length shape for stage_integrated_power and [Z, data_length] shape for radiated_power
	// # 	total_power is an array of shape data_length
	double total_power = computeRadiatedPower(impurity, experiment, iz_stage_distribution);

	cout << "Total power from all stages is "<<total_power<<" [W/m3]\n"<<endl;
}





























