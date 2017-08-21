
#include <string>
#include <vector>
#include <fstream>
#include "SD1DData.hpp"
#include "sharedFunctions.hpp"

#include "ExternalModules/json.hpp"
using json = nlohmann::json;

using namespace atomicpp;

SD1DData::SD1DData(const std::string& expt_results_json, double impurity_fraction_input){

	json data_dict = retrieveFromJSON(expt_results_json);

	impurity_fraction = impurity_fraction_input;

	// # Retrieve (normalised values)
	std::vector<std::vector<std::vector< std::vector<double> > >> extract_Ne = data_dict["Ne"];
	std::vector<std::vector<std::vector< std::vector<double> > >> extract_Nn = data_dict["Nn"];
	std::vector<std::vector<std::vector< std::vector<double> > >> extract_P  = data_dict["P"];

	// // # Retrieve normalisation factors
	double Nnorm = data_dict["Nnorm"];
	double Tnorm = data_dict["Tnorm"];
	// // # Converts N into m^-3, T into eV 

	// Attributes
	// std::vector<double> density;
	// std::vector<double> temperature;
	// std::vector<double> neutral_fraction;
	// double impurity_fraction;
	// std::vector<double> impurity_density;
	
	// Extract the 'y' dimensions, at the final time index
	// Dimensions are [t, x, y, z]
	//                [0, 1, 2, 3]
	int final_time_index = extract_Ne.size() - 1;
	int data_length = extract_Ne[0][0].size();
	
	std::vector<double> shaped_Ne(data_length);
	std::vector<double> shaped_Nn(data_length);
	std::vector<double> shaped_P(data_length);
	std::vector<double> shaped_T(data_length);

	for(int i=0; i< data_length; ++i){
  		shaped_Ne[i] = extract_Ne[final_time_index][0][i][0];
  		shaped_Nn[i] = extract_Nn[final_time_index][0][i][0];
  		shaped_P[i] = extract_P[final_time_index][0][i][0];
		// # P = 2*Ne*Te => Te = P/(2*Ne)
  		shaped_T[i] = shaped_P[i]/(2*shaped_Ne[i]);
	}
	// Verified list copy was identical for Ne
	// N.b. the loops above and below are essentially identical 
	// ways of assigning elements into a std::vector. However, since the 
	// std::vectors below are defined at classdef their length is not
	// initialised, so using push_back method instead.
	
	for(int i=0; i<data_length; ++i){
  		density.push_back(Nnorm * shaped_Ne[i]);
  		temperature.push_back(Tnorm * shaped_T[i]);
  		neutral_fraction.push_back(shaped_Nn[i]/shaped_Ne[i]);
  		impurity_density.push_back(impurity_fraction * density[i]);
	}

};
std::vector<double> SD1DData::get_temperature(){
	return temperature;
};
std::vector<double> SD1DData::get_density(){
	return density;
};
std::vector<double> SD1DData::get_neutral_fraction(){
	return neutral_fraction;
};
double SD1DData::get_impurity_fraction(){
	return impurity_fraction;
};
std::vector<double> SD1DData::get_impurity_density(){
	return impurity_density;
};

