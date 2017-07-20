
#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
#include "SD1DData.hpp"
#include "sharedFunctions.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

SD1DData::SD1DData(const string& expt_results_json, double impurity_fraction_input){

	json data_dict = retrieveFromJSON(expt_results_json);

	impurity_fraction = impurity_fraction_input;

	// # Retrieve (normalised values)
	vector<vector<vector< vector<double> > >> extract_Ne = data_dict["Ne"];
	vector<vector<vector< vector<double> > >> extract_Nn = data_dict["Nn"];
	vector<vector<vector< vector<double> > >> extract_P  = data_dict["P"];

	// // # Retrieve normalisation factors
	double Nnorm = data_dict["Nnorm"];
	double Tnorm = data_dict["Tnorm"];
	// // # Converts N into m^-3, T into eV 

	// Attributes
	// vector<double> density;
	// vector<double> temperature;
	// vector<double> neutral_fraction;
	// double impurity_fraction;
	// vector<double> impurity_density;
	
	// Extract the 'y' dimensions, at the final time index
	// Dimensions are [t, x, y, z]
	//                [0, 1, 2, 3]
	int final_time_index = extract_Ne.size() - 1;
	int data_length = extract_Ne[0][0].size();
	
	vector<double> shaped_Ne(data_length);
	vector<double> shaped_Nn(data_length);
	vector<double> shaped_P(data_length);
	vector<double> shaped_T(data_length);

	for(int i=0; i< data_length; ++i){
  		shaped_Ne[i] = extract_Ne[final_time_index][0][i][0];
  		shaped_Nn[i] = extract_Nn[final_time_index][0][i][0];
  		shaped_P[i] = extract_P[final_time_index][0][i][0];
		// # P = 2*Ne*Te => Te = P/(2*Ne)
  		shaped_T[i] = shaped_P[i]/(2*shaped_Ne[i]);
	}
	// Verified list copy was identical for Ne
	// N.b. the loops above and below are essentially identical 
	// ways of assigning elements into a vector. However, since the 
	// vectors below are defined at classdef their length is not
	// initialised, so using push_back method instead.
	
	for(int i=0; i<data_length; ++i){
  		density.push_back(Nnorm * shaped_Ne[i]);
  		temperature.push_back(Tnorm * shaped_T[i]);
  		neutral_fraction.push_back(shaped_Nn[i]/shaped_Ne[i]);
  		impurity_density.push_back(impurity_fraction * density[i]);
	}

};
vector<double> SD1DData::get_temperature(){
	return temperature;
};
vector<double> SD1DData::get_density(){
	return density;
};
vector<double> SD1DData::get_neutral_fraction(){
	return neutral_fraction;
};
double SD1DData::get_impurity_fraction(){
	return impurity_fraction;
};
vector<double> SD1DData::get_impurity_density(){
	return impurity_density;
};

