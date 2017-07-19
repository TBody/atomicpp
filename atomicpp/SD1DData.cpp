
#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
#include "SD1DData.hpp"
#include "sharedFunctions.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

SD1DData::SD1DData(const string& input_file, double impurity_fraction_input){

	json data_dict = retrieveFromJSON(input_file);

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
	
	vector<double> shaped_Ne;
	vector<double> shaped_Nn;
	vector<double> shaped_P;
	vector<double> shaped_T;

	// Extract the 'y' dimensions, at the final time index
	// Dimensions are [t, x, y, z]
	//                [0, 1, 2, 3]
	int final_time_index = extract_Ne.size() - 1;
	for(int i=0; i<extract_Ne[0][0].size(); ++i){
  		shaped_Ne.push_back(extract_Ne[final_time_index][0][i][0]);
  		shaped_Nn.push_back(extract_Nn[final_time_index][0][i][0]);
  		shaped_P.push_back(extract_P[final_time_index][0][i][0]);
		// # P = 2*Ne*Te => Te = P/(2*Ne)
  		shaped_T.push_back(shaped_P[i]/(2*shaped_Ne[i]));
	}
	// Verified list copy was identical for Ne
	
	for(int i=0; i<extract_Ne[0][0].size(); ++i){
  		density.push_back(Nnorm * shaped_Ne[i]);
  		temperature.push_back(Tnorm * shaped_T[i]);
  		neutral_fraction.push_back(shaped_Nn[i]/shaped_Ne[i]);
  		impurity_density.push_back(impurity_fraction * density[i]);
	}


};