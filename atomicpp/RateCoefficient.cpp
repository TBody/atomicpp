#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
using json = nlohmann::json;
#include <stdexcept> //For error-throwing

#include "RateCoefficient.hpp"
#include "sharedFunctions.hpp"

#include <algorithm> //for upper/lower_bound

using namespace atomicpp;
RateCoefficient::RateCoefficient(const std::string& filename){
	// # Create an instance of RateCoefficient by reading an OpenADAS JSON file

	json data_dict = retrieveFromJSON(filename);

	atomic_number   = data_dict["charge"];
	element         = data_dict["element"];
	adf11_file      = filename;

	std::vector<std::vector< std::vector<double> > > extract_log_coeff = data_dict["log_coeff"];
	std::vector<double> extract_log_temperature = data_dict["log_temperature"];
	std::vector<double> extract_log_density = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute
	log_coeff = extract_log_coeff;
	log_temperature = extract_log_temperature;
	log_density = extract_log_density;

	// std::vector<std::vector<std::vector<grid_matrix>>> alpha_coeff_init = calculate_alpha_coeff(log_temperature, log_density, log_coeff);
	// alpha_coeff = alpha_coeff_init;

};
RateCoefficient::RateCoefficient(const std::shared_ptr<RateCoefficient> source_rc){
	// # Create an instance of a blank RateCoefficient by copying from another RateCoefficient object

	atomic_number   = source_rc->get_atomic_number();
	element         = source_rc->get_element();
	adf11_file      = source_rc->get_adf11_file();

	log_temperature = source_rc->get_log_temperature();
	log_density = source_rc->get_log_density();

};
double RateCoefficient::call0D_bicubic(const int k, const double eval_Te, const double eval_Ne){
	// """Evaluate the ionisation/recombination coefficients of
	// 	k'th atomic state at a given temperature and density.
	// 	Args:
	// 		k  (int): Ionising or recombined ion stage,
	// 			between 0 and k=Z-1, where Z is atomic number.
	// 		Te (double): Temperature in [eV].
	// 		ne (double): Density in [m-3].
	// 	Returns:
	// 		c (double): Rate coefficent in [m3/s].

	// Perform a basic interpolation based on linear distance
	// values to search for
	double eval_log10_Te = log10(eval_Te);
	double eval_log10_Ne = log10(eval_Ne);
	// Look through the log_temperature and log_density attributes of RateCoefficient to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int low_Te = lower_bound(log_temperature.begin(), log_temperature.end(), eval_log10_Te) - log_temperature.begin() - 1;
	int low_Ne = lower_bound(log_density.begin(), log_density.end(), eval_log10_Ne) - log_density.begin() - 1;

	// Bounds checking -- make sure you haven't dropped off the end of the array
	if ((low_Te == (int)(log_temperature.size())-1) or (low_Te == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Te called to point off the grid for which it was defined (will give seg fault)");
	};
	if ((low_Ne == (int)(log_density.size()-1)) or (low_Ne == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Ne called to point off the grid for which it was defined (will give seg fault)");
	};

	int high_Te = low_Te + 1;
	int high_Ne = low_Ne + 1;

	double Te_norm = 1/(log_temperature[high_Te] - log_temperature[low_Te]); //Spacing between grid points
	double Ne_norm = 1/(log_density[high_Ne] - log_density[low_Ne]); //Spacing between grid points

	double x = (eval_log10_Te - log_temperature[low_Te])*Te_norm;
	double y = (eval_log10_Ne - log_density[low_Ne])*Ne_norm;

	// grid_matrix alpha_sub = alpha_coeff[k][low_Te][low_Ne];
	grid_matrix alpha_sub = {{
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}};
	double x_vector[4] = {1, x, x*x, x*x*x}; //Row vector
	double y_vector[4] = {1, y, y*y, y*y*y}; //Column vector

	double return_value = 0.0;
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			return_value += x_vector[i] * alpha_sub[i][j] * y_vector[j];
		}
	}
	std::printf("Return value: %f\n", return_value);
};
double RateCoefficient::call0D_bicubic(const int k, const std::pair<int, double> Te_interp, const std::pair<int, double> Ne_interp){
	int low_Te = Te_interp.first;
	int high_Te = low_Te+1;
	int low_Ne = Ne_interp.first;
	int high_Ne = low_Ne+1;

	double x = Te_interp.second;
	double y = Ne_interp.second;

	// grid_matrix alpha_sub = alpha_coeff[k][low_Te][low_Ne];
	grid_matrix alpha_sub = {{
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}};
	double x_vector[4] = {1, x, x*x, x*x*x}; //Row vector
	double y_vector[4] = {1, y, y*y, y*y*y}; //Column vector

	double return_value = 0.0;
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			return_value += x_vector[i] * alpha_sub[i][j] * y_vector[j];
		}
	}
	std::printf("Return value: %f\n", return_value);
};
double RateCoefficient::call0D_bilinear(const int k, const double eval_Te, const double eval_Ne){

	// """Evaluate the ionisation/recombination coefficients of
	// 	k'th atomic state at a given temperature and density.
	// 	Args:
	// 		k  (int): Ionising or recombined ion stage,
	// 			between 0 and k=Z-1, where Z is atomic number.
	// 		Te (double): Temperature in [eV].
	// 		ne (double): Density in [m-3].
	// 	Returns:
	// 		c (double): Rate coefficent in [m3/s].

	// Perform a basic interpolation based on linear distance
	// values to search for
	double eval_log10_Te = log10(eval_Te);
	double eval_log10_Ne = log10(eval_Ne);
	// Look through the log_temperature and log_density attributes of RateCoefficient to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int low_Te = lower_bound(log_temperature.begin(), log_temperature.end(), eval_log10_Te) - log_temperature.begin() - 1;
	int low_Ne = lower_bound(log_density.begin(), log_density.end(), eval_log10_Ne) - log_density.begin() - 1;

	// Bounds checking -- make sure you haven't dropped off the end of the array
	if ((low_Te == (int)(log_temperature.size())-1) or (low_Te == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Te called to point off the grid for which it was defined (will give seg fault)");
	};
	if ((low_Ne == (int)(log_density.size()-1)) or (low_Ne == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Ne called to point off the grid for which it was defined (will give seg fault)");
	};

	int high_Te = low_Te + 1;
	int high_Ne = low_Ne + 1;

	double Te_norm = 1/(log_temperature[high_Te] - log_temperature[low_Te]); //Spacing between grid points
	double Ne_norm = 1/(log_density[high_Ne] - log_density[low_Ne]); //Spacing between grid points

	double x = (eval_log10_Te - log_temperature[low_Te])*Te_norm;
	double y = (eval_log10_Ne - log_density[low_Ne])*Ne_norm;

	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    ne -> y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--Te -> x
	// //  | /     \ |      |
	// // w00 ------ w10

	double eval_log10_coeff =
	(log_coeff[k][low_Te][low_Ne]*(1-y) + log_coeff[k][low_Te][high_Ne]*y)*(1-x)
	+(log_coeff[k][high_Te][low_Ne]*(1-y) + log_coeff[k][high_Te][high_Ne]*y)*x;
	double eval_coeff = pow(10,eval_log10_coeff);

	return eval_coeff;
};
//Overloaded onto call0D - if the input is an int and two <int, double> pairs then use the SharedInterpolation method (i.e. assume that Te_interp and Ne_interp
//contain which point for which to return the coefficient - saves reevaluating)
double RateCoefficient::call0D_bilinear(const int k, const std::pair<int, double> Te_interp, const std::pair<int, double> Ne_interp){

	int low_Te = Te_interp.first;
	int high_Te = low_Te+1;
	int low_Ne = Ne_interp.first;
	int high_Ne = low_Ne+1;

	double x = Te_interp.second;
	double y = Ne_interp.second;

	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    Ne -> y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--Te -> x
	// //  | /     \ |      |
	// // w00 ------ w10

	double eval_log10_coeff =
	(log_coeff[k][low_Te][low_Ne]*(1-y) + log_coeff[k][low_Te][high_Ne]*y)*(1-x)
	+(log_coeff[k][high_Te][low_Ne]*(1-y) + log_coeff[k][high_Te][high_Ne]*y)*x;
	
	double eval_coeff = pow(10,eval_log10_coeff);

	return eval_coeff;
};
int RateCoefficient::get_atomic_number(){
	return atomic_number;
};
std::string RateCoefficient::get_element(){
	return element;
};
std::string RateCoefficient::get_adf11_file(){
	return adf11_file;
};
std::vector<std::vector< std::vector<double> > > RateCoefficient::get_log_coeff(){
	return log_coeff;
};
std::vector<double> RateCoefficient::get_log_temperature(){
	return log_temperature;
};
std::vector<double> RateCoefficient::get_log_density(){
	return log_density;
};
std::vector<std::vector<std::vector<interp_data>>> RateCoefficient::calculate_grid_coeff(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff){

	int L_k = (int)(log_coeff.size());
	int L_t = log_temperature.size();
	int L_n = log_density.size();

	std::vector<std::vector<std::vector<interp_data>>>
	grid_coeff(L_k,std::vector<std::vector<interp_data>>(L_t,std::vector<interp_data>(L_n,default_interp_data)));

	for (int k=0; k<L_k; ++k){
		for(int iT=0; iT<L_t; ++iT){
			for(int iN=0; iN<L_n; ++iN){

				// Set the function value
				grid_coeff[k][iT][iN].f = log_coeff[k][iT][iN];

				double dT_difference = 0.0;
				double dT_spacing = 0.0;
				if((iT != 0) and (iT != (int)(log_temperature.size()-1))){
					// Central difference for dT
					dT_difference = log_coeff[k][iT+1][iN] - log_coeff[k][iT-1][iN];
					dT_spacing = log_temperature[iT+1] - log_temperature[iT-1];

				} else if (iT == 0) {
					// Forward difference for dT
					dT_difference = log_coeff[k][iT+1][iN] - log_coeff[k][iT][iN];
					dT_spacing = log_temperature[iT+1] - log_temperature[iT];

				} else if (iT == (int)(log_temperature.size()-1)){
					// Backward difference for dT
					dT_difference = log_coeff[k][iT][iN] - log_coeff[k][iT-1][iN];
					dT_spacing = log_temperature[iT] - log_temperature[iT-1];
					
				}
				grid_coeff[k][iT][iN].fdT = dT_difference/dT_spacing;

				double dN_difference = 0.0;
				double dN_spacing = 0.0;
				if((iN != 0) and (iN != (int)(log_density.size()-1))){
					// Central difference for dN
					dN_difference = log_coeff[k][iT][iN+1] - log_coeff[k][iT][iN-1];
					dN_spacing = log_density[iN+1] - log_density[iN-1];

				} else if (iN == 0) {
					// Forward difference for dN
					dN_difference = log_coeff[k][iT][iN+1] - log_coeff[k][iT][iN];
					dN_spacing = log_density[iN+1] - log_density[iN];

				} else if (iN == (int)(log_density.size()-1)){
					// Backward difference for dN
					dN_difference = log_coeff[k][iT][iN] - log_coeff[k][iT][iN-1];
					dN_spacing = log_density[iN] - log_density[iN-1];
					
				}
				grid_coeff[k][iT][iN].fdN = dN_difference/dN_spacing;

			}
		}

		//Now that all axial derivatives have been calculated, use these results to calculate the mixed derivatives

		for(int iT=0; iT<L_t; ++iT){
			for(int iN=0; iN<L_n; ++iN){
				double dTN_difference = 0.0;
				double dTN_spacing = 0.0;
				if((iT != 0) and (iT != (int)(log_temperature.size()-1))){
					// Central difference for dTN
					dTN_difference = grid_coeff[k][iT+1][iN].fdN - grid_coeff[k][iT-1][iN].fdN;
					dTN_spacing = log_temperature[iT+1] - log_temperature[iT-1];

				} else if (iT == 0) {
					// Forward difference for dTN
					dTN_difference = grid_coeff[k][iT+1][iN].fdN - grid_coeff[k][iT][iN].fdN;
					dTN_spacing = log_temperature[iT+1] - log_temperature[iT];

				} else if (iT == (int)(log_temperature.size()-1)){
					// Backward difference for dTN
					dTN_difference = grid_coeff[k][iT][iN].fdN - grid_coeff[k][iT-1][iN].fdN;
					dTN_spacing = log_temperature[iT] - log_temperature[iT-1];
					
				}
				grid_coeff[k][iT][iN].fdTdN = dTN_difference/dTN_spacing;
			}
		}
	}

	return grid_coeff;
};
std::vector<std::vector<std::vector<grid_matrix>>> RateCoefficient::calculate_alpha_coeff(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff){
	// For storing the value and derivative data at each grid-point
	// Have to use vector since the array size is non-constant
	int L_k = (int)(log_coeff.size());
	int L_t = log_temperature.size();
	int L_n = log_density.size();

	std::vector<std::vector<std::vector<interp_data>>>
	grid_coeff = calculate_grid_coeff(log_temperature, log_density, log_coeff);
	
	grid_matrix default_alpha_coeff = {0.0};

	std::vector<std::vector<std::vector<grid_matrix>>>
	alpha_coeff(L_k,std::vector<std::vector<grid_matrix>>(L_t-1,std::vector<grid_matrix>(L_n-1,default_alpha_coeff)));

	const grid_matrix prematrix = {{
			{+1, +0, +0, +0},
			{+0, +0, +1, +0},
			{-3, +3, -2, -1},
			{+2, -2, +1, +1},
		}};
	const grid_matrix postmatrix = {{
			{+1, +0, -3, +2},
			{+0, +0, +3, -2},
			{+0, +1, -2, +1},
			{+0, +0, -1, +1},
		}};

	for (int k=0; k<L_k; ++k){
		for(int iT=0; iT<L_t - 1; ++iT){ //iterator over temperature dimension of grid, which is of length log_temperature.size()-1
			for(int iN=0; iN<L_n - 1; ++iN){ //iterator over density dimension of grid, which is of length log_density.size()-1

				grid_matrix f_sub = {{
					{grid_coeff[k][iT+0][iN+0].f,   grid_coeff[k][iT+0][iN+1].f,   grid_coeff[k][iT+0][iN+0].fdN,   grid_coeff[k][iT+0][iN+1].fdN},
					{grid_coeff[k][iT+1][iN+0].f,   grid_coeff[k][iT+1][iN+1].f,   grid_coeff[k][iT+1][iN+0].fdN,   grid_coeff[k][iT+1][iN+1].fdN},
					{grid_coeff[k][iT+0][iN+0].fdT, grid_coeff[k][iT+0][iN+1].fdT, grid_coeff[k][iT+0][iN+0].fdTdN, grid_coeff[k][iT+0][iN+1].fdTdN},
					{grid_coeff[k][iT+1][iN+0].fdT, grid_coeff[k][iT+1][iN+1].fdT, grid_coeff[k][iT+1][iN+0].fdTdN, grid_coeff[k][iT+1][iN+1].fdTdN},
				}};
				// grid_coeff submatrix
				grid_matrix alpha_sub = {0.0};
				
				//Matrix multiply prematrix * f_sub * postmatrix to find alpha_sub
				//As per https://en.wikipedia.org/wiki/Bicubic_interpolation
				for (int i=0; i<4; ++i){
					for(int j=0; j<4; ++j){
						for(int k=0; k<4; ++k){
							for(int l=0; l<4; ++l){
								alpha_sub[i][j] += prematrix[i][l] * f_sub[l][k] * postmatrix[k][j];
							}
						}
					}
				}
				alpha_coeff[k][iT][iN] = alpha_sub;
			}
		}
	}
	return alpha_coeff;
};