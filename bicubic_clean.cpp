// https://shiftedbits.org/2011/01/30/cubic-spline-interpolation/
// https://en.wikipedia.org/wiki/Bicubic_interpolation
// http://www.iue.tuwien.ac.at/phd/heinzl/node27.html
// https://en.wikipedia.org/wiki/Matrix_multiplication

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept> //For error-throwing
#include <array>

#include <tuple> //Not for RateCoefficient

#include <ostream>
#include <cstdio> //For print formatting (printf, fprintf, sprintf, snprintf)

#include "atomicpp/json.hpp"
using json = nlohmann::json;

#include <algorithm> //for upper/lower_bound

typedef std::array<std::array<double, 4>, 4> grid_matrix;
struct interp_data{
	// std::pair<double, double> coord; //(T,N) coordinate of point
	double f = 0.0; //Value at point
	double fdT = 0.0; //Derivative in temperature axis
	double fdN = 0.0; //Derivative in density axis
	double fdTdN = 0.0; //Cross derivative
} default_interp_data;

json retrieveFromJSON(std::string path_to_file){
	// Do not pass path_to_file by reference - results in error!
	// Reads a .json file given at path_to_file
	// Uses the json module at https://github.com/nlohmann/json/
	// This relies upon the "json.hpp" header which must be included in the same folder as the source

	// Open a file-stream at path_to_file
	std::ifstream json_file(path_to_file);
	// Initialise a json file object at j_object
	json j_object;
	json_file >> j_object;
	return j_object;
};

std::tuple< int, std::vector<double>, std::vector<double>, std::vector<std::vector< std::vector<double> > > > extract_from_json(){

	std::string filename("json_database/json_data/acd96_c.json");

	json data_dict = retrieveFromJSON(filename);

	int atomic_number   		= data_dict["charge"];
	// std::string element         = data_dict["element"];
	// std::string adf11_file      = filename;

	std::vector< std::vector< std::vector<double> > > extract_log_coeff = data_dict["log_coeff"];
	std::vector<double> extract_log_temperature = data_dict["log_temperature"];
	std::vector<double> extract_log_density = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute
	std::vector<std::vector< std::vector<double> > > log_coeff = extract_log_coeff;
	std::vector<double> log_temperature = extract_log_temperature;
	std::vector<double> log_density = extract_log_density;

	return std::make_tuple(atomic_number, log_temperature, log_density, log_coeff);
}

std::vector<std::vector<std::vector<interp_data>>> calculate_grid_coeff(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff){

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

std::vector<std::vector<std::vector<grid_matrix>>> calculate_alpha_coeff(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff){
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

int main(){

	auto json_tuple = extract_from_json();
	std::vector<double> log_temperature = std::get<1>(json_tuple);
	std::vector<double> log_density = std::get<2>(json_tuple);
	std::vector<std::vector< std::vector<double> > > log_coeff = std::get<3>(json_tuple);

	std::vector<std::vector<std::vector<grid_matrix>>> alpha_coeff = calculate_alpha_coeff(log_temperature, log_density, log_coeff);

	int k = 0;
	double eval_log10_Te = log10(50);
	double eval_log10_Ne = log10(0.8e19);

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

	grid_matrix alpha_sub = alpha_coeff[k][low_Te][low_Ne];

	double x = (eval_log10_Te - log_temperature[low_Te])/(log_temperature[low_Te + 1] - log_temperature[low_Te]);
	double x_vector[4] = {1, x, x*x, x*x*x}; //Row vector
	double y = (eval_log10_Ne - log_density[low_Ne])/(log_density[low_Ne + 1] - log_density[low_Ne]);
	double y_vector[4] = {1, y, y*y, y*y*y}; //Column vector

	double return_value = 0.0;
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			return_value += x_vector[i] * alpha_sub[i][j] * y_vector[j];
		}
	}
	std::printf("Return value: %f\n", return_value);

}



















