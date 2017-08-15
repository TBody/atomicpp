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

// typedef std::array<std::array<double, 4>, 4> grid_matrix;
typedef std::vector<std::vector<double>> grid_matrix;
grid_matrix default_grid_matrix(4, std::vector<double>(4, 0.0));
struct interp_data{
	// std::pair<double, double> coord; //(T,N) coordinate of point
	double f = 0.0; //Value at point
	double fdx = 0.0; //Derivative in temperature axis
	double fdy = 0.0; //Derivative in density axis
	double fdxdy = 0.0; //Cross derivative
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

std::tuple< int, std::vector<double>, std::vector<double>, std::vector< std::vector<double> > > extract_from_json(){

	std::string filename("json_database/json_data/acd96_c.json");

	json data_dict = retrieveFromJSON(filename);

	int atomic_number   		= data_dict["charge"];
	// std::string element         = data_dict["element"];
	// std::string adf11_file      = filename;

	std::vector< std::vector< std::vector<double> > > extract_z_values = data_dict["log_coeff"];
	std::vector<double> extract_x_values = data_dict["log_temperature"];
	std::vector<double> extract_y_values = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute
	std::vector< std::vector<double> > z_values = extract_z_values[0];
	std::vector<double> x_values = extract_x_values;
	std::vector<double> y_values = extract_y_values;

	return std::make_tuple(atomic_number, x_values, y_values, z_values);
}

std::vector<std::vector<interp_data>> calculate_grid_coeff(std::vector<double>& x_values, std::vector<double>& y_values, std::vector< std::vector<double> >& z_values){

	int Lx = x_values.size();
	int Ly = y_values.size();

	std::vector<std::vector<interp_data>>
	grid_coeff(Lx,std::vector<interp_data>(Ly,default_interp_data));

	for(int x=0; x<Lx; ++x){
		for(int y=0; y<Ly; ++y){

			// Set the function value
			grid_coeff[x][y].f = z_values[x][y];

			double dx_difference = 0.0;
			double dx_spacing = 0.0;
			if((x != 0) and (x != (int)(x_values.size()-1))){
				// Central difference for dx
				dx_difference = z_values[x+1][y] - z_values[x-1][y];
				dx_spacing = x_values[x+1] - x_values[x-1];

			} else if (x == 0) {
				// Forward difference for dx
				dx_difference = z_values[x+1][y] - z_values[x][y];
				dx_spacing = x_values[x+1] - x_values[x];

			} else if (x == (int)(x_values.size()-1)){
				// Backward difference for dx
				dx_difference = z_values[x][y] - z_values[x-1][y];
				dx_spacing = x_values[x] - x_values[x-1];
				
			}
			grid_coeff[x][y].fdx = dx_difference/dx_spacing;

			double dy_difference = 0.0;
			double dy_spacing = 0.0;
			if((y != 0) and (y != (int)(y_values.size()-1))){
				// Central difference for dy
				dy_difference = z_values[x][y+1] - z_values[x][y-1];
				dy_spacing = y_values[y+1] - y_values[y-1];

			} else if (y == 0) {
				// Forward difference for dy
				dy_difference = z_values[x][y+1] - z_values[x][y];
				dy_spacing = y_values[y+1] - y_values[y];

			} else if (y == (int)(y_values.size()-1)){
				// Backward difference for dy
				dy_difference = z_values[x][y] - z_values[x][y-1];
				dy_spacing = y_values[y] - y_values[y-1];
				
			}
			grid_coeff[x][y].fdy = dy_difference/dy_spacing;

		}
	}

	//Now that all axial derivatives have been calculated, use these results to calculate the mixed derivatives

	for(int x=0; x<Lx; ++x){
		for(int y=0; y<Ly; ++y){
			double dxy_difference = 0.0;
			double dxy_spacing = 0.0;
			if((x != 0) and (x != (int)(x_values.size()-1))){
				// Central difference for dxy
				dxy_difference = grid_coeff[x+1][y].fdy - grid_coeff[x-1][y].fdy;
				dxy_spacing = x_values[x+1] - x_values[x-1];

			} else if (x == 0) {
				// Forward difference for dxy
				dxy_difference = grid_coeff[x+1][y].fdy - grid_coeff[x][y].fdy;
				dxy_spacing = x_values[x+1] - x_values[x];

			} else if (x == (int)(x_values.size()-1)){
				// Backward difference for dxy
				dxy_difference = grid_coeff[x][y].fdy - grid_coeff[x-1][y].fdy;
				dxy_spacing = x_values[x] - x_values[x-1];
				
			}
			grid_coeff[x][y].fdxdy = dxy_difference/dxy_spacing;
		}
	}

	return grid_coeff;
};

std::vector<std::vector<grid_matrix>> calculate_alpha_coeff(std::vector<double>& x_values, std::vector<double>& y_values, std::vector< std::vector<double> >& z_values){
	// For storing the value and derivative data at each grid-point
	// Have to use vector since the array size is non-constant
	int Lx = x_values.size();
	int Ly = y_values.size();

	std::vector<std::vector<interp_data>>
	grid_coeff = calculate_grid_coeff(x_values, y_values, z_values);

	std::vector<std::vector<grid_matrix>>
	alpha_coeff(Lx-1,std::vector<grid_matrix>(Ly-1, default_grid_matrix));

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

	for(int x=0; x<Lx - 1; ++x){ //iterator over temperature dimension of grid, which is of length x_values.size()-1
		for(int y=0; y<Ly - 1; ++y){ //iterator over density dimension of grid, which is of length y_values.size()-1

			grid_matrix f_sub = {{
				{grid_coeff[x+0][y+0].f,   grid_coeff[x+0][y+1].f,   grid_coeff[x+0][y+0].fdy,   grid_coeff[x+0][y+1].fdy},
				{grid_coeff[x+1][y+0].f,   grid_coeff[x+1][y+1].f,   grid_coeff[x+1][y+0].fdy,   grid_coeff[x+1][y+1].fdy},
				{grid_coeff[x+0][y+0].fdx, grid_coeff[x+0][y+1].fdx, grid_coeff[x+0][y+0].fdxdy, grid_coeff[x+0][y+1].fdxdy},
				{grid_coeff[x+1][y+0].fdx, grid_coeff[x+1][y+1].fdx, grid_coeff[x+1][y+0].fdxdy, grid_coeff[x+1][y+1].fdxdy},
			}};
			// grid_coeff submatrix
			grid_matrix alpha_sub = default_grid_matrix;
			
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
			alpha_coeff[x][y] = alpha_sub;
		}
	}
	return alpha_coeff;
};

int main(){

	auto json_tuple = extract_from_json();
	std::vector<double> x_values = std::get<1>(json_tuple);
	std::vector<double> y_values = std::get<2>(json_tuple);
	std::vector<std::vector<double>> z_values = std::get<3>(json_tuple);
	// std::vector<double> x_values = x_values;
	// std::vector<double> y_values = y_values;
	// std::vector<std::vector<double> > z_values = z_values;


	std::vector<std::vector<grid_matrix>> alpha_coeff = calculate_alpha_coeff(x_values, y_values, z_values);

	
	double eval_log10_Te = log10(50);
	double eval_log10_Ne = log10(0.8e19);

	int low_Te = lower_bound(x_values.begin(), x_values.end(), eval_log10_Te) - x_values.begin() - 1;
	int low_Ne = lower_bound(y_values.begin(), y_values.end(), eval_log10_Ne) - y_values.begin() - 1;
	// Bounds checking -- make sure you haven't dropped off the end of the array
	if ((low_Te == (int)(x_values.size())-1) or (low_Te == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Te called to point off the grid for which it was defined (will give seg fault)");
	};
	if ((low_Ne == (int)(y_values.size()-1)) or (low_Ne == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on Ne called to point off the grid for which it was defined (will give seg fault)");
	};

	grid_matrix alpha_sub = alpha_coeff[low_Te][low_Ne];

	double x = (eval_log10_Te - x_values[low_Te])/(x_values[low_Te + 1] - x_values[low_Te]);
	double x_vector[4] = {1, x, x*x, x*x*x}; //Row vector
	double y = (eval_log10_Ne - y_values[low_Ne])/(y_values[low_Ne + 1] - y_values[low_Ne]);
	double y_vector[4] = {1, y, y*y, y*y*y}; //Column vector

	double return_value = 0.0;
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			return_value += x_vector[i] * alpha_sub[i][j] * y_vector[j];
		}
	}
	std::printf("Return value: %f\n", return_value);

}



















