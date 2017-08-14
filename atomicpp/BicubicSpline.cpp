#include <vector>
#include <array>
#include <math.h>
#include <algorithm> //for upper/lower_bound
#include <stdexcept> //For error-throwing
#include "BicubicSpline.hpp"

#include <iostream>

using namespace atomicpp;
interp_data default_interp_data;
grid_matrix default_alpha_coeff = {{
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
}};

BicubicSpline::BicubicSpline(){//Default constructor
};
BicubicSpline::BicubicSpline(
	std::vector<double>& _x_values,
	std::vector<double>& _y_values,
	std::vector< std::vector<double> > & _z_values
	){
	x_values = _x_values;
	y_values = _y_values;
	z_values = _z_values;

	alpha_coeff = calculate_alpha_coeff();	
};
double BicubicSpline::call0D(const double eval_x, const double eval_y){

	// Perform a basic interpolation based on linear distance
	// values to search for
	// Look through the x_values and y_values attributes of RateCoefficient to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int low_x = lower_bound(x_values.begin(), x_values.end(), eval_x) - x_values.begin() - 1;
	int low_y = lower_bound(y_values.begin(), y_values.end(), eval_y) - y_values.begin() - 1;

	// Bounds checking -- make sure you haven't dropped off the end of the array
	if ((low_x == (int)(x_values.size())-1) or (low_x == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on x called to point off the grid for which it was defined (will give seg fault)");
	};
	if ((low_y == (int)(y_values.size()-1)) or (low_y == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on y called to point off the grid for which it was defined (will give seg fault)");
	};

	double x_norm = 1/(x_values[low_x+1] - x_values[low_x+0]); //Spacing between grid points
	double y_norm = 1/(y_values[low_y+1] - y_values[low_y+0]); //Spacing between grid points

	double x = (eval_x - x_values[low_x+0])*x_norm;
	double y = (eval_y - y_values[low_y+0])*y_norm;

	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--x
	// //  | /     \ |      |
	// // w00 ------ w10

	double eval_coeff =
	 (z_values[low_x+0][low_y+0]*(1-y) + z_values[low_x+0][low_y+1]*y)*(1-x)
	+(z_values[low_x+1][low_y+0]*(1-y) + z_values[low_x+1][low_y+1]*y)*x;

	return eval_coeff;
};
//Overloaded onto callOD - if the input is an int and two <int, double> pairs then use the SharedInterpolation method (i.e. assume that x_interp and y_interp
//contain which point for which to return the coefficient - saves reevaluating)
double BicubicSpline::call0D_shared(const std::pair<int, double> x_interp, const std::pair<int, double> y_interp){

	int low_x = x_interp.first;
	int low_y = y_interp.first;

	double x = x_interp.second;
	double y = y_interp.second;

	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--x
	// //  | /     \ |      |
	// // w00 ------ w10

	double eval_coeff =
	 (z_values[low_x+0][low_y+0]*(1-y) + z_values[low_x+0][low_y+1]*y)*(1-x)
	+(z_values[low_x+1][low_y+0]*(1-y) + z_values[low_x+1][low_y+1]*y)*x;

	return eval_coeff;
};
std::vector< std::vector<double> > BicubicSpline::get_z_values(){
	return z_values;
};
std::vector<double> BicubicSpline::get_x_values(){
	return x_values;
};
std::vector<double> BicubicSpline::get_y_values(){
	return y_values;
};
void BicubicSpline::set_x_values(std::vector<double>& _x_values){
	x_values = _x_values;
};
void BicubicSpline::set_y_values(std::vector<double>& _y_values){
	y_values = _y_values;
};
void BicubicSpline::zero_z_values(){
	for(int i = 0; i<(int)(z_values.size()); ++i){
		for(int j = 0; j<(int)(z_values[0].size()); ++j){
			z_values[i][j] = 0.0;
		}
	}
};
std::vector<std::vector<interp_data>> BicubicSpline::calculate_grid_coeff(){

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

std::vector<std::vector<grid_matrix>> BicubicSpline::calculate_alpha_coeff(){
	// For storing the value and derivative data at each grid-point
	// Have to use vector since the array size is non-constant
	int Lx = x_values.size();
	int Ly = y_values.size();

	std::vector<std::vector<interp_data>> grid_coeff = calculate_grid_coeff();
	
	grid_matrix default_alpha_coeff = {0.0};

	std::vector<std::vector<grid_matrix>>
	alpha_coeff(Lx-1,std::vector<grid_matrix>(Ly-1,default_alpha_coeff));

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
			alpha_coeff[x][y] = alpha_sub;
		}
	}
	return alpha_coeff;
};