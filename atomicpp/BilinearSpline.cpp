#include <vector>
#include <math.h>
#include <algorithm> //for upper/lower_bound
#include <stdexcept> //For error-throwing
#include "BilinearSpline.hpp"

using namespace atomicpp;
BilinearSpline::BilinearSpline(){//Default constructor
};
BilinearSpline::BilinearSpline(
	std::vector<double>& _x_values,
	std::vector<double>& _y_values,
	std::vector< std::vector<double> > & _z_values
	){
	x_values = _x_values;
	y_values = _y_values;
	z_values = _z_values;

	//Doesn't require initialisation
};
double BilinearSpline::call0D(const double eval_x, const double eval_y){

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
double BilinearSpline::call0D_shared(const std::pair<int, double> x_interp, const std::pair<int, double> y_interp){

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
std::vector< std::vector<double> > BilinearSpline::get_z_values(){
	return z_values;
};
std::vector<double> BilinearSpline::get_x_values(){
	return x_values;
};
std::vector<double> BilinearSpline::get_y_values(){
	return y_values;
};
void BilinearSpline::set_x_values(std::vector<double>& _x_values){
	x_values = _x_values;
};
void BilinearSpline::set_y_values(std::vector<double>& _y_values){
	y_values = _y_values;
};