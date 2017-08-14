#include <vector>
#include <math.h>
#include <algorithm> //for upper/lower_bound
#include <stdexcept> //For error-throwing
#include "BilinearSpline.hpp"

using namespace atomicpp;
BilinearSpline::BilinearSpline(){//Default constructor
};
BilinearSpline::BilinearSpline(
	std::vector<double>& _temp_values,
	std::vector<double>& _dens_values,
	std::vector<std::vector< std::vector<double> > >& _coef_values
	){
	temp_values = _temp_values;
	dens_values = _dens_values;
	coef_values = _coef_values;

	//Doesn't require initialisation
};
double BilinearSpline::call0D(const int k, const double eval_temp, const double eval_dens){

	// """Evaluate the ionisation/recombination coefficients of
	// 	k'th atomic state at a given temperature and density.
	// 	Args:
	// 		k  (int): Ionising or recombined ion stage,
	// 			between 0 and k=Z-1, where Z is atomic number.
	// 		temp (double): tempmperature in [eV].
	// 		ne (double): Density in [m-3].
	// 	Returns:
	// 		c (double): Rate coefficent in [m3/s].

	// Perform a basic interpolation based on linear distance
	// values to search for
	double eval_log_temp = log10(eval_temp);
	double eval_log_dens = log10(eval_dens);
	// Look through the temp_values and dens_values attributes of RateCoefficient to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int low_temp = lower_bound(temp_values.begin(), temp_values.end(), eval_log_temp) - temp_values.begin() - 1;
	int low_dens = lower_bound(dens_values.begin(), dens_values.end(), eval_log_dens) - dens_values.begin() - 1;

	// Bounds checking -- make sure you haven't dropped off the end of the array
	if ((low_temp == (int)(temp_values.size())-1) or (low_temp == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on temp called to point off the grid for which it was defined (will give seg fault)");
	};
	if ((low_dens == (int)(dens_values.size()-1)) or (low_dens == -1)){
		// An easy error to make is supplying the function arguments already having taken the log10
		throw std::runtime_error("Interpolation on dens called to point off the grid for which it was defined (will give seg fault)");
	};

	double temp_norm = 1/(temp_values[low_temp+1] - temp_values[low_temp+0]); //Spacing between grid points
	double dens_norm = 1/(dens_values[low_dens+1] - dens_values[low_dens+0]); //Spacing between grid points

	double x = (eval_log_temp - temp_values[low_temp+0])*temp_norm;
	double y = (eval_log_dens - dens_values[low_dens+0])*dens_norm;

	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    dens -> y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--temp -> x
	// //  | /     \ |      |
	// // w00 ------ w10

	double eval_coef_values =
	 (coef_values[k][low_temp+0][low_dens+0]*(1-y) + coef_values[k][low_temp+0][low_dens+1]*y)*(1-x)
	+(coef_values[k][low_temp+1][low_dens+0]*(1-y) + coef_values[k][low_temp+1][low_dens+1]*y)*x;
	
	double eval_coeff = pow(10,eval_coef_values);

	return eval_coeff;
};
//Overloaded onto callOD - if the input is an int and two <int, double> pairs then use the SharedInterpolation method (i.e. assume that temp_interp and dens_interp
//contain which point for which to return the coefficient - saves reevaluating)
double BilinearSpline::call0D_shared(const int k, const std::pair<int, double> temp_interp, const std::pair<int, double> dens_interp){

	int low_temp = temp_interp.first;
	int low_dens = dens_interp.first;

	double x = temp_interp.second;
	double y = dens_interp.second;

	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    dens -> y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--temp -> x
	// //  | /     \ |      |
	// // w00 ------ w10

	double eval_log10_coeff =
	 (coef_values[k][low_temp+0][low_dens+0]*(1-y) + coef_values[k][low_temp+0][low_dens+1]*y)*(1-x)
	+(coef_values[k][low_temp+1][low_dens+0]*(1-y) + coef_values[k][low_temp+1][low_dens+1]*y)*x;
	
	double eval_coeff = pow(10,eval_log10_coeff);

	return eval_coeff;
};
std::vector<std::vector< std::vector<double> > > BilinearSpline::get_coef_values(){
	return coef_values;
};
std::vector<double> BilinearSpline::get_temp_values(){
	return temp_values;
};
std::vector<double> BilinearSpline::get_dens_values(){
	return dens_values;
};
void BilinearSpline::set_temp_values(std::vector<double>& _temp_values){
	temp_values = _temp_values;
};
void BilinearSpline::set_dens_values(std::vector<double>& _dens_values){
	dens_values = _dens_values;
};