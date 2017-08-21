#include <string>
#include <vector>
#include <array>
#include <fstream>
#include "json.hpp"
using json = nlohmann::json;
#include <stdexcept> //For error-throwing

#include "RateCoefficient.hpp"
#include "sharedFunctions.hpp"
#include "BivariateBSpline.hpp"

using namespace atomicpp;

RateCoefficient::RateCoefficient(const std::string& filename){
	// # Create an instance of RateCoefficient by reading an OpenADAS JSON file

	json data_dict = retrieveFromJSON(filename);

	atomic_number   = data_dict["charge"];
	element         = data_dict["element"];
	adf11_file      = filename;

	std::vector<std::vector< std::vector<double> > > log_rate = data_dict["log_coeff"];
	std::vector<double> log_temp = data_dict["log_temperature"];
	std::vector<double> log_dens = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute

	for(int k=0; k<atomic_number; ++k){
		interpolator.push_back(BivariateBSpline(log_temp, log_dens, log_rate[k]));
	}

};
RateCoefficient::RateCoefficient(const std::shared_ptr<RateCoefficient> source_RC){
	// # Create an instance of a blank RateCoefficient by copying from another RateCoefficient object

	atomic_number   = source_RC->get_atomic_number();
	element         = source_RC->get_element();
	adf11_file      = source_RC->get_adf11_file();

	interpolator = source_RC->get_interpolator();
};
double RateCoefficient::call0D(const int k, const double eval_Te, const double eval_Ne){

	double eval_log_Te = log10(eval_Te);
	double eval_log_Ne = log10(eval_Ne);

	double eval_log_rate = interpolator[k].call0D(eval_log_Te, eval_log_Ne);

	return pow(10, eval_log_rate);
};
//Overloaded onto call0D - if the input is an int and two <int, double> pairs then use the SharedInterpolation method (i.e. assume that Te_interp and Ne_interp
//contain which point for which to return the coefficient - saves reevaluating)
double RateCoefficient::call0D(const int k, const std::pair<int, double> Te_interp, const std::pair<int, double> Ne_interp){
	double eval_log_rate = interpolator[k].call0D_shared(Te_interp, Ne_interp);
	return pow(10, eval_log_rate);
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
std::vector<BivariateBSpline> RateCoefficient::get_interpolator(){
	return interpolator;
};
std::vector<double> RateCoefficient::get_log_temp(){
	return interpolator[0].get_x_values();
};
std::vector<double> RateCoefficient::get_log_dens(){
	return interpolator[0].get_y_values();
};
void RateCoefficient::zero_interpolator(){
	for(int k=0; k < (int)(interpolator.size()); ++k){
		interpolator[k].zero_z_values();
	}
};