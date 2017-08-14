#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
using json = nlohmann::json;
#include <stdexcept> //For error-throwing

#include "RateCoefficient.hpp"
#include "sharedFunctions.hpp"
#include "BilinearSpline.hpp"

using namespace atomicpp;
RateCoefficient::RateCoefficient(const std::string& filename){
	// # Create an instance of RateCoefficient by reading an OpenADAS JSON file

	json data_dict = retrieveFromJSON(filename);

	atomic_number   = data_dict["charge"];
	element         = data_dict["element"];
	adf11_file      = filename;

	std::vector<std::vector< std::vector<double> > > log_coeff = data_dict["log_coeff"];
	std::vector<double> log_temperature = data_dict["log_temperature"];
	std::vector<double> log_density = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute
	
	interpolator = BilinearSpline(log_temperature, log_density, log_coeff);

};
RateCoefficient::RateCoefficient(const std::shared_ptr<RateCoefficient> source_rc){
	// # Create an instance of a blank RateCoefficient by copying from another RateCoefficient object

	atomic_number   = source_rc->get_atomic_number();
	element         = source_rc->get_element();
	adf11_file      = source_rc->get_adf11_file();

	std::vector<double> temp_values = source_rc->get_interpolator().get_temp_values();
	std::vector<double> dens_values = source_rc->get_interpolator().get_dens_values();

	interpolator.set_temp_values(temp_values);
	interpolator.set_dens_values(dens_values);
	//Don't copy the coef_values
};
double RateCoefficient::call0D(const int k, const double Te, const double Ne){
	return interpolator.call0D(k, Te, Ne);
};
//Overloaded onto callOD - if the input is an int and two <int, double> pairs then use the SharedInterpolation method (i.e. assume that Te_interp and Ne_interp
//contain which point for which to return the coefficient - saves reevaluating)
double RateCoefficient::call0D(const int k, const std::pair<int, double> Te_interp, const std::pair<int, double> Ne_interp){
	return interpolator.call0D_shared(k, Te_interp, Ne_interp);
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
BilinearSpline RateCoefficient::get_interpolator(){
	return interpolator;
};
std::vector<std::vector< std::vector<double> > > RateCoefficient::get_log_coeff(){
	return interpolator.get_coef_values();
};
std::vector<double> RateCoefficient::get_log_temperature(){
	return interpolator.get_temp_values();
};
std::vector<double> RateCoefficient::get_log_density(){
	return interpolator.get_dens_values();
};