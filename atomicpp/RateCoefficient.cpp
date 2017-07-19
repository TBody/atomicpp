#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"

#include "RateCoefficient.hpp"
#include "sharedFunctions.hpp"
using namespace std; //saves having to prepend std:: onto common functions

#include <algorithm> //for upper/lower_bound
#include <typeinfo>

// for convenience
using json = nlohmann::json;

RateCoefficient::RateCoefficient(const string& filename){
	// # Create an instance of RateCoefficient by reading an OpenADAS JSON file
	
	json data_dict = retrieveFromJSON(filename);

	atomic_number   = data_dict["charge"];
	element         = data_dict["element"];
	adf11_file      = filename;

	vector<vector< vector<double> > > extract_log_coeff = data_dict["log_coeff"];
	vector<double> extract_log_temperature = data_dict["log_temperature"];
	vector<double> extract_log_density = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute
	log_coeff = extract_log_coeff;
	log_temperature = extract_log_temperature;
	log_density = extract_log_density;

};
ostream& operator<<(ostream& os, const RateCoefficient& RC){  
    os << "RateCoefficient object from " << RC.adf11_file << endl;
    return os;  
}
void RateCoefficient::call1D(const int k, const int data_length, const vector<double>& log10_Te, const vector<double>& log10_ne, vector<double>& log10_coeffs){
	// """Evaluate the ionisation/recombination coefficients of
	// 	k'th atomic state at a given temperature and density.
	// 	Args:
	// 		k  (int): Ionising or recombined ion stage,
	// 			between 0 and k=Z-1, where Z is atomic number.
	// 		Te (array_like): Temperature in [eV].
	// 		ne (array_like): Density in [m-3].
	// 	Returns:
	// 		c (array_like): Rate coefficent in [m3/s].

	// Iterate along the Te and ne values inputted.
	// Assume that s is the location point, Te(s) and ne(s) are temp and density values
	// as you move along the field line
	
	// Perform a basic interpolation based on linear distance
	for(int s=0; s<data_length; ++s){
		// values to search for
		double log10_Te_search = log10_Te[s];
		double log10_ne_search = log10_ne[s];
		// Look through the log_temperature and log_density attributes of RateCoefficient to find nearest (strictly lower)
		// Subtract 1 from answer to account for indexing from 0
		int low_Te = lower_bound(log_temperature.begin(), log_temperature.end(), log10_Te_search) - log_temperature.begin() - 1;
		int low_ne = lower_bound(log_density.begin(), log_density.end(), log10_ne_search) - log_density.begin() - 1;

		// log10_coeffs[s] = log_coeff[k][low_Te][low_ne];

		int high_Te = low_Te + 1;
		int high_ne = low_ne + 1;
		double Te_norm = 1/(log_temperature[high_Te] - log_temperature[low_Te]); //Spacing between grid points
		double ne_norm = 1/(log_density[high_ne] - log_density[low_ne]); //Spacing between grid points

		double x = (log10_Te_search - log_temperature[low_Te])*Te_norm;
		double y = (log10_ne_search - log_density[low_ne])*ne_norm;
		
		// // Construct the simple interpolation grid
		// // Find weightings based on linear distance
		// // w01 ------ w11    ne -> y
		// //  | \     / |      |
		// //  |  w(x,y) |    --/--Te -> x
		// //  | /     \ |      |
		// // w00 ------ w10  

		log10_coeffs[s] =
		(log_coeff[k][low_Te][low_ne]*(1-y) + log_coeff[k][low_Te][high_ne]*y)*(1-x)
		+(log_coeff[k][high_Te][low_ne]*(1-y) + log_coeff[k][high_Te][high_ne]*y)*x;
		
		// Print inspection
		// cout << endl;
		// cout << "00 - log Te: " << log_temperature[low_Te] << " log ne: " << log_density[low_ne] << " value: " << log_coeff[k][low_Te][low_ne] << endl;
		// cout << "10 - log Te: " << log_temperature[high_Te] << " log ne: " << log_density[low_ne] << " value: " << log_coeff[k][high_Te][low_ne] << endl;
		// cout << "01 - log Te: " << log_temperature[low_Te] << " log ne: " << log_density[high_ne] << " value: " << log_coeff[k][low_Te][high_ne] << endl;
		// cout << "11 - log Te: " << log_temperature[high_Te] << " log ne: " << log_density[high_ne] << " value: " << log_coeff[k][high_Te][high_ne] << endl;
		// cout << "xy - log Te: " << log10_Te_search << " log ne: " << log10_ne_search << " value: " << log10_coeffs[s] << endl;
		// cout << "x: " << x << " y: " << y << endl;
		// cout << endl;
	}

	// # Need to convert both temp and density to log-scale
	// log_temperature = np.log10(Te)
	// log_density = np.log10(ne)

	// # Find the logarithm of the rate-coefficient
	// log_coeff = self.splines[k](log_temperature, log_density, grid=False)
	// # Raise (piecewise) to the power 10 to return in m3/s
	// coeffs = np.power(10,log_coeff)

	// return coeffs
};
double RateCoefficient::call0D(const int k, const double log10_Te, const double log10_ne){
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
	// double log10_Te;
	// double log10_ne;
	// Look through the log_temperature and log_density attributes of RateCoefficient to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int low_Te = lower_bound(log_temperature.begin(), log_temperature.end(), log10_Te) - log_temperature.begin() - 1;
	int low_ne = lower_bound(log_density.begin(), log_density.end(), log10_ne) - log_density.begin() - 1;

	int high_Te = low_Te + 1;
	int high_ne = low_ne + 1;
	double Te_norm = 1/(log_temperature[high_Te] - log_temperature[low_Te]); //Spacing between grid points
	double ne_norm = 1/(log_density[high_ne] - log_density[low_ne]); //Spacing between grid points

	double x = (log10_Te - log_temperature[low_Te])*Te_norm;
	double y = (log10_ne - log_density[low_ne])*ne_norm;
	
	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    ne -> y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--Te -> x
	// //  | /     \ |      |
	// // w00 ------ w10  

	double log10_coeff =
	(log_coeff[k][low_Te][low_ne]*(1-y) + log_coeff[k][low_Te][high_ne]*y)*(1-x)
	+(log_coeff[k][high_Te][low_ne]*(1-y) + log_coeff[k][high_Te][high_ne]*y)*x;
	
	// Print inspection
	// cout << endl;
	// cout << "00 - log Te: " << log_temperature[low_Te] << " log ne: " << log_density[low_ne] << " value: " << log_coeff[k][low_Te][low_ne] << endl;
	// cout << "10 - log Te: " << log_temperature[high_Te] << " log ne: " << log_density[low_ne] << " value: " << log_coeff[k][high_Te][low_ne] << endl;
	// cout << "01 - log Te: " << log_temperature[low_Te] << " log ne: " << log_density[high_ne] << " value: " << log_coeff[k][low_Te][high_ne] << endl;
	// cout << "11 - log Te: " << log_temperature[high_Te] << " log ne: " << log_density[high_ne] << " value: " << log_coeff[k][high_Te][high_ne] << endl;
	// cout << "xy - log Te: " << log10_Te << " log ne: " << log10_ne << " value: " << log10_coeff << endl;
	// cout << "x: " << x << " y: " << y << endl;
	// cout << endl;
	
	return log10_coeff;

	// # Need to convert both temp and density to log-scale
	// log_temperature = np.log10(Te)
	// log_density = np.log10(ne)

	// # Find the logarithm of the rate-coefficient
	// log_coeff = self.splines[k](log_temperature, log_density, grid=False)
	// # Raise (piecewise) to the power 10 to return in m3/s
	// coeffs = np.power(10,log_coeff)

	// return coeffs
};
int RateCoefficient::get_atomic_number(){
	return atomic_number;
};
string RateCoefficient::get_element(){
	return element;
};
string RateCoefficient::get_adf11_file(){
	return adf11_file;
};
vector<vector< vector<double> > > RateCoefficient::get_log_coeff(){
	return log_coeff;
};
vector<double> RateCoefficient::get_log_temperature(){
	return log_temperature;
};
vector<double> RateCoefficient::get_log_density(){
	return log_density;
};