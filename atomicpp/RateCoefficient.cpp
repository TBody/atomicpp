#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
#include "RateCoefficient.hpp"
#include "sharedFunctions.hpp"
using namespace std; //saves having to prepend std:: onto common functions

#include <typeinfo>

// for convenience
using json = nlohmann::json;

RateCoefficient::RateCoefficient(const string& filename){
	// # Create an instance of RateCoefficient by reading an OpenADAS JSON file
	
	json data_dict = retrieveFromJSON(filename);

	atomic_number   = data_dict["charge"];
	element         = data_dict["element"];
	adf11_file      = filename;

	// auto log_coeff = data_dict["log_coeff"];
	// cout << typeid(log_coeff).name() << endl;

	// cout << "passed" << "\n";
	// cout << data_dict["log_coeff"];

	vector<vector< vector<double> > > log_coeff = data_dict["log_coeff"];
	vector<double> log_temperature = data_dict["log_temperature"];
	vector<double> log_density = data_dict["log_density"];

	// cout << "Temperature" << endl;
	// for (auto &T : log_temperature){
	//     cout << T << endl;
	//   }
	// cout << "Density" << endl;
	// for (auto &N : log_density){
	//     cout << N << endl;
	//   }
	
	// cout << "Coeffs" << endl;
	// cout << log_coeff[0][0][0] << endl;


	// from atomic1D import sharedFunctions

	// data_dict = sharedFunctions.retrieveFromJSON(filename)

	// self.atomic_number   = data_dict['charge']
	// self.element         = data_dict['element']
	// self.adf11_file      = filename
	// self.log_temperature = data_dict['log_temperature']
	// self.log_density     = data_dict['log_density']
	// self.log_coeff       = data_dict['log_coeff']

	// self._compute_interpolating_splines()
};
void RateCoefficient::compute_interpolating_splines(){
	// # Generate the interpolation functions for log_coeff
		



	// self.splines = []
	// for k in range(self.atomic_number):
	// 	x = self.log_temperature
	// 	y = self.log_density
	// 	z = self.log_coeff[k]
	// 	self.splines.append(RectBivariateSpline(x, y, z))

}; //Could consider fixing length, since it will always be the same shape
vector<double> RateCoefficient::call1D(int k, double Te, double ne){
	// """Evaluate the ionisation/recombination coefficients of
	// 	k'th atomic state at a given temperature and density.

	// 	For 1D case - either Te or ne vary while the other is fixed (i.e. parameter scan),
	// 	or both vary but are linked (i.e. analysis of 1D output)
	// 	>>Key difference to call2D is the grid=False argument - otherwise will return
	// 	  interpolations over the len(Te)*len(ne) 2D array

	// 	Args:
	// 		k  (int): Ionising or recombined ion stage,
	// 			between 0 and k=Z-1, where Z is atomic number.
	// 		Te (array_like): Temperature in [eV].
	// 		ne (array_like): Density in [m-3].

	// 	Returns:
	// 		c (array_like): Rate coefficent in [m3/s].
	// 	"""



	// # broadcast_arrays ensures that Te and ne are of the same shape
	// # if they are originally equal then they remain so, while if len(A) = L
	// # and len(B) = 1 then B' = L repeats of B
	// Te, ne = np.broadcast_arrays(Te, ne)

	// # Need to convert both temp and density to log-scale, since this is what the spline-interpolation is performed for
	// log_temperature = np.log10(Te)
	// log_density = np.log10(ne)

	// # Find the logarithm of the rate-coefficient
	// log_coeff = self.splines[k](log_temperature, log_density, grid=False)
	// # Raise (piecewise) to the power 10 to return in m3/s
	// coeffs = np.power(10,log_coeff)

	// return coeffs
};