#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
#include "RateCoefficient.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

RateCoefficient::RateCoefficient(string filename){
};
void RateCoefficient::compute_interpolating_splines(){
}; //Could consider fixing length, since it will always be the same shape
vector<double> RateCoefficient::call1D(int k, double Te, double ne){
};