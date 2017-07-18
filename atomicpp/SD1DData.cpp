#include <string>
#include <vector>
#include <fstream>
#include "json.hpp"
#include "SD1DData.hpp"
using namespace std; //saves having to prepend std:: onto common functions

// for convenience
using json = nlohmann::json;

SD1DData::SD1DData(string input_file){
};
void SD1DData::setImpurityFraction(float impurity_fraction){
};
void SD1DData::setImpurityDensity(float impurity_density){
};
void SD1DData::selectSingleTime(float t){
};