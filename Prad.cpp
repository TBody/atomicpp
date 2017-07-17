// Program name: atomic++/Prad.cpp
// Author: Thomas Body
// Author email: tajb500@york.ac.uk
// Date of creation: 17 July 2017
// 
// Program function: output the radiated power (Prad)
//                   by using OpenADAS rates on output JSON from SD1D run
//                   
// Based on the TBody/atomic1D code, which is in turn based on the cfe316/atomic code
// 
// Under active development: <<TODO>> indicates development goal

// Include declarations
#include <iostream>
#include <map>
#include <string>
using namespace std; //saves having to prepend std:: onto common functions

// N.b. perform all plotting routines as Python post-processing

map<string,string>datatype_abbrevs={
	{"ionisation",           "scd"},
	{"recombination",        "acd"},
	{"cx_recc",              "ccd"},
	{"continuum_power",      "prb"},
	{"line_power",           "plt"},
	{"cx_power",             "prc"},
	{"ionisation_potential", "ecd"} //N.b. ionisation_potential is not a rate-coefficient, but most of the methods are transferable
};
string input_file="sd1d-case-05.json";
string json_database_path="json_database/json_data";
string impurity_symbol="C";

class ImpuritySpecies
{
public:
	ImpuritySpecies(){
		cout << "Constructing ImpuritySpecies object for " << impurity_symbol << "\n";
	};
};

int main(){

	cout<<"Hello world\n";
	
	ImpuritySpecies();


	cout<<datatype_abbrevs["cx_power"]<<"\n";
}





























