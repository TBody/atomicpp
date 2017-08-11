# atomic++

## Contents

1. [Acknowledgements](#acknowledgements)
2. [Introduction](#introduction)
3. 

## Acknowledgements

The majority of this code is based on the [_atomic1D_](https://github.com/TBody/atomic1D), which is in turn based on the excellent OpenADAS analysis tool provided at [_cfe316/atomic_](https://github.com/cfe316/atomic). 

The _atomic1D_ code translated the _atomic_ code into Python3 and used OpenADAS JSON reads rather than directly interfacing with the Fortran helper functions. The code to generate these OpenADAS JSON files is provided as [_OpenADAS\_to\_JSON_](https://github.com/TBody/OpenADAS_to_JSON).

## Introduction

This project provides a `C++` library for the analysis of atomic processes in a SOL/divertor plasma. It is primarily intended as a library for integration into the BOUT++ [_SD1D_](https://github.com/boutproject/SD1D) project, although the testing suites provided may be generally useful.

The project consists of three main parts:

1. The [`atomicpp`](https://github.com/TBody/atomicpp/tree/master/atomicpp)library -- provides all analysis code;
2. `Prad.cpp` testing suite in `C++` -- tests basic functionality ;
3. `Prad.py` verification suite in `Python3` -- uses a Python wrapper to the `atomicpp` library to check the results from the library.

The project is built and run using a `Makefile` architecture


This code translates the [atomic1D] code into C++ so that it may be easily integrated into the . It provides a subset of the functionality of the [atomic1D] code, and therefore it is recommended that [_atomic1D_](https://github.com/TBody/atomic1D) is used in preference to this code for standalone analysis.

This program is built using a `Makefile` which controls compilation and linking of the required source files. The core development code is supplied in the main project directory as `Prad.cpp`, while the `atomicpp` folder contains contains the required header (`.hpp`) and corresponding source (`.cpp`) files which provide broadly applicable functions for storing and processing OpenADAS data. To run the `Prad.cpp` code requires

#### System requirements  

* A C++ compiler which supports C++11 standard.
* `gmake` to run the `Makefile`.

#### In the same file as the main program source (`Prad.cpp`) at time of compilation;  

* The `atomicpp` module folder containing
    - `ImpuritySpecies` (`.cpp` and `.hpp`): for storing the data of a specified plasma impurity.
    - `RateCoefficient` (`.cpp` and `.hpp`): called by `ImpuritySpecies`, for providing a callable interface to a set of OpenADAS density- and temperature-resolved rate coefficients.
    - `sharedFunctions` (`.cpp` and `.hpp`): for functions which are useful to all other modules - currently a JSON reader and a file checker.
    - `json.hpp`: header-only JSON reader from the [_nlohmann::json_](https://github.com/nlohmann/json) project. A copy of the header file is included, although it is recommended that users download an up-to-date version [here](https://github.com/nlohmann/json/blob/develop/src/json.hpp).
    - `SD1DData` (`.cpp` and `.hpp`): not required for SD1D integration. For storing and processing data from an SD1D `collect`, which is then saved as a json file. Run the `atomic1D/reference/data_dict_export.py` program in the SD1D output folder to produce the json file, and then copy this to the `atomic++` directory (where `Prad.cpp` is being executed). This data is used to train and verify the `atomic++` code.
* `impurity_user_input.json` (same directory as `Prad.cpp`): a plain-text file for providing hard-coded data relating to the impurity (such as year for OpenADAS data, etc). To modify, see the [JSON project](http://www.json.org) for nomenclature and follow the same style as the template.
* `sd1d-case-*.json` (same directory as `Prad.cpp`): a json file created by `atomic1D/reference/data_dict_export.py`, for use with `SD1DData` for training the program (not needed if integrating into SD1D).

### Understanding the code  
Effort was made to ensure reasonably comprehensive code-commenting in the source -- however, this README will not describe the code function. Instead, see the README for `atomic1D` (which has almost identical functionality) since this will be more thoroughly commented.

## Program execution  

#### Running `Prad.cpp` from terminal  

To run the main power calculation code change to the directory containing `Prad.cpp` (and `Makefile`). In a terminal window type `make`. The makefile will compile the object code (`*.o`) for the module and for the main source (`Prad.o`), link the executables and then execute `./Prad` to run the executable. Compilation options can be controlled via the `flags` make variable in `Makefile`.

#### Integrating the power function into SD1D  

*N.b. the SD1DData (.hpp and .cpp) files of the atomicpp module directory are not required*

To compute the total radiated power from a particular plasma impurity at a single location and single time the required code is
```cpp

/**Constructs the ImpuritySpecies object corresponding to the supplied impurity_symbol. N.b. this symbol corresponds to a key for the impurity_user_input.json file which records hard-coded OpenADAS parameters.
*/
string impurity_symbol="<<impurity atomic symbol, i.e. C>>";
ImpuritySpecies impurity(impurity_symbol);

/**
 * @brief Calculate the total radiated power assuming collisional-radiative equilibrium
 * @details Assumes collisional-radiative equilibrium (i.e. infinite impurity retention time,  
 no charge-exchange recombination) to provide a simple method for calculating the ionisation  
 stage distribution for the impurity. This is then used to calculate the power due to each  
 considered physics process (line power, continuum power and charge-exchange power) for each
 ionisation stage. The total power (in W/m^3) is returned.
 * 
 * @param impurity ImpuritySpecies object, which contains OpenADAS data on relevant atomic-physics rate-coefficients
 * @param Te electron temperature in eV
 * @param Ne electron density in m^-3
 * @param Nz impurity density in m^-3, summed over all ionisation stages
 * @param Nn neutral density in m^-3
 * @return Total power in W/m^3
 */
double computeRadiatedPower(ImpuritySpecies& impurity, double Te, double Ne, double Nz, double Nn){
    /**Calculates the relative distribution across ionisation stages of the impurity by  
    assuming collisional-radiative equilibrium. This is then used to calculate the density within  
    each state, allowing the total power at a point to be evaluated.
    */

    ...

    return total_power; //In W/m^3
}
```

#### Numerical evaluation of population and energy equations for SD1D  

*N.b. the SD1DData (.hpp and .cpp) files of the atomicpp module directory are not required*

To compute the derivatives of state-vector quantities required to self-consistently model the ionisation stage distribution and radiation emission from an impurity (as a function of time) the required code is

```cpp
/**
 * @brief Calculates the rate of change (input units per second) for plasma parameters due to OpenADAS atomic physics processes
 * @details Still under development
 * 
 * * @param impurity ImpuritySpecies object, which contains OpenADAS data on relevant atomic-physics rate-coefficients
 * @param Te electron temperature in eV
 * @param Ne electron density in m^-3
 * @param Nn neutral density in m^-3
 * @param Nzk impurity density in m^-3, std::vector of densities of the form [Nz^0, Nz^1+, Nz^2+, ..., Nz^Z+]
 * @param Nthres threshold density for impurity stages, below which the time evolution of this stage is ignored. Default is 1e9,
 * although it is recommended that a time-step dependance be added in the calling code.
 * return dydt;
 * //where the derivative std::vector may be unpacked as
 *   double Pcool = dydt[0]; //Electron-cooling power - rate at which energy is lost from the electron population - in W/m^3
 *   double Prad  = dydt[1]; //Radiated power - rate at which energy is dissipated as radiation (for diagnostics) - in W/m^3
 *   std::vector<double> dNzk(impurity.get_atomic_number()+1); //Density change for each ionisation stage of the impurity - in 1/(m^3 s)
 *   for(int k=0; k<=impurity.get_atomic_number(); ++k){
 *      int dydt_index = k + 2;
 *      dNzk[k] = dydt[dydt_index];
 *   }
 *   double dNe   = dydt[(impurity.get_atomic_number()+2) + 1]; //Density change for electrons due to impurity-atomic processes (perturbation) - in 1/(m^3 s)
 *   double dNn   = dydt[(impurity.get_atomic_number()+2) + 2]; //Density change for neutrals due to impurity-atomic processes (perturbation) - in 1/(m^3 s)
 */
std::vector<double> computeDerivs(ImpuritySpecies& impurity, const double Te, const double Ne, const double Nn, const std::vector<double>& Nzk, const double Nthres = 1e9){

    ...

    return dydt; //such that the following derivatives may be unpacked from the solver
}

//Electron-cooling power - rate at which energy is lost from the electron population - in W/m^3
double Pcool = dydt[0];
//Radiated power - rate at which energy is dissipated as radiation (for diagnostics) - in W/m^3
double Prad  = dydt[1];
//Density change for each ionisation stage of the impurity - in 1/(m^3 s)
std::vector<double> dNzk(impurity.get_atomic_number()+1);
for(int k=0; k<=impurity.get_atomic_number(); ++k){
    int dydt_index = k + 2;
    dNzk[k] = dydt[dydt_index];
}
//Density change for electrons due to impurity-atomic processes (perturbation) - in 1/(m^3 s)
double dNe   = dydt[(impurity.get_atomic_number()+2) + 1];
//Density change for neutrals due to charge-exchange (perturbation) - in 1/(m^3 s)
double dNn   = dydt[(impurity.get_atomic_number()+2) + 2];
```

This function relies on the following sub-functions; to calculate the scaling factors for interpolation (since, if the grids are identical, this may be shared by all the rate-coefficient calls to interpolation)
```cpp
/**
 * @brief find the lower-bound gridpoint and fraction within the grid for the given point at which to interpolate
 * @details Using bilinear interpolation, the scaling factors for interpolating the rate coefficients are the same
 * regardless of which process is called (since the underlying log_temperature and log_density grids are the same).
 * Therefore, the grid-point and fraction pair may be shared by any rate-coefficient. Have overloaded call0D such that,
 * if a <int, double> pair is supplied as an argument then the shared interpolation method will be called
 * 
 * @param log_grid grid-points for which data is given
 * @param eval point at which the interpolation should be performed
 * 
 * @return <int, double> pair where int is the lower-bound grid-point and fraction is the scaling factor (fractional distance
 * between the lower and upper-bound gridpoints)
 */
std::pair<int, double> findSharedInterpolation(const std::vector<double>& log_grid, const double eval){
    // Perform a basic interpolation based on linear distance
    // values to search for
    
    ...

    std::pair<int, double> interp_pair(interp_gridpoint, interp_fraction);
    return interp_pair;
}
```

and; to add the elements of a list with significantly varying orders of magnitude to very high precision (avoids floating point rounding error)

```cpp
/**
 * @brief Uses Neumaier algorithm to add the elements of a list
 * @details Extension on Kahan summation algorithm for an unsorted list
 * Uses a compensated sum to improve precision when summing numbers of 
 * significantly different magnitude
 * 
 * @param list_to_sum The list of numbers to sum
 * @return neumaier_pair the uncompensated sum and the compensation
 * Compensated sum is sum + correction - however, this is left external
 * in case summation is restarted
 */
std::pair<double, double> neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction = 0.0){
    
    ...

    std::pair<double, double> neumaier_pair(sum, correction);

    return neumaier_pair;
}
```







