# atomic++

## Contents

1. [Acknowledgements](#acknowledgements)
2. [Introduction](#introduction)
3. [System requirements](#system-requirements)
4. [Quick-start](#quick-start)
5. [Modifying the impurity database](#modifying-the-impurity-database)
6. [SD1D Integration](#sd1d-integration)

## Acknowledgements

The majority of this code is based on the [_atomic1D_](https://github.com/TBody/atomic1D) code, which is in turn based on the excellent OpenADAS analysis tool provided at [_cfe316/atomic_](https://github.com/cfe316/atomic). 

The data used for modelling is entirely supplied by the [OpenADAS project](http://open.adas.ac.uk/), a public-access database managed by the University of Strathclyde.

## Introduction

This project provides a `C++` library for the analysis of atomic processes in a SOL/divertor plasma. It is primarily intended as a library for integration into the BOUT++ [_SD1D_](https://github.com/boutproject/SD1D) project, although the testing suites provided may be generally useful.

The project consists of three main parts:

1. The [`atomicpp`](https://github.com/TBody/atomicpp/tree/master/atomicpp)library -- provides all analysis code as well the Cython `.pyx` wrapper;
2. `Prad.cpp` testing suite in `C++` -- tests basic functionality ;
3. `Prad.py` verification suite in `Python3` -- uses a Python wrapper to the `atomicpp` library to check the results from the library.

The project is built and run using a `Makefile` architecture. All commands are localised to the project directory -- to uninstall simply delete the project folder.

## System requirements  
It is recommended that a package manager (`apt-get` (Linux) or `MacPorts` (Mac), etc) be used to install the required packages.

For the core `atomicpp` library and the `Prad.cpp` C++ test suite
* A C++ compiler which supports C++11 standard (tested with [GNU](https://gcc.gnu.org/) compiler `gcc/g++ version 6.3.0`)
* `gmake` to run the `Makefile` (tested with `GNU Make 3.81`)

For the `Prad.py` Python3 verification suite
* A [Python3](https://www.python.org/) installation with an installed [SciPy](https://www.scipy.org/) stack and [Cython](http://cython.org/). Both the SciPy stack and Cython are available in the [Anaconda](https://www.continuum.io/anaconda-overview) distribution (tested with `Python 3.6.1 |Anaconda 4.4.0 (x86_64)`)

To extend the `json_database` of OpenADAS rate-coefficients
* The [OpenADAS_to_JSON](https://github.com/TBody/OpenADAS_to_JSON) project. See [Modifying the impurity database](#modifying-the-impurity-database) and the OpenADAS_to_JSON `README` for more details.
* A Fortran compiler (tested with [GNU](https://gcc.gnu.org/) compiler `gcc/gfortran version 6.3.0`)

#### In the same file as the main program source (`Prad.cpp`) at time of compilation;  

* The `atomicpp` module folder containing
    - `RateEquations` (`.cpp` and `.hpp`): for evaluation of the population, power and momentum balance derivatives.
    - `ImpuritySpecies` (`.cpp` and `.hpp`): Called by `RateEquations`, for storing the data of a specified plasma impurity.
    - `RateCoefficient` (`.cpp` and `.hpp`): called by `ImpuritySpecies`, for providing a callable interface to a set of OpenADAS density- and temperature-resolved rate coefficients.
    - `sharedFunctions` (`.cpp` and `.hpp`): for functions which are useful to all other modules - currently a JSON reader and a file checker.
    - `json.hpp`: header-only JSON reader from the [_nlohmann::json_](https://github.com/nlohmann/json) project. A copy of the header file is included, although it is recommended that users download an up-to-date version [here](https://github.com/nlohmann/json/blob/develop/src/json.hpp).
    - `SD1DData` (`.cpp` and `.hpp`): not required for SD1D integration. For storing and processing data from an SD1D `collect`, which is then saved as a json file. Run the `atomic1D/reference/data_dict_export.py` program in the SD1D output folder to produce the json file, and then copy this to the `atomic++` directory (where `Prad.cpp` is being executed). This data is used to train and verify the `atomic++` code.
* `impurity_user_input.json` (same directory as `Prad.cpp`): a plain-text file for providing hard-coded data relating to the impurity (such as year for OpenADAS data, etc). To modify, see the [JSON project](http://www.json.org) for nomenclature and follow the same style as the template.
* `sd1d-case-*.json` (same directory as `Prad.cpp`): a json file created by `atomic1D/reference/data_dict_export.py`, for use with `SD1DData` for training the program (not needed if integrating into SD1D).

## Program execution

#### Quick-start

The `Makefile` provides most of the desired functionality of the project. It has 5 main commands;
* `make cpp`: checks if the required source files have been built since they were last modified -- if not, compiles the `atomicpp` library and the `Prad.cpp` program into `.o` files and then links the `.o` files into an executable (`Prad`).
* `make cpp_run`: (default for `make`) performs the `make cpp` functionality and also runs the `Prad` executable via `./Prad`
* `make py`: checks if the required source files have been built since they were last modified -- if not, compiles the `atomicpp` library and generates a Python module `atomicpy` which is declared in `atomicpy.pyx` (Cython) and built with `setup.py build_ext --inplace`.
* `make py_run`: performs the `make py` functionality and also runs the `Prad.py` script via `python Prad.py`
* `make clean`: reverts the project to a fresh install state.

#### Modifying the impurity database
A separate project is supplied at [OpenADAS_to_JSON](https://github.com/TBody/OpenADAS_to_JSON). This project downloads ADF11 files for the specified impurity from [OpenADAS](http://open.adas.ac.uk/), uses Fortran helper functions to read the data and exports the data as JSON files. To change which species is being considered you'll need to do the following;
1. Modify the [`elements`](https://github.com/TBody/OpenADAS_to_JSON/blob/master/makefile#L35) tag of the `makefile`
2. Run `make json_update`
3. Copy the `OpenADAS_to_JSON/json_database` file from that project onto `atomicpp/json_database` (overwrite)
4. Update [`impurity_user_input.json`](https://github.com/TBody/atomicpp/blob/master/impurity_user_input.json) to include the impurity data required by the program
5. Update the `impurity_symbol_supplied` variable of [Python](https://github.com/TBody/atomicpp/blob/master/Prad.py#L11) and [C++](https://github.com/TBody/atomicpp/blob/master/Prad.cpp#L93) if using the supplied testing programs.


## SD1D Integration
Sub-contents
1. [Impurity rate equations](#impurity-rate-equations)
2. [Hydrogen rate equations](#hydrogen-rate-equations)
3. [The `DerivStruct` data structure](#the-derivstruct-data-structure)
4. [Print method for the data structure](#print-method-for-the-data-structure)
5. [The `computeDerivs` function](#the-computederivs-function)
6. [Shared interpolation](#shared-interpolation)
7. [Kahan-Neumaier summation](#kahan-neumaier-summation)
---
*N.b. the SD1DData (.hpp and .cpp) files of the atomicpp module directory are not required*

The principal purpose of this code is to extend the radiation model of the [SD1D](https://github.com/boutproject/SD1D) SOL/divertor plasma simulation code, which is built on the [BOUT++](https://github.com/boutproject) project. 

#### Impurity rate equations
The derivative evaluator is provided in the `RateEquations` class. This class is initialised from a `ImpuritySpecies` object, shown here for a Carbon (`c`) impurity. The derivative evaluator is supplied as `computeDerivs` -- see the `RateEquations.hpp` header for identity and units of the arguments.
```cpp
#include "atomicpp/ImpuritySpecies.hpp"
#include "atomicpp/RateEquations.hpp"

std::string impurity_symbol="c";

atomicpp::ImpuritySpecies impurity(impurity_symbol);

atomicpp::RateEquations impurity_derivatives(impurity); //Organised as a RateEquations object for cleanliness
impurity_derivatives.setThresholdDensity(1e9); //Density threshold - ignore ionisation stages which don't have at least this density
impurity_derivatives.setDominantIonMass(1.0); //Dominant ion mass in amu, for the stopping time calculation

atomicpp::DerivStruct derivative_struct = impurity_derivatives.computeDerivs(Te, Ne, Vi, Nn, Vn, Nzk, Vzk);
```

#### Hydrogen rate equations
The 'impurity' methods for evaluation of the OpenADAS rate-coefficient data can be equally applied to the dominant ion. A modified method `computeDerivsHydrogen` is supplied to avoid having to supply the same variables as the 'impurity' and 'dominant ion' values, although computationally the method is almost identical (ion-ion drag and charge-exchange are neglected). The equivalent case is
```cpp
#include "atomicpp/ImpuritySpecies.hpp"
#include "atomicpp/RateEquations.hpp"

std::string hydrogen_symbol="h";

atomicpp::ImpuritySpecies hydrogen(hydrogen_symbol);
impurity_derivatives.setThresholdDensity(1e9); //Density threshold - ignore ionisation stages which don't have at least this density

atomicpp::DerivStruct derivative_struct_H = hydrogen_derivatives.computeDerivsHydrogen(Te, Ne, Nhk, Vhk);
```

#### The `DerivStruct` data structure
The `computeDerivs` functions return a `DerivStruct` data structure. This is defined as
```cpp
struct DerivStruct{
        double Pcool;
        double Prad;
        std::vector<double> dNzk;
        std::vector<double> F_zk;
        double dNe;
        double F_i;
        double dNn;
        double F_n;
    };
```
which may be unpacked as
```cpp
double Pcool             = derivative_struct.Pcool; //Electron-cooling power, in J m^-3 s^-1 (needed for electron power balance)
double Prad              = derivative_struct.Prad;  //Radiated power, in in J m^-3 s^-1 (for comparing to diagnostic signal)
std::vector<double> dNzk = derivative_struct.dNzk;  //Change in each ionisation stage of the impurity population, in particles m^-3 s^-1
std::vector<double> F_zk = derivative_struct.F_zk;  //Force on each particle of ionisation stage k of the impurity population, in N
double dNe               = derivative_struct.dNe;   //Perturbation change in the electron density (in particles m^-3 s^-1) and
double F_i               = derivative_struct.F_i;   //  perturbation force (in N) on the electron population due to atomic processes
double dNn               = derivative_struct.dNn;   //Perturbation change in the neutral density (in particles m^-3 s^-1) and
double F_n               = derivative_struct.F_n;   //  perturbation force (in N) on the neutral population due to atomic processes
```
#### Print method for the data structure
The `RateEquations` class also includes a printing method for the `DerivStruct` data structure
```cpp
impurity_derivatives.printDerivativeStruct(derivative_struct);
```
which returned a formatted result

#### The `computeDerivs` function
The core functionality of the `RateEquations` object is provided by the `computeDerivs` function. This is given as
```cpp
DerivStruct RateEquations::computeDerivs(
    const double Te,
    const double Ne,
    const double Vi,
    const double Nn,
    const double Vn,
    const std::vector<double>& Nzk,
    const std::vector<double>& Vzk){

    resetDerivatives(); //Reset all the derivatives to zero, since the object has a memory of the previous step

    // Perform 'sharedInterpolation' - find the lower-bound gridpoint and fraction into the grid for both Te and Ne
    std::pair<int, double> Te_interp, Ne_interp;
    if (use_shared_interpolation){
        Te_interp = findSharedInterpolation(rate_coefficients["blank"]->get_log_temperature(), Te);
        Ne_interp = findSharedInterpolation(rate_coefficients["blank"]->get_log_density(), Ne);
    } else {
        throw std::runtime_error("Non-shared interpolation method requires switching of method. Declare Te_interp and Ne_interp as doubles instead of <int, double> pairs.");
        // //Pass Te_interp and Ne_interp as doubles instead of pairs and the program will auto-switch to the full interpolation method.
        // double Te_interp = Te;
        // double Ne_interp = Ne;
    }

    calculateElectronImpactPopulationEquation(Ne, Nzk, Vzk, Te_interp, Ne_interp);

    calculateChargeExchangePopulationEquation(Nn, Nzk, Vzk, Te_interp, Ne_interp);

    //Apply neumairSum corrections
    for(int k=0; k<=Z; ++k){
        dNzk[k] += dNzk_correction[k];
        F_zk[k] += F_zk_correction[k];
    }

    verifyNeumaierSummation();

    calculateIonIonDrag(Ne, Te, Vi, Nzk, Vzk);
    
    calculateElectronImpactPowerEquation(Ne, Nzk, Te_interp, Ne_interp);

    calculateChargeExchangePowerEquation(Nn, Nzk, Te_interp, Ne_interp);

    // auto derivative_tuple = makeDerivativeTuple();
    // return derivative_tuple;
    DerivStruct derivative_struct = makeDerivativeStruct();
    return derivative_struct;
};
```

#### Shared interpolation
`computeDerivs` relies 'shared interpolation' to calculate the scaling factors for interpolation, which requires that the underlying grids are identical

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

#### Kahan-Neumaier summation
To add the elements of a list with significantly varying orders of magnitude to very high precision (avoids floating point rounding error) the Kahan-Neumaier algorithm is implemented as follows;

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







