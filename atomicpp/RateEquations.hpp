#ifndef RATEEEQUATIONS_H //Preprocessor directives to prevent multiple definitions
#define RATEEEQUATIONS_H
	
	// Include declarations
	#include <ostream>
	#include <map>
	#include <string>
	#include <vector>
	#include <tuple>
	#include <fstream>
	#include <set>
	#include <stdexcept> //For error-throwing
	#include <memory> //For smart pointers
	#include <cstdio> //For print formatting (printf, fprintf, sprintf, snprintf)
	#include <cmath> //For abs

	#include <limits>
	#include <cstdint>
	#include <cinttypes>

	#include "ImpuritySpecies.hpp"
	#include "RateCoefficient.hpp"

	std::pair<int, double> findSharedInterpolation(const std::vector<double>& log_grid, const double eval);
	std::pair<double, double> neumaierSum(const std::vector<double>& list_to_sum, const double previous_correction = 0.0);
	void calculate_ElectronImpact_PopulationEquation(
		ImpuritySpecies& impurity,
		const int Z,
		const double mz,
		const double eV_to_J,
		const double Nthres,
		const double Ne,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk,
		const std::pair<int, double>& Te_interp,
		const std::pair<int, double>& Ne_interp,
		std::vector<double>& dNzk,
		std::vector<double>& F_zk,
		std::vector<double>& dNzk_correction,
		std::vector<double>& F_zk_correction,
		double dNe,
		double Pcool
		);
	std::tuple<double, double, std::vector<double>, std::vector<double>, double, double, double, double > computeDerivs(
		ImpuritySpecies& impurity,
		const double Te,
		const double Ne,
		const double Vi,
		const double Nn,
		const double Vn,
		const std::vector<double>& Nzk,
		const std::vector<double>& Vzk,
		const double Nthres = 1e9);
#endif