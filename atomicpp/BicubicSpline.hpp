#ifndef BICUBICSPLINE_H
#define BICUBICSPLINE_H
	#include <vector>
	#include <array>
	
	namespace atomicpp{
	typedef std::array<std::array<double, 4>, 4> grid_matrix;
	grid_matrix default_alpha_coeff = {{
			{0.0, 0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0, 0.0},
	}};
	
	struct interp_data{
		// std::pair<double, double> coord; //(T,N) coordinate of point
		double f = 0.0; //Value at point
		double fdT = 0.0; //Derivative in temperature axis
		double fdN = 0.0; //Derivative in density axis
		double fdTdN = 0.0; //Cross derivative
	} default_interp_data;

	class BicubicSpline{
		public:
			BicubicSpline(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff);
			/**
			 * @brief Returns the rate coefficient for a (scalar) Te and Ne supplied
			 * @details Performs a simple bivariate (multilinear) interpolation to return the rate coefficient
			 * at the supplied Te and Ne values. N.b. will throw a std::runtime_error if the supplied Te or Ne
			 * value are not on the interpolating grid (otherwise you'll get a seg fault)
			 * 
			 * @param k The charge state index process (actually k+=1 for charged-target processes, but we don't implement this here)
			 * @param eval_Te electron temperature (Te) at a point (in eV)
			 * @param eval_Ne electron density (Ne) at a point (in m^-3)
			 * @return eval_coeff evaluated rate coefficient in m^3/s
			 */
			double call0D(const int k, const double eval_Te, const double eval_Ne);
			//Overloaded onto call0D - if the input is an int and two <int, double> pairs then use the SharedInterpolation method
			//(i.e. assume that Te_interp and Ne_interp contain which point for which to return the coefficient - saves reevaluating)
			double call0D(const int k, const std::pair<int, double> Te_interp, const std::pair<int, double> Ne_interp);
			std::vector<std::vector<std::vector<interp_data>>> calculate_grid_coeff(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff);
			std::vector<std::vector<std::vector<grid_matrix>>> calculate_alpha_coeff(std::vector<double>& log_temperature, std::vector<double>& log_density, std::vector<std::vector< std::vector<double> > >& log_coeff);
		private:
			std::vector<std::vector< std::vector<double> > > log_coeff;
			std::vector<double> log_temperature;
			std::vector<double> log_density;
			std::vector<std::vector<std::vector<grid_matrix>>> alpha_coeff;
	};
	}
#endif