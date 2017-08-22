#ifndef BILINEARSPLINE_H
#define BILINEARSPLINE_H
	#include <vector>

	namespace atomicpp{
	class BilinearSpline{
	public:
		BilinearSpline(); //Default constructor
		BilinearSpline(
			std::vector<double>& _x_values,
			std::vector<double>& _y_values,
			std::vector< std::vector<double> >& _z_values
			);
		double call0D(const double eval_x, const double eval_y);
		double call0D_shared(const std::pair<int, double> x_interp, const std::pair<int, double> y_interp);
		std::vector< std::vector<double> > get_z_values();
		std::vector<double> get_x_values();
		std::vector<double> get_y_values();
		void set_x_values(std::vector<double>& _x_values);
		void set_y_values(std::vector<double>& _y_values);
		void zero_z_values();
	private:
		std::vector<double> x_values;
		std::vector<double> y_values;
		std::vector< std::vector<double> > z_values;
	};
	}
#endif