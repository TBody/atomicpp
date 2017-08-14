#ifndef BICUBICSPLINE_H
#define BICUBICSPLINE_H
	namespace atomicpp{
	
	typedef std::array<std::array<double, 4>, 4> grid_matrix;

	struct interp_data{
		// std::pair<double, double> coord; //(T,N) coordinate of point
		double f = 0.0; //Value at point
		double fdT = 0.0; //Derivative in temperature axis
		double fdN = 0.0; //Derivative in density axis
		double fdTdN = 0.0; //Cross derivative
	};

	class BicubicSpline{
	public:
		BicubicSpline(); //Default constructor
		BicubicSpline(
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