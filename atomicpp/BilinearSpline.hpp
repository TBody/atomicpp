#ifndef BILINEARSPLINE_H
#define BILINEARSPLINE_H
	namespace atomicpp{
	class BilinearSpline{
	public:
		BilinearSpline(); //Default constructor
		BilinearSpline(
			std::vector<double>& _temp_values,
			std::vector<double>& _dens_values,
			std::vector<std::vector< std::vector<double> > >& _coef_values
			);
		double call0D(const int k, const double eval_temp, const double eval_dens);
		double call0D_shared(const int k, const std::pair<int, double> temp_interp, const std::pair<int, double> dens_interp);
		std::vector<std::vector< std::vector<double> > > get_coef_values();
		std::vector<double> get_temp_values();
		std::vector<double> get_dens_values();
		void set_temp_values(std::vector<double>& _temp_values);
		void set_dens_values(std::vector<double>& _dens_values);
	private:
		std::vector<double> temp_values;
		std::vector<double> dens_values;
		std::vector<std::vector< std::vector<double> > > coef_values;
	};
	}
#endif