// This algorithm is closely based on the FITPACK routines from http://www.netlib.org/dierckx/
// by Paul Dierckx. For more information see
// 
//   dierckx p. : a fast algorithm for smoothing data on a rectangular
//                grid while using spline functions, siam j.numer.anal.
//                19 (1982) 1286-1304.
//   dierckx p. : a fast algorithm for smoothing data on a rectangular
//                grid while using spline functions, report tw53, dept.
//                computer science,k.u.leuven, 1980.
//   dierckx p. : curve and surface fitting with splines, monographs on
//                numerical analysis, oxford university press, 1993.
//
//  author:
//    p.dierckx
//    dept. computer science, k.u. leuven
//    celestijnenlaan 200a, b-3001 heverlee, belgium.
//    e-mail : Paul.Dierckx@cs.kuleuven.ac.be

#ifndef BIVARIATEBSPLINE_H
#define BIVARIATEBSPLINE_H
	
	#include <vector>
	#include <array>
	#include <utility>

	namespace atomicpp{

	struct regrid_return{
	  int nx;
	  int ny;
	  std::vector<double> tx;
	  std::vector<double> ty;
	  std::vector<double> C;
	  double fp;
	};

	class BivariateBSpline{
	public:
		BivariateBSpline(); //Default constructor
		BivariateBSpline(
			const std::vector<double>& _x_values,
			const std::vector<double>& _y_values,
			const std::vector< std::vector<double> >& _z_values,
			const bool auto_transpose  = true,
			const bool warn_transpose  = true
			);
		double call0D(const double eval_x, const double eval_y);
		std::vector< std::vector<double> > get_z_values();
		std::vector<double> get_x_values();
		std::vector<double> get_y_values();
		void set_x_values(std::vector<double>& _x_values);
		void set_y_values(std::vector<double>& _y_values);
		void zero_z_values();

		void fprota(const double cos, const double sin, double& a, double& b);
		void fpgivs(const double piv, double& ww, double& cos, double& sin);
		void fpback(const std::vector<std::vector<double>>& a, const std::vector<double>& z, const int z_start, const int n, const int k, std::vector<double>& C, const int c_start, const int nest);
		void fpback(const std::vector<std::vector<double>>& a, const std::vector<double>& z, const int n, const int k, std::vector<double>& C, const int nest);
		void fpbspl(const std::vector<double>& t,const int n, const int k, const double x, const int l, std::vector<double>& h);
		regrid_return regrid_smth(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z);
		double fpbisp(const std::vector<double>& tx, const int& nx, const std::vector<double>& ty, const int& ny, const std::vector<double>& c, const int& kx, const int& ky, const double& x, const double& y);
		double bispeu(const std::vector<double>& tx,const std::vector<double>& ty,const std::vector<double>& c,const double& x,const double& y);

	private:
		std::vector<double> x_values;
		std::vector<double> y_values;
		std::vector< std::vector<double> > z_values;
		std::vector<double> tx;
		std::vector<double> ty;
		std::vector<double> C;
	};
	}
#endif