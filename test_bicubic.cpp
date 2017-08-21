#include "atomicpp/BivariateBSpline.hpp"
#include <vector>
#include <string>

#include <stdio.h>
#include <iostream>
#include <cstring>

// g++ atomicpp/BicubicSpline.cpp test_bicubic.cpp -g -O0 -Wall -fno-inline -std=c++11 -o test_bicubic && ./test_bicubic

using namespace atomicpp;

int main(){
	std::vector<double> x_values = {-2.,-1.44444444,-0.88888889,-0.33333333,0.22222222,0.77777778,1.33333333,1.88888889,2.44444444,3.}; 
	std::vector<double> y_values = {-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.,2.5,3.}; 
	std::vector< std::vector<double>> z_values = {
		{  3.70000000e+01,   3.28641975e+01,   2.93456790e+01,   2.64444444e+01, 2.41604938e+01,   2.24938272e+01,   2.14444444e+01,   2.10123457e+01, 2.11975309e+01,   2.20000000e+01}, 
		{  3.12500000e+01,   2.73919753e+01,   2.41512346e+01,   2.15277778e+01, 1.95216049e+01,   1.81327160e+01,   1.73611111e+01,   1.72067901e+01, 1.76697531e+01,   1.87500000e+01}, 
		{  2.60000000e+01,   2.24197531e+01,   1.94567901e+01,   1.71111111e+01, 1.53827160e+01,   1.42716049e+01,   1.37777778e+01,   1.39012346e+01, 1.46419753e+01,   1.60000000e+01}, 
		{  2.12500000e+01,   1.79475309e+01,   1.52623457e+01,   1.31944444e+01, 1.17438272e+01,   1.09104938e+01,   1.06944444e+01,   1.10956790e+01, 1.21141975e+01,   1.37500000e+01}, 
		{  1.70000000e+01,   1.39753086e+01,   1.15679012e+01,   9.77777778e+00, 8.60493827e+00,   8.04938272e+00,   8.11111111e+00,   8.79012346e+00, 1.00864198e+01,   1.20000000e+01}, 
		{  1.32500000e+01,   1.05030864e+01,   8.37345679e+00,   6.86111111e+00, 5.96604938e+00,   5.68827160e+00,   6.02777778e+00,   6.98456790e+00, 8.55864198e+00,   1.07500000e+01}, 
		{  1.00000000e+01,   7.53086420e+00,   5.67901235e+00,   4.44444444e+00, 3.82716049e+00,   3.82716049e+00,   4.44444444e+00,   5.67901235e+00, 7.53086420e+00,   1.00000000e+01}, 
		{  7.25000000e+00,   5.05864198e+00,   3.48456790e+00,   2.52777778e+00, 2.18827160e+00,   2.46604938e+00,   3.36111111e+00,   4.87345679e+00, 7.00308642e+00,   9.75000000e+00}, 
		{  5.00000000e+00,   3.08641975e+00,   1.79012346e+00,   1.11111111e+00, 1.04938272e+00,   1.60493827e+00,   2.77777778e+00,   4.56790123e+00, 6.97530864e+00,   1.00000000e+01}, 
		{  3.25000000e+00,   1.61419753e+00,   5.95679012e-01,   1.94444444e-01, 4.10493827e-01,   1.24382716e+00,   2.69444444e+00,   4.76234568e+00, 7.44753086e+00,   1.07500000e+01}, 
		{  2.00000000e+00,   6.41975309e-01,  -9.87654321e-02,  -2.22222222e-01, 2.71604938e-01,   1.38271605e+00,   3.11111111e+00,   5.45679012e+00, 8.41975309e+00,   1.20000000e+01}, 
		{  1.25000000e+00,   1.69753086e-01,  -2.93209877e-01,  -1.38888889e-01, 6.32716049e-01,   2.02160494e+00,   4.02777778e+00,   6.65123457e+00, 9.89197531e+00,   1.37500000e+01}, 
		{  1.00000000e+00,   1.97530864e-01,   1.23456790e-02,   4.44444444e-01, 1.49382716e+00,   3.16049383e+00,   5.44444444e+00,   8.34567901e+00, 1.18641975e+01,   1.60000000e+01}, 
		{  1.25000000e+00,   7.25308642e-01,   8.17901235e-01,   1.52777778e+00, 2.85493827e+00,   4.79938272e+00,   7.36111111e+00,   1.05401235e+01, 1.43364198e+01,   1.87500000e+01}, 
		{  2.00000000e+00,   1.75308642e+00,   2.12345679e+00,   3.11111111e+00, 4.71604938e+00,   6.93827160e+00,   9.77777778e+00,   1.32345679e+01, 1.73086420e+01,   2.20000000e+01} 
	};
	// std::vector<double> x_values = {1, 2, 3, 4, 5};
	// std::vector<double> y_values = {1, 2, 3, 4, 5};
	// std::vector< std::vector<double>> z_values = {
	// 	{1, 2, 3, 4, 5},
	// 	{1, 2, 3, 4, 5},
	// 	{1, 2, 3, 4, 5},
	// 	{1, 2, 3, 4, 5},
	// 	{1, 2, 3, 4, 5}
	// };

	BivariateBSpline interpolator(x_values, y_values, z_values);

	int i = 1;
	int j = 2;
	double x = x_values[i];
	double y = y_values[j]+0.5;	
	std::printf("(%f,%f) = %f \n", x, y, interpolator.call0D(x, y));


};