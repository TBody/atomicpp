#include "atomicpp/BicubicSpline.hpp"
#include <vector>
#include <string>

#include <stdio.h>
#include <iostream>
#include <cstring>

// g++ atomicpp/BicubicSpline.cpp test_bicubic.cpp -g -O0 -Wall -fno-inline -std=c++11 -o test_bicubic && ./test_bicubic

using namespace atomicpp;

int main(){
	std::vector<double> x_values = {0,1,2,3};
	std::vector<double> y_values = {0,1,2,3};
	// std::vector< std::vector<double>> z_values = {
	// 	{0,1,2,3},
	// 	{4,5,6,7},
	// 	{8,9,10,11},
	// 	{12,13,14,15},
	// };
	std::vector< std::vector<double>> z_values = {
		{0,0,0,0},
		{0,0,10,0},
		{0,0,0,0},
		{0,0,0,0},
	};

	BicubicSpline interpolator(x_values, y_values, z_values);

	int i = 1;
	int j = 2;
	double x = x_values[i];
	double y = y_values[j];	
	std::printf("(%f,%f) = %f \n", x, y, interpolator.call0D(x, y));


};