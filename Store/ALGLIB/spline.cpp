// #include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

using namespace alglib;


int main(int argc, char **argv)
{
    //
    // We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled 
    // at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
    //
    real_1d_array x = "[0.0, 0.5, 1.0]";
    real_1d_array y = "[0.0, 1.0]";
    real_1d_array f = "[0.00,0.25,1.00,2.00,2.25,3.00]";
    double vx = 0.25;
    double vy = 0.50;
    double v;
    double dx;
    double dy;
    double dxy;
    spline2dinterpolant s;

    // build spline
    spline2dbuildbicubicv(x, 3, y, 2, f, 1, s);

    // calculate S(0.25,0.50)
    v = spline2dcalc(s, vx, vy);
    printf("%.4f\n", double(v)); // EXPECTED: 1.0625

    // calculate derivatives
    spline2ddiff(s, vx, vy, v, dx, dy, dxy);
    printf("%.4f\n", double(v)); // EXPECTED: 1.0625
    printf("%.4f\n", double(dx)); // EXPECTED: 0.5000
    printf("%.4f\n", double(dy)); // EXPECTED: 2.0000
    return 0;
}