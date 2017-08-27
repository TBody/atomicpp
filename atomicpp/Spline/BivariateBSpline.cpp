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

#include <vector>
#include <math.h>
#include <algorithm> //for upper/lower_bound
#include <stdexcept> //For error-throwing
#include <sstream>
#include "BivariateBSpline.hpp"

#include <stdio.h>
#include <string>
#include <iostream>

using namespace atomicpp;
BivariateBSpline::BivariateBSpline(){//Default constructor
};
BivariateBSpline::BivariateBSpline(
  const std::vector<double>& _x_values,
  const std::vector<double>& _y_values,
  const std::vector< std::vector<double> > & _z_values,
  const bool auto_transpose /* = true */,
  const bool warn_transpose /* = true */
  ){
  x_values = _x_values;
  y_values = _y_values;
  z_values = _z_values;

  int mx = x_values.size();
  int my = y_values.size();

  int zx_length = z_values.size();
  int zy_length = z_values[0].size();

  for(int i = 0; i < mx-1; ++i){
    if (x_values[i+1] <= x_values[i]){
      throw std::runtime_error("BivariateBSpline initialization error: x_values must be strictly increasing");
    }
  }
  for(int j = 0; j < my-1; ++j){
    if (y_values[j+1] <= y_values[j]){
      throw std::runtime_error("BivariateBSpline initialization error: y_values must be strictly increasing");
    }
  }

  std::vector<double> z_flattened (mx * my, 0.0); //Allocate an empty 1D vector for a flattened (i.e. ravel) of z

  if((mx == zx_length) and (my == zy_length)){
    // Flatten z and allocate it to z_flattened
    for(int i = 0; i < mx; ++i){
      for(int j = 0; j < my; ++j){
        z_flattened.at(my*i + j) = z_values.at(i).at(j);
      }
    }
  } else if((auto_transpose) and (mx == zy_length) and (my == zx_length)){
    // Easy to supply the transpose of z_values, since row and column indices swap
    if(warn_transpose){std::cout << "BivariateBSpline initialization warning: z_values has been transposed to have the correct number of elements in x and y dimensions\n";}
    // Flatten the transpose of z and allocate it to z_flattened
    for(int i = 0; i < mx; ++i){
      for(int j = 0; j < my; ++j){
        z_flattened.at(my*i + j) = z_values.at(j).at(i);
      }
    }
  // Transpose doesn't work - throw an error
  } else if(not(mx == zx_length)){
    throw std::runtime_error("BivariateBSpline initialization error: x dimension of z_values must have the same number of elements as x");
  } else if(not(my == zy_length)){
    throw std::runtime_error("BivariateBSpline initialization error: y dimension of z_values must have the same number of elements as y");
  } else{
    throw std::runtime_error("BivariateBSpline initialization error: unhandled exception");
  }

  // __INIT__ BBSpline

  regrid_return setup_spline = regrid_smth(x_values, y_values, z_flattened);

  tx = setup_spline.tx;
  ty = setup_spline.ty;
  C = setup_spline.C;
};
double BivariateBSpline::call0D(const double eval_x, const double eval_y){
  
  // Bounds checking -- make sure you haven't dropped off the end of the array
  if ((eval_x <= x_values.at(0)) or (eval_x >= x_values.at(x_values.size()-1))){
    // An easy error to make is supplying the function arguments already having taken the log10
    std::stringstream errMsg;
    errMsg << "BivariateBSpline call0D error: X value off grid - require (" << x_values.at(0) << " < " << eval_x << " < " << x_values.at(x_values.size()-1) << ")";
    throw std::runtime_error(errMsg.str());
  };
  if ((eval_y <= y_values.at(0)) or (eval_y >= y_values.at(y_values.size()-1))){
    // An easy error to make is supplying the function arguments already having taken the log10
    std::stringstream errMsg;
    errMsg << "BivariateBSpline call0D error: Y value off grid - require (" << y_values.at(0) << " < " << eval_y << " < " << y_values.at(y_values.size()-1) << ")";
    throw std::runtime_error(errMsg.str());
  };

  // __CALL__ BBSpline
  double eval_coeff = bispeu(tx,ty,C,eval_x,eval_y);

  return eval_coeff;
};
std::vector< std::vector<double> > BivariateBSpline::get_z_values(){
  return z_values;
};
std::vector<double> BivariateBSpline::get_x_values(){
  return x_values;
};
std::vector<double> BivariateBSpline::get_y_values(){
  return y_values;
};
void BivariateBSpline::set_x_values(std::vector<double>& _x_values){
  x_values = _x_values;
};
void BivariateBSpline::set_y_values(std::vector<double>& _y_values){
  y_values = _y_values;
};
void BivariateBSpline::zero_z_values(){
  for(int i = 0; i<(int)(z_values.size()); ++i){
    for(int j = 0; j<(int)(z_values.at(0).size()); ++j){
      z_values.at(i).at(j) = 0.0;
    }
  }
};

void BivariateBSpline::fprota(const double cos, const double sin, double& a, double& b){
  //  subroutine fprota applies a given rotation to a and b.
  double stor1 = a;
  double stor2 = b;

  b = cos * stor2 + sin * stor1;
  a = cos * stor1 - sin * stor2;
};

void BivariateBSpline::fpgivs(const double piv, double& ww, double& cos, double& sin){
  //  subroutine fpgivs calculates the parameters of a given
  //  transformation .

  double dd;
  if(abs(piv)>=ww){
    dd = abs(piv)*sqrt(1+(ww/piv)*(ww/piv));
  } else {
    dd = ww*sqrt(1+(piv/ww)*(piv/ww));
  }
  cos = ww/dd;
  sin = piv/dd;
  ww = dd;
};

void BivariateBSpline::fpback(const std::vector<std::vector<double>>& a, const std::vector<double>& z, const int z_start, const int n, const int k, std::vector<double>& C, const int c_start, const int nest){
  // Switch from Fortran to C++ indexing
  // Supply start as '1' to not apply a shift
  int z_shifted = z_start - 1;
  int c_shifted = c_start - 1;

  C.at(c_shifted + n-1) = z.at(z_shifted + n-1)/a.at(n-1).at(0);
  
  int i = n-1;
  double store;

  if(i!=0){
    for(int j=2; j <= n; ++j){
      store = z.at(z_shifted + i-1);
      int i1 = k-1;

      if(j <= k-1){
        i1 = j-1;
      }
      int m = i;
      for(int l=1; l <= i1; ++l){
        m = m+1;
        store = store - C.at(c_shifted + m-1)*a.at(i-1).at(l);
      }
      C.at(c_shifted + i-1) = store/a.at(i-1).at(0);
      i -= 1;
    }
  }
};

void BivariateBSpline::fpback(const std::vector<std::vector<double>>& a, const std::vector<double>& z, const int n, const int k, std::vector<double>& C, const int nest){
  //No shift version of fpback - overloaded

  // std::cout << "fpback called" << std::endl;

  C.at(n-1) = z.at(n-1)/a.at(n-1).at(0);
  
  int i = n-1;
  double store;

  if(i!=0){
    for(int j=2; j <= n; ++j){
      store = z.at(i-1);
      int i1 = k-1;

      if(j <= k-1){
        i1 = j-1;
      }
      int m = i;
      for(int l=1; l <= i1; ++l){
        m = m+1;
        store = store - C.at(m-1)*a.at(i-1).at(l);
      }
      C.at(i-1) = store/a.at(i-1).at(0);
      i -= 1;
    }
  }
};

void BivariateBSpline::fpbspl(const std::vector<double>& t,const int n, const int k, const double x, const int l, std::vector<double>& h){
  std::vector<double> hh(19,0.0);

  //  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
  //  degree k at t(l) <= x < t(l+1) using the stable recurrence
  //  relation of de boor and cox.
  //  Travis Oliphant  2007
  //    changed so that weighting of 0 is used when knots with
  //      multiplicity are present.
  //    Also, notice that l+k <= n and 1 <= l+1-k
  //      or else the routine will be accessing memory outside t
  //      Thus it is imperative that that k <= l <= n-k but this
  //      is not checked.

  h.at(1-1) = 1;

  for(int j = 1; j <= k; ++j){
    for(int i = 1; i <= j; ++i){
      hh.at(i-1) = h.at(i-1);
    }
    h.at(1-1) = 0;
    for(int i = 1; i <= j; ++i){
      int li = l+i;
      int lj = li-j;

      if(t.at(li-1) != t.at(lj-1)){
        double f = hh.at(i-1) / (t.at(li-1)-t.at(lj-1));
        h.at(i-1) = h.at(i-1)+f*(t.at(li-1)-x);
        h.at(i) = f*(x-t.at(lj-1));
      } else {
        h.at(i) = 0;
      }
    }
  }
};

/**
  given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
  i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
  iate spline approximation s(x,y) of degrees kx and ky on the rect-
  angle xb <= x <= xe, yb <= y <= ye.
  the total numbers nx and ny of these knots and their
  position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
  ally by the routine. the smoothness of s(x,y) is then achieved by
  minimalizing the discontinuity jumps in the derivatives of s(x,y)
  across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
  the amounth of smoothness is determined by the condition that f(p) =
  sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
  constant, called the smoothing factor.
  the fit is given in the b-spline representation (b-spline coefficients
  C((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
  uated by means of subroutine bispev.

  references:
  dierckx p. : a fast algorithm for smoothing data on a rectangular
  grid while using spline functions, siam j.numer.anal.
  19 (1982) 1286-1304.
  dierckx p. : a fast algorithm for smoothing data on a rectangular
  grid while using spline functions, report tw53, dept.
  computer science,k.u.leuven, 1980.
  dierckx p. : curve and surface fitting with splines, monographs on
  numerical analysis, oxford university press, 1993.

  author:
  p.dierckx
  dept. computer science, k.u. leuven
  celestijnenlaan 200a, b-3001 heverlee, belgium.
  e-mail : Paul.Dierckx@cs.kuleuven.ac.be

  creation date : may 1979
  latest update : march 1989
  **/
// nx,tx,ny,ty,C,fp,ier = regrid_smth(x,y,z,.at.at(x)b,xe,yb,ye,kx,ky,s))
regrid_return BivariateBSpline::regrid_smth(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z){
  // std::cout << "regrid_smth called" <<std::endl;

  int mx = x.size();
  int my = y.size();

  int kx = 3;
  int ky = 3;

  if (not((mx > kx) and (my > ky))){
  throw std::runtime_error("Grid too small for bicubic interpolation");
  };

  const double xb = x.at(0);
  const double xe = x.at(mx-1);
  const double yb = y.at(0);
  const double ye = y.at(my-1);

  const int nxest = mx+kx+1;
  const int nyest = my+ky+1;

  std::vector<double> tx(nxest, 0.0);
  std::vector<double> ty(nyest, 0.0);
  std::vector<double> C((nxest-kx-1)*(nyest-ky-1), 0.0);
  double fp = 0.0;

  int lwrk = 4 + nxest * (my + 2 * kx + 5) + nyest * (2 * ky + 5) + mx * (kx + 1) + my * (ky + 1) + std::max(my, nxest);
  std::vector<double> wrk(lwrk, 0.0);
  int kwrk = 3 + mx + my + nxest + nyest;
  std::vector<int> iwrk(kwrk, 0);
  
  int nminx = 2*(kx+1);
  int nminy = 2*(ky+1);
  
  int lwest = 4 + nxest * (my + 2 * (kx + 2) + 1) + nyest * (2 * (ky + 2) + 1) + mx * (kx + 1) +  my * (ky + 1)+std::max(nxest,my);
  int kwest = 3 + mx + my + nxest + nyest;
  //  before starting computations a data check is made. if the input data
  if(kx <= 0 or kx > 5) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 1");
  if(ky <= 0 or ky > 5) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 2");

  if(mx < (kx+1) or nxest < nminx) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 4");
  if(my < (ky+1) or nyest < nminy) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 5");
  if(lwrk < lwest or kwrk < kwest) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 6");
  if(xb > x.at(1-1) or xe < x.at(mx-1)) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 7");
  for(int i = 1; i< mx; ++i){
    if(x.at(i-1) >= x.at(i)) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 8");
  }
  if(yb > y.at(1-1) or ye < y.at(my-1)) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 9");
  for(int j = 1; j < my; ++j){
    if(y.at(j-1) >= y.at(j)) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 10");
  }
  if((nxest < (mx+(kx+1)) or nyest < (my+(ky+1))) ) throw std::runtime_error("Error in BivariateBSpline/regrid_smth (C++ translation) - see code 13");

  //  we partition the working space and determine the spline approximation;

  int nx     = 0;
  int ny     = 0;

  std::vector<double> fpintx(nxest, 0.0);
  std::vector<double> fpinty(nyest, 0.0);
  std::vector<int>    nrx   (mx, 0);
  std::vector<int>    nry   (my, 0);

  //  we partition the working space.
  int mm = std::max(nxest,my);
  int mynx = nxest*my;

  // std::cout << "fpregr called" <<std::endl;

  //
  // part 1: determination of the number of knots and their position.
  // ****************************************************************
  //  given a set of knots we compute the least-squares spline sinf(x,y),
  //  and the corresponding sum of squared residuals fp=f(p=inf).
  //  sinf(x,y) is the requested approximation.
  //  
  //    we have spline interpolation; in that case the number of
  //    knots equals mx+kx+1 = mx+kx+1  and  my+ky+1 = my+ky+1.
  //

  //  find mx+kx+1 and my+ky+1 which denote the number of knots in x- and y-
  //  direction in case of spline interpolation.
  //  if s = 0, s(x,y) is an interpolating spline.
  nx = mx+kx+1;
  ny = my+ky+1;
  //  find the position of the interior knots in case of interpolation.
  //  the knots in the x-direction.
  int i = (kx+1)+1;
  int j = kx/2+2;
  for(int l = 1; l <= (mx-kx-1); ++l){
    tx.at(i-1) = x.at(j-1);
    i = i+1;
    j = j+1;
  }
  //  the knots in the y-direction.
  i = ky+2;
  j = ky/2+2;

  for(int l = 1; l <= (my-ky-1); ++l){
    ty.at(i-1) = y.at(j-1);
    i = i+1;
    j = j+1;
  }

  //  main loop for the different sets of knots.mpm=mx+my is a save upper
  //  bound for the number of trials.

  //  find the position of the additional knots which are needed for the
  //  b-spline representation of s(x,y).
  i = nx;
  for(int j=1; j <= (kx+1); ++j){
    tx.at(j-1) = xb;
    tx.at(i-1) = xe;
    i = i-1;
  }

  i = ny;
  for(int j=1; j <= (ky+1); ++j){
    ty.at(j-1) = yb;
    ty.at(i-1) = ye;
    i = i-1;
  }
  
  //  find the least-squares spline sinf(x,y) and calculate for each knot
  //  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
  //  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
  //  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
  //  value belonging to that interval.
  //  fp gives the total sum of squared residuals.

  std::vector<double>  fpx = fpintx;
  std::vector<double>  fpy = fpinty;
  std::vector<std::vector<double>> spx(mx, std::vector<double>(kx+1, 0.0));
  std::vector<std::vector<double>> spy(my, std::vector<double>(ky+1, 0.0));
  std::vector<double> right(mm, 0.0);
  std::vector<double> q(mynx, 0.0);
  std::vector<std::vector<double>> ax(nx, std::vector<double>(kx+2, 0.0));
  std::vector<std::vector<double>> ay(ny, std::vector<double>(ky+2, 0.0));
  std::vector<std::vector<double>> bx(nx, std::vector<double>(kx+2, 0.0));
  std::vector<std::vector<double>> by(ny, std::vector<double>(ky+2, 0.0));


  //  ..subroutine references..
    //    fpback,fpbspl,fpgivs,fpdisc,fprota
    //  ..
    //  the b-spline coefficients of the smoothing spline are calculated as
    //  the least-squares solution of the over-determined linear system of
    //  equations  (ay) C (ax)' = q       where
    //
    //               |   (spx)    |            |   (spy)    |
    //        (ax) = | ---------- |     (ay) = | ---------- |
    //               | (1/p) (bx) |            | (1/p) (by) |
    //
    //                                | z  ' 0 |
    //                            q = | ------ |
    //                                | 0  ' 0 |
    //
    //  with C      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
    //                b-spline coefficients.
    //       z      : the my x mx matrix which contains the function values.
    //       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
    //                matrices according to the least-squares problems in
    //                the x- and y-direction.
    //       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
    //                matrices which contain the discontinuity jumps of the
    //                derivatives of the b-splines in the x- and y-direction.

  // std::cout << "fpgrre called" <<std::endl;

  std::vector<double> h(7, 0.0);

  //  calculate the non-zero elements of the matrix (spx) which is the
  //  observation matrix according to the least-squares spline approximat-
  //  ion problem in the x-direction.
  int l = (kx+1);
  int l1 = (kx+2);
  int number = 0;


  for(int it = 1; it <= mx; ++it){
    double arg = x.at(it-1);

    while(not(arg<tx.at(l1-1) or l==nx-(kx+1))){
      l = l1;
      l1 = l+1;
      number = number+1;
    }

    fpbspl(tx,nx,kx,arg,l,h);
    for(int i=1; i <= (kx+1); ++i){
      spx.at(it-1).at(i-1) = h.at(i-1);
    }
    nrx.at(it-1) = number;
  }

  //  calculate the non-zero elements of the matrix (spy) which is the
  //  observation matrix according to the least-squares spline approximat-
  //  ion problem in the y-direction.
  l = (ky+1);
  l1 = (ky+2);
  number = 0;

  for(int it = 1; it <= my; ++it){
    double arg = y.at(it-1);

    while(not(arg<ty.at(l1-1) or l==ny-(ky+1))){
      l = l1;
      l1 = l+1;
      number +=1;
    }

    fpbspl(ty,ny,ky,arg,l,h);
    for(int i=1; i <= (ky+1); ++i){
      spy.at(it-1).at(i-1) = h.at(i-1);
    }
    nry.at(it-1) = number;
  }

  //  reduce the matrix (ax) to upper triangular form (rx) using givens
  //  rotations. apply the same transformations to the rows of matrix q
  //  to obtain the my x (nx-kx-1) matrix g.
  //  store matrix (rx) into (ax) and g into q.
  l = my * (nx-(kx+1));
  //  initialization.
  for(int i = 1; i <= l; ++i){
    q.at(i-1) = 0;
  }
  for(int i = 1; i <= nx-(kx+1); ++i){
    for(int j = 1; i <= (kx+2); ++i){
      ax.at(i-1).at(j-1) = 0;
    }
  }

  l = 0;
  // int nrold = 0;
  int ibandx = (kx+1);
  //  ibandx denotes the bandwidth of the matrices (ax) and (rx).
  for(int it = 1; it <= mx; ++it){
    number = nrx.at(it-1);
    // fetch a new row of matrix (spx).
    h.at(ibandx-1) = 0.;
    for(int j = 1; j <= (kx+1); ++j){
      h.at(j-1) = spx.at(it-1).at(j-1);
    }

    // find the appropriate column of q.
    for(int j = 1; j <= my; ++j){
      l += 1;
      right.at(j-1) = z.at(l-1);
    }

    int irot = number;
    // rotate the new row of matrix (ax) into triangle.
    double sin = 0.0;
    double cos = 0.0;
    for(int i = 1; i <= ibandx; ++i){
      irot += 1;
      double piv = h.at(i-1);

      if(piv == 0){continue;}
      //calculate the parameters of the given transformation.
      fpgivs(piv, ax.at(irot-1).at(1-1), cos, sin);
      //apply that transformation to the rows of matrix q.
      int iq = (irot - 1) * my;
      for(int j = 1; j <= my; ++j){
        iq += 1;
        fprota(cos, sin, right.at(j-1), q.at(iq-1));
      }
      //apply that transformation to the columns of (ax).
      if(i == ibandx){break;} //Break out of inner loop, continue on outer loop
      int i2 = 1;
      int i3 = i+1;
      for(int j = i3; j <= ibandx; ++j){
        i2 += 1;
        fprota(cos, sin, h.at(j-1), ax.at(irot-1).at(i2-1));
      }
    }
  }

  // !  reduce the matrix (ay) to upper triangular form (ry) using givens
  // !  rotations. apply the same transformations to the columns of matrix g
  // !  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
  // !  store matrix (ry) into (ay) and h into C.
  // int ncof = (nx-(kx+1)) * (ny-(ky+1));
  // !  initialization.

  // for(int i = 1; i<= ncof; ++i){
  //   C.at(i)-1] = 0;
  // }

  // for(int i=1; i <= ny-(ky+1); ++i){
  //   for(int j=1; j <= nx-(kx+1); ++j){
  //     ay.at(i)-1].at(j)-1] = 0;
  //   }
  // }

  int ibandy = (ky+1);

  for(int it = 1; it <= my; ++it){
    number = nry.at(it-1);
    // fetch a new row of matrix (spx).
    h.at(ibandy-1) = 0.;

    for(int j = 1; j <= (ky+1); ++j){
      h.at(j-1) = spy.at(it-1).at(j-1);
    }

    l = it;
    // find the appropriate column of q.
    for(int j = 1; j <= nx-(kx+1); ++j){
      right.at(j-1) = q.at(l-1);
      l += my;
    }

    int irot = number;
    // rotate the new row of matrix (ay) into triangle.
    double sin = 0.0;
    double cos = 0.0;
    for(int i = 1; i <= ibandy; ++i){
      irot += 1;
      double piv = h.at(i-1);

      if(piv == 0){continue;}
      //calculate the parameters of the given transformation.
      fpgivs(piv, ay.at(irot-1).at(1-1), cos, sin);
      //apply that transformation to the rows of matrix q.
      int ic = irot;
      for(int j = 1; j <= nx-(kx+1); ++j){
        fprota(cos, sin, right.at(j-1), C.at(ic-1));
        ic += ny-(ky+1);
      }
      //apply that transformation to the columns of (ay).
      if(i == ibandy){break;} //Break out of inner loop, continue on outer loop
      int i2 = 1;
      int i3 = i+1;
      for(int j = i3; j <= ibandy; ++j){
        i2 += 1;
        fprota(cos, sin, h.at(j-1), ay.at(irot-1).at(i2-1));
      }
    }
  }

  int k = 1;
  for(int i=1; i <= nx-(kx+1); ++i){
    fpback(ay,C,k,ny-(ky+1),ibandy,C,k,ny);
    k = k + ny-(ky+1);
  }
  k = 0;
  for(int j=1; j <= ny-(ky+1); ++j){
    k += 1;
    l = k;
    for(int i=1; i <= nx-(kx+1); ++i){
      right.at(i-1) = C.at(l-1);
      l += ny-(ky+1);
    }

    fpback(ax,right,nx-(kx+1),ibandx,right,nx);

    l = k;
    for(int i = 1; i <= nx-(kx+1); ++i){
      C.at(l-1) = right.at(i-1);
      l += ny-(ky+1);
    }
  }
  fp = 0.;
  for(int i=1;i <= nx;++i){
    fpx.at(i-1) = 0.;
  }
  for(int i=1;i <= ny;++i){
    fpy.at(i-1) = 0.;
  }

  int iz = 0;
  int nroldx = 0;
  //  main loop for the different grid points.
  for(int i1 = 1; i1 <= mx; ++i1){
    int numx = nrx.at(i1-1);
    int numx1 = numx+1;
    int nroldy = 0;
    for(int i2 = 1; i2 <= my; ++i2){
      int numy = nry.at(i2-1);
      int numy1 = numy+1;
      iz = iz+1;
      //  evaluate s(x,y) at the current grid point by making the sum of the
      //  cross products of the non-zero b-splines at (x,y), multiplied with
      //  the appropriate b-spline coefficients.
      double term = 0.;
      int k1 = numx * (ny-(ky+1))+numy;
      for(int l1 = 1; l1 <= (kx+1); ++l1){
        int k2 = k1;
        double fac = spx.at(i1-1).at(l1-1);
        for(int l2 = 1; l2 <= (ky+1); ++l2){
          k2 = k2+1;
          term = term + fac * spy.at(i2-1).at(l2-1) * C.at(k2-1);
       }
        k1 = k1+ny-(ky+1);
      }
      //  calculate the squared residual at the current grid point.
      term = (z.at(iz-1)-term)*(z.at(iz-1)-term);
      //  adjust the different parameters.
      fp = fp+term;
      fpx.at(numx1-1) = fpx.at(numx1-1)+term;
      fpy.at(numy1-1) = fpy.at(numy1-1)+term;
      double fac = term*0.5;
      if(numy!=nroldy){
        fpy.at(numy1-1) = fpy.at(numy1-1)-fac;
        fpy.at(numy-1) = fpy.at(numy-1)+fac;
      }
      nroldy = numy;
      if(numx!=nroldx){
        fpx.at(numx1-1) = fpx.at(numx1-1)-fac;
        fpx.at(numx-1) = fpx.at(numx-1)+fac;
      }
    }
    nroldx = numx;
  }
  //  if nx=mx+kx+1 and ny=my+ky+1, sinf(x,y) is an interpolating spline.

  regrid_return setup_spline;

  // std::cout << "nx: " << nx << std::endl;
  // std::cout << "ny: " << ny << std::endl;
  // std::cout << "tx: " << tx << std::endl;
  // std::cout << "ty: " << ty << std::endl;
  // std::cout << "C: " << C << std::endl;
  // std::cout << "fp: " << fp << std::endl;

  setup_spline.nx = nx;
  setup_spline.ny = ny;
  setup_spline.tx = tx;
  setup_spline.ty = ty;
  setup_spline.C = C;
  setup_spline.fp = fp;

  return setup_spline;
  // End of translation of fitpack/regrid.f (implementation)
};

// fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
double BivariateBSpline::fpbisp(const std::vector<double>& tx, const int& nx, const std::vector<double>& ty, const int& ny, const std::vector<double>& c, const int& kx, const int& ky, const double& x, const double& y){

  std::vector<double> h(6, 0.0);
  int lx = 0;
  int ly = 0;
  std::vector<double> wx(kx+1, 0.0);
  std::vector<double> wy(ky+1, 0.0);

  // int (kx+1) = kx+1;
  int nkx1 = nx-(kx+1);
  double tb = tx.at(kx);
  double te = tx.at(nkx1);
  int l = (kx+1);
  int l1 = l+1;

  double arg = x;
  if(arg < tb){
    arg = tb;
  }
  if(arg > te){
    arg = te;
  }

  while(not(arg < tx.at(l1-1) or l==nkx1)){
    l = l1;
    l1 = l+1;
  }
  fpbspl(tx,nx,kx,arg,l,h);
  lx = l-(kx+1);

  for(int j=1; j <= (kx+1); ++j){
    wx.at(j-1) = h.at(j-1);
  }

  tb = ty.at(ky);
  te = ty.at((ny-(ky+1)));
  l = (ky+1);
  l1 = l+1;

  arg = y;
  if(arg < tb){
    arg = tb;
  }
  if(arg > te){
    arg = te;
  }
  while(not(arg<ty.at(l1-1) or l==(ny-(ky+1)))){
    l = l1;
    l1 = l+1;
  }
  fpbspl(ty,ny,ky,arg,l,h);
  ly = l-(ky+1);

  for(int j=1; j<=(ky+1); ++j){
    wy.at(j-1) = h.at(j-1);
  }

  int m = 0;

  l = lx*(ny-(ky+1));
  for(int i1=1; i1<= (kx+1); ++i1){
    h.at(i1-1) = wx.at(i1-1);
  }
  l1 = l+ly;
  double sp = 0.;
  for(int i1=1; i1<= (kx+1); ++i1){
    int l2 = l1;
    for(int j1=1; j1<= (ky+1); ++j1){
      l2 = l2+1;
      sp = sp+c.at(l2-1)*h.at(i1-1)*wy.at(j1-1);
    }
    l1 = l1+(ny-(ky+1));
  }
  m = m+1;
  double z = sp;
  return z;
}

/**
  bispeu evaluates a bivariate spline s(x,y) of degrees kx and ky, given in the b-spline representation.
  input parameters:
    tx    : real array, length nx, which contains the position of the
           knots in the x-direction.
    nx    : integer, giving the total number of knots in the x-direction
    ty    : real array, length ny, which contains the position of the
           knots in the y-direction.
    ny    : integer, giving the total number of knots in the y-direction
    C     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
           b-spline coefficients.
    kx,ky : integer values, giving the degrees of the spline.
    x     : real array of dimension (mx).
    y     : real array of dimension (my).
    m     : on entry m must specify the number points. m >= 1.
    wrk   : real array of dimension lwrk. used as workspace.
    lwrk  : integer, specifying the dimension of wrk.
           lwrk >= kx+ky+2

  output parameters:
    z     : real array of dimension m.
           on successful exit z(i) contains the value of s(x,y)
           at the point (x(i),y(i)), i=1,...,m.

  restrictions:
    m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
    tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
    ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
 **/
double BivariateBSpline::bispeu(const std::vector<double>& tx,const std::vector<double>& ty,const std::vector<double>& c,const double& x,const double& y){
  // z = bispeu(tx,ty,c,kx,ky,x,y)

  int nx = tx.size();
  int ny = ty.size();
  int kx = 3;
  int ky = 3;

  double z = fpbisp(tx, nx, ty, ny, c, kx, ky, x, y);

  return z;
}






















