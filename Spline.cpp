#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm> //for upper/lower_bound
#include <stdexcept> //For error-throwing
#include <sstream>
#include <iostream>
#include "prettyprint.hpp"

void fprota(const double cos, const double sin, double& a, double& b){
    //  subroutine fprota applies a given rotation to a and b.
    double stor1 = a;
    double stor2 = b;

    b = cos * stor2 + sin * stor1;
    a = cos * stor1 - sin * stor2;
};

void fpgivs(const double piv, double& ww, double& cos, double& sin){
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

void fpback(const std::vector<std::vector<double>>& a, const std::vector<double>& z, const int z_start, const int n, const int k, std::vector<double>& C, const int c_start, const int nest){
    // Switch from Fortran to C++ indexing
    // Supply start as '1' to not apply a shift
    int z_shifted = z_start - 1;
    int c_shifted = c_start - 1;

    C[c_shifted + n-1] = z[z_shifted + n-1]/a[n-1][0];
    
    int i = n-1;
    double store;

    if(i!=0){
      for(int j=2; j <= n; ++j){
        store = z[z_shifted + i-1];
        int i1 = k-1;

        if(j <= k-1){
          i1 = j-1;
        }
        int m = i;
        for(int l=1; l <= i1; ++l){
          m = m+1;
          store = store - C[c_shifted + m-1]*a[i-1][l];
        }
        C[c_shifted + i-1] = store/a[i-1][0];
        i -= 1;
      }
    }
};
void fpback(const std::vector<std::vector<double>>& a, const std::vector<double>& z, const int n, const int k, std::vector<double>& C, const int nest){
    //No shift version of fpback - overloaded

    // std::cout << "fpback called" << std::endl;

    C[n-1] = z[n-1]/a[n-1][0];
    
    int i = n-1;
    double store;

    if(i!=0){
      for(int j=2; j <= n; ++j){
        store = z[i-1];
        int i1 = k-1;

        if(j <= k-1){
          i1 = j-1;
        }
        int m = i;
        for(int l=1; l <= i1; ++l){
          m = m+1;
          store = store - C[m-1]*a[i-1][l];
        }
        C[i-1] = store/a[i-1][0];
        i -= 1;
      }
    }
};

void fpbspl(const std::vector<double>& t,const int n, const int k, const double x, const int l, std::vector<double>& h){
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

    h[1-1] = 1;

    for(int j = 1; j <= k; ++j){
      for(int i = 1; i <= j; ++i){
        hh[i-1] = h[i-1];
      }
      h[1-1] = 0;
      for(int i = 1; i <= j; ++i){
        int li = l+i;
        int lj = li-j;

        if(t[li-1] != t[lj-1]){
          double f = hh[i-1] / (t[li-1]-t[lj-1]);
          h[i-1] = h[i-1]+f*(t[li-1]-x);
          h[i] = f*(x-t[lj-1]);
        } else {
          h[i] = 0;
        }
      }
    }
};

void fpgrre(
  int& ifsx,
  int& ifsy,
  int& ifbx,
  int& ifby,
  const std::vector<double>& x,
  const int& mx,
  const std::vector<double>& y,
  const int& my,
  const std::vector<double>& z,
  const int& mz,
  const int& kx,
  const int& ky,
  const std::vector<double>& tx,
  const int& nx,
  const std::vector<double>& ty,
  const int& ny,
  const double p,
  std::vector<double>& C,
  const int& nc,
  double fp,
  std::vector<double>& fpx,
  std::vector<double>& fpy,
  std::vector<std::vector<double>>& spx,
  std::vector<std::vector<double>>& spy,
  std::vector<double>& right,
  std::vector<double>& q,
  std::vector<std::vector<double>>& ax,
  std::vector<std::vector<double>>& ay,
  std::vector<std::vector<double>>& bx,
  std::vector<std::vector<double>>& by,
  std::vector<int>& nrx,
  std::vector<int>& nry
  ){
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
        double arg = x[it-1];

        while(not(arg<tx[l1-1] or l==nx-(kx+1))){
          l = l1;
          l1 = l+1;
          number = number+1;
        }

        fpbspl(tx,nx,kx,arg,l,h);
        // std::cout << "Check" << h << std::endl;
        for(int i=1; i <= (kx+1); ++i){
          spx[it-1][i-1] = h[i-1];
        }
        nrx[it-1] = number;
      }

      ifsx = 1;
      //  calculate the non-zero elements of the matrix (spy) which is the
      //  observation matrix according to the least-squares spline approximat-
      //  ion problem in the y-direction.
      l = (ky+1);
      l1 = (ky+2);
      number = 0;

      for(int it = 1; it <= my; ++it){
        double arg = y[it-1];

        while(not(arg<ty[l1-1] or l==ny-(ky+1))){
          l = l1;
          l1 = l+1;
          number +=1;
        }

        fpbspl(ty,ny,ky,arg,l,h);
        // std::cout << "Check" << h << std::endl;
        for(int i=1; i <= (ky+1); ++i){
          spy[it-1][i-1] = h[i-1];
        }
        nry[it-1] = number;
      }
      ifsy = 1;

      //  reduce the matrix (ax) to upper triangular form (rx) using givens
      //  rotations. apply the same transformations to the rows of matrix q
      //  to obtain the my x (nx-kx-1) matrix g.
      //  store matrix (rx) into (ax) and g into q.
      l = my*nx-(kx+1);
      //  initialization.
      for(int i = 1; i <= l; ++i){
        q[i-1] = 0;
      }
      for(int i = 1; i <= nx-(kx+1); ++i){
        for(int j = 1; i <= (kx+2); ++i){
          ax[i-1][j-1] = 0;
        }
      }

      l = 0;
      // int nrold = 0;
      int ibandx = (kx+1);
      //  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      
      for(int it = 1; it <= mx; ++it){
        number = nrx[it-1];
        // fetch a new row of matrix (spx).
        h[ibandx-1] = 0.;

        for(int j = 1; j <= (kx+1); ++j){
          h[j-1] = spx[it-1][j-1];
        }

        // find the appropriate column of q.
        for(int j = 1; j <= my; ++j){
          l += 1;
          right[j-1] = z[l-1];
        }

        int irot = number;
        // rotate the new row of matrix (ax) into triangle.
        double sin = 0.0;
        double cos = 0.0;
        for(int i = 1; i <= ibandx; ++i){
          irot += 1;
          double piv = h[i-1];

          if(piv == 0){continue;}
          //calculate the parameters of the given transformation.
          fpgivs(piv, ax[irot-1][1-1], cos, sin);
          //apply that transformation to the rows of matrix q.
          int iq = (irot - 1) * my;
          for(int j = 1; j <= my; ++j){
            iq += 1;
            fprota(cos, sin, right[j-1], q[iq-1]);
          }
          //apply that transformation to the columns of (ax).
          if(i == ibandx){break;} //Break out of inner loop, continue on outer loop
          int i2 = 1;
          int i3 = i+1;
          for(int j = i3; j <= ibandx; ++j){
            i2 += 1;
            fprota(cos, sin, h[j-1], ax[irot-1][i2-1]);
          }
        }
      }


      // !  reduce the matrix (ay) to upper triangular form (ry) using givens
      // !  rotations. apply the same transformations to the columns of matrix g
      // !  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
      // !  store matrix (ry) into (ay) and h into C.
      int ncof = nx-(kx+1) * ny-(ky+1);
      // !  initialization.

      for(int i = 1; i<= ncof; ++i){
        C[i-1] = 0;
      }

      for(int i=1; i<=ny-(ky+1); ++i){
        for(int j=1; j<=nx-(kx+1); ++j){
          ay[i-1][j-1] = 0;
        }
      }

      int ibandy = (ky+1);

      for(int it = 1; it <= my; ++it){
        number = nry[it-1];
        // fetch a new row of matrix (spx).
        h[ibandy-1] = 0.;

        for(int j = 1; j <= (ky+1); ++j){
          h[j-1] = spy[it-1][j-1];
        }

        l = it;
        // find the appropriate column of q.
        for(int j = 1; j <= nx-(kx+1); ++j){
          right[j-1] = q[l-1];
          l += my;
        }

        int irot = number;
        // rotate the new row of matrix (ay) into triangle.
        double sin = 0.0;
        double cos = 0.0;
        for(int i = 1; i <= ibandy; ++i){
          irot += 1;
          double piv = h[i-1];

          if(piv == 0){continue;}
          //calculate the parameters of the given transformation.
          fpgivs(piv, ay[irot-1][1-1], cos, sin);
          //apply that transformation to the rows of matrix q.
          int ic = irot;
          for(int j = 1; j <= nx-(kx+1); ++j){
            fprota(cos, sin, right[j-1], C[ic-1]);
            ic += ny-(ky+1);
          }
          //apply that transformation to the columns of (ay).
          if(i == ibandy){break;} //Break out of inner loop, continue on outer loop
          int i2 = 1;
          int i3 = i+1;
          for(int j = i3; j <= ibandy; ++j){
            i2 += 1;
            fprota(cos, sin, h[j-1], ay[irot-1][i2-1]);
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
        right[i-1] = C[l-1];
        l += ny-(ky+1);
      }

      fpback(ax,right,nx-(kx+1),ibandx,right,nx);

      l = k;
      for(int i = 1; i <= nx-(kx+1); ++i){
        C[l-1] = right[i-1];
        l += ny-(ky+1);
      }
    }

    fp = 0.;
    for(int i=1;i <= nx;++i){
      fpx[i-1] = 0.;
    }
    for(int i=1;i <= ny;++i){
      fpy[i-1] = 0.;
    }

    int iz = 0;
    int nroldx = 0;

    //  main loop for the different grid points.
    for(int i1=1; i1<=mx;++i1){
      int numx = nrx[i1-1];
      int numx1 = numx+1;
      int nroldy = 0;
      for(int i2=1; i2<=my;++i2){
        int numy = nry[i2-1];
        int numy1 = numy+1;
        iz = iz+1;
        //  evaluate s(x,y) at the current grid point by making the sum of the
        //  cross products of the non-zero b-splines at (x,y), multiplied with
        //  the appropiate b-spline coefficients.
        double term = 0.;
        int k1 = numx*ny-(ky+1)+numy;
        for(int l1=1; l1<=(kx+1);++l1){
          int k2 = k1;
          double fac = spx[i1-1][l1-1];
          for(int l2=1; l2<=(ky+1);++l2){
            k2 = k2+1;
            term = term+fac*spy[i2-1][l2-1]*C[k2-1];
          }
          k1 = k1+ny-(ky+1);
        }
        //  calculate the squared residual at the current grid point.
        term = (z[iz-1]-term)*(z[iz-1]-term);
        //  adjust the different parameters.
        fp = fp+term;
        fpx[numx1-1] = fpx[numx1-1]+term;
        fpy[numy1-1] = fpy[numy1-1]+term;
        double fac = term*0.5;
        if(numy!=nroldy){
          fpy[numy1-1] = fpy[numy1-1]-fac;
          fpy[numy-1] = fpy[numy-1]+fac;
        }
        nroldy = numy;
        if(numx!=nroldx){
          fpx[numx1-1] = fpx[numx1-1]-fac;
          fpx[numx-1] = fpx[numx-1]+fac;
        }
      }
      nroldx = numx;
    }
};

int fpregr(
  const int& iopt,
  const std::vector<double>& x,
  const int& mx,
  const std::vector<double>& y,
  const int& my,
  const std::vector<double>& z,
  const int& mz,
  const double& xb,
  const double& xe,
  const double& yb,
  const double& ye,
  const int& kx,
  const int& ky,
  const double& s,
  const int& nxest,
  const int& nyest,
  const double& tol,
  const int& maxit,
  const int& nc,
  int& nx,
  std::vector<double>& tx,
  int& ny,
  std::vector<double>& ty,
  std::vector<double>& C,
  double fp,
  const std::vector<double>& fpintx,
  const std::vector<double>& fpinty,
  std::vector<int>& nrx,
  std::vector<int>& nry
  ){
  //Up to fgrre
      //  we partition the working space.
      int mm = std::max(nxest,my);
      int mynx = nxest*my;

      // std::cout << "fpregr called" <<std::endl;

      //
      // part 1: determination of the number of knots and their position.
      // ****************************************************************
      //  given a set of knots we compute the least-squares spline sinf(x,y),
      //  and the corresponding sum of squared residuals fp=f(p=inf).
      //  if iopt=-1  sinf(x,y) is the requested approximation.
      //  if iopt=0 or iopt=1 we check whether we can accept the knots:
      //    if fp <=s we will continue with the current set of knots.
      //    if fp > s we will increase the number of knots and compute the
      //       corresponding least-squares spline until finally fp<=s.
      //    the initial choice of knots depends on the value of s and iopt.
      //    if s=0 we have spline interpolation; in that case the number of
      //    knots equals mx+kx+1 = mx+kx+1  and  my+ky+1 = my+ky+1.
      //    if s>0 and
      //     *iopt=0 we first compute the least-squares polynomial of degree
      //      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.
      //     *iopt=1 we start with the knots found at the last call of the
      //      routine, except for the case that s > fp0; then we can compute
      //      the least-squares polynomial directly.
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
        tx[i-1] = x[j-1];
        i = i+1;
        j = j+1;
      }
      //  the knots in the y-direction.
      i = ky+2;
      j = ky/2+2;

      for(int l = 1; l <= (my-ky-1); ++l){
        ty[i-1] = y[j-1];
        i = i+1;
        j = j+1;
      }

      int ifsx = 0;
      int ifsy = 0;
      int ifbx = 0;
      int ifby = 0;
      double p = -1;

      //  main loop for the different sets of knots.mpm=mx+my is a save upper
      //  bound for the number of trials.

      //  find the position of the additional knots which are needed for the
      //  b-spline representation of s(x,y).
      i = nx;
      for(int j=1; j <= (kx+1); ++j){
        tx[j-1] = xb;
        tx[i-1] = xe;
        i = i-1;
      }

      i = ny;
      for(int j=1; j <= (ky+1); ++j){
        ty[j-1] = yb;
        ty[i-1] = ye;
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
      std::vector<std::vector<double>> spx(mx, std::vector<double>((kx+1), 0.0));
      std::vector<std::vector<double>> spy(mx, std::vector<double>((ky+1), 0.0));
      std::vector<double> right(mm, 0.0);
      std::vector<double> q(mynx, 0.0);
      std::vector<std::vector<double>> ax(nx, std::vector<double>((kx+2),0.0));
      std::vector<std::vector<double>> ay(ny, std::vector<double>((ky+2),0.0));
      std::vector<std::vector<double>> bx(nx, std::vector<double>((kx+2),0.0));
      std::vector<std::vector<double>> by(ny, std::vector<double>((ky+2),0.0));

    fpgrre(ifsx,
           ifsy,
           ifbx,
           ifby,
           x,
           mx,
           y,
           my,
           z,
           mz,
           kx,
           ky,
           tx,
           nx,
           ty,
           ny,
           p,
           C,
           nc,
           fp,
           fpx,
           fpy,
           spx,
           spy,
           right,
           q,
           ax,
           ay,
           bx,
           by,
           nrx,
           nry);

    //  if nx=mx+kx+1 and ny=my+ky+1, sinf(x,y) is an interpolating spline.
    if((nx==mx+kx+1) and (ny==my+ky+1)){
      int ier = -1;
      fp = 0.;
      return ier;
    } else {
      throw std::runtime_error("Not enough knots included in spline");
    }
};

struct regrid_return{
  int nx;
  int ny;
  std::vector<double> tx;
  std::vector<double> ty;
  std::vector<double> C;
  double fp;
  int ier;
};

/**
  given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
  i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
  iate spline approximation s(x,y) of degrees kx and ky on the rect-
  angle xb <= x <= xe, yb <= y <= ye.
  if iopt = -1 regrid calculates the least-squares spline according
  to a given set of knots.
  if iopt >= 0 the total numbers nx and ny of these knots and their
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

  parameters:
  iopt  : integer flag. on entry iopt must specify whether a least-
  squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
  must be determined.
  if iopt=0 the routine will start with an initial set of knots
  tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
  1,...,ky+1. if iopt=1 the routine will continue with the set
  of knots found at the last call of the routine.
  attention: a call with iopt=1 must always be immediately pre-
       ceded by another call with iopt=1 or iopt=0 and
       s != 0.
  unchanged on exit.
  mx    : integer. on entry mx must specify the number of grid points
  along the x-axis. mx > kx . unchanged on exit.
  x     : real array of dimension at least (mx). before entry, x(i)
  must be set to the x-co-ordinate of the i-th grid point
  along the x-axis, for i=1,2,...,mx. these values must be
  supplied in strictly ascending order. unchanged on exit.
  my    : integer. on entry my must specify the number of grid points
  along the y-axis. my > ky . unchanged on exit.
  y     : real array of dimension at least (my). before entry, y(j)
  must be set to the y-co-ordinate of the j-th grid point
  along the y-axis, for j=1,2,...,my. these values must be
  supplied in strictly ascending order. unchanged on exit.
  z     : real array of dimension at least (mx*my).
  before entry, z(my*(i-1)+j) must be set to the data value at
  the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
  unchanged on exit.
  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
  yb,ye   aries of the rectangular approximation domain.
  xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
  unchanged on exit.
  kx,ky : integer values. on entry kx and ky must specify the degrees
  of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
  (kx=ky=3) splines. unchanged on exit.
  s     : real. on entry (in case iopt>=0) s must specify the smoothing
  factor. s >=0. unchanged on exit.
  for advice on the choice of s see further comments
  nxest : integer. unchanged on exit.
  nyest : integer. unchanged on exit.
  on entry, nxest and nyest must specify an upper bound for the
  number of knots required in the x- and y-directions respect.
  these numbers will also determine the storage space needed by
  the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
  in most practical situation nxest = mx/2, nyest=my/2, will
  be sufficient. always large enough are nxest=mx+kx+1, nyest=
  my+ky+1, the number of knots needed for interpolation (s=0).
  see also further comments.
  nx    : integer.
  unless ier=10 (in case iopt >=0), nx will contain the total
  number of knots with respect to the x-variable, of the spline
  approximation returned. if the computation mode iopt=1 is
  used, the value of nx should be left unchanged between sub-
  sequent calls.
  in case iopt=-1, the value of nx should be specified on entry
  tx    : real array of dimension nmax.
  on succesful exit, this array will contain the knots of the
  spline with respect to the x-variable, i.e. the position of
  the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
  position of the additional knots tx(1)=...=tx(kx+1)=xb and
  tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
  if the computation mode iopt=1 is used, the values of tx(1),
  ...,tx(nx) should be left unchanged between subsequent calls.
  if the computation mode iopt=-1 is used, the values tx(kx+2),
  ...tx(nx-kx-1) must be supplied by the user, before entry.
  see also the restrictions (ier=10).
  ny    : integer.
  unless ier=10 (in case iopt >=0), ny will contain the total
  number of knots with respect to the y-variable, of the spline
  approximation returned. if the computation mode iopt=1 is
  used, the value of ny should be left unchanged between sub-
  sequent calls.
  in case iopt=-1, the value of ny should be specified on entry
  ty    : real array of dimension nmax.
  on succesful exit, this array will contain the knots of the
  spline with respect to the y-variable, i.e. the position of
  the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
  position of the additional knots ty(1)=...=ty(ky+1)=yb and
  ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
  if the computation mode iopt=1 is used, the values of ty(1),
  ...,ty(ny) should be left unchanged between subsequent calls.
  if the computation mode iopt=-1 is used, the values ty(ky+2),
  ...ty(ny-ky-1) must be supplied by the user, before entry.
  see also the restrictions (ier=10).
  C     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
  on succesful exit, C contains the coefficients of the spline
  approximation s(x,y)
  fp    : real. unless ier=10, fp contains the sum of squared
  residuals of the spline approximation returned.
  wrk   : real array of dimension (lwrk). used as workspace.
  if the computation mode iopt=1 is used the values of wrk(1),
  ...,wrk(4) should be left unchanged between subsequent calls.
  lwrk  : integer. on entry lwrk must specify the actual dimension of
  the array wrk as declared in the calling (sub)program.
  lwrk must not be too small.
  lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
  my*(ky+1) +u
  where u is the larger of my and nxest.
  iwrk  : integer array of dimension (kwrk). used as workspace.
  if the computation mode iopt=1 is used the values of iwrk(1),
  ...,iwrk(3) should be left unchanged between subsequent calls
  kwrk  : integer. on entry kwrk must specify the actual dimension of
  the array iwrk as declared in the calling (sub)program.
  kwrk >= 3+mx+my+nxest+nyest.
  ier   : integer. unless the routine detects an error, ier contains a
  non-positive value on exit, i.e.
  ier=0  : normal return. the spline returned has a residual sum of
  squares fp such that abs(fp-s)/s <= tol with tol a relat-
  ive tolerance set to 0.001 by the program.
  ier=-1 : normal return. the spline returned is an interpolating
  spline (fp=0).
  ier=-2 : normal return. the spline returned is the least-squares
  polynomial of degrees kx and ky. in this extreme case fp
  gives the upper bound for the smoothing factor s.
  ier=1  : error. the required storage space exceeds the available
  storage space, as specified by the parameters nxest and
  nyest.
  probably causes : nxest or nyest too small. if these param-
  eters are already large, it may also indicate that s is
  too small
  the approximation returned is the least-squares spline
  according to the current set of knots. the parameter fp
  gives the corresponding sum of squared residuals (fp>s).
  ier=2  : error. a theoretically impossible result was found during
  the iteration proces for finding a smoothing spline with
  fp = s. probably causes : s too small.
  there is an approximation returned but the corresponding
  sum of squared residuals does not satisfy the condition
  abs(fp-s)/s < tol.
  ier=3  : error. the maximal number of iterations maxit (set to 20
  by the program) allowed for finding a smoothing spline
  with fp=s has been reached. probably causes : s too small
  there is an approximation returned but the corresponding
  sum of squared residuals does not satisfy the condition
  abs(fp-s)/s < tol.
  ier=10 : error. on entry, the input data are controlled on validity
  the following restrictions must be satisfied.
  -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
  nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
  lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
  my*(ky+1) +max(my,nxest),
  xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
  if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
          xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
          2*ky+2<=ny<=min(nyest,my+ky+1)
          yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
      the schoenberg-whitney conditions, i.e. there must
      be subset of grid co-ordinates xx(p) and yy(q) such
      that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
             ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
  if iopt>=0: s>=0
          if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
  if one of these conditions is found to be violated,control
  is immediately repassed to the calling program. in that
  case there is no approximation returned.

  further comments:
  regrid does not allow individual weighting of the data-values.
  so, if these were determined to widely different accuracies, then
  perhaps the general data set routine surfit should rather be used
  in spite of efficiency.
  by means of the parameter s, the user can control the tradeoff
  between closeness of fit and smoothness of fit of the approximation.
  if s is too large, the spline will be too smooth and signal will be
  lost ; if s is too small the spline will pick up too much noise. in
  the extreme cases the program will return an interpolating spline if
  s=0 and the least-squares polynomial (degrees kx,ky) if s is
  very large. between these extremes, a properly chosen s will result
  in a good compromise between closeness of fit and smoothness of fit.
  to decide whether an approximation, corresponding to a certain s is
  satisfactory the user is highly recommended to inspect the fits
  graphically.
  recommended values for s depend on the accuracy of the data values.
  if the user has an idea of the statistical errors on the data, he
  can also find a proper estimate for s. for, by assuming that, if he
  specifies the right s, regrid will return a spline s(x,y) which
  exactly reproduces the function underlying the data he can evaluate
  the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
  for example, if he knows that the statistical errors on his z(i,j)-
  values is not greater than 0.1, he may expect that a good s should
  have a value not larger than mx*my*(0.1)**2.
  if nothing is known about the statistical error in z(i,j), s must
  be determined by trial and error, taking account of the comments
  above. the best is then to start with a very large value of s (to
  determine the least-squares polynomial and the corresponding upper
  bound fp0 for s) and then to progressively decrease the value of s
  ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
  and more carefully as the approximation shows more detail) to
  obtain closer fits.
  to economize the search for a good s-value the program provides with
  different modes of computation. at the first call of the routine, or
  whenever he wants to restart with the initial set of knots the user
  must set iopt=0.
  if iopt=1 the program will continue with the set of knots found at
  the last call of the routine. this will save a lot of computation
  time if regrid is called repeatedly for different values of s.
  the number of knots of the spline returned and their location will
  depend on the value of s and on the complexity of the shape of the
  function underlying the data. if the computation mode iopt=1
  is used, the knots returned may also depend on the s-values at
  previous calls (if these were smaller). therefore, if after a number
  of trials with different s-values and iopt=1, the user can finally
  accept a fit as satisfactory, it may be worthwhile for him to call
  regrid once more with the selected value for s but now with iopt=0.
  indeed, regrid may then return an approximation of the same quality
  of fit but with fewer knots and therefore better if data reduction
  is also an important objective for the user.
  the number of knots may also depend on the upper bounds nxest and
  nyest. indeed, if at a certain stage in regrid the number of knots
  in one direction (say nx) has reached the value of its upper bound
  (nxest), then from that moment on all subsequent knots are added
  in the other (y) direction. this may indicate that the value of
  nxest is too small. on the other hand, it gives the user the option
  of limiting the number of knots the routine locates in any direction
  for example, by setting nxest=2*kx+2 (the lowest allowable value for
  nxest), the user can indicate that he wants an approximation which
  is a simple polynomial of degree kx in the variable x.

  other subroutines required:
  fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
  fpknot

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
// nx,tx,ny,ty,C,fp,ier = regrid_smth(x,y,z,[xb,xe,yb,ye,kx,ky,s])
regrid_return regrid_smth(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z){
  //Up to fpregr
    // std::cout << "regrid_smth called" <<std::endl;

    // Translation of interpolate/src/fitpack.pyf (header)

    int iopt = 0;
    int mx = x.size();
    int my = y.size();

    int kx = 3;
    int ky = 3;

    if (not((mx > kx) and (my > ky))){
    throw std::runtime_error("Grid too small for bicubic interpolation");
    };

    const double s = 0.0;

    const double xb = x[0];
    const double xe = x[mx-1];
    const double yb = y[0];
    const double ye = y[my-1];

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

    int ier = 0;

    // End of translation of interpolate/src/fitpack.pyf (header)

    // call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nx,tx,ny,ty,C,fp,wrk,lwrk,iwrk,kwrk,ier)

    // Translation of fitpack/regrid.f (implementation)

    int maxit = 20;
    double tol = 0.1e-02;
    ier = 10;
    int nminx = 2*(kx+1);
    int nminy = 2*(ky+1);

    int mz = mx*my;
    int nc = (nxest-(kx+1))*(nyest-(ky+1));
    int lwest = 4+nxest*(my+2*(kx+2)+1)+nyest*(2*(ky+2)+1)+mx*(kx+1)+ my*(ky+1)+std::max(nxest,my);
    int kwest = 3+mx+my+nxest+nyest;
    //  before starting computations a data check is made. if the input data
    if(kx <= 0 or kx > 5) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 1");
    if(ky <= 0 or ky > 5) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 2");
    if(iopt < (-1) or iopt > 1) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 3");
    if(mx < (kx+1) or nxest < nminx) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 4");
    if(my < (ky+1) or nyest < nminy) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 5");
    if(lwrk < lwest or kwrk < kwest) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 6");
    if(xb > x[1-1] or xe < x[mx-1]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 7");
    for(int i = 1; i< mx; ++i){
      if(x[i-1] >= x[i]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 8");
    }
    if(yb > y[1-1] or ye < y[my-1]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 9");
    for(int j = 1; j < my; ++j){
      if(y[j-1] >= y[j]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 10");
    }
    if(not(iopt >= 0)) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 11");
    if((nxest < (mx+(kx+1)) or nyest < (my+(ky+1))) ) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 13");

    ier = 0;
    //  we partition the working space and determine the spline approximation;

    int nx     = 0; //These variables are modified by fpregr
    int ny     = 0; //

    std::vector<double> fpintx(nxest, 0.0);
    std::vector<double> fpinty(nyest, 0.0);
    std::vector<int>    nrx   (mx, 0);
    std::vector<int>    nry   (my, 0);

  ier = fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,
  nxest,nyest,tol,maxit,nc,nx,tx,ny,ty,C,fp,fpintx,fpinty,nrx,nry);

  regrid_return setup_spline;
  setup_spline.nx = nx;
  setup_spline.ny = ny;
  setup_spline.tx = tx;
  setup_spline.ty = ty;
  setup_spline.C = C;
  setup_spline.fp = fp;
  setup_spline.ier = ier;
  return setup_spline;
  // End of translation of fitpack/regrid.f (implementation)
};

// fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
double fpbisp(const std::vector<double>& tx, const int& nx, const std::vector<double>& ty, const int& ny, const std::vector<double>& c, const int& kx, const int& ky, const double& x, const double& y){

  std::vector<double> h(6, 0.0);
  int lx = 0;
  int ly = 0;
  std::vector<double> wx(kx+1, 0.0);
  std::vector<double> wy(ky+1, 0.0);

  // int (kx+1) = kx+1;
  int nkx1 = nx-(kx+1);
  double tb = tx[(kx+1)-1];
  double te = tx[nkx1];
  int l = (kx+1);
  int l1 = l+1;

  double arg = x;
  if(arg < tb){
    arg = tb;
  }
  if(arg > te){
    arg = te;
  }

  while(not(arg < tx[l1-1] or l==nkx1)){
    l = l1;
    l1 = l+1;
  }
  fpbspl(tx,nx,kx,arg,l,h);
  lx = l-(kx+1);

  for(int j=1; j <= (kx+1); ++j){
    wx[j-1] = h[j-1];
  }

  tb = ty[(ky+1)-1];
  te = ty[(ny-(ky+1))];
  l = (ky+1);
  l1 = l+1;

  arg = y;
  if(arg < tb){
    arg = tb;
  }
  if(arg > te){
    arg = te;
  }
  while(not(arg<ty[l1-1] or l==(ny-(ky+1)))){
    l = l1;
    l1 = l+1;
  }
  fpbspl(ty,ny,ky,arg,l,h);
  ly = l-(ky+1);

  for(int j=1; j<=(ky+1); ++j){
    wy[j-1] = h[j-1];
  }

  int m = 0;

  l = lx*(ny-(ky+1));
  for(int i1=1; i1<= (kx+1); ++i1){
    h[i1-1] = wx[i1-1];
  }
  l1 = l+ly;
  double sp = 0.;
  for(int i1=1; i1<= (kx+1); ++i1){
    int l2 = l1;
    for(int j1=1; j1<= (ky+1); ++j1){
      l2 = l2+1;
      sp = sp+c[l2-1]*h[i1-1]*wy[j1-1];
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
    ier   : integer error flag
    ier=0 : normal return
    ier=10: invalid input data (see restrictions)

  restrictions:
    m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
    tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
    ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
 **/
double bispeu(const std::vector<double>& tx,const std::vector<double>& ty,const std::vector<double>& c,const double& x,const double& y){
  // z = bispeu(tx,ty,c,kx,ky,x,y)

  int nx = tx.size();
  int ny = ty.size();
  int kx = 3;
  int ky = 3;

  double z = fpbisp(tx, nx, ty, ny, c, kx, ky, x, y);

  return z;
}

int main(){
  std::cout << "main called" <<std::endl;

  std::vector<double> x = {1,2,3,4,5};
  std::vector<double> y = {1,3,4,7,10};
  int mx = x.size();
  int my = y.size();
  std::vector<std::vector<double>> z(mx, std::vector<double>(my, 0.0));

  for(int i = 0; i < mx; ++i){
    for(int j = 0; j < my; ++j){
      z[i][j] = x[i] * x[i] + y[j] * y[j] + x[i] * y[j];
      // std::printf("z(%d, %d) = %f\n", i, j, z[i][j]);
    }
  }

  std::vector<double> z_flattened(begin(z[0]), end(z[0]));
  for(int i = 1; i < mx; ++i){
    z_flattened.insert(end(z_flattened), begin(z[i]), end(z[i]));
  }

  regrid_return setup_spline = regrid_smth(x, y, z_flattened);

  std::cout << "nx  = " << setup_spline.nx  << std::endl;
  std::cout << "ny  = " << setup_spline.ny  << std::endl;
  std::cout << "tx  = " << setup_spline.tx  << std::endl;
  std::cout << "ty  = " << setup_spline.ty  << std::endl;
  std::cout << "C   = " << setup_spline.C   << std::endl;
  std::cout << "fp  = " << setup_spline.fp  << std::endl;
  std::cout << "ier = " << setup_spline.ier << std::endl;

  std::vector<double> tx = setup_spline.tx;
  std::vector<double> ty = setup_spline.ty;
  std::vector<double> C = setup_spline.C;

  double eval_x = 1.5;
  double eval_y = 5;

  double z_returned = bispeu(tx,ty,C,eval_x,eval_y);

  std::cout << z_returned << std::endl;
}
























