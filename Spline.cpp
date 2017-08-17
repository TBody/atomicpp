#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm> //for upper/lower_bound
#include <stdexcept> //For error-throwing
#include <sstream>

void fprota(double cos,double sin,double a,double b){
// c  subroutine fprota applies a givens rotation to a and b.
  double stor1 = a;
  double stor2 = b;
  b = cos * stor2 + sin * stor1;
  a = cos * stor1 - sin * stor2;
};
void fpgivs(double piv, double ww,double cos,double sin){
// c  subroutine fpgivs calculates the parameters of a given
// c  transformation .
  // 1 = 0.1e+01
  double dd;
  if(abs(piv)>=ww) {dd = abs(piv)*sqrt(1+(ww/piv)*(ww/piv));}
  if(abs(piv)<ww) {dd = ww*sqrt(1+(piv/ww)*(piv/ww));}
  cos = ww/dd;
  sin = piv/dd;
  ww = dd;
};

void fpbspl(std::vector<double> t,int n,int k,double x,int l,std::vector<double> h){

};

void fpgrre(
  int ifsx,
  int ifsy,
  int ifbx,
  int ifby,
  std::vector<double> x,
  int mx,
  std::vector<double> y,
  int my,
  std::vector<double> z,
  int mz,
  int kx,
  int ky,
  std::vector<double> tx,
  int nx,
  std::vector<double> ty,
  int ny,
  double p,
  std::vector<double> c,
  int nc,
  double fp,
  std::vector<double> fpx,
  std::vector<double> fpy,
  int mm,
  int mynx,
  int kx1,
  int kx2,
  int ky1,
  int ky2,
  std::vector<std::vector<double>> spx,
  std::vector<std::vector<double>> spy,
  std::vector<double> right,
  std::vector<double> q,
  std::vector<std::vector<double>> ax,
  std::vector<std::vector<double>> ay,
  std::vector<std::vector<double>> bx,
  std::vector<std::vector<double>> by,
  std::vector<int> nrx,
  std::vector<int> nry
  ){
// c  ..subroutine references..
// c    fpback,fpbspl,fpgivs,fpdisc,fprota
// c  ..
// c  the b-spline coefficients of the smoothing spline are calculated as
// c  the least-squares solution of the over-determined linear system of
// c  equations  (ay) c (ax)' = q       where
// c
// c               |   (spx)    |            |   (spy)    |
// c        (ax) = | ---------- |     (ay) = | ---------- |
// c               | (1/p) (bx) |            | (1/p) (by) |
// c
// c                                | z  ' 0 |
// c                            q = | ------ |
// c                                | 0  ' 0 |
// c
// c  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
// c                b-spline coefficients.
// c       z      : the my x mx matrix which contains the function values.
// c       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
// c                matrices according to the least-squares problems in
// c                the x- and y-direction.
// c       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
// c                matrices which contain the discontinuity jumps of the
// c                derivatives of the b-splines in the x- and y-direction.

std::vector<double> h(7, 0.0);

// double 1 = 1;
double half = 0.5;
double nk1x = nx-kx1;
double nk1y = ny-ky1;

double pinv = 0.0;
if(p > 0.){pinv = 1/p;}
// c  calculate the non-zero elements of the matrix (spx) which is the
// c  observation matrix according to the least-squares spline approximat-
// c  ion problem in the x-direction.
int l = kx1;
int l1 = kx2;
int number = 0;

for(int it = 0; it < mx; ++it){
  double arg = x[it];
  while(not(arg<tx[l1] or l==nk1x)){
    l = l1;
    l1 = l+1;
    number = number+1;
  }
  fpbspl(tx,nx,kx,arg,l,h);
  for(int i=0; i<kx1; ++i){
    spx[it][i] = h[i];
  }
  nrx[it] = number;
}
ifsx = 1;
// c  calculate the non-zero elements of the matrix (spy) which is the
// c  observation matrix according to the least-squares spline approximat-
// c  ion problem in the y-direction.
l = ky1;
l1 = ky2;
number = 0;

for(int it = 0; it < my; ++it){
  double arg = y[it];
  while(not(arg<ty[l1] or l==nk1y)){
    l = l1;
    l1 = l+1;
    number = number+1;
  }
  fpbspl(ty,ny,ky,arg,l,h);
  for(int i=0; i<ky1; ++i){
    spy[it][i] = h[i];
  }
  nry[it] = number;
}
ifsy = 1;

// c  reduce the matrix (ax) to upper triangular form (rx) using givens
// c  rotations. apply the same transformations to the rows of matrix q
// c  to obtain the my x (nx-kx-1) matrix g.
// c  store matrix (rx) into (ax) and g into q.
l = my*nk1x;
// c  initialization.
for(int i=0; i<l-1; ++i){
  q[i] = 0;
}
for(int i=0; i<nk1x-1; ++i){
  for(int j=0; i<kx2-1; ++i){
    ax[i][j] = 0;
  }
}

l = 0;
int nrold = 0;
// c  ibandx denotes the bandwidth of the matrices (ax) and (rx).
int ibandx = kx1;
int irot, n1;
double piv, sin, cos;
int i2, i3;

for(int it = 0; it< mx; ++it){
  number = nrx[it];
  // 50
  if ( nrold==number ){
    // !  fetch a new row of matrix (spx).
    h[ibandx] = 0.;
    for(int j = 1 ; j< kx1; ++j){
      h[j] = spx[it][j];
    };
    // !  find the appropriate column of q.
    for(int j = 0; j< my; ++j){
      l = l + 1;
      right[j] = z[l];
    };
    irot = number;
  } else {
    // if ( p<=0. ) print * , "ba fpgrre 10"
    // if ( p<=0. ){ goto 150};
    if ( p<=0. ){
      nrold = nrold + 1;
      continue;
    }
    ibandx = kx2;
    // !  fetch a new row of matrix (bx).
    n1 = nrold + 1;
    for(int j = 0; j < kx2; ++j){
      h[j] = bx[n1][j]*pinv;
    }
    // !  find the appropriate column of q.
    for(int j = 0; j < my; ++j){
      right[j] = 0.;
    }
    irot = nrold;
  };
  // !  rotate the new row of matrix (ax) into triangle.
  for(int i = 0; i< ibandx; ++i){
  irot = irot + 1;
  piv = h[i];
  // if ( piv==0. ) print * , "ba fpgrre 11"
  if ( piv/=0. ){
    // !  calculate the parameters of the givens transformation.
    fpgivs(piv,ax[irot][1],cos,sin);
    // !  apply that transformation to the rows of matrix q.
    int iq = (irot-1)*my;
    for(int j = 0; j< my; ++j){
      iq = iq + 1;
      fprota(cos,sin,right[j],q[iq]);
    }
    // !  apply that transformation to the columns of (ax).
    // if ( i==ibandx ) print * , "ba fpgrre 12"
    // if ( i==ibandx ){ goto 100};
    i2 = 1;
    i3 = i + 1;
    for(int j = i3-1; j< ibandx; ++j){
      i2 = i2 + 1;
      // call fprota(cos,sin,h[j],ax[irot][i2]);
    }
  };
  }
  // if ( nrold==number ) print * , "ba fpgrre 13"
  // 100
  // if ( nrold==number ) goto 200
  // 150
  nrold = nrold + 1;
  // goto 50
  // 200
}

// double arg,fac,term;
// int i,ibandy,ic,it,iz,i1,j,k,k1,k2,l2,ncof,nroldx,nroldy,numx,numx1,numy,numy1;

};

int fpregr(int iopt,
  std::vector<double> x,
  int mx,
  std::vector<double> y,
  int my,
  std::vector<double> z,
  int mz,
  double xb,
  double xe,
  double yb,
  double ye,
  int kx,
  int ky,
  double s,
  int nxest,
  int nyest,
  double tol,
  int maxit,
  int nc,
  int nx,
  std::vector<double> tx,
  int ny,
  std::vector<double> ty,
  std::vector<double> c,
  double fp,
  double fp0,
  double fpold,
  double reducx,
  double reducy,
  std::vector<double> fpintx,
  std::vector<double> fpinty,
  int lastdi,
  int nplusx,
  int nplusy,
  std::vector<int> nrx,
  std::vector<int> nry,
  std::vector<int> nrdatx,
  std::vector<int> nrdaty,
  std::vector<double> wrk,
  int lwrk
  ){
  // c   set constants
  // double 1 = 1;
  double half = 0.5e0;
  double con1 = 0.1e0;
  double con9 = 0.9e0;
  double con4 = 0.4e-01;
  // c  we partition the working space.
  int kx1 = kx+1;
  int ky1 = ky+1;
  int kx2 = kx1+1;
  int ky2 = ky1+1;
  int lsx = 1;
  int lsy = lsx+mx*kx1;
  int lri = lsy+my*ky1;
  int mm = std::max(nxest,my);
  int lq = lri+mm;
  int mynx = nxest*my;
  int lax = lq+mynx;
  int nxk = nxest*kx2;
  int lbx = lax+nxk;
  int lay = lbx+nxk;
  int lby = lay+nyest*ky2;

  // cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  // c part 1: determination of the number of knots and their position.     c
  // c ****************************************************************     c
  // c  given a set of knots we compute the least-squares spline sinf(x,y), c
  // c  and the corresponding sum of squared residuals fp=f(p=inf).         c
  // c  if iopt=-1  sinf(x,y) is the requested approximation.               c
  // c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
  // c    if fp <=s we will continue with the current set of knots.         c
  // c    if fp > s we will increase the number of knots and compute the    c
  // c       corresponding least-squares spline until finally fp<=s.        c
  // c    the initial choice of knots depends on the value of s and iopt.   c
  // c    if s=0 we have spline interpolation; in that case the number of   c
  // c    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c
  // c    if s>0 and                                                        c
  // c     *iopt=0 we first compute the least-squares polynomial of degree  c
  // c      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c
  // c     *iopt=1 we start with the knots found at the last call of the    c
  // c      routine, except for the case that s > fp0; then we can compute  c
  // c      the least-squares polynomial directly.                          c
  // cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  // c  determine the number of knots for polynomial approximation.

  int nminx = 2*kx1;
  int nminy = 2*ky1;

  // c  acc denotes the absolute tolerance for the root of f(p)=s.
  double acc = tol*s;
  // c  find nmaxx and nmaxy which denote the number of knots in x- and y-
  // c  direction in case of spline interpolation.
  int nmaxx = mx+kx1;
  int nmaxy = my+ky1;
  // c  find nxe and nye which denote the maximum number of knots
  // c  allowed in each direction
  int nxe = std::min(nmaxx,nxest);
  int nye = std::min(nmaxy,nyest);
  // c  if s = 0, s(x,y) is an interpolating spline.
  nx = nmaxx;
  ny = nmaxy;
  // c  find the position of the interior knots in case of interpolation.
  // c  the knots in the x-direction.
  int mk1 = mx-kx1;
  int k3 = kx/2;
  int i = kx1+1;
  int j = k3+2;

  for(int l = 0; l < mk1 -1; ++l){
  tx[i] = x[j];
  i = i+1;
  j = j+1;
  }

  // c  the knots in the y-direction.
  mk1 = my-ky1;
  k3 = ky/2;
  i = ky1+1;
  j = k3+2;

  for(int l = 0; l < mk1 -1; ++l){
  ty[i] = y[j];
  i = i+1;
  j = j+1;
  }

  int mpm = mx+my;
  int ifsx = 0;
  int ifsy = 0;
  int ifbx = 0;
  int ifby = 0;
  double p = -1;

  // c  main loop for the different sets of knots.mpm=mx+my is a save upper
  // c  bound for the number of trials.

  // for(int iter = 0; iter < mpm; ++iter){
  // c  find nrintx (nrinty) which is the number of knot intervals in the
  // c  x-direction (y-direction).
  int nrintx = nx-nminx+1;
  int nrinty = ny-nminy+1;
  // c  find ncof, the number of b-spline coefficients for the current set
  // c  of knots.
  int nk1x = nx-kx1;
  int nk1y = ny-ky1;
  int ncof = nk1x*nk1y;
  // c  find the position of the additional knots which are needed for the
  // c  b-spline representation of s(x,y).
  i = nx;
  for(j=1; j < kx1; ++j){
  tx[j-1] = xb;
  tx[i-1] = xe;
  i = i-1;
  }

  i = ny;
  for(j=1; j < ky1; ++j){
  ty[j-1] = yb;
  ty[i-1] = ye;
  i = i-1;
  }
  // c  find the least-squares spline sinf(x,y) and calculate for each knot
  // c  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
  // c  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
  // c  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
  // c  value belonging to that interval.
  // c  fp gives the total sum of squared residuals.

  // ifsx = ifsx;
  // ifsy = ifsy;
  // ifbx = ifbx;
  // ifby = ifby;
  // x = x;
  // mx = mx;
  // y = y;
  // my = my;
  // z = z;
  // mz = mz;
  // kx = kx;
  // ky = ky;
  // tx = tx;
  // nx = nx;
  // ty = ty;
  // ny = ny;
  // p = p;
  // c = c;
  // nc = nc;
  // fp = fp;
  std::vector<double>  fpx = fpintx;
  std::vector<double>  fpy = fpinty;
  // mm = mm;
  // mynx = mynx;
  // kx1 = kx1;
  // kx2 = kx2;
  // ky1 = ky1;
  // ky2 = ky2;
  std::vector<std::vector<double>> spx(mx, std::vector<double>(kx1, 0.0));
  // = wrk[lsx - 1];
  std::vector<std::vector<double>> spy(mx, std::vector<double>(ky1, 0.0));
  // = wrk[lsy - 1];
  std::vector<double> right(mm, 0.0);
  // = wrk[lri - 1];
  std::vector<double> q(mynx, 0.0);
  // = wrk[lq - 1];
  std::vector<std::vector<double>> ax(nx, std::vector<double>(kx2,0.0));
  // = wrk[lax - 1];
  std::vector<std::vector<double>> ay(ny, std::vector<double>(ky2,0.0));
  // = wrk[lay - 1];
  std::vector<std::vector<double>> bx(nx, std::vector<double>(kx2,0.0));
  // = wrk[lbx - 1];
  std::vector<std::vector<double>> by(ny, std::vector<double>(ky2,0.0));
  // = wrk[lby - 1];
  // std::vector<int> nrx(mx, 0);
  // = nrx;
  // std::vector<int> nry(my, 0);
  // = nry;

  fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,ax,ay,bx,by,nrx,nry);

  double fpms = fp-s;

  // c  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
  if((nx==nmaxx) and (ny==nmaxy)){
  int ier = -1;
  fp = 0.;
  return ier;
  } else {
  throw std::runtime_error("Not enough knots included in spline");
  }
  // }


  // double f1,f2,f3,p1,p2,p3,rn;
  // int ich1,ich3,iter,l,nplx,nply,npl1;
};

struct regrid_return{
  int nx;
  int ny;
  std::vector<double> tx;
  std::vector<double> ty;
  std::vector<double> c;
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
  c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
  uated by means of subroutine bispev.

  calling sequence:
  call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
  *  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

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
  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
  on succesful exit, c contains the coefficients of the spline
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

// nx,tx,ny,ty,c,fp,ier = regrid_smth(x,y,z,[xb,xe,yb,ye,kx,ky,s])
regrid_return regrid_smth(std::vector<double> x, std::vector<double> y, std::vector<double> z){

  // Translation of interpolate/src/fitpack.pyf (header)

  int iopt = 0;
  int mx = x.size();
  int my = y.size();

  std::printf("%d, %d \n", mx, my);

  int kx = 3;
  int ky = 3;

  if (not((mx > kx) and (my > ky))){
  throw std::runtime_error("Grid too small for bicubic interpolation");
  };

  double s = 0.0;

  double xb = x[0];
  double xe = x[mx-1];
  double yb = y[0];
  double ye = y[my-1];

  std::printf("%f, %f, %f, %f \n", xb, xe, yb, ye);

  int nxest = mx+kx+1;
  int nyest = my+ky+1;

  std::vector<double> tx(nxest, 0.0);
  std::vector<double> ty(nyest, 0.0);
  std::vector<double> c((nxest-kx-1)*(nyest-ky-1), 0.0);
  double fp = 0.0;

  int lwrk = 4 + nxest * (my + 2 * kx + 5) + nyest * (2 * ky + 5) + mx * (kx + 1) + my * (ky + 1) + std::max(my, nxest);
  std::vector<double> wrk(lwrk, 0.0);
  int kwrk = 3 + mx + my + nxest + nyest;
  std::vector<int> iwrk(kwrk, 0);

  int ier = 0;

  // End of translation of interpolate/src/fitpack.pyf (header)

  // call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

  // Translation of fitpack/regrid.f (implementation)

  int maxit = 20;
  double tol = 0.1e-02;
  // c  before starting computations a data check is made. if the input data
  // c  are invalid, control is immediately repassed to the calling program.
  ier = 10;
  if(kx <= 0 or kx > 5) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 1");

  int kx1 = kx+1;
  int kx2 = kx1+1;
  if(ky <= 0 or ky > 5) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 2");

  int ky1 = ky+1;
  int ky2 = ky1+1;
  if(iopt < (-1) or iopt > 1) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 3");

  int nminx = 2*kx1;
  if(mx < kx1 or nxest < nminx) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 4");

  int nminy = 2*ky1;
  if(my < ky1 or nyest < nminy) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 5");

  int mz = mx*my;
  int nc = (nxest-kx1)*(nyest-ky1);
  int lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+ my*ky1+std::max(nxest,my);
  int kwest = 3+mx+my+nxest+nyest;
  if(lwrk < lwest or kwrk < kwest) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 6");

  if(xb > x[1-1] or xe < x[mx-1]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 7");

  for(int i = 1; i<mx-1; ++i){
  if(x[i-1] >= x[i]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 8");
  }

  if(yb > y[1-1] or ye < y[my-1]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 9");

  for(int j = 1; j<=my-1; ++j){
  if(y[j-1] >= y[j]) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 10");
  }

  if(not(iopt >= 0)){
  throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 11");
  }

  if(s < 0.) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 12");

  if(s == 0. and (nxest < (mx+kx1) or nyest < (my+ky1)) ) throw std::runtime_error("Error in BicubicSpline/regrid_smth (C++ translation) - see code 13");

  ier = 0;
  // c  we partition the working space and determine the spline approximation;
  int lfpx = 5;
  int lfpy = lfpx+nxest;
  int lww = lfpy+nyest;
  int jwrk = lwrk-4-nxest-nyest;
  int knrx = 4;
  int knry = knrx+mx;
  int kndx = knry+my;
  int kndy = kndx+nxest;

  // iopt   = iopt;
  // x      = x;
  // mx     = mx;
  // y      = y;
  // my     = my;
  // z      = z;
  // mz     = mz;
  // xb     = xb;
  // xe     = xe;
  // yb     = yb;
  // ye     = ye;
  // kx     = kx;
  // ky     = ky;
  // s      = s;
  // nxest  = nxest;
  // nyest  = nyest;
  // tol    = tol;
  // maxit  = maxit;
  // nc     = nc;
  int nx     = 0;
  // tx     = tx;
  int ny     = 0;
  // ty     = ty;
  // c      = c;
  // fp     = fp;
  double fp0    = wrk[1-1];
  double fpold  = wrk[2-1];
  double reducx = wrk[3-1];
  double reducy = wrk[4-1];
  std::vector<double> fpintx(nxest, 0.0);
  // = wrk[lfpx-1];
  std::vector<double> fpinty(nyest, 0.0);
  // = wrk[lfpy-1];
  int                 lastdi = iwrk[1-1];
  int                 nplusx = iwrk[2-1];
  int                 nplusy = iwrk[3-1];
  std::vector<int>    nrx   (mx, 0);
  // = iwrk[knrx-1];
  std::vector<int>    nry   (my, 0);
  // = iwrk[knry-1];
  std::vector<int>    nrdatx(nxest, 0);
  // = iwrk[kndx-1];
  std::vector<int>    nrdaty(nyest, 0);
  // = iwrk[kndy-1];
  int                 lwrk_fpregr   = jwrk;
  std::vector<double> wrk_fpregr   (lwrk_fpregr, 0.0);
  // = wrk[lww-1];
  // int                 ier    = ier;

  ier = fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,
  nxest,nyest,tol,maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,
  reducy,fpintx,fpinty,lastdi,nplusx,nplusy,nrx,nry,nrdatx,nrdaty,
  wrk_fpregr,lwrk_fpregr);


  // call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
  // 	tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
  // 	wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx),
  // 	iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)

  // End of translation of fitpack/regrid.f (implementation)
};

int main(){
  std::printf("Hello world\n");

  std::vector<double> x = {1,2,3,4,5};
  std::vector<double> y = {1,3,4,7,10};
  int mx = x.size();
  int my = y.size();
  std::vector<std::vector<double>> z(mx, std::vector<double>(my, 0.0));

  for(int i = 0; i < mx; ++i){
    for(int j = 0; j < my; ++j){
      z[i][j] = x[i]*x[i] + y[j]*x[j] + x[i]*y[j];
    }
  }

  std::vector<double> z_flattened(begin(z[0]), end(z[0]));
  for(int i = 1; i < mx; ++i){
    z_flattened.insert(end(z_flattened), begin(z[i]), end(z[i]));
  }



  regrid_return setup_spline = regrid_smth(x, y, z_flattened);
}























