struct return_regrid{
      int nx;
      int ny;
      vector<double> tx;
      vector<double> ty;
      vector<double> c;
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
                     s.ne.0.
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

return_regrid regrid_smth(std::vector<double> x, std::vector<double> y, std::vector<std::vector<double>> z){
      
      // Translation of interpolate/src/fitpack.pyf
      // nx,tx,ny,ty,c,fp,ier = regrid_smth(x,y,z,[xb,xe,yb,ye,kx,ky,s])

      int iopt = 0;
      int mx = x.size();
      int my = y.size();

      int kx = 3;
      int ky = 3;  

      if (not((mx > kx) and (my > ky))){
            throw runtime_error("Grid too small for bicubic interpolation")
      };

      double s = 0.0;    

      double xb = std::min(x,mx);
      double xe = std::max(x,mx);
      double yb = std::min(y,my);
      double ye = std::max(y,my);      

      int nxest = mx+kx+1;
      int nyest = my+ky+1;

      std::vector<double> tx(nxest, 0.0);
      std::vector<double> ty(nyest, 0.0);
      std::vector<double> c((nxest-kx-1)*(nyest-ky-1), 0.0);
      double fp = 0.0;

      int lwrk = 4 + nxest * (my + 2 * kx + 5) + nyest * (2 * ky + 5) + mx * (kx + 1) + my * (ky + 1) + max(my, nxest);
      std::vector<double> wrk(lwrk, 0.0);
      int kwrk = 3 + mx + my + nxest + nyest;
      std::vector<int> iwrk(kwrk, 0);

      int ier = 0;
      // End of translation of interpolate/src/fitpack.pyf

      // Start of translation of interpolate/fitpack/regrid.f

      // ..scalar arguments..
      // double xb,xe,yb,ye,s,fp;
      // int iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier;
      // ..array arguments..
      // double x(mx),y(my),z(mx*my),tx(nxest),ty(nyest), c((nxest-kx-1)*(nyest-ky-1)),wrk(lwrk);
      // int iwrk(kwrk);
      // ..local scalars..
      double tol;
      int i,j,jwrk,kndx,kndy,knrx,knry,kwest,kx1,kx2,ky1,ky2,lfpx,lfpy,lwest,lww,maxit,nc,nminx,nminy,mz;
      //  ..function references..
      int max0;

      // we set up the parameters tol and maxit.
      maxit = 20;
      tol = 0.1e-02;

      // before starting computations a data check is made. if the input data are invalid, control is immediately repassed to the calling program.
      // Try do all these tests in wrapper

      ier = 10;

      // if(kx <= 0 .or. kx > 5) throw runtime_error();
      kx1 = kx+1
      kx2 = kx1+1
      // if(ky <= 0 .or. ky > 5) throw runtime_error();
      ky1 = ky+1
      ky2 = ky1+1
      // if(iopt < (-1) .or. iopt > 1) throw runtime_error();
      nminx = 2*kx1
      // if(mx < kx1 .or. nxest < nminx) throw runtime_error();
      nminy = 2*ky1
      // if(my < ky1 .or. nyest < nminy) throw runtime_error();
      mz = mx*my
      nc = (nxest-kx1)*(nyest-ky1)
      lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+my*ky1+max0(nxest,my)
      kwest = 3+mx+my+nxest+nyest
      // if(lwrk < lwest .or. kwrk < kwest) throw runtime_error();
      // if(xb > x(1) .or. xe < x(mx)) throw runtime_error();
  //     do 10 i=2,mx
  //       if(x(i-1)>=x(i)) throw runtime_error();
  // 10  continue
      // if(yb > y(1) .or. ye < y(my)) throw runtime_error();
  //     do 20 i=2,my
  //       if(y(i-1)>=y(i)) throw runtime_error();
  // 20  continue
      // if(iopt>=0) go to 50 //If use smoothing
      // if(nx < nminx .or. nx > nxest) throw runtime_error();
      j = nx
      do 30 i=1,kx1
        tx(i) = xb
        tx(j) = xe
        j = j-1
  30  continue
      call fpchec(x,mx,tx,nx,kx,ier)
      // if(ier.ne.0) throw runtime_error();
      // if(ny < nminy .or. ny > nyest) throw runtime_error();
      j = ny
      do 40 i=1,ky1
        ty(i) = yb
        ty(j) = ye
        j = j-1
  40  continue
      call fpchec(y,my,ty,ny,ky,ier)
      if (ier.eq.0) go to 60 //Continue program if fpchec passes
      else throw runtime_error(); //otherwise crash
      // subroutine fpchec verifies the number and the position of the knots
      // t(j),j=1,2,...,n of a spline of degree k, in relation to the number
      // and the position of the data points x(i),i=1,2,...,m. if all of the
      // following conditions are fulfilled, the error parameter ier is set
      // to zero. if one of the conditions is violated ier is set to ten.
          // 1) k+1 <= n-k-1 <= m
          // 2) t(1) <= t(2) <= ... <= t(k+1)
             // t(n-k) <= t(n-k+1) <= ... <= t(n)
          // 3) t(k+1) < t(k+2) < ... < t(n-k)
          // 4) t(k+1) <= x(i) <= t(n-k)
          // 5) the conditions specified by schoenberg and whitney must hold
             // for at least one subset of data points, i.e. there must be a
             // subset of data points y(j) such that
                 // t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1

  // 50  if(s < 0.) throw runtime_error();
      // if(s.eq.0. .and. (nxest < (mx+kx1) .or. nyest < (my+ky1)) ) throw runtime_error();
      ier = 0
// c  we partition the working space and determine the spline approximation
  60  lfpx = 5
      lfpy = lfpx+nxest
      lww = lfpy+nyest
      jwrk = lwrk-4-nxest-nyest
      knrx = 4
      knry = knrx+mx
      kndx = knry+my
      kndy = kndx+nxest

      call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
            tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
            wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx),
            iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)
  // 70  return
}
















