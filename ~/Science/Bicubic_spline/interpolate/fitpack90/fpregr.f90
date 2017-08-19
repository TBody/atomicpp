!*==FPREGR.spg  processed by SPAG 6.72Dc at 11:08 on 17 Aug 2017
      SUBROUTINE FPREGR(Iopt,X,Mx,Y,My,Z,Mz,Xb,Xe,Yb,Ye,Kx,Ky,S,Nxest,  &
                      & Nyest,Tol,Maxit,Nc,Nx,Tx,Ny,Ty,C,Fp,Fp0,Fpold,  &
                      & Reducx,Reducy,Fpintx,Fpinty,Lastdi,Nplusx,      &
                      & Nplusy,Nrx,Nry,Nrdatx,Nrdaty,Wrk,Lwrk,Ier)
      IMPLICIT NONE
!*--FPREGR7
!  ..
!  ..scalar arguments..
      REAL*8 Xb , Xe , Yb , Ye , S , Tol , Fp , Fp0 , Fpold , Reducx ,  &
           & Reducy
      INTEGER Iopt , Mx , My , Mz , Kx , Ky , Nxest , Nyest , Maxit ,   &
            & Nc , Nx , Ny , Lastdi , Nplusx , Nplusy , Lwrk , Ier
!  ..array arguments..
      REAL*8 X(Mx) , Y(My) , Z(Mz) , Tx(Nxest) , Ty(Nyest) , C(Nc) ,    &
           & Fpintx(Nxest) , Fpinty(Nyest) , Wrk(Lwrk)
      INTEGER Nrdatx(Nxest) , Nrdaty(Nyest) , Nrx(Mx) , Nry(My)
!  ..local scalars
      REAL*8 acc , fpms , f1 , f2 , f3 , p , p1 , p2 , p3 , rn , one ,  &
           & half , con1 , con9 , con4
      INTEGER i , ich1 , ich3 , ifbx , ifby , ifsx , ifsy , iter , j ,  &
            & kx1 , kx2 , ky1 , ky2 , k3 , l , lax , lay , lbx , lby ,  &
            & lq , lri , lsx , lsy , mk1 , mm , mpm , mynx , ncof ,     &
            & nk1x , nk1y , nmaxx , nmaxy , nminx , nminy , nplx ,      &
            & nply , npl1 , nrintx , nrinty , nxe , nxk , nye
!  ..function references..
      REAL*8 ABS , FPRATI
      INTEGER MAX0 , MIN0
!  ..subroutine references..
!    fpgrre,fpknot
!  ..
      PRINT * , "fpregr called"
!   set constants
      one = 1
      half = 0.5E0
      con1 = 0.1E0
      con9 = 0.9E0
      con4 = 0.4E-01
!  we partition the working space.
      kx1 = Kx + 1
      ky1 = Ky + 1
      kx2 = kx1 + 1
      ky2 = ky1 + 1
      lsx = 1
      lsy = lsx + Mx*kx1
      lri = lsy + My*ky1
      mm = MAX0(Nxest,My)
      lq = lri + mm
      mynx = Nxest*My
      lax = lq + mynx
      nxk = Nxest*kx2
      lbx = lax + nxk
      lay = lbx + nxk
      lby = lay + Nyest*ky2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 1: determination of the number of knots and their position.     c
! ****************************************************************     c
!  given a set of knots we compute the least-squares spline sinf(x,y), c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1  sinf(x,y) is the requested approximation.               c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!    if fp <=s we will continue with the current set of knots.         c
!    if fp > s we will increase the number of knots and compute the    c
!       corresponding least-squares spline until finally fp<=s.        c
!    the initial choice of knots depends on the value of s and iopt.   c
!    if s=0 we have spline interpolation; in that case the number of   c
!    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c
!    if s>0 and                                                        c
!     *iopt=0 we first compute the least-squares polynomial of degree  c
!      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c
!     *iopt=1 we start with the knots found at the last call of the    c
!      routine, except for the case that s > fp0; then we can compute  c
!      the least-squares polynomial directly.                          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine the number of knots for polynomial approximation.
      nminx = 2*kx1
      nminy = 2*ky1
      IF ( Iopt>=0 ) THEN
!  acc denotes the absolute tolerance for the root of f(p)=s.
         acc = Tol*S
!  find nmaxx and nmaxy which denote the number of knots in x- and y-
!  direction in case of spline interpolation.
         nmaxx = Mx + kx1
         nmaxy = My + ky1
!  find nxe and nye which denote the maximum number of knots
!  allowed in each direction
         nxe = MIN0(nmaxx,Nxest)
         nye = MIN0(nmaxy,Nyest)
         IF ( S>0. ) THEN
!  if s > 0 our initial choice of knots depends on the value of iopt.
            IF ( Iopt/=0 ) THEN
               IF ( Fp0>S ) THEN
!  if iopt=1 and fp0 > s we start computing the least- squares spline
!  according to the set of knots found at the last call of the routine.
!  we determine the number of grid coordinates x(i) inside each knot
!  interval (tx(l),tx(l+1)).
                  l = kx2
                  j = 1
                  Nrdatx(1) = 0
                  mpm = Mx - 1
                  DO i = 2 , mpm
                     Nrdatx(j) = Nrdatx(j) + 1
                     IF ( X(i)>=Tx(l) ) THEN
                        Nrdatx(j) = Nrdatx(j) - 1
                        l = l + 1
                        j = j + 1
                        Nrdatx(j) = 0
                     ENDIF
                  ENDDO
!  we determine the number of grid coordinates y(i) inside each knot
!  interval (ty(l),ty(l+1)).
                  l = ky2
                  j = 1
                  Nrdaty(1) = 0
                  mpm = My - 1
                  DO i = 2 , mpm
                     Nrdaty(j) = Nrdaty(j) + 1
                     IF ( Y(i)>=Ty(l) ) THEN
                        Nrdaty(j) = Nrdaty(j) - 1
                        l = l + 1
                        j = j + 1
                        Nrdaty(j) = 0
                     ENDIF
                  ENDDO
                  GOTO 100
               ENDIF
            ENDIF
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial of degree kx in x and ky in y (which is a spline without
!  interior knots).
            Nx = nminx
            Ny = nminy
            Nrdatx(1) = Mx - 2
            Nrdaty(1) = My - 2
            Lastdi = 0
            Nplusx = 0
            Nplusy = 0
            Fp0 = 0.
            Fpold = 0.
            Reducx = 0.
            Reducy = 0.
         ELSE
!  if s = 0, s(x,y) is an interpolating spline.
            Nx = nmaxx
            Ny = nmaxy
!  test whether the required storage space exceeds the available one.
            IF ( Ny>Nyest .OR. Nx>Nxest ) GOTO 700
!  find the position of the interior knots in case of interpolation.
!  the knots in the x-direction.
            mk1 = Mx - kx1
            IF ( mk1/=0 ) THEN
               k3 = Kx/2
               i = kx1 + 1
               j = k3 + 2
               IF ( k3*2==Kx ) THEN
                  DO l = 1 , mk1
                     Tx(i) = (X(j)+X(j-1))*half
                     i = i + 1
                     j = j + 1
                  ENDDO
               ELSE
                  DO l = 1 , mk1
                     Tx(i) = X(j)
                     i = i + 1
                     j = j + 1
                  ENDDO
               ENDIF
            ENDIF
!  the knots in the y-direction.
            mk1 = My - ky1
            IF ( mk1/=0 ) THEN
               k3 = Ky/2
               i = ky1 + 1
               j = k3 + 2
               IF ( k3*2==Ky ) THEN
                  DO l = 1 , mk1
                     Ty(i) = (Y(j)+Y(j-1))*half
                     i = i + 1
                     j = j + 1
                  ENDDO
               ELSE
                  DO l = 1 , mk1
                     Ty(i) = Y(j)
                     i = i + 1
                     j = j + 1
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 100  mpm = Mx + My
      ifsx = 0
      ifsy = 0
      ifbx = 0
      ifby = 0
      p = -one
!  main loop for the different sets of knots.mpm=mx+my is a save upper
!  bound for the number of trials.
      DO iter = 1 , mpm
         IF ( Nx==nminx .AND. Ny==nminy ) Ier = -2
!  find nrintx (nrinty) which is the number of knot intervals in the
!  x-direction (y-direction).
         nrintx = Nx - nminx + 1
         nrinty = Ny - nminy + 1
!  find ncof, the number of b-spline coefficients for the current set
!  of knots.
         nk1x = Nx - kx1
         nk1y = Ny - ky1
         ncof = nk1x*nk1y
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(x,y).
         i = Nx
         DO j = 1 , kx1
            Tx(j) = Xb
            Tx(i) = Xe
            i = i - 1
         ENDDO
         i = Ny
         DO j = 1 , ky1
            Ty(j) = Yb
            Ty(i) = Ye
            i = i - 1
         ENDDO
!  find the least-squares spline sinf(x,y) and calculate for each knot
!  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
!  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
!  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
!  value belonging to that interval.
!  fp gives the total sum of squared residuals.
         CALL FPGRRE(ifsx,ifsy,ifbx,ifby,X,Mx,Y,My,Z,Mz,Kx,Ky,Tx,Nx,Ty, &
                   & Ny,p,C,Nc,Fp,Fpintx,Fpinty,mm,mynx,kx1,kx2,ky1,ky2,&
                   & Wrk(lsx),Wrk(lsy),Wrk(lri),Wrk(lq),Wrk(lax),       &
                   & Wrk(lay),Wrk(lbx),Wrk(lby),Nrx,Nry)
         IF ( Ier==(-2) ) Fp0 = Fp
!  test whether the least-squares spline is an acceptable solution.
         IF ( Iopt<0 ) GOTO 99999
         fpms = Fp - S
         IF ( ABS(fpms)<acc ) GOTO 99999
!  if f(p=inf) < s, we accept the choice of knots.
         IF ( fpms<0. ) GOTO 400
!  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
         IF ( Nx==nmaxx .AND. Ny==nmaxy ) GOTO 800
!  increase the number of knots.
!  if nx=nxe and ny=nye we cannot further increase the number of knots
!  because of the storage capacity limitation.
         IF ( Nx==nxe .AND. Ny==nye ) GOTO 700
         Ier = 0
!  adjust the parameter reducx or reducy according to the direction
!  in which the last added knots were located.
         IF ( Lastdi<0 ) THEN
            Reducx = Fpold - Fp
         ELSEIF ( Lastdi/=0 ) THEN
            Reducy = Fpold - Fp
         ENDIF
!  store the sum of squared residuals for the current set of knots.
         Fpold = Fp
!  find nplx, the number of knots we should add in the x-direction.
         nplx = 1
         IF ( Nx/=nminx ) THEN
            npl1 = Nplusx*2
            rn = Nplusx
            IF ( Reducx>acc ) npl1 = rn*fpms/Reducx
            nplx = MIN0(Nplusx*2,MAX0(npl1,Nplusx/2,1))
         ENDIF
!  find nply, the number of knots we should add in the y-direction.
         nply = 1
         IF ( Ny/=nminy ) THEN
            npl1 = Nplusy*2
            rn = Nplusy
            IF ( Reducy>acc ) npl1 = rn*fpms/Reducy
            nply = MIN0(Nplusy*2,MAX0(npl1,Nplusy/2,1))
         ENDIF
         IF ( nplx>=nply ) THEN
            IF ( nplx/=nply ) GOTO 200
            IF ( Lastdi<0 ) GOTO 200
         ENDIF
 150     IF ( Nx/=nxe ) THEN
!  addition in the x-direction.
            Lastdi = -1
            Nplusx = nplx
            ifsx = 0
            DO l = 1 , Nplusx
!  add a new knot in the x-direction
               CALL FPKNOT(X,Mx,Tx,Nx,Fpintx,Nrdatx,nrintx,Nxest,1)
!  test whether we cannot further increase the number of knots in the
!  x-direction.
               IF ( Nx==nxe ) GOTO 300
            ENDDO
            GOTO 300
         ENDIF
 200     IF ( Ny==nye ) GOTO 150
!  addition in the y-direction.
         Lastdi = 1
         Nplusy = nply
         ifsy = 0
         DO l = 1 , Nplusy
!  add a new knot in the y-direction.
            CALL FPKNOT(Y,My,Ty,Ny,Fpinty,Nrdaty,nrinty,Nyest,1)
!  test whether we cannot further increase the number of knots in the
!  y-direction.
            IF ( Ny==nye ) GOTO 300
         ENDDO
!  restart the computations with the new set of knots.
 300  ENDDO
!  test whether the least-squares polynomial is a solution of our
!  approximation problem.
 400  IF ( Ier/=(-2) ) THEN
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing spline sp(x,y)                c
! *****************************************************                c
!  we have determined the number of knots and their position. we now   c
!  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
!  this smoothing spline varies with the parameter p in such a way thatc
!    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c
!  is a continuous, strictly decreasing function of p. moreover the    c
!  least-squares polynomial corresponds to p=0 and the least-squares   c
!  spline to p=infinity. iteratively we then have to determine the     c
!  positive value of p such that f(p)=s. the process which is proposed c
!  here makes use of rational interpolation. f(p) is approximated by a c
!  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
!  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
!  are used to calculate the new value of p such that r(p)=s.          c
!  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  initial value for p.
         p1 = 0.
         f1 = Fp0 - S
         p3 = -one
         f3 = fpms
         p = one
         ich1 = 0
         ich3 = 0
!  iteration process to find the root of f(p)=s.
         DO iter = 1 , Maxit
!  find the smoothing spline sp(x,y) and the corresponding sum of
!  squared residuals fp.
            CALL FPGRRE(ifsx,ifsy,ifbx,ifby,X,Mx,Y,My,Z,Mz,Kx,Ky,Tx,Nx, &
                      & Ty,Ny,p,C,Nc,Fp,Fpintx,Fpinty,mm,mynx,kx1,kx2,  &
                      & ky1,ky2,Wrk(lsx),Wrk(lsy),Wrk(lri),Wrk(lq),     &
                      & Wrk(lax),Wrk(lay),Wrk(lbx),Wrk(lby),Nrx,Nry)
!  test whether the approximation sp(x,y) is an acceptable solution.
            fpms = Fp - S
            IF ( ABS(fpms)<acc ) GOTO 99999
!  test whether the maximum allowable number of iterations has been
!  reached.
            IF ( iter==Maxit ) GOTO 500
!  carry out one more step of the iteration process.
            p2 = p
            f2 = fpms
            IF ( ich3==0 ) THEN
               IF ( (f2-f3)>acc ) THEN
                  IF ( f2<0. ) ich3 = 1
               ELSE
!  our initial choice of p is too large.
                  p3 = p2
                  f3 = f2
                  p = p*con4
                  IF ( p<=p1 ) p = p1*con9 + p2*con1
                  GOTO 450
               ENDIF
            ENDIF
            IF ( ich1==0 ) THEN
               IF ( (f1-f2)>acc ) THEN
!  test whether the iteration process proceeds as theoretically
!  expected.
                  IF ( f2>0. ) ich1 = 1
               ELSE
!  our initial choice of p is too small
                  p1 = p2
                  f1 = f2
                  p = p/con4
                  IF ( p3>=0. ) THEN
                     IF ( p>=p3 ) p = p2*con1 + p3*con9
                  ENDIF
                  GOTO 450
               ENDIF
            ENDIF
            IF ( f2>=f1 .OR. f2<=f3 ) GOTO 600
!  find the new value of p.
            p = FPRATI(p1,f1,p2,f2,p3,f3)
 450     ENDDO
!  error codes and messages.
 500     Ier = 3
      ENDIF
      GOTO 99999
 600  Ier = 2
      GOTO 99999
 700  Ier = 1
      GOTO 99999
 800  Ier = -1
      Fp = 0.
99999 END
 
