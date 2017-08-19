!*==FPGRRE.spg  processed by SPAG 6.72Dc at 14:16 on 17 Aug 2017
      SUBROUTINE FPGRRE(Ifsx,Ifsy,Ifbx,Ifby,X,Mx,Y,My,Z,Mz,Kx,Ky,Tx,Nx, &
                      & Ty,Ny,P,C,Nc,Fp,Fpx,Fpy,Mm,Mynx,Kx1,Kx2,Ky1,Ky2,&
                      & Spx,Spy,Right,Q,Ax,Ay,Bx,By,Nrx,Nry)
      IMPLICIT NONE
!*--FPGRRE6
!  ..
!  ..scalar arguments..
      REAL*8 P , Fp
      INTEGER Ifsx , Ifsy , Ifbx , Ifby , Mx , My , Mz , Kx , Ky , Nx , &
            & Ny , Nc , Mm , Mynx , Kx1 , Kx2 , Ky1 , Ky2
!  ..array arguments..
      REAL*8 X(Mx) , Y(My) , Z(Mz) , Tx(Nx) , Ty(Ny) , C(Nc) ,          &
           & Spx(Mx,Kx1) , Spy(My,Ky1) , Right(Mm) , Q(Mynx) ,          &
           & Ax(Nx,Kx2) , Bx(Nx,Kx2) , Ay(Ny,Ky2) , By(Ny,Ky2) , Fpx(Nx)&
           & , Fpy(Ny)
      INTEGER Nrx(Mx) , Nry(My)
!  ..local scalars..
      REAL*8 arg , cos , fac , pinv , piv , sin , term , one , half
      INTEGER i , ibandx , ibandy , ic , iq , irot , it , iz , i1 , i2 ,&
            & i3 , j , k , k1 , k2 , l , l1 , l2 , ncof , nk1x , nk1y , &
            & nrold , nroldx , nroldy , number , numx , numx1 , numy ,  &
            & numy1 , n1
!  ..local arrays..
      REAL*8 h(7)
      PRINT * , "fpgrre called"
!  ..subroutine references..
!    fpback,fpbspl,fpgivs,fpdisc,fprota
!  ..
!  the b-spline coefficients of the smoothing spline are calculated as
!  the least-squares solution of the over-determined linear system of
!  equations  (ay) c (ax)' = q       where
!
!               |   (spx)    |            |   (spy)    |
!        (ax) = | ---------- |     (ay) = | ---------- |
!               | (1/p) (bx) |            | (1/p) (by) |
!
!                                | z  ' 0 |
!                            q = | ------ |
!                                | 0  ' 0 |
!
!  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
!                b-spline coefficients.
!       z      : the my x mx matrix which contains the function values.
!       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
!                matrices according to the least-squares problems in
!                the x- and y-direction.
!       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
!                matrices which contain the discontinuity jumps of the
!                derivatives of the b-splines in the x- and y-direction.
      one = 1
      half = 0.5
      nk1x = Nx - Kx1
      nk1y = Ny - Ky1
      IF ( P>0. ) PRINT * , "BA FPGRRE 1"
      IF ( P>0. ) pinv = one/P
!  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
!  the value of p whether the matrices (spx),(spy),(bx) and (by) still
!  must be determined.
      IF ( Ifsx/=0 ) PRINT * , "BA FPGRRE 2"
      IF ( Ifsx==0 ) THEN
!  calculate the non-zero elements of the matrix (spx) which is the
!  observation matrix according to the least-squares spline approximat-
!  ion problem in the x-direction.
         l = Kx1
         l1 = Kx2
         number = 0
         DO it = 1 , Mx
            arg = X(it)
            IF ( arg<Tx(l1) .OR. l==nk1x ) PRINT * , "BA FPGRRE 3"
 20         IF ( arg<Tx(l1) .OR. l==nk1x ) THEN
               CALL FPBSPL(Tx,Nx,Kx,arg,l,h)
               DO i = 1 , Kx1
                  Spx(it,i) = h(i)
               ENDDO
               Nrx(it) = number
            ELSE
               l = l1
               l1 = l + 1
               number = number + 1
               GOTO 20
            ENDIF
         ENDDO
         Ifsx = 1
         IF ( Ifsy/=0 ) PRINT * , "BA FPGRRE 4"
      ENDIF
      IF ( Ifsy==0 ) THEN
!  calculate the non-zero elements of the matrix (spy) which is the
!  observation matrix according to the least-squares spline approximat-
!  ion problem in the y-direction.
         l = Ky1
         l1 = Ky2
         number = 0
         DO it = 1 , My
            arg = Y(it)
            IF ( arg<Ty(l1) .OR. l==nk1y ) PRINT * , "BA FPGRRE 5"
 40         IF ( arg<Ty(l1) .OR. l==nk1y ) THEN
               CALL FPBSPL(Ty,Ny,Ky,arg,l,h)
               DO i = 1 , Ky1
                  Spy(it,i) = h(i)
               ENDDO
               Nry(it) = number
            ELSE
               l = l1
               l1 = l + 1
               number = number + 1
               GOTO 40
            ENDIF
         ENDDO
         Ifsy = 1
         IF ( P<=0. ) PRINT * , "BA FPGRRE 6"
      ENDIF
      IF ( P>0. ) THEN
!  calculate the non-zero elements of the matrix (bx).
         IF ( Ifbx/=0 .OR. Nx==2*Kx1 ) PRINT * , "BA FPGRRE 7"
         IF ( Ifbx==0 .AND. Nx/=2*Kx1 ) THEN
            CALL FPDISC(Tx,Nx,Kx2,Bx,Nx)
            Ifbx = 1
!  calculate the non-zero elements of the matrix (by).
            IF ( Ifby/=0 .OR. Ny==2*Ky1 ) PRINT * , "BA FPGRRE 8"
         ENDIF
         IF ( Ifby==0 .AND. Ny/=2*Ky1 ) THEN
            CALL FPDISC(Ty,Ny,Ky2,By,Ny)
            Ifby = 1
         ENDIF
      ENDIF
!  reduce the matrix (ax) to upper triangular form (rx) using givens
!  rotations. apply the same transformations to the rows of matrix q
!  to obtain the my x (nx-kx-1) matrix g.
!  store matrix (rx) into (ax) and g into q.
      l = My*nk1x
!  initialization.
      DO i = 1 , l
         Q(i) = 0.
      ENDDO
      DO i = 1 , nk1x
         DO j = 1 , Kx2
            Ax(i,j) = 0.
         ENDDO
      ENDDO
      l = 0
      nrold = 0
!  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = Kx1
      DO it = 1 , Mx
         number = Nrx(it)
         IF ( nrold==number ) PRINT * , "BA FPGRRE 9"
 50      IF ( nrold==number ) THEN
!  fetch a new row of matrix (spx).
            h(ibandx) = 0.
            DO j = 1 , Kx1
               h(j) = Spx(it,j)
            ENDDO
!  find the appropriate column of q.
            DO j = 1 , My
               l = l + 1
               Right(j) = Z(l)
            ENDDO
            irot = number
         ELSE
            IF ( P<=0. ) PRINT * , "BA FPGRRE 10"
            IF ( P<=0. ) GOTO 150
            ibandx = Kx2
!  fetch a new row of matrix (bx).
            n1 = nrold + 1
            DO j = 1 , Kx2
               h(j) = Bx(n1,j)*pinv
            ENDDO
!  find the appropriate column of q.
            DO j = 1 , My
               Right(j) = 0.
            ENDDO
            irot = nrold
         ENDIF
!  rotate the new row of matrix (ax) into triangle.
         DO i = 1 , ibandx
            irot = irot + 1
            piv = h(i)
            IF ( piv==0. ) PRINT * , "BA FPGRRE 11"
            IF ( piv/=0. ) THEN
!  calculate the parameters of the givens transformation.
               CALL FPGIVS(piv,Ax(irot,1),cos,sin)
!  apply that transformation to the rows of matrix q.
               iq = (irot-1)*My
               DO j = 1 , My
                  iq = iq + 1
                  CALL FPROTA(cos,sin,Right(j),Q(iq))
               ENDDO
!  apply that transformation to the columns of (ax).
               IF ( i==ibandx ) PRINT * , "BA FPGRRE 12"
               IF ( i==ibandx ) GOTO 100
               i2 = 1
               i3 = i + 1
               DO j = i3 , ibandx
                  i2 = i2 + 1
                  CALL FPROTA(cos,sin,h(j),Ax(irot,i2))
               ENDDO
            ENDIF
         ENDDO
         IF ( nrold==number ) PRINT * , "BA FPGRRE 13"
 100     IF ( nrold==number ) GOTO 200
 150     nrold = nrold + 1
         GOTO 50
 200  ENDDO
!  reduce the matrix (ay) to upper triangular form (ry) using givens
!  rotations. apply the same transformations to the columns of matrix g
!  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
!  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
!  initialization.
      DO i = 1 , ncof
         C(i) = 0.
      ENDDO
      DO i = 1 , nk1y
         DO j = 1 , Ky2
            Ay(i,j) = 0.
         ENDDO
      ENDDO
      nrold = 0
!  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = Ky1
      DO it = 1 , My
         number = Nry(it)
         IF ( nrold==number ) PRINT * , "BA FPGRRE 14"
 250     IF ( nrold==number ) THEN
!  fetch a new row of matrix (spy)
            h(ibandy) = 0.
            DO j = 1 , Ky1
               h(j) = Spy(it,j)
            ENDDO
!  find the appropiate row of g.
            l = it
            DO j = 1 , nk1x
               Right(j) = Q(l)
               l = l + My
            ENDDO
            irot = number
         ELSE
            IF ( P<=0. ) PRINT * , "BA FPGRRE 15"
            IF ( P<=0. ) GOTO 350
            ibandy = Ky2
!  fetch a new row of matrix (by).
            n1 = nrold + 1
            DO j = 1 , Ky2
               h(j) = By(n1,j)*pinv
            ENDDO
!  find the appropiate row of g.
            DO j = 1 , nk1x
               Right(j) = 0.
            ENDDO
            irot = nrold
         ENDIF
!  rotate the new row of matrix (ay) into triangle.
         DO i = 1 , ibandy
            irot = irot + 1
            piv = h(i)
            IF ( piv==0. ) PRINT * , "BA FPGRRE 16"
            IF ( piv/=0. ) THEN
!  calculate the parameters of the givens transformation.
               CALL FPGIVS(piv,Ay(irot,1),cos,sin)
!  apply that transformation to the colums of matrix g.
               ic = irot
               DO j = 1 , nk1x
                  CALL FPROTA(cos,sin,Right(j),C(ic))
                  ic = ic + nk1y
               ENDDO
!  apply that transformation to the columns of matrix (ay).
               IF ( i==ibandy ) PRINT * , "BA FPGRRE 17"
               IF ( i==ibandy ) GOTO 300
               i2 = 1
               i3 = i + 1
               DO j = i3 , ibandy
                  i2 = i2 + 1
                  CALL FPROTA(cos,sin,h(j),Ay(irot,i2))
               ENDDO
            ENDIF
         ENDDO
         IF ( nrold==number ) PRINT * , "BA FPGRRE 18"
 300     IF ( nrold==number ) GOTO 400
 350     nrold = nrold + 1
         GOTO 250
 400  ENDDO
!  backward substitution to obtain the b-spline coefficients as the
!  solution of the linear system    (ry) c (rx)' = h.
!  first step: solve the system  (ry) (c1) = h.
      k = 1
      DO i = 1 , nk1x
         CALL FPBACK(Ay,C(k),nk1y,ibandy,C(k),Ny)
         k = k + nk1y
      ENDDO
!  second step: solve the system  c (rx)' = (c1).
      k = 0
      DO j = 1 , nk1y
         k = k + 1
         l = k
         DO i = 1 , nk1x
            Right(i) = C(l)
            l = l + nk1y
         ENDDO
         CALL FPBACK(Ax,Right,nk1x,ibandx,Right,Nx)
         l = k
         DO i = 1 , nk1x
            C(l) = Right(i)
            l = l + nk1y
         ENDDO
      ENDDO
!  calculate the quantities
!    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
!    fp = sumi=1,mx(sumj=1,my(res(i,j)))
!    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
!                  tx(r+kx) <= x(i) <= tx(r+kx+1)
!    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
!                  ty(r+ky) <= y(j) <= ty(r+ky+1)
      Fp = 0.
      DO i = 1 , Nx
         Fpx(i) = 0.
      ENDDO
      DO i = 1 , Ny
         Fpy(i) = 0.
      ENDDO
      nk1y = Ny - Ky1
      iz = 0
      nroldx = 0
!  main loop for the different grid points.
      DO i1 = 1 , Mx
         numx = Nrx(i1)
         numx1 = numx + 1
         nroldy = 0
         DO i2 = 1 , My
            numy = Nry(i2)
            numy1 = numy + 1
            iz = iz + 1
!  evaluate s(x,y) at the current grid point by making the sum of the
!  cross products of the non-zero b-splines at (x,y), multiplied with
!  the appropiate b-spline coefficients.
            term = 0.
            k1 = numx*nk1y + numy
            DO l1 = 1 , Kx1
               k2 = k1
               fac = Spx(i1,l1)
               DO l2 = 1 , Ky1
                  k2 = k2 + 1
                  term = term + fac*Spy(i2,l2)*C(k2)
               ENDDO
               k1 = k1 + nk1y
            ENDDO
!  calculate the squared residual at the current grid point.
            term = (Z(iz)-term)**2
!  adjust the different parameters.
            Fp = Fp + term
            Fpx(numx1) = Fpx(numx1) + term
            Fpy(numy1) = Fpy(numy1) + term
            fac = term*half
            IF ( numy==nroldy ) PRINT * , "BA FPGRRE 19"
            IF ( numy/=nroldy ) THEN
               Fpy(numy1) = Fpy(numy1) - fac
               Fpy(numy) = Fpy(numy) + fac
            ENDIF
            nroldy = numy
            IF ( numx==nroldx ) PRINT * , "BA FPGRRE 20"
            IF ( numx/=nroldx ) THEN
               Fpx(numx1) = Fpx(numx1) - fac
               Fpx(numx) = Fpx(numx) + fac
            ENDIF
         ENDDO
         nroldx = numx
      ENDDO
      END
