!*==REGRID.spg  processed by SPAG 6.72Dc at 10:30 on 17 Aug 2017
      SUBROUTINE REGRID(Iopt,Mx,X,My,Y,Z,Xb,Xe,Yb,Ye,Kx,Ky,S,Nxest,     &
                      & Nyest,Nx,Tx,Ny,Ty,C,Fp,Wrk,Lwrk,Iwrk,Kwrk,Ier)
      IMPLICIT NONE
!*--REGRID5
 
!  ..scalar arguments..
      REAL*8 Xb , Xe , Yb , Ye , S , Fp
      INTEGER Iopt , Mx , My , Kx , Ky , Nxest , Nyest , Nx , Ny ,      &
            & Lwrk , Kwrk , Ier
!  ..array arguments..
      REAL*8 X(Mx) , Y(My) , Z(Mx*My) , Tx(Nxest) , Ty(Nyest) ,         &
           & C((Nxest-Kx-1)*(Nyest-Ky-1)) , Wrk(Lwrk)
      INTEGER Iwrk(Kwrk)
!  ..local scalars..
      REAL*8 tol
      INTEGER i , j , jwrk , kndx , kndy , knrx , knry , kwest , kx1 ,  &
            & kx2 , ky1 , ky2 , lfpx , lfpy , lwest , lww , maxit , nc ,&
            & nminx , nminy , mz
!  ..function references..
      INTEGER MAX0
!  ..subroutine references..
!    fpregr,fpchec
!  ..
      PRINT * , "regrid called"
!  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1E-02
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
      Ier = 10
      IF ( Kx>0 .AND. Kx<=5 ) THEN
         kx1 = Kx + 1
         kx2 = kx1 + 1
         IF ( Ky>0 .AND. Ky<=5 ) THEN
            ky1 = Ky + 1
            ky2 = ky1 + 1
            IF ( Iopt>=(-1) .AND. Iopt<=1 ) THEN
               nminx = 2*kx1
               IF ( Mx>=kx1 .AND. Nxest>=nminx ) THEN
                  nminy = 2*ky1
                  IF ( My>=ky1 .AND. Nyest>=nminy ) THEN
                     mz = Mx*My
                     nc = (Nxest-kx1)*(Nyest-ky1)
                     lwest = 4 + Nxest*(My+2*kx2+1) + Nyest*(2*ky2+1)   &
                           & + Mx*kx1 + My*ky1 + MAX0(Nxest,My)
                     kwest = 3 + Mx + My + Nxest + Nyest
                     IF ( Lwrk>=lwest .AND. Kwrk>=kwest ) THEN
                        IF ( Xb<=X(1) .AND. Xe>=X(Mx) ) THEN
                           DO i = 2 , Mx
                              IF ( X(i-1)>=X(i) ) GOTO 99999
                           ENDDO
                           IF ( Yb<=Y(1) .AND. Ye>=Y(My) ) THEN
                              DO i = 2 , My
                                 IF ( Y(i-1)>=Y(i) ) GOTO 99999
                              ENDDO
                              IF ( Iopt>=0 ) THEN
                                 IF ( S<0. ) GOTO 99999
                                 IF ( S==0. .AND.                       &
                                    & (Nxest<(Mx+kx1) .OR. Nyest<       &
                                    & (My+ky1)) ) GOTO 99999
                                 Ier = 0
                              ELSE
                                 IF ( Nx<nminx .OR. Nx>Nxest )          &
                                    & GOTO 99999
                                 j = Nx
                                 DO i = 1 , kx1
                                    Tx(i) = Xb
                                    Tx(j) = Xe
                                    j = j - 1
                                 ENDDO
                                 CALL FPCHEC(X,Mx,Tx,Nx,Kx,Ier)
                                 IF ( Ier/=0 ) GOTO 99999
                                 IF ( Ny<nminy .OR. Ny>Nyest )          &
                                    & GOTO 99999
                                 j = Ny
                                 DO i = 1 , ky1
                                    Ty(i) = Yb
                                    Ty(j) = Ye
                                    j = j - 1
                                 ENDDO
                                 CALL FPCHEC(Y,My,Ty,Ny,Ky,Ier)
                                 IF ( Ier/=0 ) GOTO 99999
                              ENDIF
!  we partition the working space and determine the spline approximation
                              lfpx = 5
                              lfpy = lfpx + Nxest
                              lww = lfpy + Nyest
                              jwrk = Lwrk - 4 - Nxest - Nyest
                              knrx = 4
                              knry = knrx + Mx
                              kndx = knry + My
                              kndy = kndx + Nxest
                              CALL FPREGR(Iopt,X,Mx,Y,My,Z,mz,Xb,Xe,Yb, &
                               & Ye,Kx,Ky,S,Nxest,Nyest,tol,maxit,nc,Nx,&
                               & Tx,Ny,Ty,C,Fp,Wrk(1),Wrk(2),Wrk(3),    &
                               & Wrk(4),Wrk(lfpx),Wrk(lfpy),Iwrk(1),    &
                               & Iwrk(2),Iwrk(3),Iwrk(knrx),Iwrk(knry), &
                               & Iwrk(kndx),Iwrk(kndy),Wrk(lww),jwrk,   &
                               & Ier)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
99999 END