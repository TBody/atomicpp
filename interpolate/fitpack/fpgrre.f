      subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,
     * ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,
     * ax,ay,bx,by,nrx,nry)
c  ..
c  ..scalar arguments..
      real*8 p,fp
      integer ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,
     * kx1,kx2,ky1,ky2
c  ..array arguments..
      real*8 x(mx),y(my),z(mz),tx(nx),ty(ny),c(nc),spx(mx,kx1),spy(my,ky
     *1)
     * ,right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),ay(ny,ky2),by(ny,ky2),
     * fpx(nx),fpy(ny)
      integer nrx(mx),nry(my)
c  ..local scalars..
      real*8 arg,cos,fac,pinv,piv,sin,term,one,half
      integer i,ibandx,ibandy,ic,iq,irot,it,iz,i1,i2,i3,j,k,k1,k2,l,
     * l1,l2,ncof,nk1x,nk1y,nrold,nroldx,nroldy,number,numx,numx1,
     * numy,numy1,n1
c  ..local arrays..
      real*8 h(7)
      print *, "fpgrre called"
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fprota
c  ..
c  the b-spline coefficients of the smoothing spline are calculated as
c  the least-squares solution of the over-determined linear system of
c  equations  (ay) c (ax)' = q       where
c
c               |   (spx)    |            |   (spy)    |
c        (ax) = | ---------- |     (ay) = | ---------- |
c               | (1/p) (bx) |            | (1/p) (by) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
c                b-spline coefficients.
c       z      : the my x mx matrix which contains the function values.
c       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
c                matrices according to the least-squares problems in
c                the x- and y-direction.
c       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
c                matrices which contain the discontinuity jumps of the
c                derivatives of the b-splines in the x- and y-direction.

C       print *, "fpgrre ifsx", ifsx
C       print *, "fpgrre ifsy", ifsy
C       print *, "fpgrre ifbx", ifbx
C       print *, "fpgrre ifby", ifby
C       print *, "fpgrre x", x
C       print *, "fpgrre mx", mx
C       print *, "fpgrre y", y
C       print *, "fpgrre my", my
C       print *, "fpgrre z", z
C       print *, "fpgrre mz", mz
C       print *, "fpgrre kx", kx
C       print *, "fpgrre ky", ky
C       print *, "fpgrre tx", tx
C       print *, "fpgrre nx", nx
C       print *, "fpgrre ty", ty
C       print *, "fpgrre ny", ny
C       print *, "fpgrre p", p
C       print *, "fpgrre c", c
C       print *, "fpgrre nc", nc
C       print *, "fpgrre fp", fp
C       print *, "fpgrre fpx", fpx
C       print *, "fpgrre fpy", fpy
C       print *, "fpgrre mm", mm
C       print *, "fpgrre mynx", mynx
C       print *, "fpgrre kx1", kx1
C       print *, "fpgrre kx2", kx2
C       print *, "fpgrre ky1", ky1
C       print *, "fpgrre ky2", ky2
C       print *, "fpgrre spx", spx
C       print *, "fpgrre spy", spy
C       print *, "fpgrre right", right
C       print *, "fpgrre q", q
C       print *, "fpgrre ax", ax
C       print *, "fpgrre ay", ay
C       print *, "fpgrre bx", bx
C       print *, "fpgrre by", by
C       print *, "fpgrre nrx", nrx
C       print *, "fpgrre nry", nry




      one = 1
      half = 0.5
      nk1x = nx-kx1
      nk1y = ny-ky1

C       print *, "p = ", p
C       print *, "ifsx = ", ifsx
C       print *, "ifsy = ", ifsy
C       print *, "ifbx = ", ifbx
C       print *, "ifby = ", ifby

C       if(p.gt.0.) print *, "BA FPGRRE 1"
      if(p.gt.0.) pinv = one/p
c  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
c  the value of p whether the matrices (spx),(spy),(bx) and (by) still
c  must be determined.
C       if(ifsx.ne.0) print *, "BA FPGRRE 2"
      if(ifsx.ne.0) go to 50
c  calculate the non-zero elements of the matrix (spx) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the x-direction.
      l = kx1
      l1 = kx2
      number = 0

C       print *, "tx", tx

      do 40 it=1,mx
        arg = x(it)
  10    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
C   10    if(arg.lt.tx(l1) .or. l.eq.nk1x) print *, "BA FPGRRE 3"
C         if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        do 30 i=1,kx1
          spx(it,i) = h(i)
  30    continue
        nrx(it) = number
  40  continue

C       print *, "spx", spx
C       print *, "spx(2,4)", spx(2,4)
C       print *, "nrx", nrx
C       print *, "h", h
C       stop

      ifsx = 1
  50    if(ifsy.ne.0) go to 100   
C   50  if(ifsy.ne.0) print *, "BA FPGRRE 4"
C       if(ifsy.ne.0) go to 100
c  calculate the non-zero elements of the matrix (spy) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the y-direction.
      l = ky1
      l1 = ky2
      number = 0
      do 90 it=1,my
        arg = y(it)
  60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
C   60    if(arg.lt.ty(l1) .or. l.eq.nk1y)  print *, "BA FPGRRE 5"
C         if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call fpbspl(ty,ny,ky,arg,l,h)
        do 80 i=1,ky1
          spy(it,i) = h(i)
  80    continue
        nry(it) = number
  90  continue
      ifsy = 1

C       print *, "spy", spy
C       print *, "spy(2,4)", spy(2,4)
C       print *, "nry", nry
C       print *, "h", h
C       stop


 100  if(p.le.0.) go to 120
C  100  if(p.le.0.) print *, "BA FPGRRE 6"
C       if(p.le.0.) go to 120
c  calculate the non-zero elements of the matrix (bx).
C       if(ifbx.ne.0 .or. nx.eq.2*kx1) print *, "BA FPGRRE 7"
      if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110
      call fpdisc(tx,nx,kx2,bx,nx)
      ifbx = 1
c  calculate the non-zero elements of the matrix (by).
 110  if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
C  110  if(ifby.ne.0 .or. ny.eq.2*ky1) print *, "BA FPGRRE 8"
C       if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
      call fpdisc(ty,ny,ky2,by,ny)
      ifby = 1
c  reduce the matrix (ax) to upper triangular form (rx) using givens
c  rotations. apply the same transformations to the rows of matrix q
c  to obtain the my x (nx-kx-1) matrix g.
c  store matrix (rx) into (ax) and g into q.
 120  l = my*nk1x
c  initialization.
      do 130 i=1,l
        q(i) = 0.
 130  continue
      do 140 i=1,nk1x
        do 140 j=1,kx2
          ax(i,j) = 0.
 140  continue
      l = 0
      nrold = 0
c  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1

      do 270 it=1,mx
        number = nrx(it)
  150   if(nrold.eq.number) go to 180
C  150    if(nrold.eq.number) print *, "BA FPGRRE 9"
C         if(nrold.eq.number) go to 180
C         if(p.le.0.) print *, "BA FPGRRE 10"
C         if(p.le.0.) go to 260
        nrold = number
c  fetch a new row of matrix (spx).
 180    h(ibandx) = 0.
        do 190 j=1,kx1
          h(j) = spx(it,j)
 190    continue
c  find the appropriate column of q.
        do 200 j=1,my
          l = l+1
          right(j) = z(l)
 200    continue
        irot = number
c  rotate the new row of matrix (ax) into triangle.
C         print *, "h = ", h
 210    do 240 i=1,ibandx
          
          irot = irot+1
          piv = h(i)
C           print *, ""
C           print *, "piv", piv
C           print *, "ax(irot,1)", ax(irot,1)
C           print *, "cos", cos
C           print *, "sin", sin

C           if(piv.eq.0.) print *, "BA FPGRRE 11"
          if(piv.eq.0.) go to 240
c  calculate the parameters of the given transformation.
          call fpgivs(piv,ax(irot,1),cos,sin)
c  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do 220 j=1,my
            iq = iq+1
            call fprota(cos,sin,right(j),q(iq))
 220      continue
c  apply that transformation to the columns of (ax).
C           if(i.eq.ibandx) print *, "BA FPGRRE 12"
          if(i.eq.ibandx) go to 250
          i2 = 1
          i3 = i+1
          do 230 j=i3,ibandx
            i2 = i2+1
            call fprota(cos,sin,h(j),ax(irot,i2))
 230      continue
 240    continue
 250    if(nrold.eq.number) go to 270
C  250    if(nrold.eq.number) print *, "BA FPGRRE 13"
C         if(nrold.eq.number) go to 270
C  260    nrold = nrold+1
        go to 150
 270  continue

C   19  format (' ',A6,7(f6.2))
C   18  format (' ',A6,9(f6.2))
C   17  format (' ',A6,45(f6.2))
C       print 19, "h", h
C       print 18, "right", right
C       do 275 i=1,9
C         do 276 j=1,5
C           print *, i, j, ax(i,j)
C   276   continue
C   275 continue
C       stop

c  reduce the matrix (ay) to upper triangular form (ry) using givens
c  rotations. apply the same transformations to the columns of matrix g
c  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
c  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
c  initialization.
      do 280 i=1,ncof
        c(i) = 0.
 280  continue

      do 290 i=1,nk1y
        do 290 j=1,ky2
          ay(i,j) = 0.
 290  continue

      

      nrold = 0
c  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1

      do 420 it=1,my
        number = nry(it)
C         print *, "Number = ", number
C         print *, "nrold = ", nrold
C         nrold = number
C  300    if(nrold.eq.number) go to 330
C  300    if(nrold.eq.number) print *, "BA FPGRRE 14"

c  fetch a new row of matrix (spy)
 330    h(ibandy) = 0.
        do 340 j=1,ky1
          h(j) = spy(it,j)
 340    continue
c  find the appropiate row of g.
        l = it
        do 350 j=1,nk1x
          right(j) = q(l)
          l = l+my
 350    continue
        irot = number
c  rotate the new row of matrix (ay) into triangle.
 360    do 390 i=1,ibandy
          irot = irot+1
          piv = h(i)
C           if(piv.eq.0.) print *, "BA FPGRRE 16"
          if(piv.eq.0.) go to 390
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ay(irot,1),cos,sin)
c  apply that transformation to the colums of matrix g.
          ic = irot
          do 370 j=1,nk1x
            call fprota(cos,sin,right(j),c(ic))
            ic = ic+nk1y
 370      continue
c  apply that transformation to the columns of matrix (ay).
C           if(i.eq.ibandy) print *, "BA FPGRRE 17"
          if(i.eq.ibandy) go to 400
          i2 = 1
          i3 = i+1
          do 380 j=i3,ibandy
            i2 = i2+1
            call fprota(cos,sin,h(j),ay(irot,i2))
 380      continue
 390    continue
 400    if(nrold.eq.number) go to 420
C  400    if(nrold.eq.number) print *, "BA FPGRRE 18"
C         if(nrold.eq.number) go to 420
C  410    nrold = nrold+1
C         go to 300
 420  continue

C   29  format (' ',A6,7(f6.2))
C   28  format (' ',A6,9(f6.2))
C       print 29, "h", h
C       print 28, "right", right
C       do 275 i=1,9
C         do 276 j=1,5
C           print *, i, j, ay(i,j)
C   276   continue
C   275 continue
C       stop






c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (ry) c (rx)' = h.
c  first step: solve the system  (ry) (c1) = h.
      k = 1
      do 450 i=1,nk1x
        call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
        k = k+nk1y
 450  continue

C UP TO HERE!!!

c  second step: solve the system  c (rx)' = (c1).
      k = 0
      do 480 j=1,nk1y
        k = k+1
        l = k
        print *, "l = ", l
        do 460 i=1,nk1x
          right(i) = c(l)
C           print *, 'right(i)', right(i)
          l = l+nk1y
 460    continue
        call fpback(ax,right,nk1x,ibandx,right,nx)
        l = k
        do 470 i=1,nk1x
          c(l) = right(i)
          l = l+nk1y
 470    continue
 480  continue

      do 277 i=1,25
          print *, i, c(i)
  277 continue

      stop


c  calculate the quantities
c    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
c    fp = sumi=1,mx(sumj=1,my(res(i,j)))
c    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
c                  tx(r+kx) <= x(i) <= tx(r+kx+1)
c    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
c                  ty(r+ky) <= y(j) <= ty(r+ky+1)
      fp = 0.
      do 490 i=1,nx
        fpx(i) = 0.
 490  continue
      do 500 i=1,ny
        fpy(i) = 0.
 500  continue
      nk1y = ny-ky1
      iz = 0
      nroldx = 0
c  main loop for the different grid points.
      do 550 i1=1,mx
        numx = nrx(i1)
        numx1 = numx+1
        nroldy = 0
        do 540 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
          iz = iz+1
c  evaluate s(x,y) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (x,y), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numx*nk1y+numy
          do 520 l1=1,kx1
            k2 = k1
            fac = spx(i1,l1)
            do 510 l2=1,ky1
              k2 = k2+1
              term = term+fac*spy(i2,l2)*c(k2)
 510        continue
            k1 = k1+nk1y
 520      continue
c  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
c  adjust the different parameters.
          fp = fp+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*half
C           if(numy.eq.nroldy) print *, "BA FPGRRE 19"
          if(numy.eq.nroldy) go to 530
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
 530      nroldy = numy
C           if(numx.eq.nroldx) print *, "BA FPGRRE 20"
          if(numx.eq.nroldx) go to 540
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
 540    continue
        nroldx = numx
 550  continue
      return
      end

