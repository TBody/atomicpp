      subroutine fpgivs(piv,ww,cos,sin)
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      real*8 piv,ww,cos,sin
c  ..local scalars..
      real*8 dd,one,store
c  ..function references..
      real*8 abs,sqrt
c  ..
C       print *, ""
C       print *, "fpgivs called with"
C       print *, "piv = ", piv
C       print *, " ww = ",  ww
C       print *, "cos = ", cos
C       print *, "sin = ", sin      
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd

C       print *, ""
C       print *, "fpgivs returned with"
C       print *, "piv = ", piv
C       print *, " ww = ",  ww
C       print *, "cos = ", cos
C       print *, "sin = ", sin      

      return
      end
