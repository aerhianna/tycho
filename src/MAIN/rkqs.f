      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER*4 n,NMAX
      real*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER*4 i
      real*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY
     +,PGROW,PSHRNK,ERRCON
      PARAMETER(SAFETY=0.9D0,PGROW=-0.2D0,PSHRNK=-0.25D0,ERRCON=1.89d-4)
c-------------------------------------------------------------------------
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.0D0
      do 11 i=1,n
        errmax=max(errmax,dabs(yerr(i)/yscal(i)))
11    continue
      
      errmax=errmax/eps
      if(errmax.gt.1.D0)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=dsign(max(dabs(htemp),0.1D0*dabs(h)),h)
        xnew=x+h
        if(xnew.eq.x)then
           write(*,*)xnew,x,h
           stop 'stepsize underflow in rkqs'
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
c          hnext=5.D0*h
c..added wda 8-29-07 for tighter stepsize increase control
          hnext=2.D0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
