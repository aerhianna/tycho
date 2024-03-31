      SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs,ydum)
      implicit none
      INTEGER nstep,nvar,NMAX,NSTPMX
c..   revise Numerical Recipes default dimensions
      parameter( nmax=4, nstpmx=5000)
      DOUBLE PRECISION x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX),
     1                 ydum(NMAX,NSTPMX)
      EXTERNAL derivs
      real*8 sigma2
      COMMON /path/ xx,y,sigma2
CU    USES rk4
      INTEGER i,k
      DOUBLE PRECISION h,x,dv(NMAX),v(NMAX)
c..   added for use in envel2.f, envel4.f
      real*8 rktmass
      integer*4 rkkp,rknc,rkmodeg,rki,rkj
      common/cderivs/rktmass,rkkp,rknc,rkmodeg,rkj
c---------------------------------------------------------
      do 11 i=1,nvar
        v(i)=vstart(i)
        y(i,1)=v(i)
        ydum(i,1)=v(i)
11    continue
      xx(1)=x1
      x=x1
      h=(x2-x1)/nstep
      do 13 k=1,nstep
c..   added for use in envel2.f
         rkj = rkj + 1
         
        call derivs(x,v,dv)
        call rk4(v,dv,nvar,x,h,v,derivs)
        if(x+h.eq.x)stop 'stepsize not significant in rkdumb'
        x=x+h
        xx(k+1)=x
        do 12 i=1,nvar
          y(i,k+1)=v(i)
          ydum(i,k+1)=v(i)
12      continue
13    continue

      return
      END
