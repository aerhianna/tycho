      subroutine getsurf(g,kks,nc,itfinal,accel)

      implicit none

c..   surface subroutine getsurf:
c..   solves for consistent pressure and opacity at given Teff and
c..   gravity g at the photosphere (tau = 0.7)
c..   this implies that L and R are known, so that iteration is on
c..   the remaining unknown, V(eff) = 1/rho(eff)
c..   Takes luminosity from tl(nc,kks) 
c..   Takes Teff from t(nc,kks)
c..   and puts results in arrays at (time,zone) = (nc,kks)
c..   last revised by D. Arnett, 10/11/02

      include 'dimenfile'
      include 'cconst'
      include 'comod'

      real*8 d(kdm), g(kdm)

      real*8 taus,adelv,delv,apk,apv,safe,dtv,tiny,accel

      integer*4 kks,nc,itfinal,its,itsmax

      data taus/0.7d0/,safe/0.4d0/,tiny/1.0d-8/,itsmax/300/
c      data taus/0.667d0/,safe/0.4d0/,tiny/1.0d-8/,itsmax/300/
c------------------------------------------------------------
c..   evaluate at "photosphere", so taus = tau => 0.7
c..   safe is maximum size of log extrapolation
c..   tiny is convergence criterion
c..   itsmax is maximum index of iteration loop

      if( g(kks) .le. 0.0d0 )then
         write(*,'(1pe12.3,a30)')g(kks),' ERROR in getsurf, g .le. 0 '
         stop' getsurf0'
      endif

      do its = 1,itsmax

         call state(kks,kks,nc)

         d(kks) = 1.0d0 / v(nc,kks)
c..   iterate on specific volume
         apk    = ak(kks)*p(nc,kks) - g(kks)*taus *2.0d0
c..eddington grey atmosphere, boundary at kk+1, gives factor 2
c..scaling in taus
         apv    = ak(kks)*pv(kks) + akv(kks)*p(nc,kks)
         delv   = -apk/apv

c..   limit size of change for iteration stability
         adelv = abs( delv )
         if( adelv .gt. safe*v(nc,kks) )then
            dtv = delv/adelv*safe*v(nc,kks)
         else
            dtv = delv
         endif
         if( its .gt. 100 )then
            dtv = 0.5d0*dtv
         endif
         if( adelv .lt. tiny*v(nc,kks) )then
            goto 100
         endif

         v(nc,kks) = v(nc,kks) + dtv

         if( its .eq. itsmax - 10 )then
            write(*,*)'iteration may not converge in getsurf, trial ',
     1           its

            write(*,'(a5,8a11,a16,5a11)') "its","t(nc,kks)","d(kks)",
     1           "p(nc,kks)", "ak(kks)","apk","apv","delv","dtv",
     2           "v(nc,kks)",'beta','ak*p',
     3           'taus*g','F*ak/c','g'
         elseif( its .gt. itsmax - 10 )then
            write(*,'(i5,1p8e11.3,1pe16.8,1p5e11.3)')
     1           its,t(nc,kks),d(kks),p(nc,kks),ak(kks),apk,apv,
     2           delv,dtv,v(nc,kks)
     3           ,1.0d0-arad*t(nc,kks)**4/(3.0d0*p(nc,kks) ),
     4           ak(kks)*p(nc,kks),taus*g(kks),
     5           ak(kks)*sigma*t(nc,kks)**4/crad, g(kks)
         endif

      enddo

c..   termination on nonconvergence
      write(*,*)'iteration did not converg in getsurf',its
      write(*,'(a5,8a11,a16)') "its","t(nc,kks)","d(kks)","p(nc,kks)",
     1     "ak(kks)","apk","apv","delv","dtv","v(nc,kks)"
      write(*,'(i5,1p8e11.3,1pe16.8)')
     1     its,t(nc,kks),d(kks),p(nc,kks),ak(kks),apk,apv,
     2     delv,dtv,v(nc,kks)
      write(*,'(a20,1p4e16.8)') 'ropal ',d(kks)/(1.0d-6*t(nc,kks))**3
      write(*,'(1p9e16.8)')p(nc,kks)*akv(kks),ak(kks)*pv(kks),
     1  p(nc,kks)*akv(kks)+ak(kks)*pv(kks)
      write(*,'(1p9e16.8)')p(nc,kks),akv(kks),ak(kks),pv(kks)
      write(*,'(1p9e16.8)')p(nc,kks)*ak(kks)/taus,g(kks),
     1   p(nc,kks)*ak(kks)/taus-g(kks)

      stop'getsurf1'

 100  continue
c..   success
      itfinal = its
      accel   = ak(kks)*p(nc,kks)/taus
      
      return
      end



