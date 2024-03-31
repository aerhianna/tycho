      subroutine pgplot(no,ktype,kkman)

c..   7-30-06
c..   main graphics subroutine, with multiple panels

c..   NOTE THAT REALS ARE IMPLICITLY SINGLE PRECISION
      implicit real*4 (a-h,o-z)
      implicit integer*4(i-n)

      include 'dimenfile'

      real*8 dj

      include 'compu'
      include 'comod'
      include 'cburn'
      include 'cpgplot'
      include 'cnabla'
      include 'cconst'
      include 'conline'
      include 'cbug'
      include 'cadvec'
      include 'caeps'
      include 'crate'
      include 'cenv'
      include 'csurface'
      include 'ctmstp'
      include 'cenrchk'
      include 'czone'

      integer*4 nzpg
      parameter (nzpg = 29)

      integer*4 ixlogflg

      real*4 zmet,xmet,ymet
      real*4 xx(kdme),yy(kdme),zz(kdme,nzpg),xh(kdme),xb(kdme)
      real*4 setem(kdm),serho(kdm),ser(kdm),sem(kdm),sek(kdm)
      real*4 dummy, dxx

c..gravitational energy release and luminosity
      real*8 lgrav(kdm), egrav(kdm)


      character*20 xlabel, ylabel, zlabel(nzpg)
      character*10 ctime, cdt
      character*5 cnum
      character*15 cdummy
      character cpgdtlab(4),cdtdum
      character*1 cmode,cmodec,cmodes
      character*15 cmodeall

      real*4 tabul(12),tabhi(8)
      real*4 tabc1(6),tabc2(5)
      real*4 fontsize
      real*4 dmhmax,dmhmin
      data dmhmax/0.0/, dmhmin/-8.0/

c..   device dependent parameters, set for 1280x1024 screen (5-20-01) 

c      data tabul/ -0.05, 0.13, 0.21, 0.34, 0.42, 0.55,
c     1     0.63, 0.76, 0.84, 0.94/
      data tabul/ -0.04, 0.09, 0.17, 0.30, 0.38, 0.51,
     1     0.59, 0.72, 0.80, 0.90, 0.98, 1.04/

      data tabhi/ 0.00,-1.40,-2.80,-4.20,
     1     -8.40,-5.60,-7.00,-9.8/

      data tabc2/ 0.05, 0.25, 0.50, 0.75, 1.00/
      data tabc1/ 2.0, 2.95, 3.90, 4.85, 5.80, 6.75/

      data fjust/1.0/

      data zlabel/" log rho  ", " log T   "," log P  ",
     1     " Luminosity", " v(con) ",
     1     " Ye ", " H ", " n ", " He4 ", " C12 ",
     2     " O16 "," Ne20 "," Mg24 "," SiCa "," Ni56 ",
     3     " Co56 ", "Fe56 ",
     4     " dT ", " dV "," dR "," dL ",
     5     " Radius   "," Mass ",
     6     " hp ", " u ",
     7     " nab ", "nad", "nrad", "dm"/

      data ctime/" time(sec)"/, cdt/"dtime(sec)"/

      data cpgdtlab/"R","T","V","L"/

      integer*4 pgopen

c..pgscal scales the size of the pgplot, hrplt, cvplt screens
c     standard (1028-768 screen)
      data fontsize/1.2/,pgscal/1.0/
c     higher resolution (more points, so closer points)
c      data fontsize/1.2/,pgscal/1.20/

c..colors:
c..0=background (black)
c..1=white
c..2=red 
c..3=green
c..4=blue
c..5=turquoise
c..6=magenta(red+blue)
c..7=yellow
c..8=gold
c..9=green+yellow
c..10=blue+cyan
c..11=medium blue
c..12=violet
c..13=deep pink
c..14=dark grey
c..15=grey
c..15<i<256 (16to255) undefined

c--------------------------------------------------------
      if( igraf .ne. 0 )return

      if( l .eq. 1 )then

c..   scratch file for character conversion
         open(10,file='scratch.pg')

         idpg1 = pgopen(device1)

         IF ( idpg1 .LE. 0 ) STOP' pgopen error'

c..   10 works fine on chandra..................................
         call pgpap(8.0*pgscal,0.8)

         call pgsubp(1,5)

c..   no query for device
         call PGASK (.FALSE.)
c..   roman font = 2
         call PGSCF(1)
c..   font scaling (size)
         call PGSCH(3.5*fontsize)
c..   new page
         call PGPAGE
      endif


c..   select device1
      call pgslct(idpg1)

c..   setup buffering
      call PGBBUF

      sptime = time
      sdt    = dth(2)

      nx  = kk
      nxb = kk + 1
c      nx  = kk-1
c      nxb = kk

      if( modes .eq. 2 )then
         nxe = nx + 1
      else
         nxe = nx
      endif
      do i = 1, nxe
         if( v(1,i+1) .gt. 0.0d0 )then
            zz(i,1) = 1.0d0/v(1,i+1)
         else
            zz(i,1) = 0.0
         endif

         if( ixflag .eq. 0 )then
            xh(i)   = 0.5 * ( r(1,i+1) + r(1,i) )
         elseif( ixflag .eq. 1 )then
            xh(i)   = dlog10( 0.5 * ( r(1,i+1) + r(1,i) ) )
         elseif( ixflag .eq. 2 )then
            xh(i)   = (xm(i) + 0.5*dmh(i+1) )/sol
         elseif( ixflag .eq. 3 )then
            xh(i)   = float(i) - 0.5
         else
            write(*,*)'pgplot: ixflag error'
            stop
         endif

         zz(i,2) =  t(1,i+1)
         zz(i,3) =  p(1,i+1)
         zz(i,18)= dt(2,i+1)
         zz(i,19)= dv(2,i+1)
         zz(i,29)= dmh(i+1)/sol
      enddo
      if( ixflag .eq. 0 )then
         xh(nxe)   = r(1,nxe)
      elseif( ixflag .eq. 1 )then
         xh(nxe)   = dlog10( r(1,nxe) )
      elseif( ixflag .eq. 2 )then
         xh(nxe)   = xm(nxe) /sol
      elseif( ixflag .eq. 3 )then
         xh(nxe)   = nxe
      else
         write(*,*)'pgplot: ixflag error'
         stop
      endif
      
c..   zone boundary quantities
      do i = 1, nxb
         zz(i,4) = tl(1,i)/sollum
         zz(i,5) =   hp(i)
         zz(i,20)=   dr(i)
         zz(i,21)=  dtl(i)/sollum
         zz(i,22)=  r(1,i)
         zz(i,23)=   xm(i)
         zz(i,24)=   hp(i)
         zz(i,25)=  u(1,i)

         if( ixflag .eq. 0 )then
            xb(i)   = r(1,i)
         elseif( ixflag .eq. 1 )then
            if( i .eq. 1 )then
               xb(i) = dlog10( 0.5d0*r(1,2) )
            else
               xb(i) = dlog10( r(1,i) )
            endif
         elseif( ixflag .eq. 2 )then
            xb(i)   = xm(i)/sol
         elseif( ixflag .eq. 3 )then
            xb(i)   = float(i)
         endif
      enddo


c..   reset y array                            (25=fluid velocity)
      umin = u(1,1)
      umax = u(1,1)
      kun = 1
      kux = 1
      do i=1, nx
         if( u(1,i) .gt. umax) then
            umax = u(1,i)
            kux  = i
         endif
         if( u(1,i) .lt. umin) then
            umin = u(1,i)
            kun  = i
         endif
      enddo
c..   reset label                            (nuclear burning)
c..   find max and min, array indices
      epmin = 0
      epmax = 0
      kepn = 1
      kepx = 1
      do i=1, nx
         if( s(5,i) .gt. dble(epmax) ) then
            epmax = s(5,i)
            kepx  = i
         endif
         if( s(4,i) .lt. dble(epmin) ) then
            epmin = s(4,i)
            kepn  = i
         endif
c..   negative energy generation allowed
         if( s(5,i) .lt. dble(epmin) ) then
            epmin = s(5,i)
            kepn  = i
         endif
      enddo

c..   gravitational energy generation rate
      egrav(1) = 0.0d0
      lgrav(1) = 0.0d0
      do i = 2, nx
         egrav(i) = -(       e(2,i) - e(1,i)
     1        + 0.5d0*( p(2,i) + p(1,i) )
     2        *( v(2,i) - v(1,i) ) )/dth(2)
         lgrav(i) = lgrav(i-1) + dmh(i) * egrav(i)
      enddo

cccccc      call minmax(egrav,1,nx,egmin,egmax)

c..   first panel.....................DIAGNOSTIC SCALARS
      call PGPANL(1, 1)
c..   erase screen
      call PGERAS
c..   viewplane
      call PGSVP(0.1, 0.9, 0.1, 0.9)
c..   window
      call PGSWIN(0.0, 1.0, 0.0, 1.0)

c..   font scaling (size) for first line

      call PGSCH(3.0*fontsize)
      call  PGSCI(1)
c..   TYCHO version
      call PGMTXT('T', tabhi(1), -0.1, 0.0, ctstamp)
c..   ll value for this run
      call modflg(ll,cnum)
      call PGMTXT('T', tabhi(1), tabul(2), fjust, cnum)
c..   hydrodynamic mode
      if( mode .eq. 0 )then
         call PGMTXT('T', tabhi(1), tabul(3), fjust, 'dynam')
      elseif( mode .eq. 1 )then
         call PGMTXT('T', tabhi(1), tabul(3), fjust, 'relax')
      else
         call PGMTXT('T', tabhi(1), tabul(3), fjust, 'damped')
      endif
c..   burn option and netsize
      if( nburn .eq. -1 )then
         call PGMTXT('T', tabhi(1), tabul(4), fjust, 'noburn')
      else
         call modflg(netsize,cnum)
         cdummy = cnum//' nuc'
         call PGMTXT('T', tabhi(1), tabul(4), fjust, cdummy)
      endif
c..   envelope option
      if( modes .eq. 2 )then
         call PGMTXT('T', tabhi(1), tabul(5), fjust, 'envel')
      else
         call PGMTXT('T', tabhi(1), tabul(5), fjust, '0-surf')
      endif
c..   opacity option
      if( nopac .eq. -1 )then
         call PGMTXT('T', tabhi(1), tabul(6), fjust, 'nopac')
      elseif( nopac .eq. 0 )then
         call PGMTXT('T', tabhi(1), tabul(6), fjust, 'OPAL2')
      elseif( nopac .eq. 1 )then
         call PGMTXT('T', tabhi(1), tabul(6), fjust, 'OPAL1')
c         call PGMTXT('T', tabhi(1), tabul(7), fjust, copal)
      endif
c..   EOS option
      if( nopaleos .eq. 0   )then
         call PGMTXT('T', tabhi(1), tabul(7), fjust, 'HELM')
      else
         call PGMTXT('T', tabhi(1), tabul(7), fjust, 'OPALEOS')
      endif
c..   gravity option
      if( newt .eq. 0   )then
         call PGMTXT('T', tabhi(1), tabul(8), fjust, 'GR')
      else
         call PGMTXT('T', tabhi(1), tabul(8), fjust, 'Newton')
      endif
c.. 
c..   rezoning option
      if( ktot .eq. 0   )then
         call PGMTXT('T', tabhi(1), tabul(9), fjust, 'nozone')
      else
         call PGMTXT('T', tabhi(1), tabul(9), fjust, 'Rezone')
      endif
c.. 
c      call itoa(mode ,cmode)
c      call itoa(modec,cmodec)
c      call itoa(modes,cmodes)
c      cmodeall = 'mode( cs)'//' '//cmode//' '//cmodec//' '//cmodes
c..   convection option
      if( modec .eq. -1 )then
         call PGMTXT('T', tabhi(1), tabul(10), fjust, 'no cv')
      elseif( modec .eq. 0 )then
         call PGMTXT('T', tabhi(1), tabul(10), fjust, 'dS')
      elseif( modec .eq. 1 )then
         call PGMTXT('T', tabhi(1), tabul(10), fjust, 'Schw')
      elseif( modec .eq. 2 )then
         call PGMTXT('T', tabhi(1), tabul(10), fjust, 'L+mix')
      endif

c..end of first line in first panel

c..   font scaling (size) for body
      call PGSCH(3*fontsize)
c..   annotate
      if( it .lt. iter .and. no .eq. 0 )then
         call  PGSCI(3)
      else
         call  PGSCI(2)
      endif

c..........................
      call PGMTXT('T', tabhi(2), tabul(1), fjust, 'time(s)')
      call ftoc(sptime,ctime)
      call PGMTXT('T', tabhi(2), tabul(2), fjust, ctime)

      call PGMTXT('T', tabhi(2), tabul(3), fjust, 'dt(s)')
      call ftoc(sdt,cdt)
      call PGMTXT('T', tabhi(2), tabul(4), fjust, cdt)

      sing = telog
      call PGMTXT('T', tabhi(2), tabul(5), fjust, 'log Te')
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(2), tabul(6), fjust, cdt)

      sing = xlol
      call PGMTXT('T', tabhi(2), tabul(7), fjust, 'log L')
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(2), tabul(8), fjust, cdt)

      sing = r(1,kk+1)*omeg(kk+1)
      call PGMTXT('T', tabhi(2), tabul(9), fjust, 'vrot')
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(2), tabul(10), fjust, cdt)

      sing = alphaml
      call PGMTXT('T', tabhi(2), tabul(11), fjust, 'alphaml')
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(2), tabul(12), fjust, cdt)

c..........................
      call PGMTXT('T', tabhi(3), tabul(1), fjust, 'Model')
      call modflg(model,cnum)
      call PGMTXT('T', tabhi(3), tabul(2), fjust, cnum)

      call PGMTXT('T', tabhi(3), tabul(3), fjust, 'steps')
      call modflg(l,cnum)
      call PGMTXT('T', tabhi(3), tabul(4), fjust, cnum)

      call PGMTXT('T', tabhi(3), tabul(5), fjust, 'kk')
      call modflg(kk,cnum)
      call PGMTXT('T', tabhi(3), tabul(6), fjust, cnum)

      call PGMTXT('T', tabhi(3), tabul(7), fjust, 'it')
      call modflg(it,cnum)
      call PGMTXT('T', tabhi(3), tabul(8), fjust, cnum)

      call PGMTXT('T', tabhi(3), tabul(9), fjust, 'J53')
      dj = 0.0d0
      do j=2, kk+1
         dj = dj + dmi(j)*r(1,j)**2*omeg(j)
      enddo
      sing = dj*1.0d-53
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(3), tabul(10), fjust, cdt)

      i = 0
      do j = 2, kk
         i = i + iadd(j)
      enddo
      call PGMTXT('T', tabhi(3), tabul(11), fjust, 'nadd')
      call modflg(i,cnum)
      call PGMTXT('T', tabhi(3), tabul(12), fjust, cnum)


c..........................

      call PGMTXT('T', tabhi(4), tabul(9), fjust, 'eta(2)')
      sing = 1.0-2.0*x(ndim,2)
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(4), tabul(10), fjust, cdt)

      call PGMTXT('T', tabhi(4), tabul(7), fjust, 'echk')
      sing = echeck
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(4), tabul(8), fjust, cdt)


      call PGMTXT('T', tabhi(4), tabul(5), fjust, 'wtest')
      if( nbug .eq. 0 )then
         sing = ditmax
      else
         ditmax = wtest(1)
         nitest = 1
         do ii = 2,4
            if( abs( wtest(ii) ) .gt. abs(ditmax) )then
               ditmax = wtest(ii)
               kitmax = nwtest(ii)
               sing   = wtest(ii)
               nitest = ii
            endif
         enddo
      endif
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(4), tabul(6), fjust, cdt)

      call modflg(kitmax,cnum)
      call PGMTXT('T', tabhi(4), tabul(4), fjust, cnum)

      if( nitest .gt. 0 .and. nitest .le. 4 )then
         cdtdum = cpgdtlab( nitest )
         call PGMTXT('T', tabhi(4), tabul(3), fjust, cdtdum)
      endif

      if( ktype+1 .ge. 1 .and. ktype+1 .le. condim )then

         call PGMTXT('T', tabhi(4), tabul(1), fjust,
     1        dthflag(ktype+1))
         call modflg(kkman,cnum)
         call PGMTXT('T', tabhi(4), tabul(2), fjust, cnum)

      endif

      i = 0
      do j = 2, kk
         i = i + idel(j)
      enddo
      call PGMTXT('T', tabhi(4), tabul(11), fjust, 'ndel')
      call modflg(i,cnum)
      call PGMTXT('T', tabhi(4), tabul(12), fjust, cnum)


c.......   
      call PGMTXT('T', tabhi(5), tabul(1), fjust, 'knu')
      call modflg(kepn,cnum)
      call PGMTXT('T', tabhi(5), tabul(2), fjust, cnum)
      call PGMTXT('T', tabhi(5), tabul(3), fjust, 'epnu')
      sing = abs( epmin )
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(5), tabul(4), fjust, cdt)

      call PGMTXT('T', tabhi(5), tabul(5), fjust, 'knuc')
      call modflg(kepx,cnum)
      call PGMTXT('T', tabhi(5), tabul(6), fjust, cnum)

      call PGMTXT('T', tabhi(5), tabul(7), fjust, 'epnuc')
      sing = epmax
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(5), tabul(8), fjust, cdt)


c..   xlnukk = (nuclear + neutrino) luminosity 
      xlnukk = 0.0
      do i = 2, nxb
         fact  = ( s(5,i) + s(4,i) )*dmh(i)/sollum
         xlnukk = xlnukk + fact
      enddo
c      call PGMTXT('T', tabhi(5), tabul(9), fjust, 'dm(2)')
c      sing = dmh(2)/sol
      call PGMTXT('T', tabhi(5), tabul(9), fjust, 'Lnuc')
      sing = xlnukk
      call ftoc(sing,cdt)
      call PGMTXT('T', tabhi(5), tabul(10), fjust, cdt)


    
      call PGMTXT('T', tabhi(5), tabul(11), fjust, 'ncytot')
      call modflg(ncytot,cnum)
      call PGMTXT('T', tabhi(5), tabul(12), fjust, cnum)


c..........................

      call PGMTXT('T', tabhi(6), tabul(1), fjust, 'M/sol')
      if( modes .eq. 2 )then
         sing = ( xm(kk) + dmh(kk+1) )/sol
      else
         sing = xm(kk)/sol
      endif

      call ftoc(sing,cdt)

      call PGMTXT('T',  tabhi(6), tabul(2), fjust, cdt)

      call PGMTXT('T',  tabhi(6), tabul(3), fjust, 'Men/sol')
      sing = dmh(kk+1)/sol
      call ftoc(sing,cdt)
      call PGMTXT('T',  tabhi(6), tabul(4), fjust, cdt)

      if( mloss .ne. 0 )then
         call PGMTXT('T', tabhi(6), tabul(5), fjust, 'dM/dt')
         sing = peryear
         call ftoc(sing,cdt)
         call PGMTXT('T',  tabhi(6), tabul(6), fjust, cdt)

c        call PGMTXT('T', tabhi(6), tabul(7), fjust, 'dM(sol)')
c         sing = peryear/secpy*dth(2)
c         call ftoc(sing,cdt)
c         call PGMTXT('T',  tabhi(6), tabul(8), fjust, cdt)
      endif

      if( omeg(kk+1) .ne. 0.0d0 )then
        call PGMTXT('T', tabhi(6), tabul(7), fjust, 'rot/crit')
         sing = omeg(kk+1)*r(1,kk+1)*sqrt( r(1,kk+1)/grav/xm(kk+1) )
         call ftoc(sing,cdt)
         call PGMTXT('T',  tabhi(6), tabul(8), fjust, cdt)
      endif

c..metals all but H and He by this definition
      if( modes .eq. 2 )then
         xmet = 0.0d0
         ymet = 0.0d0
         zmet = 0.0d0
         do n = 1, netsize
            if( lz(n) .eq. 1 )then
               xmet = xmet + x(n,kk+1)
            elseif( lz(n) .eq. 2 )then
               ymet = ymet + x(n,kk+1)
            else
               zmet = zmet + x(n,kk+1)
            endif
         enddo
      else
         xmet = 0.0d0
         ymet = 0.0d0
         zmet = 0.0d0
         do n = 1, netsize
            if( lz(n) .eq. 1 )then
               xmet = xmet + x(n,kk)
            elseif( lz(n) .eq. 2 )then
               ymet = ymet + x(n,kk)
            else
               zmet = zmet + x(n,kk)
            endif
         enddo

      endif

      call PGMTXT('T', tabhi(6), tabul(9), fjust, 'He')
      call ftoc(ymet,cdt)
      call PGMTXT('T',  tabhi(6), tabul(10), fjust, cdt)


    
      call PGMTXT('T', tabhi(6), tabul(11), fjust, 'netadd')
      call modflg(netadd,cnum)
      call PGMTXT('T', tabhi(6), tabul(12), fjust, cnum)


c..........................

      call PGMTXT('T', tabhi(7), tabul(9), fjust, 'metals')
      call ftoc(zmet,cdt)
      call PGMTXT('T',  tabhi(7), tabul(10), fjust, cdt)

      call PGMTXT('T',  tabhi(7), tabul(1), fjust, 'R(cm)')
      if( modes .eq. 2 )then
         sing = r(1,kk+1)
      else
         sing = r(1,kk)
      endif
      call ftoc(sing,cdt)
      call PGMTXT('T',  tabhi(7), tabul(2), fjust, cdt)

      call PGMTXT('T',  tabhi(7), tabul(3), fjust, 'T(c)')
      sing = t(1,2)
      call ftoc(sing,cdt)
      call PGMTXT('T',  tabhi(7), tabul(4), fjust, cdt)

      fact = 1.0d0/v(1,2)
      call PGMTXT('T',  tabhi(7), tabul(5), fjust, 'rho(c)')
      sing = fact
      call ftoc(sing,cdt)
      call PGMTXT('T',  tabhi(7), tabul(6), fjust, cdt)


      call PGMTXT('T',  tabhi(7), tabul(7), fjust, 'T(kk)')
      sing = t(1,kk)
      call ftoc(sing,cdt)
      call PGMTXT('T',  tabhi(7), tabul(8), fjust, cdt)

    
      call PGMTXT('T', tabhi(7), tabul(11), fjust, 'mixmode')
      call modflg(mixmode,cnum)
      call PGMTXT('T', tabhi(7), tabul(12), fjust, cnum)


c..........................


c..   reset color
      call PGSCI(1)
c..   font scaling (size)
      call PGSCH(3.0*fontsize)


c..   second panel.....................LOG T, RHO
c..   reset x label                            (1 = density)
      xlabel = zlabel(1)

      nxe = nx+1

c..   set x array
      do i=1, nx
         if( zz(i,1) .gt. 0.0 )then
c..   log10 density
            xx(i) = alog10( zz(i,1) )
         else
            xx(i) = 0.0
         endif
      enddo

      if( v(1,kk+1) .gt. 0.0d0 )then
         xx(nx+1) = dlog10( 1.0d0/v(1,kk+1) )
      else
         xx(nx+1) = 0.0d0
      endif
      if( v(1,kk+2) .gt. 0.0d0 )then
         xx(nx+2) = dlog10( 1.0d0/v(1,kk+2) )
      else
         xx(nx+2) = 0.0d0
      endif
      if( modes .eq. 2 )then
         call minmax(xx,1,nxe+1,xmin,xmax)
         xmax = amax1( xmax, 6.0 )
      elseif( mode .eq. 0 .and. modes .le. 1 )then
         call minmax(xx,1,nx-1,xmin,xmax)
      else
         call minmax(xx,1,nxe,xmin,xmax)
      endif
      xmin0 = xmin
      xmax0 = xmax
      xmin  = xmin0 - 0.05*(xmax0 - xmin0)
      xmax  = xmax0 + 0.05*(xmax0 - xmin0)

      call PGPANL(1, 2)
c..   erase screen
      call PGERAS
c..   viewplane
      call PGSVP(0.1, 0.9, 0.18, 1.0)

c..   reset color
      call PGSCI(1)
c..   font scaling (size)
      call PGSCH(3.0*fontsize)

c..   reset y label                        ( 2 = temperature )
      ylabel = zlabel(2)
c..   reset y array
      do i=1, nxe
         if( t(1,i+1) .gt. 0.0d0 )then
            yy(i) = dlog10( t(1,i+1) )
         else
            yy(i) = 0
         endif
      enddo
      if( t(1,kk+1) .gt. 0.0d0 )then
         yy(nx+1) = dlog10( t(1,kk+1) )
      else
         yy(nx+1) = 0.0d0
      endif
      if( t(1,kk+2) .gt. 0.0d0 )then
         yy(nx+2) = dlog10( t(1,kk+2) )
      else
         yy(nx+2) = 0.0d0
      endif
      if( modes .eq. 2 )then
         call minmax(yy,1,nxe+1,ymin,ymax)
         ymax = amax1( ymax, 9.0 )
      elseif( mode .eq. 0 .and. modes .le. 1 )then
         call minmax(yy,1,nx-1  ,ymin,ymax)
      else
c     call minmax(yy,1,nxe  ,ymin,ymax)
         call minmax(yy,1,nx-1,ymin,ymax)

      endif
c..keep T minimum low enough to see envelope up to photosphere
      ymin = amin1( ymin, 3.4 )

c..use photospheric values and extra space to avoid the axis
      if( modes .eq. 2 )then
         if( vrho(1,3) .gt. 0.0d0 )then
            xmin = dlog10( vrho(1,3) ) -0.5d0
         endif
         if( vtem(1,3) .gt. 0.0d0 )then
            ymin = dlog10( vtem(1,3) ) -0.5d0
         endif
      endif

c..   window
      call PGSWIN(xmin,xmax,ymin,ymax)

c..   set color
      call PGSCI(1)
c..   tickmarks on x
      call PGBOX('BCNST',0.0,0,'BNST',0.0,0)
c..   tickmarks on y(left) and right
      call PGBOX('CNST',0.0,0,'BCNST',0.0,0)

c..   Left label
      call PGMTXT('L',2.0,0.5,0.5,ylabel)
      call PGMTXT('R',2.0,0.5,0.5,ylabel)
c..   Bottom label
      call PGMTXT('B',2.0,0.5,0.5,xlabel)

c..   set color (blue)
      call PGSCI(4)
c..   draw curve
      if( modes .eq. 2 )then
         call PGLINE(nxe,xx,yy)
      else
        call PGLINE(nx-1,xx,yy)
      endif

      if( nxflg .eq. 0 )then
c..   colors indicate fractions of x-axis interval
         do i = 1, nxe
            modxx = ( xb(i)/xb(nx) * 10.0 )
            if( modxx .gt. 10 .or.  modxx .lt. 0  )then
               write(*,*)modxx,i,nxe,xb(i),xb(nx)
               stop'pgplot 1'
            else
               mcolor = 2 + modxx
               call PGSCI(mcolor)
            endif
            if( ic(i) .eq. 0 )then
               call PGPNTS(1,xx(i),yy(i),2,1)
            else
               call PGPNTS(1,xx(i),yy(i),7,1)
            endif
         enddo

      else
c..   colors indicate nucleus of maximum abundance
c..   use abundances to color points (mass fractions)
         do i = 1, nxe
c..   find largest abundance in zone i+1
            xxmx = 0.0
            nnxx = 0
            do j = 1, nucpg
               n = nucp(j)
               if( n .gt. 0 )then
                  xxxx = x(n,i+1)
                  if( xxxx .gt. xxmx )then
                     xxmx = xxxx
                     nnxx = j
                  endif
               endif
            enddo
            if( xxmx .gt. 0.0 )then
               modxx = nnxx
            else
               do n = 1, nnuc
                  write(*,*)n,lz(n),ln(n),x(n,i+1)
               enddo
               write(*,*)'xxmx ',xxmx,i+1
               stop' pgplot maximum X <0'
            endif
            if( modxx .gt. nnuc .or.  modxx .lt. 1  )then
               write(*,*)modxx,i,nx,xb(i),xb(nx)
               write(*,*)nnuc,xxmx
               stop'pgplot colors'
            else
               mcolor = mod(modxx-1,15)+1
               call PGSCI(mcolor)
            endif

c..   use shape of point to indicate convection
            if( ic(i) .eq. 0 )then
               call PGPNTS(1,xx(i),yy(i),2,1)
            else
               call PGPNTS(1,xx(i),yy(i),7,1)
            endif
         enddo
      endif

c..   mark interesting zones......................
c..   most iterations
c..   reset color
      call PGSCI(2)
      if( kitmax .ge. 1 .and. kitmax .le. nxb )then
         call PGMOVE(xx(kitmax),ymax)
         call PGDRAW(xx(kitmax),ymin)
         call PGPTXT(xx(kitmax),ymin,0.0,0.0,'it')
      endif
c..   most restrictive time step
c..   reset color
      call PGSCI(1)
      if( kkman .ge. 1 .and. kkman .le. nx )then
         call PGMOVE(xx(kkman),ymax)
         call PGDRAW(xx(kkman),ymin)
         call PGPTXT(xx(kkman),ymin,0.0,0.0,'dth')
      endif

c..   nuclear network cycles
      call PGSCI(3)
      if( kcycmax .ge. 1 .and. kcycmax .le. nx )then
         call PGMOVE(xx(kcycmax),ymax)
         call PGDRAW(xx(kcycmax),ymin)
         if( nxid(kcycmax) .ge. 1 .and. nxid(kcycmax) .le. nnuc )then
            call PGPTXT(xx(kcycmax),ymax,0.0,0.0,xidk(kcycmax))
         else
            call PGPTXT(xx(kcycmax),ymax,0.0,0.0,'net')
         endif
      endif

c..   limits of loburn in log10 T(K)
      call PGSCI(3)
      call PGMOVE(xmin,5.7)
      call PGDRAW(xmax,5.7)

      call PGMOVE(xmin,7.0)
      call PGDRAW(xmax,7.0)


c..   make box showing eos table limits (Timmes table)
c..   set color to grey
      call PGSCI(14)
      call PGMOVE(-10.0,11.0)
      call PGDRAW(-10.0,4.0)
      call PGDRAW(11.0,4.0)
      call PGDRAW(11.0,11.0)
      call PGDRAW(-10.0,11.0)

c..   make trapezoids showing opal table limits
c..   opal21
      call PGMOVE(-14.75,3.75)
      call PGDRAW( -5.75,3.75)
      call PGDRAW(  9.1,8.7)
      call PGDRAW(  0.1,8.7)
      call PGDRAW(-14.75,3.75)
c..   boundary of table region without values
      call PGMOVE(4.3,7.1)
      call PGDRAW(4.1,7.2)
      call PGDRAW(3.9,7.3)
      call PGDRAW(4.2,7.4)
      call PGDRAW(4.0,7.5)
      call PGDRAW(4.9,7.8)
      call PGDRAW(4.7,7.9)
      call PGDRAW(6.5,8.5)
      call PGDRAW(6.6,8.7)

c..   reset color
      call PGSCI(1)


      if( modes .eq. 2 )then
c..   plot envelope using middle value of stencil (see cenvel)
c..   reset y array
         if( nvmax(3) .gt. kdme .or. nvmax(3) .lt. 1 )then
            write(*,*)nvmax(3),kdme
            stop'pgplot: envel array size'
         endif
         do i=1, nvmax(3)
            if( vtem(i,3) .gt. 0.0d0 .and. vrho(i,3) .gt. 0.0d0 )then
               yy(i) = dlog10( vtem(i,3))
               xx(i) = dlog10( vrho(i,3))
            else
               write(*,*)'pgplot: neg. envelope T,rho '
               write(*,'(i5,1p8e12.3)')i,vtem(i,3),vrho(i,3)
               stop'pgplot neg. envelope'
            endif
            setem(i) = yy(i)
            serho(i) = xx(i)
            ser(i) = vr(i,3)
            sem(i) = (xm(kk+1) + vm(i,3))/sol
c..   envelope mass are measured from the surface inward
         enddo

c..   set color (2=red)
         call PGSCI(2)
c..   draw curve
         call PGLINE(nvmax(3),serho,setem)
      endif

c..   reset color
      call PGSCI(1)


c..   third panel...................Luminosity AND VELOCITIES

c..   go to next panel
      call PGPANL(1, 3)
c..   erase screen
      call PGERAS

      if( ixflag .eq. 0 )then
c..   Radius
         xlabel = zlabel(22)
      elseif( ixflag .eq. 1 )then
         xlabel = 'log r'
      elseif( ixflag .eq. 2 )then
c..   Mass
         xlabel = zlabel(23)
      else
         xlabel = 'zone k'
      endif
c..   zone boundary quantities
      do i = 1, nxb
         xx(i) = xb(i)
      enddo
      call minmax(xx,1,nxb,xmin,xmax)
      if( ixflag .eq. 3 )then
c..revise to see envelope zones with 10 zone extra border
         if( modes .eq. 2 )then
            xmax = kk + nvmax(3) + 10
         else
            xmax = kk + 10
         endif
      endif
c..   directly read in gen.f from params.d
      if( fbot .ne. ftop .and. ftop .gt. 0. )then
         xmin = fbot
         xmax = ftop
      endif

c..   reset y array                         (4 = Luminosity )

      do i=1, nxb
         yy(i) = tl(2,i)/sollum
      enddo
c..   set envelope luminosity to interpolated value tl(2,kk), not 
c..   to previous base tl(1,kk+1) or tl(2,kk+1)
      yy(nxb) = tl(2,kk)/sollum

      call minmax(yy,1,nxb,ymin,ymax)
c..force L=0 to center of y axis
      ymin = -ymax
c..   reset color (green=3)
      call PGSCI(3)
c..   viewport
      call PGSVP(0.1, 0.9, 0.01, 1.0)
c..   window
      call PGSWIN(xmin,xmax,ymin,ymax)
c..   ticks
      call PGBOX('BCNST',0.0,0,'BNST',0.0,0)
c..   Left label
      ylabel ='Luminosity/Sol'
      call PGMTXT('L',2.0,0.5,0.5,ylabel)
c..   draw curve
      call PGLINE(nxb,xx,yy)

c..   reset y array              (nuclear + neutrino luminosity )
      yy(1) = 0.0
      xlnumin = yy(1)
      xlnumax = xlnumin
      do i = 2, nxb
         fact  = ( s(5,i) + s(4,i) )*dmh(i)/sollum
         yy(i) = yy(i-1) + fact
         xlnumin = amin1(xlnumin,yy(i))
         xlnumax = amax1(xlnumax,yy(i))
      enddo
      xlnukk = yy(nxb)
c..   window
c      call PGSWIN(xmin,xmax,xlnumin,xlnumax)
c..   reset color (15=grey)
      call PGSCI(15)
c..   Left label
      ylabel = 'L(n+n)'
      call PGMTXT('L',3.0,0.4,0.5,ylabel)
c..   draw curve
      call PGLINE(nxb,xx,yy)

c..   reset y array              (convective luminosity )
      yy(1) = 0.0
      do i = 1, nxb-1
         yy(i) = b(i)*a(i)/sollum
      enddo
c..   reset color (13=1,0,0.5=deep pink)
      call PGSCI(13)
c..   draw curve
      call PGLINE(nxb-1,xx,yy)
c..   Left label
      ylabel = 'Lcnv'
      call PGMTXT('L',3.0,0.8,0.5,ylabel)

c..   reset y variable to gravitational energy release
         do i = 1, nx
            yy(i) = lgrav(i)/sollum
         enddo
c..   reset color(8=gold)
c         call PGSCI(8)
c..   Left label
c         ylabel = 'L(grav)'
c         call PGMTXT('L',6.0,1.0,1.0,ylabel)
c..   draw points
c         call PGLINE(nx-1,xh,yy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c..   reset y array     (total nuclear heating/neutrino cooling )
      do i=1, nx
         yy(i) = ss(i+1) + snu(i+1)
      enddo

c     epmin = amin1( epmin, egmin )
c     epmax = amax1( epmax, egmax )
      if( epmin .ne. epmax )then
c..   reset color (5=turquoise)
         call PGSCI(5)
c..   window
         call PGSWIN(xmin,xmax,epmin,epmax)
c..   tick marks
c         call PGBOX('A',0.0,0,'T',0.0,0)
c..   Left label
         ylabel = 'eps (norm)'
         call PGMTXT('L',5.0,0.5,0.5,ylabel)
c..   draw points
         call PGPT(nx-1,xh,yy,2)

c..   reset y variable to gravitational energy release
         do i = 1, nx
            yy(i) = egrav(i)
         enddo
      endif

c..number of leqs cycles per zone
      do i = 1, nx
         yy(i) = float( ncyc(i) )
      enddo
      call minmax(yy,1,nx,ymin,ymax)
      if( ymin .ne. ymax )then
c..   reset color (14=darl grey)
         call PGSCI(14)
c..   window
         call PGSWIN(xmin,xmax,ymin,ymax)
c..   draw line
         call PGLINE(nx-1,xh,yy)
      endif

c..id of solven algorithm for zone
c      do i = 1, nx
c         yy(i) = float( itbrn(i) )
c      enddo
c      ymin = 0.
c      ymax = 3.
c..   reset color (14=dark grey)
c         call PGSCI(14)
c..   window
c         call PGSWIN(xmin,xmax,ymin,ymax)
c..   draw line
c         call PGLINE(nx-1,xh,yy)



c     (sound speed)  1.0e-5 to convert to km/sec (sound speed)
c..   reset y array
      do i=1, nx
         yy(i) = sound(i+1) * 1.0e-5
      enddo
      call minmax(yy,1,nx,ymin2,ymax2)

c..   get largest velocity for window
c..   reset y array                          (25=fluid velocity)
      do i=1, nxb
         yy(i) = zz(i,25) * 1.0e-5
      enddo
      call minmax(yy,1,nx,ymin1,ymax1)

c..   reset y array                           (5=conv velocity)
      do i=1, nxb
         yy(i) = zz(i,5) * 1.0e-5
      enddo
      call minmax(yy,1,nx,ymin,ymax)

c..   reset label                            (rotational velocity)
c..   reset y array
      do i=1, nxb
         yy(i) = omeg(i)*r(1,i) * 1.0e-5
      enddo

      call minmax(yy,1,nxb,ymin3,ymax3)

c      ymin = amin1(ymin3,ymin2,ymin1,ymin)
c      ymax = amax1(ymax3,ymax2,ymax1,ymax)
c.. rotation, fluid, convection
      ymin = amin1(ymin3,ymin1,ymin)
      ymax = amax1(ymax3,ymax1,ymax)

      if( cvsc .ne. 0.0 )then
         if( mode .ne. 0 )then
c..   allows read of params.d to set velocity scale for hydrostatic
c..   and damped (implicit) hydrodynamic modes
            ymin = -cvsc * 1.0e-5
            ymax =  cvsc * 1.0e-5
         else
c..   explicit hydro mode uses fluid velocity and input limits
            ymin = -cvsc * 1.0e-5
            ymax =  cvsc * 1.0e-5
            ymim = amin1(ymin,ymin1)
            ymax = amax1(ymax,ymax1)
         endif
      endif

c..   reset color
      call PGSCI(1)

      if( ymin .lt. ymax )then
c..   window
         call PGSWIN(xmin,xmax,ymin,ymax)
c..   tick marks
         call PGBOX('A',0.0,0,'CMST',0.0,0)

c..   reset label                            (25=fluid velocity)
         ylabel = zlabel(25)
c..   reset y array
         do i=1, nxb
            yy(i) = u(1,i) * 1.0d-5
         enddo

c..   color 8=gold
         call PGSCI(8)
c..   Right label
         ylabel = 'u(r)'
         call PGMTXT('R',3.0,0.5,0.5,ylabel)
c..   draw curve
         if( modes .eq. 2 )then
            call PGLINE(nxb,xb,yy)
         else
            call PGLINE(nx,xb,yy)
         endif

c..   reset label                            (sound speed)
         ylabel = 'sound'
c..   reset y array
         do i=1, kk-1
            yy(i) = sound(i+1) * 1.0d-5
         enddo
c..   reset color (11=medium blue)
         call PGSCI(11)
         call PGMTXT('R',4.0,0.5,0.5,ylabel)
c..   draw curve
         call PGLINE(kk-1,xh,yy)


c..   set color (11=medium blue)
         call PGSCI(11)
c..   plot envelope ---------------------------------------
         if( modes .eq. 2 )then
c..   reset y array
            do i=1, nvmax(3)
               if( ixflag .eq. 0 )then
                  xx(i) = vr(i,3)
               elseif( ixflag .eq. 1 )then
                  xx(i) = dlog10( vr(i,3) )
               elseif( ixflag .eq. 2 )then
                  xx(i) = sem(i)
               else
                  xx(i) = float( kk -i + nvmax(3) -2 )
               endif
               if( vrho(i,3) .gt. 0.0d0 .and. vp(i,3) .gt. 0.0d0 )then
                  yy(i) = sqrt( 1.5d0 * vp(i,3)/vrho(i,3) ) * 1.0d-5
               else
                  yy(i) = 0.0d0
               endif
            enddo

c..   draw curve
            call PGLINE(nvmax(3),xx,yy)
         endif

c..set velocity scale
c..   window
c         call PGSWIN(xmin,xmax,ymin3,ymax3)
c..   reset label                            (rotational speed)
c..   reset y array
         do i=1, nxb
            yy(i) = omeg(i) * r(1,i) * 1.0d-5 
         enddo
c..   reset color (12=violet)
         call PGSCI(12)
         ylabel = 'v(rot) km/s'
         call PGMTXT('R',5.0,0.5,0.5,ylabel)
c..   draw curve
         call PGLINE(nxb,xb,yy)
      endif

c..   reset velocity scale
c..   window
c         call PGSWIN(xmin,xmax,ymin,ymax)

c..   reset color (7=yellow)
      call PGSCI(7)
c..   reset label                            (convective speed)
      ylabel = zlabel(5)
c..   reset y array
      do i=1, kk-1
         yy(i) = h(i) * 1.0d-5
      enddo
      call PGMTXT('R',2.0,0.5,0.5,ylabel)
c..   draw curve
      call PGLINE(kk-1,xb,yy)
c..   draw points colored according to convection algorithm
      do i = 1, kk-1
         if( ic(i) .eq. 1 )then
c..   reset color (7=yellow) for full convection
            call PGSCI(7)
            call PGPT(1,xb(i),yy(i),2)
         elseif( ic(i) .eq. 2 )then
c..   reset color (3=green) for semiconvection
            call PGSCI(9)
            call PGPT(1,xb(i),yy(i),2)
c         elseif( ic(i) .eq. 3 )then
c         else
c..   reset color (12=violet) for time dependent convection
c            call PGSCI(12)
c            call PGPT(1,xb(i),yy(i),2)
         endif
      enddo

c..   set color (2=red)
      call PGSCI(2)

c..   plot envelope ----------------------------------------
      if( modes .eq. 2 )then
c..   reset y array
         do i=1, nvmax(3)
            if( ixflag .eq. 0 )then
               xx(i) = vr(i,3)
            elseif( ixflag .eq. 1 )then
               xx(i) = dlog10( vr(i,3) )
            elseif( ixflag .eq. 2 )then
               xx(i) = sem(i)
            else
               xx(i) = float( kk + nvmax(3) - i -2 )
            endif
            
            yy(i) = vvel(i,3) * 1.0d-5
         enddo

         do i=1, nvmax(3)
            if( yy(i) .ne. 0.0 )then
               call PGPT(1,xx(i),yy(i),2)
            endif
         enddo
c..   draw curve
         call PGLINE(nvmax(3),xx,yy)
      endif

c..   mark interesting zones......................
c..   most iterations
c..   reset color (2=red)
      call PGSCI(2)
      if( kitmax .ge. 2 .and. kitmax .le. nxb )then
c..fix 1 unit offset in r(1,k) and xh(k)
         call PGMOVE(xb(kitmax-1),ymin)
         call PGDRAW(xb(kitmax-1),ymax)
         call PGPTXT(xb(kitmax-1),ymax,0.0,0.0,'it')
      endif
c..   most restrictive time step
c..   reset color (1=white)
      call PGSCI(1)
      if( kkman .ge. 2 .and. kkman .le. nx )then
c..fix 1 unit offset in r(1,k) and xh(k)
         call PGMOVE(xh(kkman-1),ymin)
         call PGDRAW(xh(kkman-1),ymax)
         call PGPTXT(xh(kkman-1),ymax,0.0,0.0,'dth')
      endif

c..   nuclear network cycles
c..6=red+blue=magenta
c..maximum cycles: destruction of nucleus xid at zone k
      call PGSCI(6)
      if( kcycmax .ge. 1 .and. kcycmax .le. nx )then
         call PGMOVE(xh(kcycmax),ymax)
         call PGDRAW(xh(kcycmax),ymin)
         if( nxid(kcycmax) .ge. 1 .and. nxid(kcycmax) .le. nnuc )then
            call PGPTXT(xh(kcycmax),ymax,0.0,0.0,xidk(kcycmax))
         else
            call PGPTXT(xh(kcycmax),ymax,0.0,0.0,'net')
         endif
      endif

      lobu  = 2
      lobud = 0
      do j = 2,kk
         if( t(1,j) .gt. 1.0d7 )then
            lobu = j
         endif
         if(  t(1,j) .gt. 5.0d5 )then
            lobud = j
         endif
      enddo
c..11=blue+cyan
      call PGSCI(11)
      if( lobu .gt. 2 .and. lobu .lt. nx )then
         call PGMOVE(xh(lobu),ymax)
         call PGDRAW(xh(lobu),ymin)
         call PGPTXT(xh(lobu),ymin,0.0,0.0,'1e7')
      endif
      if( lobud .gt. 2 .and. lobud .lt. nx )then
         call PGMOVE(xh(lobud),ymax)
         call PGDRAW(xh(lobud),ymin)
         call PGPTXT(xh(lobud),ymin,0.0,0.0,'5e5')
      endif

c..   reset y label                        ( temperature = 2)
      ylabel = 'log T (norm)'
c..   reset y array
      do i=1, kk-1
         if( t(2,i+1) .gt. 0.0d0 )then
            yy(i) = dlog10( t(2,i+1) )
         else
            yy(i) = 0.0
         endif
      enddo
      call minmax(yy,1,kk-1,ymin,ymax)
      if( modes .eq. 2 )then
         ymin = dlog10( vtem(1,3) )
         ymin = amin1(ymin,3.4)
      endif

c..   window
      call PGSWIN(xmin,xmax,ymin,ymax)
c..   set color (1=white)
      call PGSCI(1)
c..   tickmarks on x(bottom=B,top=C) and y(left)
c     call PGBOX('A',0.0,0,'CMST',0.0,0)
c..   Left label
      call PGMTXT('L',4.0,0.5,0.5,ylabel)

      call PGLINE(kk-1,xh,yy)
c..   font scaling (size)
      call PGSCH(3.0*fontsize)

c..   set color for envelope value (2=red)
      call PGSCI(2)
      if( ixflag .eq. 0 )then
         call PGLINE(nvmax(3),ser,setem)
      elseif( ixflag .eq. 1 )then
         do j = 1, nvmax(3)
            ser(j) = alog10( ser(j) )
         enddo
         call PGLINE(nvmax(3),ser,setem)
      elseif( ixflag .eq. 2 )then
         call PGLINE(nvmax(3),sem,setem)
      else
         do j = 1, nvmax(3)
            sek(j) = nvmax(3) - j + 1 + kk
         enddo
         call PGLINE(nvmax(3),sek,setem)
      endif

c..   reset color (1=white)
      call PGSCI(1)


c..   fourth panel................ABUNDANCES.........................

c..   zone center quantities

c..   go to next panel
      call PGPANL(1, 4)
c..   erase screen
      call PGERAS
c..   reset label
c      ylabel = cnuc(nucp(1))
c..   reset y array
c      do i=1, nxe
c         yy(i) = x(nucp(1),i+1)*float( lz(nucp(1)) + ln(nucp(1)) )
c      enddo

c..   log abundances
c..   xlogmin and xlogmax are from gen.f read of params.d
      if( xlogmax .ge. 0.0 .and. xlogmin .ge. -2.0 )then
         ixlogflg = 1
         ymin = -0.05
         ymax =  1.05
      else
         ixlogflg = 0
         ymin = xlogmin - 0.1
         ymax = xlogmax + 0.1
      endif

c      if( ixlogflg .eq. 0 )call vlog(1,nxe,ymin,ymax,yy)

c..   viewport
      call PGSVP(0.1, 0.9, 0.01, 0.99)
c..   window
      call PGSWIN(xmin,xmax,ymin,ymax)
c..   ticks
      call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
c..   reset color
      call PGSCI(1)

c..   font scaling (size) to avoid clipping
      call PGSCH(2.8*fontsize)
      jjmin = 0
c..   only graph and id nuclei in nucp array; nucp set in abinit.f
      do j = 1, nucpg
         if( nucp(j) .gt. 0 )then
            jjmin = jjmin + 1
         endif
      enddo
      do j = 1, jjmin
c..   5 per line, 2 sides, 6 levels, take H1, He4, and lowest values

c..   reset color
c..   j starts at 1, sets initial color index to 1 (WHITE) for 
c..   consistency with panel 2 above, runs to 15, and repeats
         call PGSCI( mod(j-1,15) + 1 )
c..   reset y array
         do i=1, nxe
            yy(i) = x(nucp(j),i+1)
         enddo
         if( ixlogflg .eq. 0 )call vlog(1,nxe,ymin,ymax,yy)
         ylabel = cnuc(nucp(j))

c..   draw only those nuclei which have a gradient, to reduce clutter
c..   first, find largest relative difference
         dummy = 0.0
         do i=1,nxe
            dummy = amax1( abs(yy(i)-yy(1)), dummy )
         enddo

c..   force deuterium plot
         if( cnuc(nucp(j)) .eq. '    d'  .and. 
     1        yy(1) .gt. ymin+(ymax-ymin)*1.0e-6 )then
            dummy = 1.0d0
         endif

         if( dummy .gt. 1.0d-2 )then

            if( (j-1)/5 + 1 .le. 6 )then
c..   tabc1 is displacement in character height away from axis
c..   (L=left side)
c..   tabc2 is coordinate along this axis in fractions of length
               call PGMTXT('L',tabc1((j-1)/5+1),tabc2(mod(j-1,5)+1),
     1              0.5,ylabel)
            elseif( (j-1)/5 + 1 .le. 12 )then
c..   cycle through table tabc1 again
c               call PGMTXT('R',tabc1((j-1)/5+1-6)+1.0,
c     1              tabc2(mod(j-1,5)+1),0.5,ylabel)
c..adjusts x position of nuclei symbols on right hand side of panel

               call PGMTXT('R',tabc1((j-1)/5+1-6)+0.0,
     1              tabc2(mod(j-1,5)+1),0.5,ylabel)

            else
               stop' pgplot: error in symbol (He,etc) location'
            endif
c..   draw curve
            call PGLINE(nxe,xh,yy)
         endif
      enddo

c..   reset color
      call PGSCI(1)

c..   fifth panel.....................................(options)
      if( mode .ne. 0 )then

         if( nabflg .eq. 0 .or. modec .lt. 0 )then
c.....................iterated variables

c..   zone center quantities
            do i = 1, nx
               xx(i) = xh(i)
            enddo

            ymin0 = - cdelt*2.0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            ymax0 = - ymin0

c..   go to next panel
            call PGPANL(1, 5)
c..   erase screen
            call PGERAS
c..   reset label
            ylabel = zlabel(18)
c..   reset y array                         (18 = dT )
            do i=1, nx
               if( zz(i,2) .ne. 0.0 )then
                  yy(i) = zz(i,18)/zz(i,2)
               else
                  yy(i) = 0.0
               endif
            enddo
            call minmax(yy,1,nx,ymin,ymax)
c..   viewport
            call PGSVP(0.1, 0.9, 0.2, 0.99)
c..   window
            call PGSWIN(xmin,xmax,ymin0,ymax0)

c..   ticks
            call PGBOX('BCNST',0.0,0,'BNST',0.0,0)

c..   Bottom label
            call PGMTXT('B',2.0,0.5,0.5,xlabel)


c..   reset color
            call PGSCI(3)
c..   and character string
            call PGMTXT('L',3.0,0.25,0.5,ylabel)
c..   draw curve
            call PGLINE(nx,xx,yy)
c..   reset color
            call PGSCI(1)

c..   reset label                            (19= dV)
            ylabel = 'dV/3'
c..   reset y array
            do i=1, nx
               yy(i) = zz(i,19)*zz(i,1)/3.0
            enddo
            call minmax(yy,1,nx,ymin,ymax)
            if( ymin .lt. ymax )then
c..   window
               call PGSWIN(xmin,xmax,ymin0,ymax0)

c..   tick marks
               call PGBOX('BCNST',0.0,0,'BNST',0.0,0)
c..   color
               call PGSCI(7)
c..   Right(Left) label
               call PGMTXT('L',3.0,0.75,0.5,ylabel)
c..   draw curve
               call PGLINE(nx,xx,yy)
            endif
c..   reset color
            call PGSCI(1)

c..most rapidly depleted nucleus
            if( nucdep .gt. 0 .and. nucdep .le. netsize+1 )then
c..   reset label                            (dHe4)
               ylabel = cnuc(nucdep)
c..   reset y array
               do i=1, nx
                  yy(i) = (x(nucdep,i+1)-xold(nucdep,i+1))*xa(nucdep)
               enddo
               call minmax(yy,1,nx,ymin,ymax)
               if( ymin .lt. ymax )then
c..   window
                  call PGSWIN(xmin,xmax,ymin0,ymax0)
c..   tick marks
                  call PGBOX('BCNST',0.0,0,'BNST',0.0,0)
c..   color
                  call PGSCI(8)
c..   Right(Left) label
                  call PGMTXT('R',3.0,0.05,0.5,ylabel)
c..   draw curve
                  call PGLINE(nx,xx,yy)
               endif
c..   reset color
               call PGSCI(1)
            endif

c..   zone boundary quantities
            do i = 1, nxb
               xx(i) = xb(i)
            enddo


c..   reset label
            ylabel = zlabel(20)
c..   reset y array                         (20 = dR )
            do i=1, nxb
               if( zz(i,22) .ne. 0.0 )then
                  yy(i) = zz(i,20)/ zz(i,22)
               else
                  yy(i) = 0.0
               endif
            enddo

c..   color index
            call PGSCI(2)
c..   label
            call PGMTXT('L',2.0,0.25,0.0,ylabel)
c..   draw curve
            call PGLINE(nx,xx,yy)
c..   reset color
            call PGSCI(1)

c..   reset label
            ylabel = zlabel(21)
c..   reset y array                         (21 = dL )
c..   find maximum and minimum luminosity
            do i=1, nxb
               yy(i) = zz(i,4)
            enddo
            call minmax(yy,1,nxb,ymin,ymax)
c..   normalize dL to this range
            if( ymin .lt. ymax )then
               do i=1, nxb
                  yy(i) = zz(i,21)/(ymax - ymin) * ymax0
               enddo
c..   color
               call PGSCI(12)
c..   Left label
               call PGMTXT('L',2.0,0.75,0.0,ylabel)
c..   draw curve
               call PGLINE(nx,xx,yy)
            endif

c..   reset color
            call PGSCI(5)
c..   Bottom label
            call PGMTXT('B',2.0,0.5,0.5,xlabel)

c..   reset label
            ylabel = 'log dm'
c..   reset y array                         (29 = dmh )
            do i=1,nx-1
               if( dmh(i+1) .gt. 0.0d0 )then
                  yy(i) = dlog10(  dmh(i+1)/sol )
               else
                  write(*,'(2a5,8a12)')'kk','i+1','dmh(i+1)'
                  write(*,'(2i5,1p8e12.3)')kk,i+1,dmh(i+1)
                  write(*,'(a5,8a12)')'k','dmh','xm','V'
                  do j = 1, kk
                     write(*,'(2i5,1p8e12.3)')i,j,dmh(j),xm(j),v(2,j)
                  enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  stop'pgplot1: dmh<=0'
               endif
            enddo
c..   uses fixed scale for dmh
            ymax = dmhmax
            ymin = dmhmin

c..   size of window determined by zones on grid, not dmh(envelope)
c..   window
            call PGSWIN(xmin,xmax,ymin,ymax)
c..   ticks
            call PGBOX('A',0.0,0,'CMST',0.0,0)
c..   right label
            call PGMTXT('R',3.0,0.5,0.0,ylabel)
c..   draw curve
            call PGLINE(nx-1,xh,yy)

c..   plot envelope
            do i = 1, nvmax(3)-2
               yy(i) = (vm(i,3)-vm(i+1,3))/sol
               if( ixflag .eq. 0 )then
                  xx(i) = vr(i,3)
               elseif( ixflag .eq. 1 )then
                  if( vr(i,3) .gt. 0.0d0 )then
                     xx(i) = dlog10( vr(i,3) )
                  else
                     xx(i) = 0.0
                  endif
               elseif( ixflag .eq. 2 )then
                  xx(i) = 0.5*( sem(i) + sem(i+1) )
               elseif( ixflag .eq. 3 )then
                  xx(i) = float(kk - i + nvmax(3)-2 )
               endif
               if( yy(i) .gt. 0.0 )then
                  yy(i) = alog10( yy(i) )
               else
                  yy(i) = 0.0
               endif
            enddo

            call PGSCI(5)
c..   draw curve
            call PGLINE(nvmax(3)-2,xx,yy)

c..   reset color
            call PGSCI(1)


         elseif( nabflg .eq. 1 )then
c.....................convection variables..................

c..   zone boundary quantities
            nxm = nxb - 1
            do i = 1, nxb
               xx(i) = xb(i)
            enddo

c..   go to next panel
            call PGPANL(1, 5)
c..   erase screen
            call PGERAS
c..   reset label
            ylabel = zlabel(28)
c..   reset y array                         (28 = dnrad )
            do i=2, nxm
c               yy(i) = dnab(i)
               yy(i) = dnrad(i)
            enddo
            yy(1) = 0.0

            ymax = 0.8
            ymin = 0.0

c..   viewport
            call PGSVP(0.1, 0.9, 0.3, 0.99)
c..   window
            call PGSWIN(xmin,xmax,ymin,ymax)
c..   ticks
            call PGBOX('BCNST',0.0,0,'BNST',0.0,0)
c..   Bottom label
            call PGMTXT('B',2.0,0.5,0.5,xlabel)
c..   reset color (red=2)
            call PGSCI(2)

c..   left label
            call PGMTXT('L',3.0,0.25,0.5,ylabel)
c..   draw curve
            call PGLINE(nxm,xx,yy)

c..   draw points at each zone
c            do i = 2,nxm
c               if( ic(i) .eq. 0 )then
c..   nonconvective
c                  call PGSCI(1)
c                  call PGPNTS(1,xx(i),yy(i),7,1)
c               elseif( ic(i) .eq. 1 )then
c..   convective
c                  call PGSCI(2)
c                  call PGPNTS(1,xx(i),yy(i),7,1)
c               else
c..   semiconvective
c                  call PGSCI(11)
c                  call PGPNTS(1,xx(i),yy(i),-4,1)
c               endif
c            enddo


c..   reset color
            call PGSCI(1)

c..   reset label                            (27 = dnad)
            ylabel = zlabel(27)
c..   reset y array
            do i=2, nxm
               yy(i) = dnad(i)
            enddo
c..   uses previous symmetry at origin (or at inner zone)
            yy(1) = yy(2)

c..   color 7 = yellow
            call PGSCI(7)
            call PGMTXT('L',3.0,0.75,0.5,ylabel)
c..   draw curve
            call PGLINE(nxm,xx,yy)

c..   reset color 11=medium blue
            call PGSCI(11)
c..   reset label
            ylabel = 'rich'
            call PGMTXT('L',2.0,0.75,0.0,ylabel)
c..   richardson stiffness = v(rms)**2/g*H_P = doux+dnad-dnab
            yy(1) = 0.
            do i = 2, nxm
               if( p(1,i)*v(1,i) .gt. 0.0d0 )then
                  yy(i) = h(i)**2/( p(1,i)*v(1,i) )
               else
                  yy(i) = 0.
               endif
            enddo
            call PGLINE(nxm,xx,yy)

c..   reset label
            ylabel = zlabel(26)
c..   reset y array                         (26 = dnab )
            do i=2, nxm
               yy(i) = dnab(i)
            enddo

c..   reset color (green=3)
            call PGSCI(3)
            call PGMTXT('L',4.0,0.5,0.0,ylabel)
c..   draw curve
            call PGLINE(nxm,xx,yy)

c..   reset label
            ylabel = 'doux'
c..   reset y array                         (doux )
            do i=2, nxm
               yy(i) = doux(i)
            enddo
c..   reset color
            call PGSCI(8)
            call PGMTXT('L',4.0,0.0,0.0,ylabel)
c..   draw curve
            call PGLINE(nxm,xx,yy)

c..   set color (3=green=dnab)
            call PGSCI(3)
c..   plot envelope -------------------
c..   reset y array
            do i=1, nvmax(3)
               if( ixflag .eq. 0 )then
                  xx(i) = vr(i,3)
               elseif( ixflag .eq. 1 )then
                  xx(i) = dlog10( vr(i,3) )
               elseif( ixflag .eq. 2 )then
                  xx(i) = xb(nxb) + vm(i,3)/sol
               else
                  xx(i) = float( kk + nvmax(3) - i +1 )
               endif
               yy(i) = vnab(i,3)
            enddo

c..   draw curve
            call PGLINE(nvmax(3),xx,yy)

c..   reset color (7=yellow=dnad)
            call PGSCI(7)
c..   reset y array
            do i=1, nvmax(3)
               yy(i) = vnad(i,3)
            enddo
c..   draw curve
            call PGLINE(nvmax(3),xx,yy)

c..   reset color (2-red=dnrad )
            call PGSCI(2)
c..   reset y array
            do i=1, nvmax(3)
               yy(i) = vnrad(i,3)
            enddo
c..   draw curve
            call PGLINE(nvmax(3),xx,yy)

c..   reset color
            call PGSCI(1)
c..   reset label
            ylabel = 'log dm'
c..   reset y array                         (29 = dmh )
            do i=1,nx-1
               if( dmh(i+1) .gt. 0.0d0 )then
                  yy(i) = dlog10(  dmh(i+1)/sol )
               else
                  write(*,'(2a5,8a12)')'kk','i+1','dmh(i+1)'
                  write(*,'(2i5,1p8e12.3)')kk,i+1,dmh(i+1)
                  write(*,'(a5,8a12)')'k','dmh','xm','V'
                  do j = 1, kk
                     write(*,'(2i5,1p8e12.3)')i,j,dmh(j),xm(j),v(2,j)
                  enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  stop'pgplot2: dmh<=0'
               endif
            enddo
c..   uses fixed scale for dmh
            ymax = dmhmax
            ymin = dmhmin

c..   window
            call PGSWIN(xmin,xmax,ymin,ymax)
c
c..   ticks
            call PGBOX('A',0.0,0,'CMST',0.0,0)
c..   reset color
            call PGSCI(5)
c..   right label
            call PGMTXT('R',3.0,0.25,0.0,ylabel)
c..   draw curve
            call PGLINE(nx-1,xh,yy)

c..   plot envelope
            do i = 1, nvmax(3)-2
               yy(i) = (vm(i,3)-vm(i+1,3))/sol
               if( ixflag .eq. 0 )then
                  xx(i) = vr(i,3)
               elseif( ixflag .eq. 1 )then
                  if( vr(i,3) .gt. 0.0d0 )then
                     xx(i) = dlog10( vr(i,3) )
                  else
                     xx(i) = 0.0
                  endif
               elseif( ixflag .eq. 2 )then
                  xx(i) = 0.5*( sem(i) + sem(i+1) )
               elseif( ixflag .eq. 3 )then
                  xx(i) = float(kk) - i + nvmax(3)-2
               endif
               if( yy(i) .gt. 0.0 )then
                  yy(i) = alog10( yy(i) )
               else
                  yy(i) = 0.0
               endif

            enddo

c..   draw curve
            call PGLINE(nvmax(3)-2,xx,yy)


c..   reset color
            call PGSCI(6)
            ylabel = 'add zone'
c..   right label
            call PGMTXT('R',4.0,0.25,0.0,ylabel)


         elseif( nabflg .eq. 2 )then


c.....................monitor equation error

c..   zone center quantities
            do i = 1, nx-1
               xx(i) = xh(i)
            enddo

c..   go to next panel
            call PGPANL(1, 5)
c..   erase screen
            call PGERAS
c..   reset label
            ylabel = 'dM'
c..   reset y array                         (1=mass eqn.)
            do i=1, nx-1
                  yy(i) = aerr(1,i+1)/dmh(i+1)
            enddo
c..find max and min
            call minmax(yy,1,nx-1,ymin0,ymax0)
            do i=1, nx-1
                  yy(i) = aerr(2,i+1)/e(1,i+1)
            enddo
            call minmax(yy,1,nx-1,ymin,ymax)
            ymin0 = amin1(ymin0,ymin)
            ymax0 = amax1(ymax0,ymax)
            do i=2, nxb-1
                  yy(i) = aerr(3,i)/g(i)
            enddo
            yy(1) = 0.0
            call minmax(yy,1,nxb-1,ymin,ymax)
            ymin0 = amin1(ymin0,ymin)
            ymax0 = amax1(ymax0,ymax)
            do i=2, nxb-1
               yy(i) = aerr(4,i)/dnad(i)
            enddo
            yy(1) = 0.0
            call minmax(yy,1,nxb-1,ymin,ymax)
            ymin0 = amin1(ymin0,ymin)
            ymax0 = amax1(ymax0,ymax)
c..   limits found, so proceed
cccccccccccccccccccccccccccccccccccccccccccccc
c..over-ride; resid is iteration tolerance
            ymin0 = -2.0d0*resid
            ymax0 = -ymin0

c..   reset y array                         (1=mass eqn.)
            do i=1, nx-1
                  yy(i) = aerr(1,i+1)/dmh(i+1)
            enddo

c..   viewport
            call PGSVP(0.1, 0.9, 0.2, 0.99)
c..   window
            call PGSWIN(xmin,xmax,ymin0,ymax0)
c..   ticks
            call PGBOX('BCNST',0.0,0,'BNST',0.0,0)
c..   Bottom label
            call PGMTXT('B',2.0,0.5,0.5,xlabel)
c..   reset color
            call PGSCI(3)
c..   and character string
            call PGMTXT('L',3.0,0.25,0.5,ylabel)
c..   draw curve
            call PGLINE(nx-1,xx,yy)
c..   reset color
            call PGSCI(1)
c..   reset label                            (2=energy eqn.)
            ylabel = 'dE'
c..   reset y array
            do i=1, nx-1
               yy(i) = aerr(2,i+1)/e(1,i+1)
            enddo
c..   color
               call PGSCI(7)
c..   Right(Left) label
               call PGMTXT('L',3.0,0.75,0.5,ylabel)
c..   draw curve
               call PGLINE(nx-1,xx,yy)

c..   reset color
            call PGSCI(1)
c..   zone boundary quantities
            do i = 1, nxb-1
               xx(i) = xb(i)
            enddo
c..   reset label
            ylabel = 'du/dt'
c..   reset y array                         (3=accel. eqn. )
            do i=2, nxb-1
                  yy(i) = aerr(3,i)/g(i)
            enddo
            yy(1) = 0.0
c..   color index
            call PGSCI(2)
c..   label
            call PGMTXT('L',2.0,0.25,0.0,ylabel)
c..   draw curve
            call PGLINE(nxb-1,xx,yy)
c..   reset color
            call PGSCI(1)

c..   reset label
            ylabel = 'dnab'
c..   reset y array              (4 = gradient/L definition )
            do i=2, nxb-1
               yy(i) = aerr(4,i)/dnad(i)
            enddo
            yy(1) = 0.0
ccccccccccccccccccccccccccccccccccccc NOT GENERAL! ccccccccccccccccccc

c..   color
               call PGSCI(12)
c..   Left label
               call PGMTXT('L',2.0,0.75,0.0,ylabel)
c..   draw curve
               call PGLINE(nxb-1,xx,yy)

c..   reset color
            call PGSCI(5)
c..   Bottom label
            call PGMTXT('B',2.0,0.5,0.5,xlabel)

c..   reset label
c            ylabel = 'log dm'
c..   reset y array                         (29 = dmh )
c            do i=1,nx-1
c               if( dmh(i+1) .gt. 0.0d0 )then
c                  yy(i) = dlog10(  dmh(i+1)/sol )
c               else
c                  write(*,*)i,i+1,dmh(i+1),nx
c                  stop'pgplot: dmh<0'
c               endif
c            enddo
c..   uses fixed scale for dmh
c            ymax = dmhmax
c            ymin = dmhmin

c..   size of window determined by zones on grid, not dmh(envelope)
c..   window
c            call PGSWIN(xmin,xmax,ymin,ymax)
c..   ticks
c            call PGBOX('A',0.0,0,'CMST',0.0,0)
c..   right label
c            call PGMTXT('R',3.0,0.5,0.0,ylabel)
c..   draw curve
c            call PGLINE(nx-1,xh,yy)

c..   plot envelope
c            sollog = dlog10( sol )
c            do i = 1, jmaxz-2
c               yy(i) = 0.5*( sem(i) - sem(i+1) )
c               xx(i) = 0.5*( sem(i) + sem(i+1) )
c               if( yy(i) .gt. 0.0 )then
c                  yy(i) = alog10( yy(i) ) + sollog
c               else
c                  yy(i) = 0.0
c               endif
c            enddo
c            call PGSCI(5)
c..   draw curve
c            call PGLINE(jmaxz-2,xx,yy)

         else
c......................nabflg error
            write(*,*)'pgplot error in nabflg ',nabflg
            stop'pgplot: nabflg'
         endif


c..   reset color
         call PGSCI(1)



      else
c.....................hydro variables

c..   zone center quantities
         do i = 1, nx
            xx(i) = xh(i)
         enddo

         ymin0 = -0.01
         ymax0 = - ymin0

c..   go to next panel
         call PGPANL(1, 5)
c..   erase screen
         call PGERAS

c..   reset y label                        ( temperature = 2)
         ylabel = zlabel(2)
c..   reset y array
         do i=1, kk-1
            if( t(2,i+1) .gt. 0.0d0 )then
               yy(i) = dlog10( t(2,i+1) )
            else
               yy(i) = 0.0
            endif
         enddo
         if( modes .eq. 1 .or. modes .eq. 0 )then
            call minmax(yy,1,nx-1,ymin,ymax)
         else
            call minmax(yy,1,nx,ymin,ymax)
         endif

c..   viewport
         call PGSVP(0.1, 0.9, 0.2, 0.99)
c..   window
         call PGSWIN(xmin,xmax,ymin,ymax)

c..   ticks
         call PGBOX('BCNST',0.0,0,'BNST',0.0,0)

c..   Bottom label
         call PGMTXT('B',2.0,0.5,0.5,xlabel)

c..   reset color
         call PGSCI(3)
c..   and character string
         call PGMTXT('L',3.0,0.25,0.5,ylabel)
c..   draw curve
         call PGLINE(nx,xx,yy)
c..   reset color
         call PGSCI(1)


c..   reset y label                        ( pressure = 3)
         ylabel = zlabel(3)
c..   reset y array
         do i=1, kk-1
            if( p(2,i+1) .gt. 0.0d0 )then
               yy(i) = dlog10( p(2,i+1) )
            else
               yy(i) = 0.0
            endif
         enddo
         if( modes .eq. 1 .or. modes .eq. 0  )then
            call minmax(yy,1,nx-1,ymin,ymax)
         else
            call minmax(yy,1,nx,ymin,ymax)
         endif
c..   window
         call PGSWIN(xmin,xmax,ymin,ymax)
c..   reset color
         call PGSCI(7)
c..   ticks
         call PGBOX('A',0.0,0,'CMST',0.0,0)
c..   Bottom label
         call PGMTXT('R',3.0,0.25,0.5,ylabel)
c..   draw curve
         call PGLINE(nx,xx,yy)
c..   reset color
         call PGSCI(1)


c..   reset y label                        ( density = 1)
         ylabel = zlabel(1)
c..   reset y array
         do i=1, kk-1
            if( v(2,i+1) .gt. 0.0d0 )then
               yy(i) = -dlog10( v(2,i+1) )
            else
               yy(i) = 0.0
            endif
         enddo
         if( modes .eq. 1 .or. modes .eq. 0 )then
            call minmax(yy,1,nx-1,ymin,ymax)
         else
            call minmax(yy,1,nx,ymin,ymax)
         endif
c..   window
         call PGSWIN(xmin,xmax,ymin,ymax)
c..   reset color
         call PGSCI(8)
c..   Bottom label
         call PGMTXT('R',3.0,0.75,0.5,ylabel)
c..   draw curve
         call PGLINE(nx,xx,yy)
c..   reset color
         call PGSCI(1)

      endif

c..   rezoning 
      call PGSCH(4.0*fontsize)
      do i = 1,nx
         if( iadd(i) .ne. 0 )then
c..   reset color
            call PGSCI(6)
            call PGPT(1,xx(i),ymax-0.1*(ymax-ymin),30) 
         elseif( idel(i) .ne. 0 )then
c..   reset color
            call PGSCI(5)
            call PGPT(1,xx(i),ymin+0.1*(ymax-ymin),31) 
         endif
      enddo

c..   flush buffer
      call PGEBUF


      if( device1 .eq. '/ps' .or. device1 .eq. '/cps' )then
c..   terminate for hardcopy option
         call PGEND
         write(*,'(a20,2x,a5)')"pgplot using ",device
         write(*,*)'This option is "device" in input (see gen.f)'
         write(*,*)'graphical output in "pgplot.ps" file'
         stop' pgplot 2'
      endif

      return
      end

