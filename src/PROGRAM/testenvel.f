      program testenvel

c..   written darnett (wda) 12-18-99
c..   revised wda 9-24-05
c..   produces a set of envelopes using TYCHO subroutines
c..   steps through Teff holding L constant (jvar=0)
c..   steps through L holding Teff constant (jvar=1)
c..   steps through Teff,L holding R constant (jvar=2)

c..   pgplot output 
c..   calls getsurf.f which calls state.f
c..   calls abinit.f which reads ssystem.dat
c..   calls masses.f, kappinit.f, envel4.f

      implicit none

      include 'dimenfile'
      include 'cconst'
      include 'comod'
      include 'cenv'
      include 'cburn'
      include 'crate'

      real*8 yy(ndim),xx(ndim)
      real*8 fact,fact1,fact2,xsum,zhyd,zpop,xxc,xxo,ropal,ropals
      real*8 xlol0,telog0,xmass0,xmenv0,dlxlol,dtelog,drlog,aklog
      real*8 radlog0

      integer*4 ibu,idummy
      integer*4 kks,kp,nuc
      integer*4 nc12,nn14,no16

      logical tobe

c      integer*4 nc,itfinal,n,l,jvar,k,m
      integer*4 nc,itfinal,n,jvar,k,m,j
      real*8 g(kdm), a(kdm), accel
      real*8 teff,radius,pb,db,xmass,xhe,telog,xlol
      real*8 sgamma,tb,akb,xmenv,xminner
      real*8 xtlum,stl,sp,sr,st,sm,sv0,sd,stau,sent,sv
      real*8 tledd
      real*8 eps

      integer*4 modeg

      integer*4 lpldim,lprdim,nvardim,loopl,loopr
      parameter( lpldim=120, lprdim = 20, nvardim=12)

      real*4    garay(nvardim,lprdim,lpldim)

      real*4    xxg(kdm),yyg(kdm),xmin,xmax,ymin,ymax
      real*4    xmin0,xmax0,ymin0,ymax0
      real*4    xminf,xmaxf,yminf,ymaxf

      integer*4 pgbeg, linesty
      integer*4 ixvar, iyvar

      character*20 cxlabel, cylabel
      character*5  device
      character*72 text
      character*44 txt
      character*8  labl, labl1


      data modeg/0/,alphaml/2.0d0/,uuml/1.0d0/
c------------------------------------------------------

c      mixmode  = 0
      nopaleos = 0
      modes    = 2
      modec    = 2
      nopac    = 1
      eps      = 3.0d-7

      open(2,file='testenvel.in')

c..   input parameters

      read (2,*)text
      read (2,*)text
      read (2,*)text

c..   vary Te, L fixed
c     jvar = 0
c..   vary L, Te fixed
c     jvar = 1
c..   vary L, R fixed
c     jvar = 2

c..   variation flag
      read (2,*)txt,labl,jvar
      labl1 = 'jvar'
      if( labl .ne. labl1 )goto 2000

c..   loop on luminosity
      read (2,*)txt,labl,loopl
      labl1 = 'loopl'
      if( labl .ne. labl1 )goto 2000

c..   loop on radius
      read (2,*)txt,labl,loopr
      labl1 = 'loopr'
      if( labl .ne. labl1 )goto 2000

c..   initial luminosity in solar units
      read (2,*)txt,labl,xlol0
      labl1 = 'xlol'
      if( labl .ne. labl1 )goto 2000

c..   initial log Teff
      read (2,*)txt,labl,telog0
      labl1 = 'telog'
      if( labl .ne. labl1 )goto 2000

c..   initial log Radius
      read (2,*)txt,labl,radlog0
      labl1 = 'radlog'
      if( labl .ne. labl1 )goto 2000

c..   Total mass/sol
      read (2,*)txt,labl,xmass0
      labl1 = 'xmass'
      if( labl .ne. labl1 )goto 2000

c..   envelope mass/sol
      read (2,*)txt,labl,xmenv0
      labl1 = 'xmenv'
      if( labl .ne. labl1 )goto 2000


c..   change interval for log luminosity
      read (2,*)txt,labl,dlxlol
      labl1 = 'dlxlol'
      if( labl .ne. labl1 )goto 2000

c..   change interval for log Teff
      read (2,*)txt,labl,dtelog
      labl1 = 'dtelog'
      if( labl .ne. labl1 )goto 2000

c..   change interval for log radius
      read (2,*)txt,labl,drlog
      labl1 = 'drlog'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text

c..   pgplot graphics output device (eg, /xwin, /cps)
      read (2,*)txt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
c..   x variable choice
      read (2,*)txt,labl,ixvar
      labl1 = 'ixvar'
      if( labl .ne. labl1 )goto 2000

c..   y variable choice
      read (2,*)txt,labl,iyvar
      labl1 = 'iyvar'
      if( labl .ne. labl1 )goto 2000

c..   lines or points?
      read (2,*)txt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

c..   force graphical limits?
      read (2,*)txt,labl,xminf
      labl1 = 'xminf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmaxf
      labl1 = 'xmaxf'
      if( labl .ne. labl1 )goto 2000

      if( xminf .ne. xmaxf )then
         write(*,*)'      x override'
      endif

      read (2,*)txt,labl,yminf
      labl1 = 'yminf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymaxf
      labl1 = 'ymaxf'
      if( labl .ne. labl1 )goto 2000

      if( yminf .ne. ymaxf )then
         write(*,*)'      y override'
      endif

      read (2,*)text

      close(2)
c........................................................

c..   kk is last internal zone, 
c..   kk+1 is envelope, 
c..   kp is photospheric zone
      nc  = 1
      kk = 2
      kp = kk+2
      kks = 2

      dmh(kk) = 1.0d-1 * xmenv0 * sol

c..   input composition for state.f
      do k = 1,kp
         do n = 1, ndim
            x(n,k) = 0.0d0
         enddo
      enddo
c....................................................
c..   resource file for analysis programs
      inquire(file='net.rc',exist=tobe)
      if( .not. tobe )then
         write(*,*)'testenv: no net.rc'
         stop'testenv 1'
      else
         write(*,*)'testenv: net.rc exists'
      endif
      write(*,*)'testenv: checking netrc'

      open(30,file='net.rc',status='old')
      ibu = 1
 100  read(30,'(3i5,a5,0pf10.4,1pe12.4)',end=101)
     1     idummy,lz(ibu),ln(ibu),cnuc(ibu),qq(ibu),xx(ibu)

      if( lz(ibu) .eq. 2 .and. ln(ibu) .eq. 2 )goto 101
      ibu = ibu + 1
      goto 100
 101  continue
      netsize = ibu
      read(30,'(10i5)')nuc
      if( ibu .ne. netsize )then
         write(*,*)'conflict in dimenfile and net.rc ',ibu,netsize
         stop'testenv 2'
      endif
c..   net.rc is consistent
      write(*,*)'testenv: net.rc is left unchanged ',ibu,' nuclei'
      close(30)

ccccccccccccccccccccccccccccccccccccccc
c..fix for testenvel
      do n = 1, ndim
         nz(n) = lz(n)
         nn(n) = ln(n)
         xid(n) =cnuc(n)
      enddo

c..   input nuclear data
      newnet = 0

c..   this may be redundant with the logic above
c..   the logic above could be edited for brevity?
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call abinit

      call masses(qq,w,nnuc)

      do n = 1, nnuc
         qex(n) = qq(n)
         xa(n) = dble( lz(n) + ln(n) ) + qex(n)/931.5d0
      enddo

c      write(*,*)'fudge so that ag88-->gs98 : not validated!!!'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   fudge so that ag88-->gs98 approximately, and OPAL table
c..   values are reproduced
c..   adjust H-He
c      fact = 0.0183
c..   adjust metals - H+He
c      fact1= 0.031
c      xx(nnuc-1) = xx(nnuc-1) - fact -fact1
c      xx(nnuc)   = xx(nnuc)   + fact -fact1

c..   renormalize and reset Ye
      xsum = 0.0d0
      do n = 1, nnuc
         xsum = xsum + xx(n)
      enddo
      write(*,*)'testenv: sumX-1 = :',xsum-1.0d0
      do n = 1, nnuc
         xx(n) = xx(n)/xsum
         yy(n) = xx(n)/dble( lz(n) + ln(n) )
cccccccccccccccccccccccccccccccccccccccccccc
      enddo
      yy(ndim) = 0.0d0
      do n = 1, nnuc
         yy(ndim) = yy(ndim) + yy(n)*dble( lz(n) )
      enddo
c..   check renormalization
      xsum = 0.0d0
      do n = 1, nnuc
         xsum = xsum + xx(n)
      enddo
      write(*,*)'testenv: sumX-1 = :',xsum-1.0d0,' renormalized'

c..   modify C12 and He4
      fact2 = 0.0d0
      write(*,*)'revising C12 by ',fact2,' and He4 by',-fact2
      xx(nnuc) = xx(nnuc) - fact2
      nc12 = 0
      nn14 = 0
      no16 = 0
      do n = 1, nnuc
         if( lz(n) .eq. 6 .and. ln(n) .eq. 6 )then
            nc12 = n
         endif
         if( lz(n) .eq. 7 .and. ln(n) .eq. 7 )then
            nn14 = n
         endif
         if( lz(n) .eq. 8 .and. ln(n) .eq. 8 )then
            no16 = n
         endif
      enddo
      if( nc12 .lt. 1 )then
         write(*,*)nc12,' nc12 not found'
         stop'testenv 3'
      endif
      xx(nc12) = xx(nc12) + fact2

c..   renormalize and reset Ye
      xsum = 0.0d0
      do n = 1, nnuc
         xsum = xsum + xx(n)
      enddo
      write(*,*)'testenv 4: sumX-1 = :',xsum-1.0d0
      do n = 1, nnuc
         xx(n) = xx(n)/xsum
         yy(n) = xx(n)/dble( lz(n) + ln(n) )
cccccccccccccccccccccccccccccc
      enddo
      yy(ndim) = 0.0d0
      do n = 1, nnuc
         yy(ndim) = yy(ndim) + yy(n)*dble( lz(n) )
      enddo
c..   check renormalization
      xsum = 0.0d0
      do n = 1, nnuc
         xsum = xsum + xx(n)
      enddo
      write(*,*)'testenv 5: sumX-1 = :',xsum-1.0d0,' renormalized'


      zhyd  = yy(ndim-2)
      xhe   = xx(ndim-1)
      zpop  = 1.0d0 - zhyd - xhe
c..for state.f
      zpop0 = zpop
      write(*,'(3(a12,1pe11.3))')'hyd',zhyd,'he4',xhe,'z',zpop

c..   enchanced alpha elements
      xxc = 0.0d0
      xxo = 0.0d0
      do n = 1, nnuc
         if( lz(n) .eq. 6 .and. ln(n) .eq. 6 )then
            xxc = xxc + xx(n)
         endif
         if( lz(n) .eq. 7 .and. ln(n) .eq. 7 )then
            xxc = xxc + xx(n)*0.5d0
            xxo = xxo + xx(n)*0.5d0
         endif
         if( lz(n) .eq. 8 .and. ln(n) .eq. 8 )then
            xxo = xxo + xx(n)
         endif
      enddo
      write(*,'(2(a5,1pe12.3))')'Xc',xxc,'Xo',xxo
c..   define xxc and xxo as EXCESS cno relative to solar
c     xxc = xxc - 3.823d-3
c     xxo = xxo - 1.082d-2
c     write(*,'(2(a5,1pe12.3))')'dXc',xxc,'dXo',xxo
      do k = 1,kp
c..   define the x array for state.f
         x(ndim,k)   = yy(ndim)
         do n = 1, nnuc
c           x(n,k) = yy(n)
            x(n,k) = xx(n)
         enddo
      enddo

c..   abundance arrays filled...............................

      call kappinit

      write(*,'(4(a10,0pf10.3))')'mass',xmass0,'Menv/M',xmenv0

      do m = 1, loopr
         do l = 1, loopl

c..   initial values
            if( jvar .eq. 0 )then
               xlol   = xlol0*1.0d1**( dlxlol*dble(m-1) )
               xtlum  = xlol * sollum
c..   march through teff(l) for fixed luminosity(m)
               telog  = telog0 + dtelog*dble(l-1)
               teff   = 10.0d0**telog
               radius = dsqrt( xtlum /(sigma*teff**4)/pi4 )

            elseif( jvar .eq. 1 )then
c..   march through L for fixed Teff, mass
               xlol   = xlol0*1.0d1**( dlxlol*dble(l-1) )
               xtlum  = xlol * sollum
               telog  = telog0
               teff   = 10.0d0**telog
               radius = dsqrt( xtlum /(sigma*teff**4)/pi4 )

            elseif( jvar .eq. 2 )then
c..   march through L(l) for fixed radius(m)
               radius = 10.0d0**( radlog0 + drlog*dble(m-1) )
               telog  = telog0 + dtelog*dble(l-1)
               teff = 10.0d0**telog
               xtlum = pi4*radius**2 * sigma*teff**4
               xlol  = xtlum/sollum

            else
               write(*,*)jvar
               stop'jvar error in testenvel'
            endif
            xmass  = xmass0*sol
            xmenv  = xmenv0*xmass
            xminner = xmass - xmenv

c..   estimate surface pressure, rho for first guess in iteration
            pb  = 1.0d2
            db   = pb/rgas/teff

            v(nc,kp) = 1.0d0/db
            t(nc,kp) = teff
            g(kp)    = grav*xmass/radius**2
            a(kp)    = pi4 * radius**2
            r(nc,kp) = radius
            xm(kp)   = xmass

            call getsurf(g,kp,nc,itfinal,accel)

c..   pressure and opacity at boundary
            db       = 1.0d0/v(nc,kp)
            pb       = p(nc,kp)
            tb       = t(nc,kp)
            akb      = ak(kp)

            sgamma = 4.0d0*arad*teff**3/(3.0d0*db*rgas)

            sm  = 0.0d0
            sr  = radius
            stl = xtlum

            tledd = pi4*crad*grav*xmass/akb
            if( db .gt. 0.0d0 .and. tb .gt. 0.0d0 )then
               ropals = dlog10( db/(1.0d-6*tb)**3 )
            else
               write(*,*)'ropal error: ',db,tb
               stop'testenvel 6'
            endif

            call envel4(kp,nc,db,tb,pb,xmass,xmenv,
     1           stl,sp,sr,st,sent,sm,sv,sd,stau,modeg,eps)

c            write(*,'(i5,1p8e12.3)')jmaxz,xmass,xmenv,sp,sr,st,pb,tb

             if( sv .gt. 0.0d0 )then
                sd = 1.0d0/sv
             else
                write(*,*)'ERROR(testenvel): sv < 0 ',sv
                stop'testenvel'
             endif


c           call envel(kp,nc,db,tb,xmass,xmenv,
c    1           stl,sp,sr,st0,sm,sv0,sd,stau,modeg)

            if( sd .gt. 0.0d0 .and. st .gt. 0.0d0 )then
               ropal = dlog10( sd/(1.0d-6*st)**3 )
            else
               write(*,*)'ropal error: ',sd,st
               stop'testenvel 7'
            endif
            if( zak(jmaxz) .gt. 0.0d0 )then
               aklog = dlog10( zak(jmaxz)  )
            else
               write(*,*)'ropal error: ',jmaxz,zak(jmaxz)
               stop'testenvel 8'
            endif

            if( xlol .gt. 0.0d0 )then
               garay(1,m,l) = dlog10( xlol )
            else
               garay(1,m,l) = -20.
            endif
            garay(2,m,l) = telog
            garay(3,m,l) = radius
            garay(4,m,l) = pb
            garay(5,m,l) = db
            garay(6,m,l) = sr/radius
            if( sp .gt. 0.0d0 )then
               garay(7,m,l) = dlog10( sp )
            else
               garay(7,m,l) = -20.
            endif
            garay(8,m,l) = sd
            garay(9,m,l) = st
            garay(10,m,l) = ropal
            garay(11,m,l) = ropals
            garay(12,m,l) = aklog

            if( l .le. 1 .and. loopr .ge. 2 )then
               write(*,'(2a5,2a8,3a10,2a8,2a10,3a8)')
     1              'l','nvmax','L/Lsol','log Te','Radius','P(surf)',
     2              'rho(s)','Rj/R','logPj',
     3              'rho(j)','T(j)','lropalj','lropal','lg(ak)'
            endif

            if( loopr .le. 2 )then
               write(*,'(2i5,0p2f8.3,1p3e10.2,0p2f8.4,1p2e10.2,
     1         0p3f8.3)') l,nvmax(3), (garay(k,m,l),k=1,nvardim)
            endif

            if( st .gt. 1.0d9 )goto 1000
         enddo
c        if( m .eq. 1 .and. loopr .gt. 1)then
c           write(*,'(a5,8a12)')'m','R','L/Lsol','log Te','Pj','Tj','rj'
c        endif
c        if( m .eq. loopr .and. m .gt. 1)then
c           write(*,'(a5,8a12)')'m','R','L/Lsol','log Te','Pj','Tj','rj'
c        endif
      enddo

 1000 continue

      close(7)

c..   graphics section................................................
c..   scratch file for character conversion
      open(10,file='scratch.pg')

      IF ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
      call PGASK (.FALSE.)
c..   roman font
      call PGSCF(2)
c..   font scaling (size)
      call PGSCH(1.3)
c..   line width
      call PGSLW(4)
c..   new page
      call PGPAGE

      xmin = garay(ixvar,1,1)
      xmax = garay(ixvar,1,1)
      ymin = garay(iyvar,1,1)
      ymax = garay(iyvar,1,1)

      do m = 1, loopr
         do l = 1, loopl
            xmin = amin1( garay(ixvar,m,l) , xmin)
            xmax = amax1( garay(ixvar,m,l) , xmax)

            ymin = amin1( garay(iyvar,m,l) , ymin)
            ymax = amax1( garay(iyvar,m,l) , ymax)
         enddo
      enddo

      if( ixvar .eq. 1 )then
         cxlabel = 'log L/sol'
      elseif( ixvar .eq. 2 )then
         cxlabel = 'log Te'
c..invert order as in HR diagram
         fact = xmin
         xmin = xmax
         xmax = fact
      else
         cxlabel = 'x?'
      endif
      if( iyvar .eq. 7 )then
         cylabel = 'log Pj'
      else
         cylabel = 'y?'
      endif

c..   override x
      if( xminf .ne. xmaxf )then
         xmin = xminf
         xmax = xmaxf
      endif
c..   override y
      if( yminf .ne. ymaxf )then
         ymin = yminf
         ymax = ymaxf
      endif

c..   adjust for more aesthetic margins
      xmin0 = xmin - (xmax-xmin)*0.05
      xmax0 = xmax + (xmax-xmin)*0.05
      ymin0 = ymin - (ymax-ymin)*0.05
      ymax0 = ymax + (ymax-ymin)*0.05

c..   window
      call PGSWIN(xmin0,xmax0,ymin0,ymax0)
c..   set color
      call PGSCI(1)
c..   tickmarks on x(bottom=B,top=C) and y(left)
      call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

c..   Left label
      call PGMTXT('L',2.0,0.5,0.5,cylabel)
c..   bottom label
      call PGMTXT('B',2.0,0.5,0.5,cxlabel)

      do m = 1, loopr
         do l = 1, loopl
            xxg(l) = garay(ixvar,m,l)
         enddo
         do l = 1, loopl
            yyg(l) = garay(iyvar,m,l)
         enddo

         if( m .lt. 14 )then
            call PGSCI(m+1)
         else
            call PGSCI(m-13)
         endif


c..   draw curve
         if( linesty .eq. 0 )then
            call PGLINE(loopl,xxg,yyg)
         elseif( linesty .eq. 1 )then
            call PGPT(loopl,xxg,yyg,20)
         else
            call PGLINE(loopl,xxg,yyg)
            call PGPT(loopl,xxg,yyg,20)
            if( jvar .eq. 2 .and. iyvar .eq. 7 )then
               do l = 1, loopl-1
                  if( yyg(l+1) .le. yyg(l) )then
c..put white dot if Pj gradient is negative in L
                     call PGSCI(1)
                     call PGPT(1,xxg(l),yyg(l),22)
                  endif
               enddo
            endif
         endif

      enddo

c..   window
      call PGSWIN(0.0,1.0,0.0,1.0)
c..   set color
      call PGSCI(1)
      if( jvar .eq. 2 )then
         call PGPTXT(0.4,1.05,0.0,0.5,'Lines of constant Radius')
      else

      endif
      pause

      call PGEND

      stop'successful termination'

 2000 continue

      write(*,*)'testenvel: error in  input file: testenvel.in'
      stop'testenvel error'

      end
