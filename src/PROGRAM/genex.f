      program genex
c..   7-29-06
      implicit none
c----------------------------------------------------------------------
c..   reads new gen format, writes in new format
c..   after revising composition (does not read old composition)
c..   uses mass fractions in x(n,k)
c     
c     reads old.model, writes new.model in local directory
c..   solid body rotation an option
c----------------------------------------------------------------------

c..   comod common block must follow dimenfile
      include 'dimenfile'
      include 'comod'

      integer*4 izbu(ndim),inbu(ndim),ibu,idummy
      character*5 cbu

      include 'compu'
      include 'cconst'
      include 'czone'
      include 'cgen'
      include 'cburn'
      include 'ceoset'
      include 'cnabla'

c...  used in eos subroutine
      integer*4 ikd(kdm),ikt(kdm), intab(kdm)
      common/ineos/ikd,ikt,intab

      real*8    scr(kdm), hetoz
      real*8    xmsol, sum, yions, yee, ytot, delta, fact
      real*8    fm, fr, xhe, xc12, xo16, sumx, sumz
c..   array for new abundances
      real*8    rx(ndim)

c..   angular velocity of rotation, critical angular velocity
      real*8    omegbar,omegcrit,ajaysum

      real*8    tlold(kdm),aki,aom,dtry,dtmin

      integer*4 ifm, k, nc, n

      integer*4 nh1, nhe4, j

      integer*4 lumin

      character*9  filein
      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2  txtxt
      character*1  clause

c..   solar system abundances taken from net.rc file
      real*8 zscale

      logical tobe

      data txtxt/'  '/
c..   initialize rx abundance vector
      data rx/ ndim*0.0d0 /
c-----------------------------------------------------------------
      inquire(file='genex.in',exist=tobe)
      if( tobe )then
         open(2,file='genex.in')
      else
         stop'genex: no genex.in file'
      endif

      inquire(file='old.model',exist=tobe)

c      if(  tobe )then
c         open(8,file='old.model')
c      else

c      if(  tobe )then
c         open(8,file='old.model')
c      else
      if( .not. tobe )then

         stop'genex: no old.model file'
      endif

c..   reset the time
      time = 0.0d0
      write(*,'(a30,1pe12.3)')'time reset to ',time

c...............................................................
c..   input parameters for model modification from genex.in

      read (2,*)text
      write(*,*)text
c..   dummy read
      read (2,*)text
      write(*,*)text

      read (2,*)text
      write(*,*)text

      read (2,*)txt,labl,ifm
      write(*,*)txt,txtxt,labl,ifm
      labl1 = 'ifm'
      if( labl .ne. labl1 )goto 2000


      if( ifm .lt. 0 .or. ifm .gt. 2 )then
         write(*,*)'genex: read ifm',ifm
         write(*,*)'this option outside valid range'
         stop'genex'
      endif

      read (2,*)txt,labl,lumin
      write(*,*)txt,txtxt,labl,lumin
      labl1 = 'lumin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nopac
      write(*,*)txt,txtxt,labl,nopac
      labl1 = 'nopac'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nopaleos
      write(*,*)txt,txtxt,labl,nopaleos
      labl1 = 'nopaleos'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
      write(*,*)text

      read (2,*)txt,labl,fm
      write(*,*)txt,txtxt,labl,fm
      labl1 = 'fm'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,fr
      write(*,*)txt,txtxt,labl,fr
      labl1 = 'fr'
      if( labl .ne. labl1 )goto 2000

      if( ifm .eq. 0 )then
         if( fm .le. 0.0d0 )then
            write(*,*)'Gen: fm error, ',fm
            stop
         endif
         if( fr .le. 0.0d0 )then
            write(*,*)'Gen: fr error, ',fr
            stop
         endif
      endif

      read (2,*)text
      write(*,*)text
c...............................................

c..   zscale=1 for population I
      read (2,*)txt,labl,zscale
      write(*,*)txt,txtxt,labl,zscale
      labl1 = 'zscale'
      if( labl .ne. labl1 )goto 2000

c..   hetoz=1 gives as much new He4 as heavy elements z (by mass)
      read (2,*)txt,labl,hetoz
      write(*,*)txt,txtxt,labl,hetoz
      labl1 = 'hetoz'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
      write(*,*)text

c...............................................

c..   zscale
      read (2,*)txt,labl,zpop
      write(*,*)txt,txtxt,labl,zpop
      labl1 = 'zpop'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xhe
      write(*,*)txt,txtxt,labl,xhe
      labl1 = 'xhe'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xc12
      write(*,*)txt,txtxt,labl,xc12
      labl1 = 'xc12'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xo16
      write(*,*)txt,txtxt,labl,xo16
      labl1 = 'xo16'
      if( labl .ne. labl1 )goto 2000

      read (2,*)note
      write(*,*)note

      read (2,*)txt,labl,omegbar
      write(*,*)txt,txtxt,labl,omegbar
      labl1 = 'omegbar'
      if( labl .ne. labl1 )goto 2000
      if( omegbar .lt. 0.0d0 )then
         write(*,*)'genex input error, omegbar < 0 ',omegbar
         goto 2000
      elseif( omegbar .eq. 0.0d0 )then
         write(*,*)'no rotation, omegbar = ',omegbar
      else
         write(*,*)'solid body rotation with angular velocity = ',
     1        omegbar,' radian/second'
      endif

      read (2,*)text
      write(*,*)text
c.............................................
      read (2,*)txt,labl,mixmode
      write(*,*)txt,txtxt,labl,mixmode
      labl1 = 'mixmode'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
      write(*,*)text

      close(2)

c..   get nuclear network info
c..   nuclear and solar abundance resource file for analysis programs
      inquire(file='net.rc',exist=tobe)
      if( .not. tobe )then
         write(*,*)'genex: no net.rc'
         stop'genex'
      endif

      write(*,*)'genex: checking netrc'

      open(30,file='net.rc',status='old')
      ibu = 1
c..initializing qex and solarx
 100  read(30,'(3i5,a5,0pf10.4,1pe12.4)',end=101)
     1     idummy,izbu(ibu),inbu(ibu),cbu,qex(ibu),solarx(ibu)
      lz(ibu)   = izbu(ibu)
      ln(ibu)   = inbu(ibu)
      cnuc(ibu) = cbu
      xa(ibu) = dble( lz(ibu) + ln(ibu) ) + qex(ibu)/931.5d0
      if( izbu(ibu) .eq. 2 .and. inbu(ibu) .eq. 2 )goto 101
      ibu = ibu + 1
      goto 100
 101  continue

c..   network size for this net.rc
      netsize = ibu
c..   net.rc is consistent
      read(30,'(10i5)')nucp

      write(*,*)'net.rc network has netsize = ',netsize

c..   determine nonzero entries
      mnucpg = 0
      do j = 1, nucpg
         if( nucp(j) .ne. 0 )then
            mnucpg = mnucpg + 1
         endif
      enddo

      write(*,*)'genex: net.rc exists and is left unchanged'
      close(30)

c..   locate special nuclei
      do n = 1, nnuc
         if( lz(n) .eq. 6 .and. ln(n) .eq. 6 )lc12 = n
         if( lz(n) .eq. 7 .and. ln(n) .eq. 7 )ln14 = n
         if( lz(n) .eq. 8 .and. ln(n) .eq. 8 )lo16 = n
      enddo

c..   special nuclei (enhanced C12 and O16)
      do n = 1, ndim
         if(     lz(n) .eq. 6 .and. ln(n) .eq. 6 )then
            rx(n) = rx(n) + xc12
            lc12 = n
         elseif( lz(n) .eq. 8 .and. ln(n) .eq. 8 )then
            rx(n) = rx(n) + xo16
            lo16 = n
         elseif( lz(n) .eq. 2 .and. ln(n) .eq. 2 )then
            rx(n) = xhe
         endif
      enddo
      nh1  = nnuc -1
      nhe4 = nnuc

c..   file for i(nitial) model (formatted data)
      filein = 'old.model'

      write(*,*)'filein is ',filein
cccccccccccccccccccc

      inquire(file=filein,exist=tobe)
      if( tobe )then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         call system('file imodel >dummy')
c         call system('file imodel')

c..get first character (blank for formatted)
         close(13)
         open(13,file=filein,status='OLD')
         rewind(13)
         read(13,'(a72)')note
         write(*,'(a72)')note
         close(13)
c         call system('rm dummy')
c..flush clean
         clause = ' '
         clause = note(1:1)
         write(*,*)note(1:1)

c..note: does not give correct answer!!!!!!!!!!!!!!!!!!!!!!!!!!!
c         inquire(file=filein,FORMATTED=clause)
c         write(*,*)clause
c

         if( clause .eq. 'T' )then

c..   formatted
            write(*,*)'reading formatted file ',filein
            open(8,file=filein,form='formatted',status='old')

            rewind 8
            nc = 1

            
            read(8,'(a72)')note
            write(*,*)note
            rewind 8

c            if( note(7:7) .eq. '7' .or. note(7:7) .eq. '8' )then
               write(*,*)'CALLING READF'
               call readf(note,nc)
c            else
c               write(*,*)'CALLING NREADF'
c               call nreadf(note,nc)
c            endif

         elseif( clause .eq. 'H' )then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   workaround: inquire not identifying unformatted files
c..   binary
            write(*,*)'reading binary file ',filein
            open(8,file=filein,form='unformatted',status='old')

            rewind 8
            nc = 1

            call readb(note,nc)

         else
            stop'gen format read error'
         endif

      else
         write(*,*)'gen: no file ',filein,' in directory'
         stop'gen'
      endif

      close(8)


c..............................................................
      if( ifm .eq. 0 )then
         write(*,*)'ifm = ',ifm,' M,R scaling'
      elseif( ifm .ge. 1 .and. ifm .le. 4 )then
         write(*,*)'ifm = ',ifm,' abundance changing'
         fm = 1.0d0
         fr = 1.0d0
      else
         write(*,*)'ifm= ',ifm
         stop'genex error: ifm'
      endif
c..   mass, radius scaling
      call homol(fm,fr)
c..............................................................

c..   define mass coordinate
      dmi(1)   = dmh(2)
      dmi(kk+1) = dmh(kk+1)
      do k = 2, kk
         xm(k)  = xm(k-1) + dmh(k)
         dmi(k) = 0.5d0*( dmh(k+1) + dmh(k) )
      enddo
      xm(kk+1) = xm(kk) + dmh(kk+1)

      xmsol = xm(kk)/sol
      write(*,*)' grid MASS =',xmsol,' solar masses,',' model =',model
      write(*,*)' Envelope Mass =',dmh(kk+1)/sol,', Total Mass =',
     1     xmsol + dmh(kk+1)/sol

c...............................................................
      write(*,*)'reaction network has nnuc =',nnuc,' nuclei'

      if( ibu .ne. nnuc )then
         write(*,*)'net.rc inconsistent with dimenfile'
         write(*,*)ibu,' entries in net.rc'
         write(*,*)nnuc,' nuclei of network in dimenfile'
         stop'net.rc or dimenfile needs to be reset'
      else
         write(*,*)'net.rc consistent with dimenfile: ibu=',ibu
      endif


c..   revise abundances...................................
c..   these use mass fractions
      if( ifm .eq. 1 )then
c..   solar scaling of z, he4 back to big bang
         call subifm1(rx,zscale,hetoz,ifm)
      elseif( ifm .eq. 2 )then
c..   opal abundances (fixed z if type 2)
         call chkopal(rx,zpop,xhe,xc12,xo16)
      elseif( ifm .eq. 0 )then
c..   default abundances from net.rc
         do j = 1,netsize
            rx(j) = solarx(j)
         enddo
      endif

c..   map revised abundances into model abundance arrays
c..   use mass fraction x(n,k)
      do k = 1, kk+2
         do j = 1, netsize
            x(j,k) = rx(j)
         enddo
      enddo
      write(*,*)netsize,' revised abundances scattered over ',
     1        kk+2,' zones'

c...........................................................

      if( ifm .ne. 0 )then
         sumx = 0.0d0
         sumz = 0.0d0
         do n = 1, nnuc
            sumx     = sumx     + rx(n)
            if( lz(n) .gt. 2 )then
               sumz     = sumz     + rx(n)
            endif
         enddo
         write(*,*)'renormalize: err =',sumx-1.0d0
c..   get Ye and Yion
         rx(ndim) = 0.0d0
         sum      = 0.0d0
         do n = 1, nnuc
            rx(ndim) = rx(ndim) + rx(n) * dble( lz(n) )/xa(n)
            sum      = sum      + rx(n)/xa(n)
         enddo
         yee   = rx(ndim)
         yions = sum
         ytot = yee + yions
         write(*,'(5(a6,1pe12.3))')'Ye',yee,'Yion',yions,'Ytot',ytot
         write(*,'(5(a6,1pe12.3))')'sumz',sumz,'H1',rx(nh1),
     1        'He4',rx(nhe4),'err',sumx-1.0d0
         write(*,'(5(a6,1pe12.3))')'c12',rx(lc12),'n14',rx(ln14),
     1        'o16',rx(lo16)
ccccccccccccccccccccc

c..   put mass fractions and ye into x(n,k) array for tycho and state
         do k = 2,kk+1
            do n = 1, ndim
               x(n,k) = rx(n)
            enddo
         enddo
      endif

c..   iterate for consistent temperature, keeping pressure
c..   and density fixed; hydrodynamic time scales are shorter
c..   than thermal times scales for stars, so this is reasonable
c      write(*,*)'Iterating T at constant P and V using state.f'
      do k = 1, kk+1
         scr(k) = p(2,k)
      enddo

      do n = 1, 30

         call state(2,kk+1,2)

         iflag = 0
         do k = 2, kk + 1
            delta  = ( scr(k) - p(2,k) )/pt(k)
            if( n .gt. 8 )then
c..   use bisection for poor convergence cases
               delta = 0.5d0*delta
            endif
            t(2,k) = t(2,k) + delta
            fact   = dabs( delta / t(2,k) )
            if( fact .gt. 1.0d-7 )then
               iflag = 1
            endif
         enddo
         if( iflag .eq. 0) goto 1010
      enddo

      write(*,*)'temperature iteration not converved in',n,' tries'
      stop'genex error 3'

 1010 continue
      write(*,*)'temperature iteration converged after', n,' tries'

      if( lumin .ne. 0 )then
         write(*,*)'reevaluating luminosity, lumin ',lumin
         do k = 1, kk+2
            tlold(k) = tl(2,k)
         enddo
         do k = 2, kk
            aki  = 0.5d0*( ak(k) + ak(k+1) )
            a(k) = 4.0d0 * pi * r(2,k)**2
            aom  = a(k) / dmi(k)
            f(k) = aom*arad*crad/3.0d0
     1           *( t(2,k)**4 - t(2,k+1)**4 )/aki
            tl(2,k) = a(k) * f(k)
         enddo

         tl(2,kk)   = tl(1,kk-1)
         tl(2,kk+1) = tl(1,kk)
         tl(2,kk+2) = tl(1,kk+1)
         do k = kk-5, kk+2
            write(*,'(i5,1p11e10.2)')k,tl(2,k),tlold(k),h(k),
     1           t(1,k),t(2,k),a(k),f(k),ak(k),
     2           x(ndim-1,k),x(ndim-2,k),
     3           1.0d0 -x(ndim-1,k)-x(ndim-2,k)
         enddo
         write(*,'(a5,11a10)')'k','L(new)','L(old)','vconv',
     1        'Tnew','Told','Area','flux','opac','Xhe','XH','z'
         write(*,*)'where kk =',kk
      endif

c..   find thermal time scale
      dtmin = 1.0d100
      do k = 2, kk
         s(1,k) = - ( tl(2,k) - tl(2,k-1))/dmh(k)
         if( s(1,k) .ne. 0.0d0 )then
            dtry = abs( e(2,k)/s(1,k) )
         else
            dtry = dtmin
         endif
         dtmin = dmin1( dtmin, dtry )
      enddo
      if( dtmin .ge. 1.0d100 )then 
         stop'GENEX: dt minimun error'
      else
         dth(1) = dtmin * 1.0d-2
      endif
      dth(2) = dth(1)
      write(*,'(a20,1pe12.3)')'new minimum timestep is ',dth(1)

      write(*,*)'new total mass ',xmsol+dmh(kk+1)/sol
      if(     ifm .eq. 0 )then
         write(*,*)'mass and radius scaling ',fm,fr
      elseif( ifm .eq. 1 )then
         write(*,*)'abundances: solar values scaled, zscale ',zscale
      elseif( ifm .eq. 2 )then
         write(*,'(a12,3(a6,1pe15.7))') 'abundances:',
     1        'H1',x(nnuc-1,kk),'He4',x(nnuc,kk),'zpop', zpop
         write(*,'(12x,4(a6,1pe15.7))') 'z/H',zpop/x(nnuc-1,kk),
     1        'c12',x(lc12,kk),'n14',x(ln14,kk),'o16',x(lo16,kk)
      elseif( ifm .eq. 3 )then
         write(*,*)'make helium star '
      else
         write(*,*)'ifm error',ifm
      endif

      write(*,'(a40,1p8e12.3)')'new central temperature(K)',t(1,2)

c..   maximum zone extent for rezoning
      ktot = kdm - 2
      iflag = 0
      ifm   = 0
      dlnv  = 1.0d-1
      xmmax = 1.0d-5
      xmmin = 1.0d-3
      dmmax = 1.0d-2
      drmax = 1.0d-1

c..   solid body rotation
      if( omegbar .gt. 0.0d0 )then
c..   critical angular velocity
         omegcrit = grav*xm(kk)/r(1,kk)**3
         if( omegcrit .gt. 0.0d0 )then
            omegcrit = sqrt( omegcrit )
         else
            write(*,*)' OMEGA CRITICAL ERROR ',omegcrit
            stop'genex omega critical'
         endif

         write(*,'(3(a20,1pe12.3))')'omega crit',omegcrit,
     1        'omega',omegbar,'omega/omcrit',omegbar/omegcrit

         if( omegbar .gt. omegcrit )then
            write(*,*)'spinning too fast'
            write(*,*)'omega must be less than omega crit'
            stop'genex spin error'
         endif
c..   set angular velocities
         do k = 1, kk+1
            omeg(k) = omegbar
         enddo
c..   set angylar momenta per unit mass
         ajaysum = 0.0d0
         do k = 1, kk+1
            ajay(k) = omeg(k)*r(1,k)**2
            ajaysum = ajaysum + ajay(k)*dmi(k)
         enddo
         write(*,'(3(a20,1pe12.3))')'total ang.mom.',ajaysum,
     1        'surf. speed(km/s)',omeg(kk+1)*r(1,kk+1)*1.0d-5,
     2        'free fall(km/s)',sqrt( 2.0d0*grav*xm(kk+1)/r(1,kk+1) )
     3        *1.0d-5
      endif

c..   reset note to flag new model, will be redefined in gen.f
      note(11:13) = 'new'
      write(*,*)note

      write(*,*)'using more complete format'
c..   adjusted model in nc=2 arrays
      open(8,file='new.model')
      model = 0
c      call nritef(note,2)
      call ritef(note,2)
      close(8)

      write(*,*)'leaving genex, new.model written '

      stop'successful termination'


 2000 continue
      write(*,*)'genex: error in  input file: genex.in'
      stop'input file error'

      end




      subroutine chkopal(rx,zpop,xhe,xxc,xxo)
      implicit none

c..reads opal abundances
c..writes summary
c..returns opal abundances (for all isotopes) in array rx
c..scale metallicity for type 1
c..force metallicity of OPAL table for type 2

      include 'dimenfile'
      include 'comod'
      include 'cburn'

c..   array for new abundances
      real*8    rx(ndim),zpop,xhe,xxc,xxo,xh1
      real*8    sumx,sumz,sumhe,sumh
      real*8    gndum(28),xael(28),xx(28)
      real*8    sumc,sumn,sumo

      integer*4 i,j

c..opal force variables
      real*8 xnum(28),xfra(28),xelem(28),varz
      integer*4 inum
      character*2 celem(28)
      character*11 ch11
      character*10 cdummy
      character*21 cdum21

c..   amu elemental mass from OPAL website where available
c..   this insures that the number densities used in OPAL opacities
c..   (and EOS?) are consistent with the values used in nuclear burn
      data xael/  1.00790d0, 4.00260d0, 6.941d0, 9.012d0, 10.811d0,
     1     12.0110d0,  14.0067d0,  15.9994d0,  18.99840d0, 20.1790d0,
     2     22.98977d0, 24.3050d0, 26.98154d0,  28.0855d0,  30.97376d0,
     3     32.0600d0,  35.453d0,  39.948d0,    39.0983d0,  40.080d0, 
     4     44.9559d0,  47.900d0,  50.9415d0,   51.9960d0,  54.9380d0, 
     5     55.847d0,   58.9332d0, 58.7000d0 /     

      logical tobe
c------------------------------------------------------------

c      write(*,'(//a20)')'entering chkopal'
 
      if( nopac .eq. 1 )then

         write(*,*)'OPAL type 1'
c..   check abundances of type 1 opacity tables
         inquire(file='GN93hz',exist=tobe) 
         if( .not. tobe )then
            write(*,*)'No file codataa found'
            stop'state.f GN93hz'
         endif
         open(20,file='GN93hz',status='old')
         do j = 1, 62
            if( j .ge. 32 .and. j .le. 59 )then
               inum = j-31
c..   inum is charge number Z
c..   celem is symbol for element
c..   xnum is fraction of metals by number
c..   xfra is fraction of metals by mass
               read(20,'(3x,a2,a11,f6.3,a21,0pf8.6,8x,0pf8.6)'),
     1              celem(inum),ch11,gndum(inum),
     2              cdum21,xnum(inum),xfra(inum)

               write(*,'(3x,a2,a11,f6.3,a21,0pf8.6,8x,0pf8.6)'),
     1              celem(inum),ch11,gndum(inum),
     2              cdum21,xnum(inum),xfra(inum)

               if( inum .eq. 3 )then
c..   use meteoritic value (12-12-06)
                  gndum(inum) = gndum(inum) + 2.15d0
                  write(*,*)' Updating ',celem(inum),
     1                 ' to meteoritic value log(A)= ',gndum(inum)
               endif
            else
               read(20,'(a10)')cdummy
            endif
         enddo

         read(20,'(14x,a12,25x,0pf6.5)')copal,zmetal
         write(*,'(a12,a14,a20,0pf6.5)')'OPAL file:',copal
         close(20)
         write(*,*)'Setting zmetal to ',zpop
         zmetal = zpop
c..   set H abundance (xhe and zmetal are defined)
         xh1 = 1.0d0 -xhe-zmetal
c..   GN93 version of solar mass fractions
c..   define elemental mass fraction of each element

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   fake bsp by reducing log10(c)
c         gndum(6) = gndum(6) - 0.04
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         sumz = 0.0d0
         do i = 3, 28
c..   gndum(i) are log10(Ni/NH) for log10 NH = 12.
c..   mass fractions (unnormalized) relative to H
c..   this preserves number ratios among metals
            xx(i) = 10.0d0 ** (gndum(i) - 12.0d0) *xael(i)
     1           *xh1/xael(1)
            sumz = sumz + xx(i)
         enddo

c..to have the cno as low as bsp do requires that either their
c..ratios to H are wrong, or their sum over metals is
c         write(*,'(8a12)')'sumz','c','n','o','fe','sum-zpop'
c         write(*,'(1p8e12.4)')sumz,xx(6),xx(7),xx(8),xx(26),
c     1        sumz-zpop
c         stop'a'
ccccccccccc

c..   normalize to zmetal imposed
         do i = 3, 28
            xx(i) = xx(i)/sumz * zmetal
         enddo
         sumz = 0
         do i = 3, 28
            sumz = sumz + xx(i)
            write(*,'(a10,1p8e12.3)')celem(i),xx(i)
ccccccccccccc
         enddo

c..   set H and He abundance (xhe and zmetal are defined)
         xh1 = 1.0d0 -xhe-zmetal
         xx(1) = xh1
         xx(2) = xhe
c..check normalization
         sumx = 0
         do i = 1, 28
c..   uses elemental atomic weights which is badly wrong for argon
c..   but this gives consistent number densities with OPAL
            sumx = sumx + xx(i)
         enddo
         write(*,*)'normalization error ',sumx-1.0d0
         if( abs(sumx - 1.0d0) .gt. 1.0d-14 )then
            write(*,*)'too large: in genex'
            stop'genex: chkopal'
         endif
c..write summary, with check of metal consistency
         sumz = 0
         do i = 3,28
            sumz = sumz + xx(i)
         enddo
         write(*,*)' H ', xx(1), ' He ',xx(2),
     1        ' metals ',sumz, ' diff z ',zmetal-sumz

      elseif( nopac .eq. 0 )then

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..check for mass fraction consistency
         write(*,*)'OPAL type 2'
c..   check abundances of type 2 opacity tables
         inquire(file='codataa',exist=tobe) 
         if( .not. tobe )then
            write(*,*)'No file codataa found'
            stop'state.f codataa'
         endif
         open(20,file='codataa',status='old')
c..read abundance data
         do j = 1, 62
            if( j .ge. 32 .and. j .le. 59 )then
               inum = j-31
c..   inum is charge number Z
c..   celem is symbol for element
c..   xnum is fraction of metals by number
c..   xfra is fraction of metals by mass
               read(20,'(3x,a2,a11,27x,0pf8.6,8x,0pf8.6)'),
     1              celem(inum),ch11,xnum(inum),xfra(inum)

               write(*,'(3x,a2,a11,27x,0pf8.6,8x,0pf8.6)'),
     1              celem(inum),ch11,xnum(inum),xfra(inum)
cccccccc
            else
               read(20,'(a10)')cdummy
            endif
         enddo
         read(20,'(14x,a12,25x,0pf6.5)')copal,zmetal
         write(*,'(a12,a14,a20,0pf6.5)')'OPAL file:',copal,
     1        'metallicity =',zmetal
cccccccc
         close(20)
         zpop = zmetal
         write(*,*)'Overwriting zpop<--zmetal ',zpop,
     1        ', the value in codataa'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
         write(*,*)'OPAL type error'
      endif



c..   fill working array
c..   scale the reference nuclear mass fractions to have the
c..   isotopic fractions of solarx but the elemental values of xx
c..   which has been taken from OPAL
      do i = 1, 28
         sumz = 0.0d0
         do j = 1, nnuc
c..   find the relevant isotopes of this element i 
c..   and sum their mass fractions
            if( lz(j) .eq. i )then
               sumz = sumz + solarx(j)
            endif
         enddo
c..   scale these nuclei to get the desired elemental abundance
         do j = 1, nnuc
            if( lz(j) .eq. i .and. sumz .gt. 0.0d0 )then
c..ratio of isotopic to elemental abundance * desired elemental abundance
               rx(j) = solarx(j)/sumz * xx(i)
            endif
         enddo
      enddo

c..   compare opal abundances with solar
c..   sum for elemental mass fractions xelem
         do inum = 1, 28
            xelem(inum) = 0.0d0
            do j = 1, nnuc
               if( inum .eq. lz(j) )then
                  xelem(inum) = xelem(inum) + solarx(j)
               endif
            enddo
         enddo
         write(*,*)'Opal abundances relative to solar (initial values)' 
         write(*,'(2(a5,a5,2a12,a15))')
     1        'Z','El','Xopal','Xsol','opal/sol',
     2        'Z','El','Xopal','Xsol','opal/sol'
         do inum = 1, 28, 2
            write(*,'(2(i5,a5,1p2e12.3,0pf15.10))')inum,celem(inum),
     1           xfra(inum)*zmetal,xelem(inum),
     2           xfra(inum)*zmetal/xelem(inum),
     3           inum+1,celem(inum+1),
     1           xfra(inum+1)*zmetal,xelem(inum+1),
     2           xfra(inum+1)*zmetal/xelem(inum+1)
         enddo

c..   compare opal abundances with solar
c..   sum for elemental mass fractions xelem
         do inum = 1, 28
            xelem(inum) = 0.0d0
            do j = 1, nnuc
               if( inum .eq. lz(j) )then
                  xelem(inum) = xelem(inum) + rx(j)
               endif
            enddo
         enddo

         sumx = 0
         do inum = 1, 28
            sumx = sumx + xelem(inum)
         enddo
         write(*,*)'check of xelements(OPAL) ',sumx-1.0d0
         sumz = 0
         do inum = 3, 28
            sumz = sumz + xelem(inum)
         enddo
         write(*,*)'revised metalicity: xelements(OPAL) ',sumz


         write(*,*)'Opal abundances relative to solar (revised values)' 
         write(*,'(2(a5,a5,2a12,a15))')
     1        'Z','El','Xopal','Xsol','opal/sol',
     2        'Z','El','Xopal','Xsol','opal/sol'
         do inum = 1, 28, 2
            write(*,'(2(i5,a5,1p2e12.3,0pf15.10))')inum,celem(inum),
     1           xfra(inum)*zmetal,xelem(inum),
     2           xfra(inum)*zmetal/xelem(inum),
     3           inum+1,celem(inum+1),
     1           xfra(inum+1)*zmetal,xelem(inum+1),
     2           xfra(inum+1)*zmetal/xelem(inum+1)
         enddo
c..   save opal values for use in defining excess C and O
         opalc =  xfra(6)*zmetal
         opaln =  xfra(7)*zmetal
         opalo =  xfra(8)*zmetal
         write(*,'(a10,4(a12,1pe12.4))')'OPAL','zmetal',zmetal,
     1        'opalC',opalc,'opalN',opaln,'opalO',opalo
c..   check nucleon conservation
         sumx = 0.0d0
c..   mass fraction of nuclei above helium (lithium and up)
         sumz = 0.0d0
         do j = 1, netsize
            sumx = sumx + rx(j)
            if( lz(j) .gt. 2 )then
               sumz = sumz + rx(j)
            endif
            if( lz(j) .eq. 6 )then
               sumc = sumc + rx(j)
            endif
            if( lz(j) .eq. 7 )then
               sumn = sumn + rx(j)
            endif
            if( lz(j) .eq. 8 )then
               sumo = sumo + rx(j)
            endif
         enddo
         write(*,'(a10,4(a12,1pe12.4))')'TYCHO','zsum(ty)',sumz,
     1        'tychoC',sumc,
     2        'tychoN',sumn,
     3        'tychoO',sumo
         write(*,'(a10,4(a12,1pe12.4))')'TYC-OP',
     1        'z',sumz-zmetal,
     1        'C',sumc-opalc,
     2        'N',sumn-opaln,
     3        'O',sumo-opalo
         if( abs(1.0d0-sumx) .gt. 1.0d-2 )then
            write(*,*)'state init error: sumx-1 ',sumx-1.0d0,kk
            stop'state sumx'
         endif

c..   metallicity from abundance variable (burn)
c..   test against OPAL labels
c..   metallicity .lt. 1.0d-10 is essentially zero
         varz = (zmetal-sumz) /
     1        (zmetal+sumz+1.0d-10) * 2.0d0
         if( abs( varz ) .gt. 1.0d-2 )then
            write(*,'(a26,1pe10.2,a14,1pe10.2,a8,1pe10.2,a4,i5)')
     1           'WARNING IN STATE: OPAL z =',zmetal,
     2           ' .ne. BURN z =', sumz,
     3           ' error =',sumz-zmetal,' kk=',kk
         else
            write(*,'(a48,1pe12.3)')
     1           'State: OPAL z agrees with actual metallicity to',
     2           sumz-zmetal
         endif

c..elemental H and He for scaling
         sumh = 0.0d0
         sumhe = 0.0d0
         do j = 1, nnuc
            if( lz(j) .eq. 1 )then
               sumh = sumh + rx(j)
            endif
            if( lz(j) .eq. 2 )then
               sumhe = sumhe + rx(j)
            endif
         enddo
            
         xh1 = 1.0d0 - sumz - xhe -xxc -xxo
         write(*,'(6(a5,1pe12.3))')'zpop',zpop,'xhe',xhe,
     1        'xc12',xxc,'xo16',xxo,'xh',xh1

         do j = 1, nnuc
            if( lz(j) .eq. 1 )then
c..scale H isotopes
               rx(j) = rx(j)*xh1/sumh
            endif
            if( lz(j) .eq. 2 )then
c..scale He isotopes
               rx(j) = rx(j)*xhe/sumhe
            endif
         enddo
         sumz = 0.0d0
         sumx = 0.0d0
         sumhe = 0.0d0
         sumh = 0.0d0
         do j = 1, nnuc
            sumx = sumx + rx(j)
            if( lz(j) .gt. 2 )then
               sumz = sumz + rx(j)
            elseif( lz(j) .eq. 2 )then
               sumhe = sumhe + rx(j)
            elseif( lz(j) .eq. 1 )then
               sumh = sumh + rx(j)
            endif
         enddo

         write(*,'(6(a7,1pe13.5))')'sumz',sumz,
     1        'xhe',xhe,'sumhe',sumhe,'sumh',sumh
         write(*,'(6(a7,1pe13.5))')'sumx-1',sumx-1.0d0,
     1        'he-she',xhe-sumhe,'z-sz',zmetal-sumz
         write(*,'(6(a7,1pe13.5))')'C12',rx(lc12),'N14',rx(ln14),
     1        'O16',rx(lo16)

         write(*,*)'leaving chkopal'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end






      subroutine subifm1(rx,zscale,hetoz,ifm)
      implicit none
c..   scales abundances from Big Bang through solar to super-metal-rich

      include 'dimenfile'
      include 'comod'
      include 'cburn'

c..   array for new abundances
      real*8    rx(ndim),hetoz,sumx,sumz,solscale,zsol

      integer*4 ifm,nh1,nhe4,n,j,i

c..   big bang abundances
      integer*4 nucbb
      parameter( nucbb = 13 )
      integer*4 ibb(nucbb),ibbz(nucbb),ibbn(nucbb)
      real*8    xbb(nucbb), fint, omfint
      character*5 xidbb(nucbb)

c..   scaling of metallicity
      real*8 fact1,sum1,zscale,bz

c..   big bang abundances by mass and nuclear z and n
c..   h3 decayed to he3, c14 decayed to n14
      data ibbz/ 1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 1, 2 /
      data ibbn/ 1, 1, 3, 4, 3, 5, 5, 6, 6, 7, 7, 0, 2 /
      data xidbb/'d','he3','li6','li7',
     1     'be7','be9','b10','b11','c12',
     2     'c13','n14','p','he4'    /
      data xbb/
     1     1.0d-4,  5.0d-5,  1.0d-13, 5.0d-10,
     2     5.0d-10, 1.0d-16, 2.0d-18, 2.0d-16, 1.0d-14,
     3     1.0d-15, 1.2d-16, 0.76d0,  2.4d-1   /
      save
c-------------------------------------------------------------
      write(*,'(/a25,i5,a10,i5)')'entering subifm1, ifm= ',ifm,
     1     ' nopac= ',nopac
      nh1  = nnuc-1
      nhe4 = nnuc
c..   define new abundances (rx = nucleon fraction)

c..   solar scaling
      if( zscale .lt. 0.0d0 )then
         write(*,*)'genex error: zscale',zscale,' ifm ',ifm
         stop'subifm1 genex'
      endif
c..   solar abundances as defined in net.rc
c..   fill working array
      do j = 1, nnuc
         rx(j) = solarx(j)
      enddo
c..   renormalize to machine accuracy
      sumx = 0.0d0
      do j = 1, nnuc
         sumx = sumx + rx(j)
      enddo
      write(*,*)'initial normalization error ', sumx-1.0d0
      do j = 1, nnuc
         rx(j) = rx(j)/sumx
      enddo
      sumx = 0.0d0
      sumz = 0.0d0
      do j = 1, nnuc
         sumx = sumx + rx(j)
         if( lz(j) .gt. 2 )then
            sumz = sumz + rx(j)
         endif
      enddo
      write(*,*)'renormalization error ',sumx-1.0d0

c..   adjust heavy elements (not H1 and He4)
c..   what factor will keep X(H1) and X(He4) fixed?
      bz = 1.0d0 - rx(nnuc) - rx(nnuc-1)
      sum1 = 0.0d0
      do n = 1, nnuc-2
         sum1 = sum1 + rx(n)
      enddo
      solscale = bz/sum1
      do n = 1, nnuc-2
         rx(n) = rx(n)*solscale
      enddo
      fact1 = -1.0d0
      do n = 1, nnuc
         fact1 = fact1 + rx(n)
      enddo
      if( abs(fact1) .gt. 1.0d-8 )then
         write(*,*)'bad normalization ',fact1
         stop'rscale error'
      else
         write(*,*)'rscale normalization ',fact1
      endif

      if( zscale .le. 0.02d0 )then
c..   zscale=0.02 is extreme pop II 
c..   interpolate back to big bang abundances
c..   normalize big bang abundances for precision
         sumx = 0.0d0
         do i = 1, nucbb
            sumx = sumx + xbb(i)
         enddo
         do i = 1, nucbb
            xbb(i) = xbb(i)/sumx
         enddo
         sumx = 0.0d0
         do i = 1, nucbb
            sumx = sumx + xbb(i)
         enddo
c..   identify big bang nuclei
         do i = 1, nucbb
            ibb(i) = 0
         enddo
         do n = 1, nnuc
            do i = 1, nucbb
               if( ibbz(i) .eq. lz(n) .and. ibbn(i) .eq. ln(n) )then
                  ibb(i) = n
c     xidbb(i) = cnuc(n)
               endif
            enddo
         enddo
         do i = 1, nucbb
            if( ibb(i) .eq. 0 )then
               write(*,'(35x,a5,i5,a5,i3,a6,i3,a15)')xidbb(i),
     1              i,' Z =',ibbz(i),' N =',ibbn(i),' not in net.rc'
            endif
         enddo
c..   replace with interpolated value, linear interpolation in mass
         fint   = zscale
         fint   = dmax1(fint,0.0d0)
         fint   = dmin1(fint,1.0d0)
         omfint = 1.0d0 - fint
c..   scaled solar fraction
         do j = 1, nnuc
            rx(j) = fint*rx(j)
         enddo
         if( omfint .gt. 0.0d0 )then
c..   add rest as big bang
            do i = 1, nucbb
               if( ibb(i) .gt. 0 )then
                  j = ibb(i)
                  rx(j) = rx(j) + omfint*xbb(i)
               endif
            enddo
         endif

c..   test normalization and shift He for nucleosynthesis in stars
         sumx = 0.0d0
         do n = 1, nnuc
            sumx = sumx + rx(n)
         enddo
         write(*,'(a20,1p2e12.3)')" normalization (*) ",
     1        sumx,sumx-1.0d0

c..   test normalization and shift He for nucleosynthesis in stars
         sumx = 0.0d0
         do i = 1, nucbb
            sumx = sumx + xbb(i)
         enddo
         write(*,'(a20,1p2e12.3)')" normalization (bb)",
     1        sumx,sumx-1.0d0

         do n = 1, nnuc
            rx(n) = rx(n)/sumx
         enddo
c..   test normalization and shift He for nucleosynthesis in stars
         sumx = 0.0d0
         do n = 1, nnuc
            sumx = sumx + rx(n)
         enddo
         write(*,'(a20,1p2e12.3)')" normalization ",sumx,sumx-1.0d0
c..   add in extra hydrogen (dHe/dz = 2) for lower metallicity
         write(*,*)'reset abundances:'
         write(*,'(3(3i5,a5,1pe12.3))')
     1        (n,lz(n),ln(n),cnuc(n),rx(n), n = 1,nnuc)
         write(*,'( 3(a10,1pe12.3) )')
     1        ' he4',rx(nnuc),' h1',rx(nnuc-1),
     2        ' other',1.0d0-rx(nnuc)-rx(nnuc-1)

      else

c..   pop II to pop I to super-metal-rich(SMR); scale solar values
         write(*,*)'Super-metal-rich: ifm = ',ifm,' zscale ',zscale
         zsol = 1.0d0-rx(nh1)-rx(nhe4)

         write(*,'(a20,1pe12.3)')'solar metallicity',
     1        zsol
         write(*,'(a20,1pe12.3)')'solar He4',
     1        rx(nhe4)
         write(*,'(a20,1pe12.3)')'solar H1',
     1        rx(nh1)
         write(*,'(a20,1pe12.3)')'He/z production',
     1        hetoz
         write(*,'(a20,1pe12.3)')'new metallicity',
     1        zsol*zscale

c..   revise H1, He4 with production of He with metals
         rx(nh1)  = solarx(nh1) - zsol*(zscale-1.0d0)*(1.0d0 + hetoz)
         rx(nhe4) = solarx(nhe4) + zsol*(zscale-1.0d0)*hetoz

         write(*,'(a20,1pe12.3)')'new He4',
     1        rx(nhe4)
         write(*,'(a20,1pe12.3)')'new H1',
     1        rx(nh1)

c..   scale up Li and higher, not D and He3
         do n = 1, nh1-1
            if( lz(n) .gt. 2 )then
               rx(n) = rx(n)*zscale
            endif
         enddo

c..   test normalization and shift He for nucleosynthesis in stars
         sumx = 0.0d0
         do n = 1, nnuc
            sumx = sumx + rx(n)
         enddo
         write(*,'(a20,1p2e15.7)')" normalization ",sumx,
     1        sumx-1.0d0
         do n = 1, nnuc
            rx(n) = rx(n)/sumx
         enddo
c..   test normalization and shift He for nucleosynthesis in stars
         sumx = 0.0d0
         do n = 1, nnuc
            sumx = sumx + rx(n)
         enddo
         write(*,'(a20,1p2e15.7)')" renormalization ",sumx,
     1        sumx-1.0d0

         write(*,*)'reset abundances:'
         write(*,'(3(3i5,a5,1pe12.3))')
     1        (n,lz(n),ln(n),cnuc(n),rx(n), n = 1,nnuc)
         write(*,'( 4(a11,1pe12.3) )')
     1        ' zed',1.0d0-rx(nnuc)-rx(nnuc-1),' he4',rx(nnuc),
     2        ' hydrogen',rx(nnuc-1),'Prod:He/z',hetoz
      endif

      return
      end
