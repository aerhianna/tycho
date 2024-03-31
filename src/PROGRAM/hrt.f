      program hrt

c..   reads hr.? file from standard input and
c..   uses pgplot to plot 
c..   modified for extended network
c..   detailed abundances
c..   readf uses unit 5

      implicit none

      include "dimenfile"
      include "cconst"

      integer*4 icolor

      integer*4 l, ll, m, ixvar, iyvar, lstyle, n, nfiles, numfiles
      integer*4 pgbeg

      integer*4 ndm
      parameter( ndm = 50000, numfiles = 40 )

      integer*4 model1(ndm), it(ndm)
      real*8 xej(ndim),xxej(5,ndm)
      real*8 time(ndm),dth(ndm),rphot(ndm),rhoph(ndm),uphot(ndm),
     1     xlol(ndm),telog(ndm),xmkk(ndm),
     2     peryear(ndm),speryear(ndm),omegkk(ndm),vinf(ndm),
     3     omegkk1(ndm),dmhkk1(ndm)
      real*8 smass
      real*8 feh(ndm), glog(ndm), inner(7,ndm), outer(4,ndm)
      real*8 teff, Mbol
      real*8 color(5)
      real*8 V(ndm), BV(ndm)
      real*8 bound(5,ndm), unlogl(ndm), seff(5,ndm),chi2

      real*8 Tnot(ndm),lx(ndm),omegadot(ndm),
     1        Br_star(ndm),Mdot_star(ndm),rossby(ndm),
     2        pnot(ndm),bo_omega(ndm),bo_omega2(ndm),
     3        tauc(ndm),taug(ndm),vconv(ndm),
     4        xbosurfomeg(ndm)

      integer*4 i,j
      integer*4 lw
      real*4 fsize

      real*4 xx(ndm),yy(ndm)
      real*4 xmin,xmax,ymin,ymax,xtest,ytest
      real*4 xmin0,xmax0,ymin0,ymax0

      character*5 hrfile(numfiles)
      character*40 toplab
      character*5 device

      integer*4 nuclei,nucleus(5)


      integer*4 idum
      real*4    aa(14)
      real*4    fehhyad

      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt

      logical tobe

      data l/1/
      data lw/3/,fsize/1.2/
c..   ag88 value relative to opal normalization = 0.02
      data fehhyad/0.0/
      data txtxt/'  '/
c--------------------------------------------------------------------
c..   read graphics options
      inquire(file='hrt.in',exist=tobe)
      if( tobe )then
         open(9,file='hrt.in')
      else
         stop'hrt: no hrt.in input file error'
      endif

c..   dummy reads for labelling header lines
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text

      read (9,*)txt,labl,device
      write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000


      read (9,*)txt,labl,nfiles
      write(*,*)txt,txtxt,labl,nfiles
      labl1 = 'nfiles'
      if( labl .ne. labl1 )goto 2000


      if( nfiles .gt. numfiles )then
         write(*,*)'too many file:',nfiles,numfiles
         stop'hrplot'
      endif

      do n = 1 , nfiles
         read (9,*)txt,labl,hrfile(n)
         write(*,*)txt,txtxt,labl,hrfile(n)
         labl1 = 'hrfile'
         if( labl .ne. labl1 )goto 2000
      enddo

      write(*,*)n-1,' files have been read'

      read (9,*)txt,labl,ixvar
      write(*,*)txt,txtxt,labl,ixvar
      labl1 = 'ixvar'
      if( labl .ne. labl1 )goto 2000

      if( ixvar .eq. 0 )then
         write(*,*)'x axis is time in seconds'
      elseif( ixvar .eq. 1 )then
         write(*,*)'x axis is model number'
      elseif( ixvar .eq. 2 )then
         write(*,*)'x axis is log10( Teff(K) )'
      elseif( ixvar .eq. 3 )then
         write(*,*)'x axis is time in years'
      elseif( ixvar .eq. 4 )then
         write(*,*)'x axis is B-V color'
      elseif( ixvar .eq. 5 )then
         write(*,*)'x axis is time in log10(year)'
      else
         write(*,*)'error in hrt.in, ixvar = ',ixvar
         stop
      endif

      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text

      read (9,*)txt,labl,iyvar
      write(*,*)txt,txtxt,labl,iyvar
      labl1 = 'iyvar'
      if( labl .ne. labl1 )goto 2000

      if( iyvar .eq. 0 )then
         write(*,*)'y axis is log10( L/Lsun )'
      elseif( iyvar .eq. 1 )then
         write(*,*)'y axis is -mass loss in gm/sec'
      elseif( iyvar .eq. 2 )then
         write(*,*)'y axis is log10( Teff(K) )'
      elseif( iyvar .eq. 3 )then
         write(*,*)'y axis is mass loss in Msun/year'
      elseif( iyvar .eq. 4 )then
         write(*,*)'y axis is abundance Y(i) loss'
      elseif( iyvar .eq. 5 )then
         write(*,*)'y axis is V magnitude'
      elseif( iyvar .eq. 6 )then
         write(*,*)'y axis is number of iterations'
      elseif( iyvar .eq. 7 )then
         write(*,*)'y axis is mass'
      elseif( iyvar .eq. 8 )then
         write(*,*)'y axis is Radius'
      elseif( iyvar .eq. 9 )then
         write(*,*)'y axis is log timestep'
      elseif( iyvar .eq. 10 )then
         write(*,*)'y axis is log time'
      elseif( iyvar .eq. 11 )then
         write(*,*)'y axis is habitable zone radius'
      elseif( iyvar .eq. 12 )then
         write(*,*)'y axis is log10(coronal temperature)'
      elseif( iyvar .eq. 13 )then
         write(*,*)'y axis is l_x [L_x/L_x_sun]'
      elseif( iyvar .eq. 14 )then
         write(*,*)'y axis is omega_dot [d_omega/d_t]'
      elseif( iyvar .eq. 15 )then
         write(*,*)'y axis is Br_star [radial magnetic field]'
      elseif( iyvar .eq. 16 )then
         write(*,*)'y axis is -Mdot_star [d_M/d_t]'
      elseif( iyvar .eq. 17 )then
         write(*,*)'y axis is Rossby number'
      elseif( iyvar .eq. 18 )then
         write(*,*)'y axis is pnot [pressure at corona base]'
      elseif( iyvar .eq. 19 )then
         write(*,*)'y axis is omega(kk)'
      elseif( iyvar .eq. 20 )then
         write(*,*)'y axis is omega(kk+1)'
      elseif( iyvar .eq. 21 )then
         write(*,*)'y axis is B.O. omega(kk)'
      elseif( iyvar .eq. 22 )then
         write(*,*)'y axis is B.O. omega(kk+1)'
      elseif( iyvar .eq. 23 )then
         write(*,*)'y axis is tau_c'
      elseif( iyvar .eq. 24 )then
         write(*,*)'y axis is tau_g'
      elseif( iyvar .eq. 25 )then
         write(*,*)'y axis is v_conv'
      else
         write(*,*)'error in hrt.in, iyvar = ',iyvar
         stop
      endif
      if( ixvar .eq. 4 .and. iyvar .ne. 5 )then
         write(*,*)'at present ixvar=4 requires iyvar=5'
         stop'input CMD error in y'
      endif
      if( iyvar .eq. 5 .and. ixvar .ne. 4 )then
         write(*,*)'at present iyvar=5 requires ixvar=4'
         stop'input CMD error in x'
      endif

      read (9,*)text
      write(*,*)text

c      read (9,*)text
c      write(*,*)text
c      read (9,*)text
c      write(*,*)text
c      read (9,*)text
c      write(*,*)text
c      read (9,*)text
c      write(*,*)text

      read(9,*)nuclei

      if( nuclei .eq. 0 )then
         write(*,*)'no nuclei used; nuclei =',nuclei
      elseif( nuclei .gt. 0 .and. nuclei .le. 5 )then
         read(9,*)(nucleus(i),i=1,nuclei)
         write(*,'(a35,5i5)')'number of nuclei used here:',nuclei
         write(*,'(a35,5i5)')'indices of nuclei used here are:',
     1        (nucleus(i),i=1,nuclei)
      else
         write(*,*)'nuclei ',nuclei
         stop'hrt: error in nuclei'
      endif

      read (9,*)text
      write(*,*)text

      read (9,*)txt,labl,xmin0
      write(*,*)txt,txtxt,labl,xmin0
      labl1 = 'xmin'
      if( labl .ne. labl1 )goto 2000

      read (9,*)txt,labl,xmax0
      write(*,*)txt,txtxt,labl,xmax0
      labl1 = 'xmax'
      if( labl .ne. labl1 )goto 2000

      read (9,*)txt,labl,ymin0
      write(*,*)txt,txtxt,labl,ymin0
      labl1 = 'ymin'
      if( labl .ne. labl1 )goto 2000

      read (9,*)txt,labl,ymax0
      write(*,*)txt,txtxt,labl,ymax0
      labl1 = 'ymax'
      if( labl .ne. labl1 )goto 2000

      read (9,*)txt,labl,lstyle
      write(*,*)txt,txtxt,labl,lstyle
      labl1 = 'lstyle'
      if( labl .ne. labl1 )goto 2000

      if( lstyle .eq. 0 )then
         write(*,*)'line style is 0 (dots)'
      elseif(  lstyle .eq. 1 )then
         write(*,*)'line style is 1 (line)'
      elseif(  lstyle .eq. 2 )then
         write(*,*)'line style is 2 (line+dots)'
      else
         write(*,*)'hrt.in: lstyle error ',lstyle
         write(*,*)'0 = dots, 1 = line, 2 = both'
         stop
      endif

      read (9,*)text
      write(*,*)text
      close(9)

      if( nfiles .le. 5 )then
         toplab = hrfile(1)
         do n = 2, nfiles
            toplab = hrfile(n)//' '//toplab
         enddo
      endif
      toplab = ' '
      write(*,*)'top label is:'
      write(*,*)toplab

c..   THE FOLLOWING VARIABLES WERE WRITTEN IN EDIT.F IN TYCHO
c..   model1   corresponding model number
c..   it       iterations used 
c..   time     seconds elapsed in whole sequence
c..   dth(2)   time step
c..   rphot    photospheric radius
c..   rhoph    photospheric density (g/cc)
c..   uphot    photospheric velocity (cm/s)
c..   xlol     log( L / Lsun )              log => base 10
c..   telog    log( Teff )                  ln  => base e
c..   xm(kk)   mass (gm) interior to outer fitting radius r(1,kk)
c..   peryear  mass loss(-) or gain(+) in solar masses per year
c..   speryear solar masses ejected in this time step
c..   omeg(kk) angular rotational velocity (1/s)
c..   xej      ejected nucleon fraction array of size ndim

c..   read first model to get limits of data
c..   hr.?? data file
      inquire(file=hrfile(1),exist=tobe)
      if( tobe )then
         open(8,file=hrfile(1))


 100     continue
         if( l .gt. ndm )goto 101

         read(8,70,end=101) model1(l), it(l), time(l), dth(l),
     1        rphot(l), rhoph(l), uphot(l), xlol(l), telog(l), 
     2        xmkk(l), peryear(l), vinf(l), omegkk(l),
     3        dmhkk1(l),omegkk1(l)
c         chi2 = ((rphot(l) - 1.19d11)/6.959d10)**2.0d0
c         chi2 = ((rphot(l) - 1.19d11)/(0.013*6.959d10))**2.0d0 + 
c     1          ((xlol(l) - 1.405d0)/0.022)**2.0d0
c         chi2 = chi2/6.53d-3
c	 write(*,*)chi2,time(l),rphot(l),xlol(l),telog(l)
         read(8,71) xej
         read(8,73) Tnot(l),lx(l),omegadot(l),Br_star(l),
     1            Mdot_star(l),rossby(l),pnot(l),bo_omega(l),
     2           bo_omega2(l),tauc(l),taug(l),vconv(l),
     3           xbosurfomeg(l)
         Mdot_star(l) = -Mdot_star(l)/sol*secpy
         if( iyvar .eq. 4 )then
            do j = 1, nuclei
               xxej(j,l) = xej(nucleus(j))
            enddo
         endif
         l = l+1
         goto 100
 101     continue
         ll = l-1
c..   echo last data read in
         write(*,*)'last data read in is:'
         write(*,'(2i6,i3,1p11e9.1)')ll, model1(ll), it(ll), time(ll), 
     1        dth(ll),
     1        rphot(ll), rhoph(ll), uphot(ll), xlol(ll), telog(ll), 
     2        xmkk(ll), peryear(ll), speryear(ll), omegkk(ll)
         if( iyvar .eq. 4 )then
            write(*,'(2a5,a12)')'j','nuc','Y(j)'
            do j = 1, nuclei
               write(*,'(2i5,1pe12.3)') j,nucleus(j), xxej(j,ll)
            enddo
         endif
      else
         write(*,'(a10,a5)')'no file ',hrfile(1)
         stop'hrt: input error'
      endif

      close(8)

      if( ixvar .eq. 0 )then
         xmin = time(1)
         xmax = time(ll)
      elseif( ixvar .eq. 1 )then
         xmin = model1(1)
         xmax = model1(ll)
      elseif( ixvar .eq. 2 )then
         xmin = telog(1)
         xmax = xmin
c..   invert to get astronomers log Te for HR diagram
         do l = 1, ll
            xmin = amax1(xmin,sngl(telog(l)))
            xmax = amin1(xmax,sngl(telog(l)))
         enddo
      elseif( ixvar .eq. 3 )then
         xmin = time(1)  /secpy
         xmax = time(ll) /secpy        
      elseif( ixvar .eq. 4 )then
         do m = 1, ll
            glog(m) = log10( grav*xmkk(m)/rphot(m)**2 )
         enddo  

         do m = 1, ll
            feh(m) = fehhyad
         enddo

         inquire(file='kurcoljkl.tbl',exist=tobe)
         if( tobe )then
c..   subroutine kurintpt.f opens file kurcoljkl.tbl on unit 30
         else
            stop'hrt: no kurcoljkl.tbl input file error'
         endif

         do m = 1, ll
            teff = 10**telog(m)
            call kurintpt(feh(m), teff, glog(m), color)
            BV(m) = color(3)
            Mbol = 4.75 - 2.5 * xlol(m)
            V(m) = Mbol - color(1)
         enddo
         xmin = BV(1)
         xmax = BV(ll)
      elseif( ixvar .eq. 5 )then
         xmin = log10(time(1)/secpy)
         xmax = log10(time(ll)/secpy)
                  
      endif

      if( iyvar .eq. 0 )then
         ymin = xlol(1)
      elseif( iyvar .eq. 1 )then
         ymin = -peryear(1)
      elseif( iyvar .eq. 2 )then
         ymin = telog(1)
      elseif( iyvar .eq. 3 )then
         ymin = 0.0d0
      elseif( iyvar .eq. 4 )then
         ymin = xxej(1,1)
         do l = 1, ll
            do j = 1, nuclei
               ymin = dmin1(xxej(j,l),dble(ymin))
            enddo
         enddo
      elseif( iyvar .eq. 5 )then
         ymin = V(1)
      elseif( iyvar .eq. 8)then
         ymin = rphot(1)
         do l = 1, ll
            ymin = dmin1(rphot(l),dble(ymin))
         enddo
      elseif( iyvar .eq. 9 )then
         ymin = dth(1) 
         do l = 1, ll
            ymin = dmin1(dth(l),dble(ymin))
         enddo
      elseif( iyvar .eq. 10 )then
         ymin = time(1)
      elseif( iyvar .eq. 11 )then
         ymin = 0.0e0
      elseif( iyvar .eq. 12 )then
         ymin = LOG10(Tnot(1))
      elseif( iyvar .eq. 13 )then
         ymin = lx(1)
      elseif( iyvar .eq. 14 )then
         ymin = omegadot(1)
      elseif( iyvar .eq. 15 )then
         ymin = Br_star(1)
      elseif( iyvar .eq. 16 )then
         ymin = Mdot_star(1)
      elseif( iyvar .eq. 17 )then
         ymin = rossby(1)
      elseif( iyvar .eq. 18 )then
         ymin = pnot(1)
      elseif( iyvar .eq. 19 )then
         ymin = omegkk(1)
      elseif( iyvar .eq. 20 )then
         ymin = omegkk1(1)
      elseif( iyvar .eq. 21 )then
         ymin = bo_omega2(1)
      elseif( iyvar .eq. 22 )then
         ymin = bo_omega(1)
      elseif( iyvar .eq. 23 )then
         ymin = tauc(1)
      elseif( iyvar .eq. 24 )then
         ymin = taug(1)
      elseif( iyvar .eq. 25 )then
         ymin = vconv(1)
      endif
      ymax = ymin

      if( ixvar .eq. 4 )then
         do l = 2, ll
            xtest = BV(l)
            xmin = amin1(xmin,xtest )
            xmax = amax1(xmax,xtest )
         enddo
      endif

      if( iyvar .eq. 0 )then
         do l = 2, ll
            ytest = xlol(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 1 )then
         do l = 2, ll
            ytest = -peryear(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 2 )then
         do l = 2, ll
            ytest = telog(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 3 )then
         do l = 2, ll
            ytest = speryear(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
         do l = 1, ll
            smass = smass + speryear(l)
         enddo
      elseif( iyvar .eq. 5 )then
         do l = 2, ll
c..   astronomers inverted scale; negative is brighter
            ytest = V(l)
            ymin = amax1(ymin,ytest )
            ymax = amin1(ymax,ytest )
            
         enddo
      elseif( iyvar .eq. 8 )then
         do l = 2, ll
            ymax = dmax1( rphot(l), dble(ymax) )
         enddo
      elseif( iyvar .eq. 9 )then
         do l = 2, ll
            ymax = dmax1( dth(l), dble(ymax) )
         enddo
c..   set up for log10
         ymin = alog10( ymin )
         ymax = alog10( ymax )
      elseif( iyvar .eq. 4 )then
         do l = 1, ll
            do j = 1, nuclei
               ymax = dmax1( xxej(j,l), dble(ymax) )
            enddo
         enddo
c..   set up for log10 
         if( ymin .gt. 0.0d0 )then
            ymin = amax1( alog10(ymin), -12.0 )
         else
            ymin = -12.0
         endif
         if( ymax .gt. 0.0d0 )then
            ymax = amin1( alog10(ymax), 0.0 )
         else
            ymax = 0.0
         endif
      elseif( iyvar .eq. 6 )then
         ymin = 0
         ymax = ymin
         do l = 1, ll
            ymax = amax1( float( it(l) ), ymax )
         enddo
      elseif( iyvar .eq. 7 )then
         ymin = 0
         ymax = ymin
         do l = 1, ll
            ymax = amax1( sngl( xmkk(l)/sol ), ymax )
         enddo
      elseif( iyvar .eq. 10 )then
         ymax = time(ll)
c..   set up for log10 
         if( ymin .gt. 0.0d0 )then
            ymin = amax1( alog10(ymin), -12.0 )
         else
            ymin = -12.0
         endif
c..   set up for log10 
         if( ymax .gt. 0.0d0 )then
            ymax = amax1( alog10(ymax), -12.0 )
         else
            ymax = -12.0
         endif
      elseif( iyvar .eq. 11 )then
         ymax = ymin
         do l = 2, ll
            ytest = (10.0e0**(xlol(l))*3.9e33/(16.0e0* pi * sigma * 
     1              5.555e9))**0.5d0/1.5e13
            ymax = amax1(ymax, ytest)
         enddo
      elseif( iyvar .eq. 12 )then
         do l = 2, ll
            ytest = LOG10(Tnot(l))
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 13 )then
         do l = 2, ll
            ytest = lx(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 14 )then
         do l = 2, ll
            ytest = omegadot(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 15 )then
         do l = 2, ll
            ytest = Br_star(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 16 )then
         do l = 2, ll
            ytest = Mdot_star(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 17 )then
         do l = 2, ll
            ytest = rossby(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 18 )then
         do l = 2, ll
            ytest = pnot(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 19 )then
         do l = 2, ll
            ytest = omegkk(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 20 )then
         do l = 2, ll
            ytest = omegkk1(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 21 )then
         do l = 2, ll
            ytest = bo_omega2(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 22 )then
         do l = 2, ll
            ytest = bo_omega(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 23 )then
         do l = 2, ll
            ytest = tauc(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 24 )then
         do l = 2, ll
            ytest = taug(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      elseif( iyvar .eq. 25 )then
         do l = 2, ll
            ytest = vconv(l)
            ymin = amin1(ymin,ytest )
            ymax = amax1(ymax,ytest )
         enddo
      endif

c..   allow input values in hrt.in to override if they are
c..   more restrictive
      if( xmin0 .ne. xmax0 )then
         if( xmax0 .lt. xmin .and. ixvar .ne. 2 )then
c..   ixvar is log Te, astronomers invert axis for HR diagram
            write(*,*)'hrt warning: xmax0 < xmin, reset your xmax0?'
         elseif( xmin0 .gt. xmax )then
            write(*,*)'hrt warning: xmin0 > xmin, reset your xmin0?'
         endif
      endif
      if( ymin0 .ne. ymax0 )then
         if( ymax0 .lt. ymin .and. iyvar .ne. 5 )then
            write(*,*)'hrt warning: ymax0 < ymin, reset your ymax0?'
         elseif( ymin0 .gt. ymax )then
            write(*,*)'hrt warning: ymin0 > ymax, reset your ymin0?'
         endif
      endif

c..   force use of read-in limits
      if( xmin0 .ne. xmax0 )then
         xmin = xmin0
         xmax = xmax0
      endif
      if( ymin0 .ne. ymax0 )then
         ymin = ymin0
         ymax = ymax0
      endif

c..   add border region to keep graph off axes
      xtest = xmax - xmin
      ytest = ymax - ymin
      xmin  = xmin -0.05*xtest
      xmax  = xmax +0.05*xtest
      ymin  = ymin -0.05*ytest
      ymax  = ymax +0.05*ytest

c..   initialize graphics...................................

c..   scratch file for character conversion
      open(10,file='scratch.pg')

      if ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
      call PGASK (.FALSE.)
c..   roman font
      call PGSCF(2)
c..   font scaling (size)
      call PGSCH(fsize)
c..   new page
      call PGPAGE

c..   window
      call PGSWIN(xmin,xmax,ymin,ymax)
c..   reset color
      call PGSCI(1)
c..   line width
      call PGSLW(lw)
c..   tickmarks on x(bottom=B,top=C) and y(left)
      call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

c..   Left label, dispacement outside in character heights,
c..   location in fractions of edge, justification(0.5=centered)
c..   and character string

      if( iyvar .eq. 0 )then
         call PGMTXT('L',2.0,0.5,0.5,'log L/L(sol)')
      elseif( iyvar .eq. 1 )then
         call PGMTXT('L',2.0,0.5,0.5,'M(sol)/year')
      elseif( iyvar .eq. 2 )then
         call PGMTXT('L',2.0,0.5,0.5,'log Teff(K)')
      elseif( iyvar .eq. 3 )then
         call PGMTXT('L',2.0,0.5,0.5,' ')
      elseif( iyvar .eq. 4 )then
         call PGMTXT('L',2.0,0.5,0.5,'log Y(i)')
      elseif( iyvar .eq. 5 )then
         call PGMTXT('L',2.0,0.5,0.5,'V')
      elseif( iyvar .eq. 6 )then
         call PGMTXT('L',2.0,0.5,0.5,'it')
      elseif( iyvar .eq. 7 )then
         call PGMTXT('L',2.0,0.5,0.5,'M/M(sun)')
      elseif( iyvar .eq. 8 )then
         call PGMTXT('L',2.0,0.5,0.5,'R(cm)')
      elseif( iyvar .eq. 9 )then
         call PGMTXT('L',2.0,0.5,0.5,'log dth(sec)')
      elseif( iyvar .eq.10 )then
         call PGMTXT('L',2.0,0.5,0.5,'log time(s)')
      elseif( iyvar .eq.11 )then
         call PGMTXT('L',2.0,0.5,0.5,'Habitable zone radius (AU)')
      elseif( iyvar .eq.12 )then
         call PGMTXT('L',2.0,0.5,0.5,'log Coronal Temperature (K)')
      elseif( iyvar .eq.13 )then
         call PGMTXT('L',2.0,0.5,0.5,'Converged l_x [L_x/L_x_sun]')
      elseif( iyvar .eq.14 )then
         call PGMTXT('L',2.0,0.5,0.5,'d_omega/d_t (rad/sec^2)')
      elseif( iyvar .eq.15 )then
         call PGMTXT('L',2.0,0.5,0.5,'B_r (B-field units?)')
      elseif( iyvar .eq.16 )then
         call PGMTXT('L',2.0,0.5,0.5,'-d_M/d_t (sol/yr)')
      elseif( iyvar .eq.17 )then
         call PGMTXT('L',2.0,0.5,0.5,'Rossby number')
      elseif( iyvar .eq.18 )then
         call PGMTXT('L',2.0,0.5,0.5,'pnot [pressure at corona base]')
      elseif( iyvar .eq.19 )then
         call PGMTXT('L',2.0,0.5,0.5,'omega(kk) [rad/sec]')
      elseif( iyvar .eq.20 )then
         call PGMTXT('L',2.0,0.5,0.5,'omega(kk+1) [rad/sec]')
      elseif( iyvar .eq.21 )then
         call PGMTXT('L',2.0,0.5,0.5,'B.O. omega(kk) [rad/sec]')
      elseif( iyvar .eq.22 )then
         call PGMTXT('L',2.0,0.5,0.5,'B.O. omega(kk+1) [rad/sec]')
      elseif( iyvar .eq.23 )then
         call PGMTXT('L',2.0,0.5,0.5,'tauc [sec]')
      elseif( iyvar .eq.24 )then
         call PGMTXT('L',2.0,0.5,0.5,'taug [sec]')
      elseif( iyvar .eq.25 )then
         call PGMTXT('L',2.0,0.5,0.5,'v_conv [cm/sec]')
      endif
      if( ixvar .eq. 0 )then
         call PGMTXT('B',2.0,0.5,0.5,'time(s)')
      elseif( ixvar .eq. 1 )then
         call PGMTXT('B',2.0,0.5,0.5,'model')
      elseif( ixvar .eq. 2 )then
         call PGMTXT('B',2.0,0.5,0.5,'log Teff(K)')
      elseif( ixvar .eq. 3 )then
         call PGMTXT('B',2.0,0.5,0.5,'time(y)')
      elseif( ixvar .eq. 4 )then
         call PGMTXT('B',2.0,0.5,0.5,'B-V')
      elseif( ixvar .eq. 5 )then
         call PGMTXT('B',2.0,0.5,0.5,'log(t) [y]')   
      endif

      call PGMTXT('T',2.0,0.5,0.5,toplab)

c..   reset color
      call PGSCI(1)
c..   line width
      call PGSLW(lw)

c..   pgplot
      if( ixvar .eq. 0 )then
         do i = 1, ll
            xx(i) = time(i)
         enddo
      elseif( ixvar .eq. 1 )then
         do i = 1, ll
            xx(i) = model1(i)
         enddo
      elseif( ixvar .eq. 2 )then
         do i = 1, ll
            xx(i) = telog(i)
         enddo
      elseif( ixvar .eq. 3 )then
         do i = 1, ll
            xx(i) = time(i)/secpy
         enddo
      elseif( ixvar .eq. 4 )then
         do i = 1, ll
            xx(i) = BV(i)
         enddo
      elseif( ixvar .eq. 5 )then
         do i = 1, ll
            xx(i) = log10(time(i)/secpy)
         enddo   
      endif

      if( iyvar .eq. 0 )then
         do i = 1, ll
            yy(i) = xlol(i)
         enddo
      elseif( iyvar .eq. 1 )then
         do i = 1, ll
            yy(i) = -peryear(i)
         enddo
      elseif( iyvar .eq. 2 )then
         do i = 1, ll
            yy(i) = telog(i)
         enddo
      elseif( iyvar .eq. 3 )then
         do i = 1, ll
            yy(i) = speryear(i)
         enddo
      elseif( iyvar .eq. 5 )then
         do i = 1, ll
            yy(i) = V(i)
         enddo
      elseif( iyvar .eq. 6 )then
         do i = 1, ll
            yy(i) = it(i)
         enddo
      elseif( iyvar .eq. 7 )then
         do i = 1, ll
            yy(i) = xmkk(i)/sol
         enddo
      elseif( iyvar .eq. 8 )then
         do i = 1, ll
            yy(i) = rphot(i)
         enddo
      elseif( iyvar .eq. 9 )then
         do i = 1, ll
            yy(i) = dlog10( dth(i) )
         enddo
      elseif( iyvar .eq. 10 )then
         do i = 1, ll
c..   set up for log10 
            if( time(i) .gt. 0.0d0 )then
               yy(i) =  dlog10( time(i) )
               yy(i) = amax1( yy(i), -12.0 )
            else
               yy(i) = -12.0
            endif
         enddo
      elseif( iyvar .eq. 11 )then
         open(111, file="hzonerange.dat")
         do i = 1, ll
            yy(i) = (10.0e0**(xlol(i))*3.9e33/(16.0e0 * pi * sigma * 
     1              5.555e9))**0.5d0/1.5e13
            teff = 10.0d0**telog(i)-5780.0d0
            unlogl(i) = 10.0d0**xlol(i)
            seff(1,i) = 1.7753 + 1.4316d-4*teff + 2.9875d-9*teff**2.0d0
     1                - 7.5702d-12*teff**3.0d0 - 1.1635d-15*teff**4.0d0
            seff(2,i) = 1.0512 + 1.3242d-4*teff + 1.5418d-8*teff**2.0d0
     1                - 7.9895d-12*teff**3.0d0 - 1.8328d-15*teff**4.0d0
            seff(3,i) = 1.0140 + 8.1774d-5*teff + 1.7063d-9*teff**2.0d0
     1                - 4.3241d-12*teff**3.0d0 - 6.6462d-16*teff**4.0d0
            seff(4,i) = 0.3438 + 5.8942d-5*teff + 1.6558d-9*teff**2.0d0
     1                - 3.0045d-12*teff**3.0d0 - 5.2983d-16*teff**4.0d0
            seff(5,i) = 0.3179 + 5.4513d-5*teff + 1.5313d-9*teff**2.0d0
     1                - 2.7786d-12*teff**3.0d0 - 4.8997d-16*teff**4.0d0
            bound(1,i)= (unlogl(i)/seff(1,i))**0.5d0
            bound(2,i)= (unlogl(i)/seff(2,i))**0.5d0
            bound(3,i)= (unlogl(i)/seff(3,i))**0.5d0
            bound(4,i)= (unlogl(i)/seff(4,i))**0.5d0
            bound(5,i)= (unlogl(i)/seff(5,i))**0.5d0
            write(*,*)teff,(seff(n,i), n=1,5)
c            inner(1,i) = (0.72d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            inner(2,i) = (0.84d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            inner(3,i) = (0.95d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            inner(4,i) = (0.68d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            inner(5,i) = (0.76d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            inner(6,i) = (0.46d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            inner(7,i) = (0.51d0 - 2.7619d-5*(10.0d0**telog(i)-5700) -
c     1           3.8095d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0 
c            outer(1,i) = (1.77d0 - 1.3786d-4*(10.0d0**telog(i)-5700) -
c     1           1.4286d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0
c            outer(2,i) = (1.67d0 - 1.3786d-4*(10.0d0**telog(i)-5700) -
c     1           1.4286d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0
c            outer(3,i) = (1.95d0 - 1.3786d-4*(10.0d0**telog(i)-5700) -
c     1           1.4286d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0
c            outer(4,i) = (2.4d0 - 1.3786d-4*(10.0d0**telog(i)-5700) -
c     1           1.4286d-9*(10.0d0**telog(i)-5700)**2.0d0)*
c     1           (10.0d0**xlol(i))**0.5d0
c            write(111,'(1p12e13.6)')time(i), (inner(n,i), n=1,7), 
c     1           (outer(n,i), n=1,4)
            write(111, '(1p6e13.6)')time(i), (bound(n,i), n=1,5)
         enddo
         close(111)
      elseif( iyvar .eq. 4 )then
c..   done separately below
      elseif( iyvar .eq. 12 )then
         do i = 1, ll
            yy(i) = LOG10(Tnot(i))
         enddo
      elseif( iyvar .eq. 13 )then
         do i = 1, ll
            yy(i) = lx(i)
         enddo
      elseif( iyvar .eq. 14 )then
         do i = 1, ll
            yy(i) = omegadot(i)
         enddo
      elseif( iyvar .eq. 15 )then
         do i = 1, ll
            yy(i) = Br_star(i)
         enddo
      elseif( iyvar .eq. 16 )then
         do i = 1, ll
            yy(i) = Mdot_star(i)
         enddo
      elseif( iyvar .eq. 17 )then
         do i = 1, ll
            yy(i) = rossby(i)
         enddo
      elseif( iyvar .eq. 18 )then
         do i = 1, ll
            yy(i) = pnot(i)
         enddo
      elseif( iyvar .eq. 19 )then
         do i = 1, ll
            yy(i) = omegkk(i)
         enddo
      elseif( iyvar .eq. 20 )then
         do i = 1, ll
            yy(i) = omegkk1(i)
         enddo
      elseif( iyvar .eq. 21 )then
         do i = 1, ll
            yy(i) = bo_omega2(i)
         enddo
      elseif( iyvar .eq. 22 )then
         do i = 1, ll
            yy(i) = bo_omega(i)
         enddo
      elseif( iyvar .eq. 23 )then
         do i = 1, ll
            yy(i) = tauc(i)
         enddo
      elseif( iyvar .eq. 24 )then
         do i = 1, ll
            yy(i) = taug(i)
         enddo
      elseif( iyvar .eq. 25 )then
         do i = 1, ll
            yy(i) = vconv(i)
         enddo
      endif

      if( iyvar .ne. 4 )then
c..   draw curve
         if( lstyle .eq. 0 )then
            call PGPT(ll,xx,yy,20)
         elseif( lstyle .eq. 1 )then
            call PGLINE(ll,xx,yy)
         else
            call PGPT(ll,xx,yy,20)
            call PGLINE(ll,xx,yy)
         endif
         if( iyvar .eq. 5 .or. iyvar .eq. 0 )then
            call PGPT(1,xx(ll),yy(ll),-4)
         endif
      else

c..   draw remaining curves
         do j = 1, nuclei
c..   reset color
            call PGSCI(j+2)
            do i = 1, ll
               if( xxej(j,i) .gt. 1.0d-30 )then
                  yy(i) = dlog10( xxej(j,i) )
               else
                  yy(i) = -30.0
               endif
            enddo
            if( lstyle .eq. 0 )then
               call PGPT(ll,xx,yy,20)
            elseif( lstyle .eq. 1 )then
               call PGLINE(ll,xx,yy)
            else
               call PGPT(ll,xx,yy,20)
               call PGLINE(ll,xx,yy)
            endif
         enddo
      endif

c..............................................................

      do n = 2, nfiles

         l = 1
c..   hr.? data file
         open(8,file=hrfile(n))

 200     continue
         if( l .gt. ndm )goto 201

         read(8,70,end=201) model1(l), it(l), time(l), dth(l),
     1        rphot(l), rhoph(l), uphot(l), xlol(l), telog(l), 
     2        xmkk(l), peryear(l), vinf(l), omegkk(l),
     3        dmhkk1(l),omegkk1(l)
         read(8,71) xej
         read(8,73) Tnot(l),lx(l),omegadot(l),Br_star(l),
     1            Mdot_star(l),rossby(l),pnot(l),bo_omega(l),
     2           bo_omega2(l),tauc(l),taug(l),vconv(l),
     3           xbosurfomeg(l)
         Mdot_star(l) = -Mdot_star(l)/sol*secpy
         if( iyvar .eq. 4 )then
            do j = 1, nuclei
               xxej(j,l) = xej(nucleus(j))
            enddo
         endif
         l = l+1
         goto 200
 201     continue
         ll = l-1

         close(8)

c..   construct B-V, V
         if( ixvar .eq. 4 .and. iyvar .eq. 5 )then
            do m = 1, ll
               glog(m) = log10( grav*xmkk(m)/rphot(m)**2 )
c..   write(*,'(1p3e9.1)')xmkk(m),rphot(m),glog(m)
            enddo  
            
            do m = 1, ll
               feh(m) = fehhyad
            enddo
            do m = 1, ll
               teff = 10**telog(m)
               call kurintpt(feh(m), teff, glog(m), color)
               BV(m) = color(3)
               Mbol = 4.75 - 2.5 * xlol(m)
               V(m) = Mbol - color(1)
c     write(10,'(i5,1p6e11.3)')m,V(m),BV(m),telog(m),glog(m)
            enddo
         else
c     stop'second or greater model: error'
         endif

c..   begin buffer
         call PGBBUF

c..   reset color(cycle of 16 omitting 0=background)
         if( n .lt. 16 )then
            icolor = n
         elseif( n .lt. 32 )then
            icolor = n-15
         else
            icolor = n-30
         endif
         call PGSCI(icolor) 

         write(*,'(a20,2i5)')hrfile(n),n,icolor

c..   line width
         call PGSLW(lw)
         
c..   pgplot
         if( ixvar .eq. 0 )then
            do i = 1, ll
               xx(i) = time(i)
            enddo
         elseif( ixvar .eq. 1 )then
            do i = 1, ll
               xx(i) = model1(i)
            enddo
         elseif( ixvar .eq. 2 )then
            do i = 1, ll
               xx(i) = telog(i)
            enddo
         elseif( ixvar .eq. 3 )then
            do i = 1, ll
               xx(i) = time(i)/secpy
            enddo
         elseif( ixvar .eq. 4 )then
            do i = 1, ll
               xx(i) = BV(i)
            enddo
         endif

         if( iyvar .eq. 0 )then
            do i = 1, ll
               yy(i) = xlol(i)
            enddo
         elseif( iyvar .eq. 1 )then
            do i = 1, ll
               yy(i) = -peryear(i)
            enddo
         elseif( iyvar .eq. 2 )then
            do i = 1, ll
               yy(i) = telog(i)
            enddo
         elseif( iyvar .eq. 3 )then
            do i = 1, ll
               yy(i) = speryear(i)
            enddo
         elseif( iyvar .eq. 5 )then
            do i = 1, ll
               yy(i) = V(i)
            enddo
         elseif( iyvar .eq. 6 )then
            do i = 1, ll
               yy(i) = it(i)
            enddo
         elseif( iyvar .eq. 7 )then
            do i = 1, ll
               yy(i) = xmkk(i)/sol
            enddo
         elseif( iyvar .eq. 8 )then
            do i = 1, ll
               yy(i) = rphot(i)
            enddo
         elseif( iyvar .eq. 9 )then
            do i = 1, ll
               yy(i) = dlog10( dth(i) )
            enddo
         elseif( iyvar .eq. 10 )then
            do i = 1, ll
               yy(i) = dlog10( time(i) )
            enddo
         elseif( iyvar .eq. 11 )then
            do i = 1, ll
               yy(i) = (10.0e0**(xlol(i))*3.9e33/(16.0e0 * pi * sigma * 
     1              5.555e9))**0.5d0/1.5e13
            enddo
         elseif( iyvar .eq. 4 )then
c..   done separately below
         elseif( iyvar .eq. 12 )then
            do i = 1, ll
               yy(i) = LOG10(Tnot(i))
            enddo
         elseif( iyvar .eq. 13 )then
            do i = 1, ll
               yy(i) = lx(i)
            enddo
         elseif( iyvar .eq. 14 )then
            do i = 1, ll
               yy(i) = omegadot(i)
            enddo
         elseif( iyvar .eq. 15 )then
            do i = 1, ll
               yy(i) = Br_star(i)
            enddo
         elseif( iyvar .eq. 16 )then
            do i = 1, ll
               yy(i) = Mdot_star(i)
            enddo
         elseif( iyvar .eq. 17 )then
            do i = 1, ll
               yy(i) = rossby(i)
            enddo
         elseif( iyvar .eq. 18 )then
            do i = 1, ll
               yy(i) = pnot(i)
            enddo
         elseif( iyvar .eq. 19 )then
            do i = 1, ll
               yy(i) = omegkk(i)
            enddo
         elseif( iyvar .eq. 20 )then
            do i = 1, ll
               yy(i) = omegkk1(i)
            enddo
         elseif( iyvar .eq. 21 )then
            do i = 1, ll
               yy(i) = bo_omega2(i)
            enddo
         elseif( iyvar .eq. 22 )then
            do i = 1, ll
               yy(i) = bo_omega(i)
            enddo
         elseif( iyvar .eq. 23 )then
            do i = 1, ll
               yy(i) = tauc(i)
            enddo
         elseif( iyvar .eq. 24 )then
            do i = 1, ll
               yy(i) = taug(i)
            enddo
         elseif( iyvar .eq. 25 )then
            do i = 1, ll
               yy(i) = vconv(i)
            enddo
         endif

         if( iyvar .ne. 4 )then
c..   draw curve
            if( lstyle .eq. 0 )then
               call PGPT(ll,xx,yy,20)
            elseif( lstyle .eq. 1 )then
               call PGLINE(ll,xx,yy)
            else
               call PGPT(ll,xx,yy,20)
               call PGLINE(ll,xx,yy)
            endif
            if( iyvar .eq. 5 .or. iyvar .eq. 0 )then
               call PGPT(1,xx(ll),yy(ll),-4)
            endif
         else
c..   draw remaining curves
            do j = 1, nuclei
c..   reset color
               call PGSCI(j+2)
c..   only works for j+2 < 16
               do i = 1, ll
                  if( xxej(j,i) .gt. 1.0d-30 )then
                     yy(i) = dlog10( xxej(j,i) )
                  else
                     yy(i) = -30.0
                  endif
               enddo
               if( lstyle .eq. 0 )then
                  call PGPT(ll,xx,yy,20)
               elseif( lstyle .eq. 1 )then
                  call PGLINE(ll,xx,yy)
               else
                  call PGPT(ll,xx,yy,20)
                  call PGLINE(ll,xx,yy)
               endif
            enddo
         endif

c..   get buffer
         call PGEBUF
         call PGSCI(1)

      enddo

      if( ixvar .eq. 4 .and. iyvar .eq. 5 )then
         open(22,file='hyades2.dat')
 300     continue
         read(22,'(i6,0pf8.4,0p13f9.4)',end=301)idum,aa
         call pgsci(1)
c..   small error bars obscured by points
c     call pgpt(1,aa(1),aa(3),16)
         call pgerr1(5,aa(1),aa(3),aa(2),2)
         call pgerr1(6,aa(1),aa(3),aa(4),2)
         go to 300
 301     continue
         close(22)

      elseif( ixvar .eq. 2 .and. iyvar .eq. 0 )then
         open(22,file='hyades2.dat')
 310     continue
         read(22,'(i6,0pf8.4,0p13f9.4)',end=311)idum,aa
         call pgsci(1)
c..   small error bars obscured by points
         call pgerr1(5,aa(7),aa(11),aa(8),2)
         call pgerr1(6,aa(7),aa(11),aa(12),2)

         go to 310
 311     continue
         close(22)
      endif

      pause
      write(*,*)'normal termination'

      close(10)
      call PGEND
      stop

 70   format(2i6,1pe14.6,1p8e12.4,1p4e11.3)
 71   format(1p10e11.3) 
 73   format(1p13e11.3)


 2000 continue

      write(*,*)ixvar,labl

      write(*,'(//a40)')'hrt: error in input file: hrt.in'
      stop'hrt error'


      end 

