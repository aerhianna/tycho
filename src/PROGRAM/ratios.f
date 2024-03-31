      program ratios

c..   reads hr.? file from standard input and
c..   uses pgplot to plot 
c..   modified for extended network
c..   detailed abundance ratios
c..   readf uses unit 5

      implicit none

      include "dimenfile"
      include "cconst"
      include "cburn"

      integer*4 icolor

      integer*4 l, ll, m, ixvar, iyvar, lstyle, n, nfiles, numfiles
      integer*4 pgbeg

      integer*4 ndm
      parameter( ndm = 100000, numfiles = 20 )

      integer*4 model1(ndm), it(ndm)
      real*8 xej(ndim),xxej(5,ndm)
      real*8 time(ndm),dth(ndm),rphot(ndm),rhoph(ndm),uphot(ndm),
     1     xlol(ndm),telog(ndm),xmkk(ndm),
     2     peryear(ndm),speryear(ndm),omegkk(ndm)
      real*8 smass
      real*8 feh(ndm), glog(ndm)
      real*8 teff, Mbol
      real*8 color(5)
      real*8 V(ndm), BV(ndm)

      integer*4 i,j
      integer*4 lw
      real*4 fsize

      real*4 xx(ndm),yy(ndm)
      real*4 xmin,xmax,ymin,ymax,xtest,ytest
      real*4 xmin0,xmax0,ymin0,ymax0

      character*5 hrfile(numfiles)
      character*40 toplab
      character*5 device

      integer*4 nuclei,nucleus(2)
      integer*4 ipz(2),ipn(2)

      integer*4 idum
      real*4    aa(14)
      real*4    fehhyad

      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt
      character*12 cylabel

      logical tobe

      data l/1/
      data lw/3/,fsize/1.2/
c..   ag88 value relative to opal normalization = 0.02
      data fehhyad/0.0/
      data txtxt/'  '/
c--------------------------------------------------------------------

c..   check network and model consistency (for abundance sanity)
      inquire(file='net.rc',exist=tobe)
      if( tobe )then
         open(9,file='net.rc',form='formatted',status='old')
 1001    continue
         read(9,'(3i5,a5,0pf10.4,1pe12.4)',end=1000)
     1        n,lz(n),ln(n),cnuc(n),qex(n),solarx(n)
         if( lz(n) .eq. 2 .and. ln(n) .eq. 2 )then
            goto 1000
         endif
         goto 1001
      else
         write(*,*)'no file net.rc in current directory'
         stop' ratios error 1'
      endif
 1000 continue
      close(9)


c..   read graphics options
      inquire(file='ratios.in',exist=tobe)
      if( tobe )then
         open(9,file='ratios.in')
      else
         stop'ratios: no ratios.in input file error'
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
      else
         write(*,*)'error in ratios.in, ixvar = ',ixvar
         stop
      endif

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
         write(*,*)'y axis is abundance ratio'
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
         write(*,*)'y axis is photospheric velocity'
      else
         write(*,*)'error in ratios.in, iyvar = ',iyvar
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

      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text
      read (9,*)text
      write(*,*)text

      read(9,*)nuclei

      if( nuclei .eq. 0 )then
         write(*,*)'no nuclei used; nuclei =',nuclei
      elseif( nuclei .gt. 0 .and. nuclei .le. 2 )then
c         read(9,*)(nucleus(i),i=1,nuclei)
         write(*,'(a35,5i5)')'number of nuclei used here:',nuclei
c         write(*,'(a35,5i5)')'indices of nuclei used here are:',
c     1        (nucleus(i),i=1,nuclei)
         do i=1,nuclei
            read(9,*)ipz(i),ipn(i)
c            write(*,'(10i5)')i,ipz(i),ipn(i)
         enddo
c..identify nuclei
c         write(*,*)nnuc
         do i = 1,nuclei
            do j = 1, nnuc
               if( ipz(i) .eq. lz(j) .and. ipn(i) .eq. ln(j) )then
                  nucleus(i) = j
c                  write(*,*)lz(j),ln(j),cnuc(j),nucleus(i)
               endif
            enddo
         enddo
         do i = 1,nuclei,2
            write(*,'(i5,2a12)')i,' pair: ', 
     1           cnuc(nucleus(i))//'/'//cnuc(nucleus(i+1))
         enddo
      else
         write(*,*)'nuclei ',nuclei
         stop'ratios: error in nuclei'
      endif

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
         write(*,*)'ratios.in: lstyle error ',lstyle
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
     2        xmkk(l), peryear(l), speryear(l), omegkk(l)
         read(8,71) xej

c         write(*,*)nuclei,nucleus(1),nucleus(2),
c     1        xej(nucleus(1)),xej(nucleus(2))

         if( iyvar .eq. 4 )then
            do j = 1, nuclei,2
               xxej(j,l) = xej(nucleus(j))/xej(nucleus(j+1))
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
c         stop'ccc'
ccccccccccc
      else
         write(*,'(a10,a5)')'no file ',hrfile(1)
         stop'ratios: input error'
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
            stop'ratios: no kurcoljkl.tbl input file error'
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
         j = 1
         do l = 1, ll
            ymin = dmin1(xxej(j,l),dble(ymin))
         enddo
         write(*,*)ymin,' ymin'
cccccccccc
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
         ymin = 0.
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
         j = 1
         ymax = xxej(j,1)
         do l = 1, ll
            ymax = dmax1( xxej(j,l), dble(ymax) )
         enddo
         write(*,*)ymax,' ymax'
ccccccccccc
c..   set up for log10 
c         if( ymin .gt. 0.0d0 )then
c            ymin = amax1( alog10(ymin), -12.0 )
c         else
c            ymin = -12.0
c         endif
c         write(*,*)ymin,' ymin'
cccccccccc
c         if( ymax .gt. 0.0d0 )then
c            ymax = amin1( alog10(ymax), 0.0 )
c         else
c            ymax = 0.0
c         endif
c         write(*,*)ymax,' ymax'
cccccccccccccccc
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
         ymin = 0
         ymax = ymin
         do l=1, ll
            ymax = dmax1( uphot(l), dble(ymax) )
            ymin = dmin1( uphot(l), dble(ymin) )
         enddo
      endif

c..   allow input values in ratios.in to override if they are 
c..   more restrictive
      if( xmin0 .ne. xmax0 )then
         if( xmax0 .lt. xmin .and. ixvar .ne. 2 )then
c..   ixvar is log Te, astronomers invert axis for HR diagram
            write(*,*)'ratios warning: xmax0 < xmin, reset your xmax0?'
         elseif( xmin0 .gt. xmax )then
            write(*,*)'ratios warning: xmin0 > xmin, reset your xmin0?'
         endif
      endif
      if( ymin0 .ne. ymax0 )then
         if( ymax0 .lt. ymin .and. iyvar .ne. 5 )then
            write(*,*)'ratios warning: ymax0 < ymin, reset your ymax0?'
         elseif( ymin0 .gt. ymax )then
            write(*,*)'ratios warning: ymin0 > ymax, reset your ymin0?'
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


      write(*,*)xmin,xmax,ymin,ymax
cccccccccccc

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
         call PGMTXT('L',2.0,0.5,0.5,'log L/sol')
      elseif( iyvar .eq. 1 )then
         call PGMTXT('L',2.0,0.5,0.5,'M(sol)/year')
      elseif( iyvar .eq. 2 )then
         call PGMTXT('L',2.0,0.5,0.5,'log Te')
      elseif( iyvar .eq. 3 )then
         call PGMTXT('L',2.0,0.5,0.5,' ')
      elseif( iyvar .eq. 4 )then
         cylabel = cnuc(nucleus(1))//'/'//cnuc(nucleus(2))
         write(*,*)cylabel
         call PGMTXT('L',2.0,0.5,0.5,cylabel)
ccccccccccccccccc
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
      elseif( iyvar .eq. 11 )then
         call PGMTXT('L',2.0,0.5,0.5,'u(phot)')
      endif
      if( ixvar .eq. 0 )then
         call PGMTXT('B',2.0,0.5,0.5,'time(s)')
      elseif( ixvar .eq. 1 )then
         call PGMTXT('B',2.0,0.5,0.5,'model')
      elseif( ixvar .eq. 2 )then
         call PGMTXT('B',2.0,0.5,0.5,'log Te')
      elseif( ixvar .eq. 3 )then
         call PGMTXT('B',2.0,0.5,0.5,'time(y)')
      elseif( ixvar .eq. 4 )then
         call PGMTXT('B',2.0,0.5,0.5,'B-V')
      endif

      call PGMTXT('T',2.0,0.5,0.5,toplab)

c..   reset color
      call PGSCI(2)
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
c..   set up for log10 
            if( time(i) .gt. 0.0d0 )then
               yy(i) =  dlog10( time(i) )
               yy(i) = amax1( yy(i), -12.0 )
            else
               yy(i) = -12.0
            endif
         enddo      
      elseif( iyvar .eq. 11 )then
         do i = 1, ll
            yy(i) = uphot(i)
         enddo
      elseif( iyvar .eq. 4 )then
c..   done separately below
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
         do j = 1, nuclei,2
c..   reset color
            call PGSCI(j+2)
            do i = 1, ll
                  yy(i) = xxej(j,i) 
c               if( xxej(j,i) .gt. 1.0d-30 )then
c                  yy(i) = dlog10( xxej(j,i) )
c               else
c                  yy(i) = -30.0
c               endif
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
     2        xmkk(l), peryear(l), speryear(l), omegkk(l)
         read(8,71) xej
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
         elseif( iyvar .eq. 4 )then
c..   done separately below
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

 70   format(2i6,1pe14.6,1p8e12.4,1p2e11.3)
 71   format(1p10e11.3) 


 2000 continue

      write(*,*)ixvar,labl

      write(*,'(//a40)')'ratios: error in input file: ratios.in'
      stop'ratios error'


      end 

