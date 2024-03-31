      program gennuc

c..   reads model for analysis of nuclear burning
c..   uses new gen format, pgplot
c..   last revised 2-6-2008 Dave Arnett

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'czone'
      include 'cgen'
      include 'cconst'
      include 'cburn'

      integer*4 mods, imods, i, k, nc, kworst, linesty, kmax
      integer*4 istyle, j
      integer*4 pgbeg

      integer*4 n,nz(ndim), nn(ndim), netrc

      parameter( mods = 3 )

      real*8    xmsol(mods), sum, yion, ye, ytot
      real*8    sume, diffold, difzold,xmet,ymet,zmet,deltay

      real*4    xx(kdm),yy(kdm),xmin,xmax,ymin,ymax
      real*4    xmin0,xmax0,ymin0,ymax0
      real*4    ax(mods,kdm),ay(mods,ndim,kdm)
      real*4    axmin(mods),axmax(mods),aymin(mods),aymax(mods)

      real*4    elemen(ndim),solarel(ndim),massinv,xlabel,ylabel

      real*4 pgfontsz
      integer*4 ifont,linewd

      integer*4 kkold(mods)

      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt
      character*5 cvarx

      character*7 cmod(mods)
      character*12 cxlabel, cylabel
      character*20 ctlabel
      character*5 device

      character*2 celem(32)

      data celem/ ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     1     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     2     'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3     'Ga','Ge'/

      logical tobe
      data txtxt/'  '/
c..defaults
c..   pgfontsz=1.4, ifont=1, linewd=5 work well (for transparencies)
c..   on HP 5mp laser printer (postscript)
      data pgfontsz /1.4/, ifont/1/, linewd/5/
c-----------------------------------------------------------------

      write(*,*)'starting gennuc'

      inquire(file='gennuc.in',exist=tobe)
      if( tobe )then
         open(2,file='gennuc.in',form='formatted',status='old')
      else
         write(*,*)'gennuc: no file gennuc.in in directory'
         stop'gennuc: no gennuc.in input file error'
      endif
      write(*,*)'reading file: gennuc.in'
c..   input parameters for model modification from gentran.in

c..   dummy reads for labelling header lines
      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text


      read (2,*)txt,labl,device
      write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,istyle
      write(*,*)txt,txtxt,labl,istyle
      labl1 = 'istyle'

      if( labl .ne. labl1 )goto 2000
      read (2,*)txt,labl,linesty
      write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,deltay
      write(*,*)txt,txtxt,labl,deltay
      labl1 = 'deltay'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,imods
      write(*,*)txt,txtxt,labl,imods
      labl1 = 'imods'
      if( labl .ne. labl1 )goto 2000

      if( imods .le. 0 .or. imods .gt. mods )goto 2000

c..   read cycle for model identities
      do i = 1, imods
         read (2,*)txt,labl,cmod(i)
         write(*,*)txt,txtxt,labl,cmod(i)
         labl1 = 'cmod'
         if( labl .ne. labl1 )goto 2000
      enddo

      read (2,*)txt,labl,cvarx
      write(*,*)txt,txtxt,labl,cvarx
      labl1 = 'cvarx'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmin0
      write(*,*)txt,txtxt,labl,xmin0
      labl1 = 'xmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmax0
      write(*,*)txt,txtxt,labl,xmax0
      labl1 = 'xmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymin0
      write(*,*)txt,txtxt,labl,ymin0
      labl1 = 'ymin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymax0
      write(*,*)txt,txtxt,labl,ymax0
      labl1 = 'ymax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)note
      write(*,*)note

      read (2,*)text
      write(*,*)text


c..   check network and model consistency (for abundance sanity)
      inquire(file='net.rc',exist=tobe)
      if( tobe )then
         open(9,file='net.rc',form='formatted',status='old')
 1001    continue
         read(9,'(3i5,a5,0pf10.4,1pe12.4)',end=1000)
     1        n,nz(n),nn(n),cnuc(n),qex(n),solarx(n)
         if( nz(n) .eq. 2 .and. nn(n) .eq. 2 )then
            read(9,'(10i5)')nucp
            goto 1000
         endif
         goto 1001
      else
         write(*,*)'no file net.rc in current directory'
         stop' gennuc error 1'
      endif
 1000 continue
      netrc   = n
      netsize = n
      write(*,*)'read net.rc, netsize = ',netsize
      close(9)

c..   loop to read models and extract xx, yy
      do i = 1, imods
         inquire(file=cmod(i),exist=tobe)
         if( tobe )then
c..   formatted read in new format
            open(8,file=cmod(i))
            rewind 8
            nc = 1
            call readf(note,nc)
            close(8)
         else
            write(*,*)'gennuc: no file ',cmod(i),' in directory '
            write(*,*)'line ',i,' in list'
            stop'gennuc file list error'
         endif

         if( modes .eq. 2 )then
            kmax = kk+1
         else
            kmax = kk
         endif

c..   test model and net.rc consistency
         if( netrc .ne. netsize )then
            write(*,*)'net.rc inconsistent with model'
            write(*,*)netrc,' entries in net.rc'
            write(*,*)netsize,' nuclei of network in model'
            stop'gennuc: netsize'
         endif

c..   define mass coordinate
         dmi(1)    = dmh(2)
         dmi(kmax+1) = dmh(kmax+1)
         open(32,file=cmod(i)//".nuc")
         open(33,file=cmod(i)//".vconv")
         do k = 2, kmax
            xm(k)  = xm(k-1) + dmh(k)
            dmi(k) = 0.5*( dmh(k+1) + dmh(k) )
            write(32,'(177(1pe12.4))')(x(n,k),n=1,nnuc)
            write(33,'(e12.4)')h(k)*100.0d0/6.96d10
	    write(33,'(e12.4)')h(k)*100.0d0/6.96d10
         enddo
         close(32)
         close(33)
         xmsol(i) = xm(kmax)/sol
         write(*,'(a10,a10,0pf10.4,a15,a8,i5,a8)')
     1        cmod(i),' MASS =',xmsol(i),' solar masses,',
     2        ' nnuc =',nnuc,' nuclei'

         
            
         
c..   test normalization
         diffold = 0.0d0
         difzold = 0.0d0
         kworst  = 0
         do k = 2, kmax
            sum = 0.0d0
            do n = 1, nnuc
               sum = sum + x(n,k)
            enddo

            sume      = 0.0d0
            do n = 1, nnuc
               sume  = sume + x(n,k) * dble( nz(n) )/xa(n)
            enddo
            ye   = x(ndim,k)
            yion = sume
            ytot = ye + yion
            if( abs( sum-1.0d0 ) .gt. abs(diffold) )then
               kworst = k
               diffold = sum  - 1.0d0
               difzold = sume - ye
            endif
         enddo

         if( abs(diffold) .gt. 1.0d-6 .or.
     1        abs(difzold) .gt. 1.0d-6 )then
            write(*,'(a10,i5,a10,1pe12.3,a10,1pe12.3)')
     1           'kworst',kworst,'diff X',diffold,'diff Z',difzold
         endif

         xmet = x(ndim-2,kmax)
         ymet = x(ndim-1,kmax)
         zmet = 1.0d0 - xmet - ymet
         write(*,'(3(a5,0pf10.6),a10,1pe12.3)')
     1        'X', xmet, 'Y', ymet, 'Z', zmet,
     2        'z/zsun',zmet/0.0153d0
c..   0.0191 for 80 specie network
c..   0.0148 for Lodders
         kkold(i) = kmax

         if( istyle .eq. 0 )then
c..   lineout of mass fractions

c..   defines x-coordinate

            call getvec(xx,xmin,xmax,cvarx,cxlabel)

            do k = 1, kmax
               ax(i,k) = xx(k)
            enddo

c..   defines y-coordinate values
            cylabel = 'log X(Z,N)'
c..Ye
            do k = 2, kmax
               x(nnuc+1,k) = 0.0d0
               do n = 1, nnuc
                  x(nnuc+1,k) = x(nnuc+1,k) + x(n,k)*dble(nz(n))/xa(n)
               enddo
            enddo
c..Xi
            do k = 2, kmax
               do n = 1, nnuc
                  ay(i,n,k) = x(n,k)
               enddo
               ay(i,nnuc+1,k) = x(nnuc+1,k)
            enddo
c..set inner ghost zone
            do n = 1, nnuc+1
               ay(i,n,1) = ay(i,n,2)
            enddo
            ymin = -3.0
            ymax = 0.0

         else

c..   elemental abundances versus Z
            xmin = 0.0
            xmax = xmin
            do n = 1, nnuc
               xmax = amax1( float( nz(n) ) , xmax )
            enddo
            do j = 1, int(xmax)
               xx(j) = float( j )
               ax(i,j) = xx(j)
            enddo

c..   get solar elemental abundances in solarel(nz)
            do j = 1, int( xmax)
               solarel(j) = 0.0
               do n = 1, nnuc
                  if( nz(n) .eq. j )then
                     solarel(j) = solarel(j) 
     1                    + solarx(n)
                  endif
               enddo
            enddo
c..   get abundances for all zones
            massinv = 1.0d0/(xm(kk) - xm(1))
c..   each element
            do j = 1, int( xmax)
               elemen(j) = 0.0
c..   each zone
               do k = 2, kk
                  do n = 1, nnuc
                     if( nz(n) .eq. j )then
                        elemen(j) = elemen(j) 
     1                       + x(n,k)*dmh(k)*massinv
                     endif
                  enddo
c..   use log of sum relative to solar
                  if( k.eq. kk)then
                     if( elemen(j)/solarel(n) .gt. 0.0 )then
                        yy(j) = alog10( elemen(j)/solarel(j) )
                     endif
                  endif
               enddo
            enddo
c..store in first element; array compacted over zone index k
            do j = 1, int(xmax)
               ay(i,j,1) = yy(j)
            enddo
            ymin = 0.0
            ymax = ymin
            do j = 1, int(xmax)
               ymin = amin1(ymin,yy(j))
               ymax = amax1(ymax,yy(j))
            enddo

            cxlabel = 'Z'
            cylabel = 'log(X/Xsol)'
            kkold(i) = int(xmax)
         endif

         axmin(i) = xmin
         axmax(i) = xmax
         aymin(i) = ymin
         aymax(i) = ymax
      enddo

c..   find limits over all models
      do i = 1, imods-1
         xmin = amin1( xmin, axmin(i) )
         xmax = amax1( xmax, axmax(i) )
         ymin = amin1( ymin, aymin(i) )
         ymax = amax1( ymax, aymax(i) )
      enddo

c..   override if xmin0, etc. are defined
      if( xmin0 .ne. xmax0 )then
         xmin = xmin0
         xmax = xmax0
      endif
      if( ymin0 .ne. ymax0 )then
         ymin = ymin0
         ymax = ymax0
      endif

      xmin0 = xmin - (xmax-xmin)*0.05
      xmax0 = xmax + (xmax-xmin)*0.05
      ymin0 = ymin - (ymax-ymin)*0.05
      ymax0 = ymax + (ymax-ymin)*0.05


c..   scratch file for character conversion
      open(10,file='scratch.pg')

      IF ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
      call PGASK (.FALSE.)
c..   roman font=2, clean=1
      call PGSCF(ifont)
c..   font scaling (size)
      call PGSCH(pgfontsz)
c..   line width
      call PGSLW(linewd)

c..   new page
      call PGPAGE

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

      if( istyle .eq. 0 )then
         do i = 1, imods
            kmax = kkold(i)
            do k = 1, kmax
               xx(k) = ax(i,k)
            enddo

            do n = 1, ndim
               do k = 1, kmax
                  if( ay(i,n,k) .gt. 10.0**(ymin-1.0d0) )then
                     yy(k) = log10( ay(i,n,k) )
                  else
                     yy(k) = ymin -1
                  endif
               enddo

c..   plot only nuclei which change, to minimize congestion
               if( abs(yy(2)-yy(kmax) ) 
     1             .ge. abs(deltay*yy(kmax)) )then
c..   set color
                  call PGSCI(mod(n-1,15)+1)

                  if( n .eq. nnuc+1 )then
c..   Ye
                     call PGSCI(1)
                  endif
c..   draw curve
                  if( linesty .eq. 0 )then
                     call PGLINE(kmax,xx,yy)
                  else
                     call PGPT(kmax,xx,yy,20)
                  endif
c..   label
                  ctlabel = cnuc(n)
                  call PGSCH(0.7)
c..   is value significant?
                  if( yy(kmax) .gt. ymin )then
c..   is value increasing?
                     if( yy(kmax) .gt. yy(1) )then
c..   labels on right hand side
                        call PGTEXT(xmax,yy(kmax),ctlabel)
                     endif
                  endif
c..   is value significant?
                  if( yy(1) .gt. ymin )then
c..   is value decreasing or constant?
                     if( yy(1) .ge. yy(kmax) )then
c..   labels on left hand side
                        call PGTEXT(xmin0,yy(1),ctlabel)
                        if(aint(real(n/2))-n/2 .lt. 0)xmin0 = xmin0+0.5
                     endif
                  endif
               endif
            enddo
         enddo

c..   identity of model if only one plotted
         if( imods .eq. 1 )then
            call PGSCI(1)
            call PGTEXT(0.03*(xmin0+xmax0),
     1           ymax0+0.03*(ymax0-ymin0),cmod(1))
         endif
      else
         
c..   elemental abundances
         do i = 1, imods
            kmax = kkold(i)
            do k = 1, kmax
               xx(k) = ax(i,k)
               yy(k) = ay(i,k,1)
            enddo

            do k = 1, kmax
               yy(k) = ay(i,k,1) 
c..   set color
               call PGSCI(mod(k-1,15)+1)
c..   draw points
               call PGPT1(xx(k),yy(k),16)
c..   label elements
               xlabel = k
               ylabel = yy(k)
               call PGTEXT(xlabel,ylabel,celem(k))
            enddo
         enddo

c..   add zero line
         call PGSCI(15)
         call PGSLS(4)
         call PGMOVE(xmin0,0.0)
         call PGDRAW(xmax0,0.0)
         call PGSLS(1)
         call PGSCI(1)

c..   identity of model if only one plotted
         if( imods .eq. 1 )then
            call PGSCI(1)
            call PGTEXT(0.03*(xmin0+xmax0),
     1           ymax0+0.03*(ymax0-ymin0),cmod(1))
         endif
      endif

      if( device .eq. '/cps' .or. device .eq. '/ps' )then
c..   avoid pause, make hardcopy
      else
c..   xwindow
         pause
      endif

      close(10)

      call PGEND

      stop'successful termination'

 2000 continue

      write(*,*)'gennuc: error in  input file: gennuc.in'
      stop'gennuc error'

      end


