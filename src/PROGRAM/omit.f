      program omit

c..   reads  model for mantle removal

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'czone'
      include 'cgen'
      include 'cconst'
      include 'cburn'
      include 'cnabla'

c      integer*4 n, nz(ndim), nn(ndim)
c      character*5 cnuc(ndim)
c      common /cnetrc/ nz,nn,cnuc

      integer*4 i, k, n, nc, kworst, linesty
c..new kk after slicing is ksl
      integer*4 ksl,inside
      integer*4 pgbeg


      real*8    xmsol, sum, yion, ye, ytot
      real*8    sume, diffold, difzold,xmet,ymet,zmet,slice

c..   explosion energy (ergs) = explode
c..   < 0 or = 0 gives no explosion
c..   explosion mass into which it is put =massex
c..   (at bottom of grid by default)
      real*8    explode,massex
c..   index of zone boundary outside explosion
      integer*4 kexplode

      real*4    xx(kdm),yy(kdm),xmin,xmax,ymin,ymax
      real*4    xmin0,xmax0,ymin0,ymax0

      real*4    xminf,xmaxf,yminf,ymaxf
      real*4    xlabel,ylabel

      real*4    tmass

      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt
      character*5 cvarx, cvary

      character*12 cxlabel, cylabel
      character*20 ctlabel
      character*5  cnum
      character*5 device
      character*9 cmod

      logical tobe
      data txtxt/'  '/

c-----------------------------------------------------------------
      inquire(file='omit.in',exist=tobe)
      if( tobe )then
         open(2,file='omit.in',form='formatted',status='old')
      else
         write(*,*)'omit: no file omit.in in directory'
         stop'omit: no .in file error'
      endif

c..   dummy read
      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text

c..   input parameters
      read (2,*)txt,labl,device
      write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,linesty
      write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nopac
      write(*,*)txt,txtxt,labl,nopac
      labl1 = 'nopac'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nkscale
      write(*,*)txt,txtxt,labl,nkscale
      labl1 = 'nkscale'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,fthoul
      write(*,*)txt,txtxt,labl,fthoul
      labl1 = 'fthoul'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,noion
      write(*,*)txt,txtxt,labl,noion
      labl1 = 'noion'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,norad
      write(*,*)txt,txtxt,labl,norad
      labl1 = 'norad'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nomole
      write(*,*)txt,txtxt,labl,nomole
      labl1 = 'nomole'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,nocoul
      write(*,*)txt,txtxt,labl,nocoul
      labl1 = 'nocoul'
      if( labl .ne. labl1 )goto 2000

c      read (2,*)txt,labl,noburn
c      write(*,*)txt,txtxt,labl,noburn
c      labl1 = 'noburn'
c      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,neutro
      write(*,*)txt,txtxt,labl,neutro
      labl1 = 'neutro'
      if( labl .ne. labl1 )goto 2000

c      read (2,*)txt,labl,jnb
c      write(*,*)txt,txtxt,labl,jnb
c      labl1 = 'jnb'
c      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nopaleos
      write(*,*)txt,txtxt,labl,nopaleos
      labl1 = 'nopaleos'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,mixmode
      write(*,*)txt,txtxt,labl,mixmode
      labl1 = 'mixmode'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,cmod
      write(*,*)txt,txtxt,labl,cmod
      labl1 = 'cmod'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,cvarx
      write(*,*)txt,txtxt,labl,cvarx
      labl1 = 'cvarx'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,cvary
      write(*,*)txt,txtxt,labl,cvary
      labl1 = 'cvary'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xminf
      write(*,*)txt,txtxt,labl,xminf
      labl1 = 'xminf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmaxf
      write(*,*)txt,txtxt,labl,xmaxf
      labl1 = 'xmaxf'
      if( labl .ne. labl1 )goto 2000

      if( xminf .eq. xmaxf )then
         write(*,*)'     no x override'
      endif

      read (2,*)txt,labl,yminf
      write(*,*)txt,txtxt,labl,yminf
      labl1 = 'yminf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymaxf
      write(*,*)txt,txtxt,labl,ymaxf
      labl1 = 'ymaxf'
      if( labl .ne. labl1 )goto 2000

      if( yminf .eq. ymaxf )then
         write(*,*)'     no y override'
      endif

      read (2,*)txt,labl,inside
      write(*,*)txt,txtxt,labl,inside
      labl1 = 'inside'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,slice
      write(*,*)txt,txtxt,labl,slice
      labl1 = 'slice'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,explode
      write(*,*)txt,txtxt,labl,explode
      labl1 = 'explode'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,massex
      write(*,*)txt,txtxt,labl,massex
      labl1 = 'massex'
      if( labl .ne. labl1 )goto 2000



c..   dummy read
      read (2,*)text
      write(*,*)text

      close(2)

      call kappinit

c..   read network file
      inquire(file='net.rc',exist=tobe)
      if( tobe )then
         open(9,file='net.rc',form='formatted',status='old')
 1001    continue
         read(9,'(3i5,a5,0pf10.4,1pe12.4)',end=1000)n,lz(n),ln(n),
     1        cnuc(n),qex(n),solarx(n)
         xa(n) = dble( lz(n) + ln(n) ) + qex(n)/931.5d0

         if( lz(n) .eq. 2 .and. ln(n) .eq. 2 )then
            goto 1000
         endif
         goto 1001
      else
         write(*,*)'no file net.rc in current directory'
         stop' omit error 1'
      endif
 1000 continue
      close(9)
      netsize = n


c..   formatted read of model files
      inquire(file=cmod,exist=tobe)
      if( tobe )then
         open(8,file=cmod)

         rewind 8
         nc = 1
         call readf(note,nc)
         close(8)

c..   define mass coordinate
         dmi(1)    = dmh(2)
         dmi(kk+1) = dmh(kk+1)
         do k = 2, kk
            xm(k)  = xm(k-1) + dmh(k)
            dmi(k) = 0.5*( dmh(k+1) + dmh(k) )
         enddo
         xm(kk+1) = xm(kk) + dmh(kk+1)
         xmsol    = xm(kk+1)/sol

         write(*,*)cmod,' TOTAL MASS =',xmsol,
     1        ' solar masses,',' nnuc =',nnuc,' nuclei'
         write(*,*)xm(kk)/sol,' solar masses on grid'
         write(*,*)dmh(kk+1)/sol,' solar masses in envelope'

c..   check network and model consistency (for abundance sanity)


         if( n .ne. nnuc )then
            write(*,*)'net.rc inconsistent with dimenfile'
            write(*,*)n,' entries in net.rc'
            write(*,*)nnuc,' nuclei of network in dimenfile'
            write(*,*)'resetting to net.rc abundances'
         endif

c..   test normalization
         diffold = 0.0d0
         difzold = 0.0d0
         kworst  = 0
         do k = 2, kk
            sum = 0.0d0
            do n = 1, nnuc
               sum = sum + x(n,k)*xa(n)
            enddo

            sume      = 0.0d0
            do n = 1, nnuc
               sume  = sume + x(n,k) * dble( lz(n) )
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

         xmet = x(ndim-2,kk)
         ymet = x(ndim-1,kk)
         zmet = 1.0d0 - xmet - ymet
         write(*,'(3(a5,0pf10.6),a10,1pe12.3)')
     1        'X', xmet, 'Y', ymet, 'Z', zmet,
     2        'z/zsun',zmet/0.0153d0

c..   0.017605 is for 37 nucleus network,
c..   full network would be 0.0191

c..   fill kk eos arrays
         call state(2,kk,nc)

c..   x-coordinate
         call getvec(xx,xmin,xmax,cvarx,cxlabel)

c..   y-coordinate
         call getvec(yy,ymin,ymax,cvary,cylabel)

      else
         write(*,*)'omit: no file ',cmod
         stop'omit: no model file error'
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

      xmin0 = xmin - (xmax-xmin)*0.05
      xmax0 = xmax + (xmax-xmin)*0.05
      ymin0 = ymin - (ymax-ymin)*0.05
      ymax0 = ymax + (ymax-ymin)*0.05


c..   scratch file for character conversion
      open(10,file='scratch.pg')

      IF ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
      call PGASK (.FALSE.)
c..   roman font
      call PGSCF(2)
c..   font scaling (size)
      call PGSCH(1.0)
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

c..   set color
      call PGSCI(2)
c..   draw curve
      if( linesty .eq. 0 )then
         call PGLINE(kk,xx,yy)
      elseif( linesty .eq. 1 )then
         call PGPT(kk,xx,yy,20)
      else
         call PGLINE(kk,xx,yy)
         call PGPT(kk,xx,yy,20)
      endif
c..   label
      tmass   = xmsol
      call ftoc(tmass,cnum)
      ctlabel = cnum//' '//cmod
      xlabel = xmin
      ylabel = ymax - float(i-1)*(ymax-ymin)*0.05

      call PGTEXT(xlabel,ylabel,ctlabel)

      ksl = kk
      do k = 1, kk
         if( xm(k)/sol .gt. slice )then
            write(*,*)k,xx(k),slice
            ksl = k
            goto 100
         endif
      enddo
 100  continue

      xx(1) = slice
      xx(2) = xx(1)
      yy(1) = ymin
      yy(2) = ymax
      call PGSCI(14)
      call PGLINE(2,xx,yy)

      pause
cccccccccccccccccccccccccccccccccc

      close(10)

      call PGEND

      open(8,file='newomit')

      rewind 8
      nc = 1

c      write(*,*)modes,kk
c      do k = kk-3, kk+3
c         write(*,'(i5,1p8e12.4)')k,dmh(k),xm(k),r(nc,k),t(nc,k),v(nc,k),
c     1        p(nc,k),x(netsize,k)
c      enddo
ccccccccccccccccccccccccccccccc

c..do sliceing

      if( inside .ge. 0 )then
c..slice off the outside
         kk = ksl
      else
c..slice off the inside
c..move k to k-ksl+1

c..   map outer zones, overwriting inner zones which are discarded
         do i = 1, kk+1-ksl+2
            k = i+ksl-1
            r(nc,i) = r(nc,k)
            u(nc,i) = u(nc,k)
            omeg(i) = omeg(k)
            ajay(i) = ajay(k)
            h(i)    = h(k)
            xm(i)  = xm(k)
            tl(nc,i) = tl(nc,k)

            t(nc,i+1) = t(nc,k+1)
            v(nc,i+1) = v(nc,k+1)
            p(nc,i+1) = p(nc,k+1)
            dmh(i+1)  = dmh(k+1)

            do n = 1, netsize+1
               x(n,i+1) = x(n,k+1)
            enddo

         enddo
         kk = kk - ksl 
      endif

      if( modes .eq. 2 )then
c..   envelope is outer zone kk+1
         kk = kk+1
         write(*,*)kk,ksl
c..   zero for safety (these were nonzero j504547)
ccccccccccccccccccccccccccc this may need to be extended to other variables
         xm(kk+2) = 0.0d0
         xm(kk+3) = 0.0d0
         r(nc,kk+3) = 0.0d0
      endif

c      write(*,*)modes,kk
c      do k = kk-3, kk+3
c         write(*,'(i5,1p8e12.4)')k,dmh(k),xm(k),r(nc,k),t(nc,k),v(nc,k),
c     1        p(nc,k),x(netsize,k)
c      enddo
ccccccccccccccccccccccccccc

c..add explosion
      if( explode .gt. 0.0d0 .and. massex .gt. 0.0d0 )then
         kexplode = 0
         do k = 2, kk
            if( xm(k)-xm(1) .gt. massex )then
               kexplode = k
               go to 900
            endif
         enddo
         write(*,*)'OMIT error: massex too large? ',massex,xm(kk)-xm(1)
         stop'OMIT: massex error'
 900     continue
         write(*,*)'ADDING EXPLOSION INSIDE'
         write(*,*)kexplode,xm(kexplode)-xm(1)
         u(1,2) = sqrt( 2.0d0 * explode /massex )
         write(*,'(a20,1p8e12.3)')'explosion velocity',u(1,2)
c..   put in explosion as kinetic energy
         do k = 3, kexplode
            u(1,k) = u(1,2)
         enddo
c..   inner radius is static
         u(1,1) = 0.0d0
      endif

c..write revised model
      call ritef(note,nc)
      close(8)

      stop'successful termination, newomit written'

 2000 continue
      write(*,*)'omit: error in  input file: omit.in'
      stop'omit error'

      end


