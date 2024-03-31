      program apsides

c..   evaluates apsidal motion structure constant k2
c..   uses new gen format
c..   uses zeta = 3 - eta variable for better asymptotic behavior
c..   eta = d ln c2 / d ln r
c..   where c2 is the quadrapole coef. in the distorting potential
c..   m. schwarzschild, structure and evolution of stars,
c..   princeton university press, 1957, p.151

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'czone'
      include 'cgen'
      include 'cconst'

      include 'compu'
      include 'compua'
      include 'ceoset.h'
      include 'ceos'
      include 'cburn'

      include 'cenv'
      include 'csurface'
      include 'cnabla'


      real*8 xmsol,zeta2(kdm),apside2,uu(kdm),rr(kdm)
      real*8 vv(kdm),r5k1(kdm)
      real*8 rhobar(kdm),pbar(kdm),dri(kdm)
      real*8 top,bottom
      real*8 sumx
      real*8 ue(kdme),ve(kdme),dre(kdme),re(kdme),zetae(kdme),
     1     r5ke(kdme)

      real*4    xx(kdm),yy(kdm),xmin,xmax,ymin,ymax
      real*4    xmin0,xmax0,ymin0,ymax0
      real*4    xminf,xmaxf,yminf,ymaxf
      real*4    tmass

      integer*4 pgbeg, linesty, journal

      integer*4 k, nc, j, igraph, mm,m1,m2

      character*72 text
      character*44 txt
      character*8  labl, labl1

      character*20 cxlabel, cylabel
      character*80 ctlabel
      character*5  cnum,cnum1
      character*1  cmm(3)

      character*7  cmod
      character*5 device

      data linesty/1/
      data alphaml/2.0d0/, epsilon/0.4d0/
c-----------------------------------------------------------------

      open(2,file='apsides.in')

c..   input parameters for model modification from gentran.in

      read (2,*)text
c     write(*,*)text
c..   dummy read
      read (2,*)text
c     write(*,*)text

      read (2,*)text
c     write(*,*)text

c..   pgplot graphics output device (eg, /xwin, /cps)
      read (2,*)txt,labl,device
c     write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

c..   model identity
      read (2,*)txt,labl,cmod
c     write(*,*)txt,txtxt,labl,cmod
      labl1 = 'cmod'
      if( labl .ne. labl1 )goto 2000

c..   lines or points?
      read (2,*)txt,labl,linesty
c     write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

c..   id header? no=1
      read (2,*)txt,labl,journal
c     write(*,*)txt,txtxt,labl,journal
      labl1 = 'journal'
      if( labl .ne. labl1 )goto 2000

c..   force graphical limits?
      read (2,*)txt,labl,xminf
c     write(*,*)txt,txtxt,labl,xminf
      labl1 = 'xminf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmaxf
c     write(*,*)txt,txtxt,labl,xmaxf
      labl1 = 'xmaxf'
      if( labl .ne. labl1 )goto 2000

      if( xminf .ne. xmaxf )then
         write(*,*)'      x override'
      endif

      read (2,*)txt,labl,yminf
c     write(*,*)txt,txtxt,labl,yminf
      labl1 = 'yminf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymaxf
c     write(*,*)txt,txtxt,labl,ymaxf
      labl1 = 'ymaxf'
      if( labl .ne. labl1 )goto 2000

      if( yminf .ne. ymaxf )then
         write(*,*)'      y override'
      endif

      read (2,*)txt,labl,igraph
c     write(*,*)txt,labl,igraph
      labl1 = 'igraph'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
c     write(*,*)text
      close(2)

c..   get the model
c..   formatted read in new format
      open(8,file=cmod)
      rewind 8
      nc = 1
      call nreadf(note,nc)
      close(8)

c..   define mass coordinate
      dmi(1)    = dmh(2)
      dmi(kk+1) = dmh(kk+1)
      do k = 2, kk
         xm(k)  = xm(k-1) + dmh(k)
         dmi(k) = 0.5*( dmh(k+1) + dmh(k) )
      enddo
      xm(kk+1) = xm(kk) + dmh(kk+1)
      xmsol = xm(kk+1)/sol

      write(*,'(a8,a7,a10,0pf12.6,a10,i5,a8)')'Model  ',cmod,
     1     'MASS', xmsol, 'nnuc', nnuc,'nuclei'

      write(*,'( 3(a12,0pf12.6) )')
     1     'envelop mass',dmh(kk+1)/sol,'r(env)/R',
     1     r(1,kk)/r(1,kk+1),'R/Rsun',r(1,kk+1)/solrad

c..   get microphysics
c..   set values at start of time step interval
      nc   = 1
c..   check nucleon conservation
      do k = 2, kk
         sumx = 0.0d0
         do j = 1, ndim-1
            sumx = sumx + x(j,k)*dble(lz(j)+ln(j))
         enddo
         if( abs(sumx-1.0d0) .gt. 1.0d-4 )then
c..   large error
            write(*,*)'bad abundance error ',sumx,sumx-1.0d0
            stop'apsides: abnormal'
         elseif(  abs(sumx-1.0d0) .gt. 1.0d-10 )then
c..   renormalize if small error
            do j = 1, ndim-1
               x(j,k) = x(j,k)/sumx
            enddo
c..   check new scaled values
            sumx = 0.0d0
            do j = 1, ndim-1
               sumx = sumx + x(j,k)*dble(lz(j)+ln(j))
            enddo
            if( abs(1.0d0-sumx) .gt. 1.0d-10 )then
               write(*,*)'apsides init error: sumx-1 ',sumx-1.0d0,k
               stop'apsides sumx'
            endif
c..   new Ye
            x(ndim,k) = 0.0d0
            do j = 1, ndim-1
               x(ndim,k) = x(ndim,k) + x(j,k)*dble(lz(j))
            enddo
         endif
      enddo

c..   fill kk eos arrays
      call state(kk,kk,nc)

c..   initialize additional variables
      do k = 1, kk+1
         a(k) = pi4*r(1,k)**2
         g(k) = grav*xm(k)/r(1,k)**2
      enddo
c..   see data statement for alphaml, epsilon

c..   get envelope

      call fitenv

c..............................................................
c..   apsidal constant computation
c..   interior
      vv(1) = 0.0d0
      rr(1) = 0.0d0
      uu(1) = 3.0d0
      do k = 2, kk+1
         rhobar(k) = (1.0d0/v(1,k) + 1.0d0/v(1,k+1) )*0.5d0
         pbar(k)   = (p(1,k) + p(1,k+1))*0.5d0
         uu(k)     = pi4 * r(1,k)**3 * rhobar(k) / xm(k)
         vv(k)     = grav*xm(k)*rhobar(k)/( r(1,k) * pbar(k))
         rr(k)     = r(1,k)/r(1,kk+1)
      enddo

      do k = 2, kk
         dri(k) = 2.0d0*( r(1,k) - r(1,k-1) )/(r(1,k)+r(1,k-1))
      enddo
      dri(1) = dri(2)

      zeta2(1) = 3.0d0

      do k = 2, kk
c..   second order differences as used in henyey solution in hstat.f
         bottom  = dri(k)*( 2.5d0-zeta2(k-1) + (uu(k)+uu(k-1)/2.0d0) )
         top     = dri(k)*( zeta2(k-1)**2 -5.0d0*zeta2(k-1)
     1        + (uu(k)+uu(k-1))*( 4.0d0 - zeta2(k-1) )  )
         zeta2(k) = zeta2(k-1) + top/( 1.0d0 + bottom )
         r5k1(k) = zeta2(k)/(10.0d0-2.0d0*zeta2(k))*rr(k)**5
      enddo

      write(*,*)'interior contribution alone: '

      apside2 = zeta2(kk)/(10.0d0 - 2.0d0*zeta2(kk))*rr(kk)**5

      write(*,'(3(a14,1pe12.3))')'k2',apside2,'log k2',dlog10(apside2), 
     1     'zeta2(kk)',zeta2(kk)

      write(*,'(3(a14,1pe12.3))')'R/Rsun',r(1,kk+1)/solrad,
     1     'M/Msol',xmsol,
     2     'k2*R**5',apside2*(r(1,kk+1)/solrad)**5
      write(*,*)' '

c..   envelope
c..   map envelope to ascending index order with increasing radius
      do j = 1, jmaxz
         k = jmaxz - j + 1
         ue(k) = pi4 * zr(j)**3 * zrho(j) /(xm(kk+1)+zm(j))
         ve(k) = grav*( xm(kk+1) + zm(j) )*zrho(j)/( zr(j) * zp(j) )
         re(k) = zr(j)/r(1,kk+1)
      enddo
      do j = 2, jmaxz
         dre(j) = 2.0d0*( re(j) - re(j-1) )/( re(j) + re(j-1))
      enddo
      dre(1) = dre(2)

      zetae(1) = zeta2(kk)
      r5ke(1) = r5k1(kk)
      do k = 2, jmaxz
         bottom = dre(k)*( 2.5d0 - zetae(k-1) + ue(k) )
         top    = dre(k)*(zetae(k-1)**2 -5.0d0*zetae(k-1)
     1        + ue(k)*2.0d0*(4.0d0 - zetae(k-1)) )
         zetae(k) = zetae(k-1) + top/(1.0d0 + bottom)
         r5ke(k) = zetae(k)/(10.0d0-2.0d0*zetae(k))*re(k)**5
      enddo

      write(*,*)'envelope contribution added:'
      apside2 = zetae(jmaxz)/(10.0d0 - 2.0d0*zetae(jmaxz))
     1     *re(jmaxz)**5

      write(*,'(3(a14,1pe12.3))')'k2',apside2,'log k2',dlog10(apside2), 
     1     'zetae(jmaxz)',zetae(jmaxz)

      write(*,'(3(a14,1pe12.3))')'R/Rsun',r(1,kk+1)/solrad,
     1     'M/Msol',xmsol,
     2     'k2*R**5',apside2*(r(1,kk+1)/solrad)**5


c..   graphics section................................................
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

      if( igraph .eq. 0 )then
         xmin = rr(1)
         xmax = re(jmaxz)

         ymin = r5k1(1)
         ymax = r5ke(jmaxz)
         do k = 1, jmaxz
            ymin = dmin1( r5ke(k) , dble(ymin) )
            ymax = dmax1( r5ke(k) , dble(ymax) )
         enddo
         cylabel = 'k\\d2\\u(r/R)\\u5'

         cxlabel = 'r/R'
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

         do k = 1, kk
            xx(k) = rr(k)
            yy(k) = r5k1(k)
         enddo

c..   set color
         call PGSCI(2)
c..   draw curve
         if( linesty .eq. 0 )then
            call PGLINE(kk,xx,yy)
         else
            call PGPT(kk,xx,yy,20)
         endif


         do k = 1, jmaxz
            xx(k) = re(k)
            yy(k) = r5ke(k)
         enddo

c..   set color
         call PGSCI(4)
c..   draw curve for envelope as points
         call PGPT(jmaxz,xx,yy,20)

         if( journal .eq. 0 )then
c..   label
c..   set color
         call PGSCI(1)
         tmass   = xm(kk+1)/sol
         call ftoc(tmass,cnum)
         cxlabel = 'M= '//cnum
         call PGTEXT(xmin+0.1*(xmax-xmin),ymax-0.1*(ymax-ymin),
     1     cxlabel)
         write(*,*)cxlabel

         tmass   = r5ke(jmaxz)*(r(1,kk+1)/solrad)**5
c..convert to character string
         mm = alog10(tmass)
         tmass = tmass/10.0**(mm)
c..get rescaled number
         call ftoc(tmass,cnum)
c..get exponent in form +12, -00, etc.
         m2 = abs(mm)/10
         m1 = abs(mm)-m2*10
         call itoa(m1,cmm(3))
         call itoa(m2,cmm(2))
         if( mm .gt. 0 )then
            cmm(1) = ' '
         else
            cmm(1) = '-'
         endif


         cxlabel = 'kR\\u5\\d= '
     1   //cnum//'e'//cmm(1)//cmm(2)//cmm(3)
         call PGTEXT(xmin+0.1*(xmax-xmin),ymax-0.15*(ymax-ymin),
     1     cxlabel)

         tmass = r5ke(jmaxz)
         tmass = -alog10( tmass )
         call ftoc(tmass,cnum1)

         ctlabel = cmod//'  '//'log k2 = -'//cnum1

           call PGMTXT('T',2.0,0.5,0.5,ctlabel)
         endif

      else
c..   UV plane
         xmin = uu(1)
         xmax = xmin
         ymin = vv(1)
         ymax = ymin
         do k = 1, kk
            xmin = dmin1( uu(k) , dble(xmin) )
            xmax = dmax1( uu(k) , dble(xmax) )
            ymin = dmin1( vv(k) , dble(ymin) )
            ymax = dmax1( vv(k) , dble(ymax) )
         enddo
         do k = 1, jmaxz
            xmin = dmin1( ue(k) , dble(xmin) )
            xmax = dmax1( ue(k) , dble(xmax) )
            ymin = dmin1( ve(k) , dble(ymin) )
            ymax = dmax1( ve(k) , dble(ymax) )
         enddo
         ymax = amin1( ymax, 10. )

         cylabel = 'V'
         cxlabel = 'U'

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

         do k = 1, kk
            xx(k) = uu(k)
            yy(k) = vv(k)
         enddo

c..   set color
         call PGSCI(2)
c..   draw curve
         if( linesty .eq. 0 )then
            call PGLINE(kk,xx,yy)
         else
            call PGPT(kk,xx,yy,20)
         endif


         do k = 1, jmaxz
            xx(k) = ue(k)
            yy(k) = ve(k)
         enddo

c..   set color
         call PGSCI(4)
c..   draw curve
         if( linesty .eq. 0 )then
            call PGLINE(jmaxz,xx,yy)
         else
            call PGPT(jmaxz,xx,yy,20)
         endif

c..   label
c..   set color
         call PGSCI(1)
         tmass = xm(kk+1)/sol
         call ftoc(tmass,cnum1)

         ctlabel = cmod//'  '//'M = '//cnum1
         call PGMTXT('T',2.0,0.5,0.5,ctlabel)

      endif

      pause

      close(10)
      
      call PGEND

      stop'successful termination'

 2000 continue

      write(*,*)'apsides: error in  input file: apsides.in'
      stop'apsides error'

      end



