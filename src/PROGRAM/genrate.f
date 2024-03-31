      program genrate

c..   reads "mods" different models for comparison
c..   uses pgplot

      implicit none

c..   dimenfile include definition of dimensions for arrays,
c..   and must come first
      include 'dimenfile'
      include 'comod'
      include 'compu'
      include 'czone'
      include 'cgen'
      include 'cconst'
      include 'crate'
      include 'cburn'
      include 'caeps'
      include 'comlink'
      include 'cdtnuc'
      include 'ceoset'

      integer*4 i, k, n, kz, nc, kworst, linesty
      integer*4 pgbeg
      integer*4 idebug, mbflag
      integer*4 kkold, imax
      integer*4 iabflag, iequil
      integer*4 forcex,forcey

      real*8    xmsol, sum, sumx, yions, yee, ytot
      real*8    sume, diffold, difzold,xmet,ymet,zmet,epsilon
      real*8    sumall(kdm)

      real*8    yeq(ndim),enc,t9,rho
      real*8    fact
      real*8    epp1,epp2,epp3,ecn,e3al,sum1,dlog,tlog
      real*8    xhyd,xhe4,xhe3,xc12,xc13,xn14,xo16,xo18
      real*8    jnbrate(13),jnbfkt(13)
      real*8    bigdx, dtover
      integer*4 nbigdx
cccccccccccccccccccccccccc

      real*4    xx(kdm),yy(kdm),xmin,xmax,ymin,ymax
      real*4    xmin0,xmax0,ymin0,ymax0
      real*4    xxx(2),yyy(2)

      real*4    tmass,fscale
      integer*4 lwidth

      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt
      character*3 cvarx, cvary

      character*9 cmod
      character*12 cxlabel, cylabel
      character*11 cfnum
      character*20 ctlabel
      character*5  cnum
      character*5 device

      logical tobe

c..   note that weakscr gets mapped into weakscreening in exportenergy.f
      real*8 weakscr, Flux
      Common/Fluxes/Flux(10),weakscr

      data txtxt/'  '/

c..define for solven, dtnuc
      data delchi/0.2d0/, chimin/1.0d-6/, fdtn/1.414d0/,
     1 fdysum/2.0d-2/
      data fscale /1.4/,lwidth/6/

c..bahcall limit for weak screening
      data weakscr/ 0.03d0/

      save
c-----------------------------------------------------------------
      it      = 0
      mbflag  = 0
c..explicit nucleosynthesis
      nburn   = 1
      ncymax  = 300
      ncytest = ncymax - 10

      open(2,file='genrate.in')
      open(3,file='genrate.out')

c..   input parameters for model modification from gentran.in

      read (2,*)text
      write(*,*)text
c..   dummy read
      read (2,*)text
      write(*,*)text

      read (2,*)text
      write(*,*)text


      read (2,*)txt,labl,device
      write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,linesty
      write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,lwidth
      write(*,*)txt,txtxt,labl,lwidth
      labl1 = 'lwidth'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,fscale
      write(*,*)txt,txtxt,labl,fscale
      labl1 = 'fscale'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,cmod
      write(*,*)txt,txtxt,labl,cmod
      labl1 = 'cmod'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,kz
      write(*,*)txt,txtxt,labl,kz
      labl1 = 'kz'
      if( labl .ne. labl1 )goto 2000

c..   smallest abundance allowed (epsilon=xmin/xmax)
      read (2,*)txt,labl,epsilon
      write(*,*)txt,txtxt,labl,epsilon
      labl1 = 'epsilon'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,iabflag
      write(*,*)txt,txtxt,labl,iabflag
      labl1 = 'iabflag'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,iequil
      write(*,*)txt,txtxt,labl,iequil
      labl1 = 'iequil'
      if( labl .ne. labl1 )goto 2000

      if( iequil .lt. 0 .or. iequil .gt. 2 )then
         write(*,*)'GENRATE: error in iequil, see genrate.in'
         write(*,*)'iequil = ',iequil
         stop'GENRATE: iequil error'
      endif

      read (2,*)txt,labl,jnb
      write(*,*)txt,txtxt,labl,jnb
      labl1 = 'jnb'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,weakscr
      write(*,*)txt,txtxt,labl,weakscr
      labl1 = 'weakscr'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,dtover
      write(*,*)txt,txtxt,labl,dtover
      labl1 = 'dtover'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,mbflag
      write(*,*)txt,txtxt,labl,mbflag
      labl1 = 'mbflag'
      if( labl .ne. labl1 )goto 2000
      if( mbflag .lt. 0 .or. mbflag .gt. 1 )then
         write(*,*)'GENRATE: mbflag out of range ',mbflag
         stop'mbflag'
      endif

      read (2,*)txt,labl,forcex
      write(*,*)txt,txtxt,labl,forcex
      labl1 = 'forcex'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,forcey
      write(*,*)txt,txtxt,labl,forcey
      labl1 = 'forcey'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
      write(*,*)text

c..   BUILD NETWORK
      call build(0)

c..   check network and model consistency (for abundance sanity)
      inquire(file='net.rc',exist=tobe)
      if( tobe )then
         open(9,file='net.rc',form='formatted',status='old')
 1001    continue
         read(9,'(3i5,a5,0pf10.4,1pe12.4)',end=1000)n,nz(n),nn(n),
     1        cnuc(n),qex(n),solarx(n)
         lz(n)  = nz(n)
         ln(n)  = nn(n)
         xid(n) = cnuc(n)
         if( nz(n) .eq. 2 .and. nn(n) .eq. 2 )then
            read(9,'(10i5)')nucp
            goto 1000
         endif
         goto 1001
      else
         write(*,*)'no file net.rc in current directory'
         stop' genrate error 1'
      endif
 1000 continue
      close(9)

      if( n .ne. netsize )then
         write(*,*)'net.rc inconsistent with dimenfile'
         write(*,*)n,' entries in net.rc'
         write(*,*)nnuc,' nuclei of network in dimenfile'
         write(*,*)'resetting to net.rc abundances'
         stop'genrate: network data'
      endif


c..   formatted read in new format
      inquire(file=cmod,exist=tobe)
      if( tobe )then
         open(8,file=cmod)
         rewind 8
         nc = 1
c         call nreadf(note,nc)
         call readf(note,nc)
         close(8)
      else
         write(*,*)'NO FILE ',cmod,' IN THIS DIRECTORY'
         stop'no model'
      endif

c..   override time step
      if( dtover .gt. 0.0d0 )then
         write(*,*)'Overriding timestep ',dth(2),' with ',dtover
         dth(2) = dtover
      endif

c..   define mass coordinate
      dmi(1)    = dmh(2)
      dmi(kk+1) = dmh(kk+1)
      do k = 2, kk
         xm(k)  = xm(k-1) + dmh(k)
         dmi(k) = 0.5d0*( dmh(k+1) + dmh(k) )
      enddo
      xm(kk+1) = xm(kk) + dmh(kk+1)
      xmsol    = xm(kk+1)/sol

      write(*,'(a8,a10,a8,0pf10.5,a13,a10,i5,a8)')
     1     'Model =',cmod,'MASS =',xmsol,'solar masses',
     2     'netsize =',netsize,'nuclei'


c..   test normalization
      diffold = 0.0d0
      difzold = 0.0d0
      kworst  = 0
      do k = 2, kk
         sum = 0.0d0
         do n = 1, nnuc
c            sum = sum + x(n,k)*xa(n)
            sum = sum + x(n,k)
         enddo
         sumall(k) = sum
         sume      = 0.0d0
         do n = 1, nnuc
c            sume  = sume + x(n,k) * dble( nz(n) )
            sume  = sume + x(n,k) * dble( nz(n) )/xa(n)
         enddo
         yee   = x(nnuc+1,k)
cccccccccccccccc
         yions = sume
         ytot = yee + yions
         if( abs( sum-1.0d0 ) .gt. abs(diffold) )then
            kworst = k
            diffold = sum  - 1.0d0
            difzold = sume - yee
         endif
      enddo

      if( abs(diffold) .gt. 1.0d-10 .or.
     1     abs(difzold) .gt. 1.0d-10 )then
         write(*,'(a10,i5,a10,1pe12.3,a10,1pe12.3)')
     1        'kworst',kworst,'diff X',diffold,'diff Z',difzold
      endif

c..   renormalize
      do k = 2, kk
         do n = i, nnuc
            x(n,k) = x(n,k)/sumall(k)
         enddo
      enddo
c..   check
      do k = 2, kk
         sum = 0.0d0
         do n = 1, nnuc
            sum = sum + x(n,k)
         enddo
         if( abs( sum - 1.0d0 ) .gt. 1.0d-7 )then
            write(*,'(i6,1p8e12.3)')k,sum-1.0d0
            stop'GENRATE: bad renormalization'
         endif
      enddo
      write(*,*)'abundances renormalized in all zones '

c..   define Ye
      do k = 2, kk
         x(nnuc+1,k) = 0.0d0
         do n = 1, nnuc
            x(nnuc+1,k) = x(nnuc+1,k) + x(n,k)*dble(nz(n))/xa(n)
         enddo
      enddo
c..   define metallicity from outer zone kk
      zmet = 0.0d0
      xmet = 0.0d0
      ymet = 0.0d0
      do n = 1, nnuc
         if( nz(n) .eq. 1 )then
            xmet = xmet + x(n,kk)
         elseif( nz(n) .eq. 2 )then
            ymet = ymet + x(n,kk)
         else
            zmet = zmet + x(n,kk)
         endif
      enddo

      write(*,'(/3(a5,0pf10.6),a5,1pe10.2)')
     1     'X', xmet, 'Y', ymet, 'Z', zmet,'err', xmet+ymet+zmet-1.0d0

      do n = 1, nnuc
         if( nz(n) .eq. 6 .and. nn(n) .eq. 6 )then
            lc12 = n
         elseif( nz(n) .eq. 6 .and. nn(n) .eq. 7 )then
            lc13 = n
         elseif( nz(n) .eq. 7 .and. nn(n) .eq. 7 )then
            ln14 = n
         elseif( nz(n) .eq. 8 .and. nn(n) .eq. 8 )then
            lo16 = n
         elseif( nz(n) .eq. 8 .and. nn(n) .eq. 10 )then
            lo18 = n
         elseif( nz(n) .eq. 2 .and. nn(n) .eq. 1 )then
            lhe3 = n
         endif
      enddo
      write(*,'(4(a5,0pf10.6))')
     1     cnuc(nnuc-2),x(nnuc-2,kz),
     2     cnuc(nnuc-1),x(nnuc-1,kz),
     3     cnuc(nnuc),x(nnuc,kz),
     4     cnuc(lhe3),x(lhe3,kz)
      write(*,'(4(a5,0pf10.6))')
     1     cnuc(lc12),x(lc12,kz),
     2     cnuc(lc13),x(lc13,kz),
     3     cnuc(ln14),x(ln14,kz),
     4     cnuc(lo16),x(lo16,kz),
     5     cnuc(lo18),x(lo18,kz)



c..   for bahcall, testing consistency
      if( jnb .ne. 0 .and. t(1,kz) .lt. 5.0d7 )then
c..   bahcall values for T< 5e Kelvin
         sumx = 0
         sume = 0
         sum  = 0
         do n=1, nnuc
            sumx = sumx + x(n,kz)/xa(n)
            sume = sume + x(n,kz)*dble( nz(n) )/xa(n)
            sum  = sum  + x(n,kz)*dble( nz(n)*nz(n) )/xa(n)
         enddo
         write(*,*)'averages of Yi, Yi*Zi, Yi*Zi**2'
         write(*,'(3(a5,1pe15.6))')'Y',sumx,'YZ',sume,'YZZ',sum

      else
         jnb = 0
         write(*,'(a30,1pe12.3,a24,i2)')
     1        'Too hot for Bahcall rates, T =',
     1        t(1,kz),'Kelvin, setting jnb to',jnb
      endif

c..   0.017605 is for 37 nucleus network,
c..   full network would be 0.0191 for anders-grevesse
c..   full network would be 0.0148 for lodders

      write(*,'(1x,a12,1pe12.3,a8)')'Stellar age:',time/secpy,' years'

      if( iabflag .ne. 0 )then
c..   abundance table
         write(*,'(a5,a5,2a12,a14)')'n','nuc','X','sol X','X/sol X'
         do i = 1, nnuc
            if( solarx(i) .gt. 0.0d0 )then
               write(*,'(i5,a5,1p2e12.3,0pf14.6)')
     1              i,cnuc(i),x(i,kz),solarx(i),
     2              x(i,kz)/solarx(i)
            endif
         enddo
         write(*,'(a5,a5,2a12,a10)')'n','nuc','X','sol X','X/sol X'
      endif

      call state(kz,kz,1)

c..   define original abundance for burn.f
      do k = 2, kk
         do i = 1, ndim
            xold(i,k) = x(i,k)
         enddo
      enddo

c..require only one network solution (no derivatives d/dT, d/dV)
c      mbflag = 1

      call burn(1,mbflag,0,kz,kz)

      call select(x,t,v,epsilon,jnbfkt,kz,1,iabflag,iequil,jnb)

      write(*,'(/2(a10,1pe12.4,a10))')
     1     's(5,k) ',s(5,kz), ' ergs/g-s', 'dth(2)' , dth(2), ' seconds'

      write(*,'(3(a10,1pe12.4))')
     1     'ss(k)',ss(kz), 'dlns/dlnT', sa(kz),'dlns/dlnV',sb(kz)
      write(*,'(3(a10,1pe12.4))')
     1     'snu(k)',snu(kz), 'dlns/dlnT', snua(kz),'dlns/dlnV',snub(kz)
      write(*,'(a15,i5,a30,i5)')'burn flag',mbflag,
     1     'Number of leqs solutions ',ncyc(kz)
cccccccccccccccc

c..   implicit estimate
      write(*,'(3(a10,1pe12.4,a10))')
     1     'aeps ',   aeps(kz),    ' ergs/g-s',
     2     'aenu ',   aenu(kz),    ' ergs/g-s',
     3     'aenucnu ',aenucnu(kz), ' ergs/g-s',
     4     'nuclear', aeps(kz)+aenucnu(kz),  ' ergs/g-s',
     5     'total',   aeps(kz)+aenu(kz)+aenucnu(kz),  ' ergs/g-s'


c      if( mbflag .eq. 0 )then
c..   derivatives
c         if( aeps(kz) .eq. 0.0d0 )then
c            write(*,'(2(a18,1pe12.4,a10))')
c     1           'aepst ',aepst(kz)
c            write(*,'(2(a18,1pe12.4,a10))')
c     1           'aepsv ',aepsv(kz)
c         else
c            write(*,'(2(a18,1pe12.4))')
c     1           'aepst ',aepst(kz),'dlneps/dlnT',
c     2           aepst(kz)/aeps(kz)*t(1,kz)
c            write(*,'(2(a18,1pe12.4))')
c     1           'aepsv ',aepsv(kz),'dlneps/dlnV',
c     2           aepsv(kz)/aeps(kz)*v(1,kz)
c         endif
c      endif


      if( jnb .ne. 0 )then
c..   john bahcall's exportenergy
c     call Energy(EPP1,EPP2,EPP3,ECN,E3AL,SUM1,DL,TL,
c     $     X,Y,XHe3,XC12,XC13,XN14,
c     $     XO16,xo18,IU)
         dlog = log10( rhoz(kz) )
         tlog = log10( t(1,kz) )

         xhyd = x(nnuc-1,kz)*xa(nnuc-1)
         xhe4 = x(nnuc,kz)*xa(nnuc)
         do n = 1, nnuc-2
            if( nz(n) .eq. 2 .and. nn(n) .eq. 1 )then
               xhe3 = x(n,kz)*xa(n)
            elseif(  nz(n) .eq. 6 .and. nn(n) .eq. 6 )then
               xc12 = x(n,kz)*xa(n)
            elseif(  nz(n) .eq. 6 .and. nn(n) .eq. 7 )then
               xc13 = x(n,kz)*xa(n)
            elseif(  nz(n) .eq. 7 .and. nn(n) .eq. 7 )then
               xn14 = x(n,kz)*xa(n)
            elseif(  nz(n) .eq. 8 .and. nn(n) .eq. 8 )then
               xo16 = x(n,kz)*xa(n)
            elseif(  nz(n) .eq. 8 .and. nn(n) .eq. 10 )then
               xo18 = x(n,kz)*xa(n)
            endif
         enddo

         call Energy(EPP1,EPP2,EPP3,ECN,E3AL,SUM1,DLog,TLog,
     $        Xhyd,xhe4,XHe3,XC12,XC13,XN14,
     $        XO16,jnbfkt,kz)

      endif


      write(*,*)'NUCLEON FRACTION CHANGES'
      write(*,'(5(a5,1pe11.3))')
     1     ( cnuc(i),xd(i,kz)*dth(2), i=1,netsize)
      sum = 0.0d0
      sumx = -1.0d0
      sume = 0.0d0
      do i = 1, netsize
         sum = sum + xd(i,kz)*dth(2)
         sumx = sumx + x(i,kz)
         sume = sume + xd(i,kz)*dble( nz(i) )/xa(i)
      enddo

      write(*,'(4(a14,1pe12.4))')'sum del x',sum,'sumx-1',sumx,
     1     'sumdz',sume,'sumdz*dth',sume*dth(2),
     2     'yeold',xold(nnuc+1,kz), 'ye',x(nnuc+1,kz)

      bigdx  = 0.0d0
      nbigdx = 0
      do i=1,netsize
         if( abs(xd(i,kz)*dth(2) ) .gt. abs(bigdx) )then
            bigdx = xd(i,kz)*dth(2)
            nbigdx = i
         endif
      enddo
      if( nbigdx .gt. 0 .and. nbigdx .le. netsize )then
         write(*,'(a25,a5,1p8e12.3)')'fastest changing nucleus',
     1        cnuc(nbigdx),bigdx
      endif

      write(*,*)
     1     'BIG FRACTIONAL CHANGES BY NUMBER (Y>1e-20, dlnY>0.01)'
      write(*,'(a5,a5,8a12)')'k','cnuc','dlnY','Yold','Ynew',
     1     'dY','dY/dt'
      do i =1, netsize
         if( x(i,kz) .gt. 1.0d-20 )then
            fact = dth(2)*xd(i,kz)/x(i,kz)
            if( abs( fact) .gt. 1.0d-2 )then
               write(*,'(i5,a5,1p8e12.3)')i,cnuc(i),fact,
     1              x(i,kz),x(i,kz)+dth(2)*xd(i,kz),
     2              dth(2)*xd(i,kz),xd(i,kz)
            endif
         endif
      enddo

      write(*,*)'printing nuclei (array nucp)'
      write(*,'(10(i5,a5))')(nucp(i),cnuc(nucp(i)),i=1,nucpg)

c      t9 = t(nc,kz)*1.0d-9
c      rho = rhoz(kz)
c      do i = 1, netsize+1
c         yeq(i) = x(i,kz)
c      enddo
c      call sqse(t9,rho,yeq,enc,kz)
c      write(*,'(5(a5,1p2e10.2,2x))')(cnuc(i),yeq(i),x(i,kz),i=1,netsize)
c      write(*,'(1p10e10.2)')enc,t9,rho,ss(kz)
c
c      t9 = t(nc,kz)*1.0d-9
c      rho = rhoz(kz)
c      do i = 1, netsize+1
c         yeq(i) = x(i,kz)
c      enddo
c      call nse(t9,rho,yeq,enc,kz)
c      write(*,'(5(a5,1p2e10.2,2x))')(cnuc(i),yeq(i),x(i,kz),i=1,netsize)
c      write(*,'(1p10e10.2)')enc,t9,rho,ss(kz)

      kkold = kk

c..   x-coordinate
      xmin = 0
      xmax = max0( 24, nn(nnuc-3) )

      do k = 1, nnuc
         xx(k) = float( nn(k) )
      enddo
      cxlabel = 'N'
c..   y-coordinate
      ymin = 0
      ymax = xmax*0.75
c..   preserve aspect ratio
      if( ymax .lt. float( nz(nnuc-3) ) )then
         ymax = nz(nnuc-3)
         xmax = ymax/0.75
      endif
      do k = 1, nnuc
         yy(k) = float( nz(k) )
      enddo
      cylabel = 'Z'
      kk = nnuc

c..override axis size
      if( forcex .gt. 0 .and. forcey .gt. 0 )then
         xmax = forcex
         ymax = forcey
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
c..   roman font=2
      call PGSCF(1)
c..   font scaling (size)
      call PGSCH(fscale)
c..   set line width (1-201)
      call PGSLW(lwidth)
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

c..   set color(white=1)
      call PGSCI(1)

c..   labels
      tmass   = xmsol
      call ftoc(tmass,cnum)
      ctlabel = 'M='//cnum//' '//cmod
      call PGTEXT(xmin,ymax,ctlabel)

      call modflg(kz,cnum)
      ctlabel = 'zone'//cnum
      call PGTEXT(xmin+(xmax-xmin)*0.35,ymax,ctlabel)

      tmass = t(nc,kz)
      rewind 10
      write(10,'(1pe11.3)')tmass
      rewind 10
      read(10,'(a11)')cfnum
      ctlabel = 'T='//cfnum//' K'
      call PGTEXT(xmin,ymax-0.05*(ymax-ymin),ctlabel)

      tmass = 1.0d0/v(nc,kz)
      rewind 10
      write(10,'(1pe11.3)')tmass
      rewind 10
      read(10,'(a11)')cfnum
      ctlabel = 'rho='//cfnum//' g/cc'
      call PGTEXT(xmin,ymax-0.1*(ymax-ymin),ctlabel)

      tmass = s(5,kz)
      rewind 10
      write(10,'(1pe11.3)')tmass
      rewind 10
      read(10,'(a11)')cfnum
      ctlabel = 's='//cfnum//' erg/gs'
      call PGTEXT(xmin,ymax-0.15*(ymax-ymin),ctlabel)
      

      do k = 1, nreac
         if( in(k) .gt. 0 )then
            xxx(1) = x1(1,k)
            xxx(2) = x1(2,k)
            yyy(1) = y1(1,k)
            yyy(2) = y1(2,k)
            call PGSLW(ltyp(k)+1)
            call PGSCI(ltyp(k)+2 )
            call PGSLS(ldot(k))
            call PGLINE(2,xxx,yyy)
         endif
      enddo

c..   reset to default
      call PGSLW(3)
      call PGSCI(1)
      call PGSLS(1)

c..   draw points for nuclei
      call PGPT(kk,xx,yy,-4)

      if( device .eq. '/ps' .or. device .eq. '/cps' )then
c..   direct hardcopy, so skip pause
      else
c..   xwindow
         pause
      endif

      close(10)

      call PGEND

      stop'successful termination'

 2000 continue

      write(*,*)'genrate: error in  input file: genrate.in'
      stop'genrate error'

      end


