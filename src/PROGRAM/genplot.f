      program genplot

c..   8-2-06
c..   reads and plots "mods" different models for comparison
c..   uses binary or formatted models, pgplot

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'czone'
      include 'cgen'
      include 'cconst'
      include 'cburn'

      include 'cenv'
      include 'compu'
      include 'csurface'
      include 'cnabla'
      include 'ceoset'
      include 'cdtnuc'

      integer*4 n, j

      integer*4 mods, imods, i, k, nc, kworst, linesty, lmods
      integer*4 pgbeg

      parameter( mods = 15 )

      real*8    xmsol(mods), sum, ytot,timm(mods)
      integer*4 nopeos(mods)
      real*8    sume, diffold, difzold,xmet,ymet,zmet

      real*8    xmass
      real*8    fint(4),xint(4),xxint,xmeint
      integer*4 ji(4)

      real*8    ff
      real*8 scr(kdm)

      real*4    xx(kdm),yy(kdm),xmin,xmax,ymin,ymax
      real*4    xmin0,xmax0,ymin0,ymax0
      real*4    ax(mods,kdm),ay(mods,kdm)
      real*4    axmin(mods),axmax(mods),aymin(mods),aymax(mods)
      real*4    xminf,xmaxf,yminf,ymaxf
      real*4    xlabel,ylabel

      real*4  xxa(kdm),yya(kdm)

c..   arrays to hold kk from model file (kkold) and effective kk
c..   if envelope is included (kkeff)
      integer*4 kkold(mods),kkeff(mods)

      real*4    tmass

      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt
      character*5 cvarx, cvary

      character*9  cmod(mods),cmodd
      character*12 cxlabel, cylabel
      character*20 ctlabel
      character*5  cnum
      character*5 device
      character*7 cdevice
      character*10 axistag

      logical tobe
      data txtxt/'  '/
      data delchi/0.2d0/, chimin/1.0d-4/, fdtn/1.414d0/,
     1     fdysum/2.0d-2/
c..   params.d variable used in fitenv.f
      data resid/1.0d-5/
c-----------------------------------------------------------------

      inquire(file='genplot.in',exist=tobe)
      if( tobe )then
         open(2,file='genplot.in',form='formatted',status='old')
      else
         write(*,*)'genplot: no file genplot2.in in directory'
         stop'genplot: no .in file error'
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
c..   use pgplot window 10 to avoid conflicts with TYCHO which 
c..   uses 1,2,3,4
      if( device .eq. '/xwin' )then
         cdevice = '10'//device
      else
c..   uses default
         cdevice = device
      endif

      read (2,*)txt,labl,linesty
      write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nopaleos
      write(*,*)txt,txtxt,labl,nopaleos
      labl1 = 'nopaleos'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,nopac
      write(*,*)txt,txtxt,labl,nopac
      labl1 = 'nopac'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nkscale
      write(*,*)txt,txtxt,labl,nkscale
      labl1 = 'nkscale'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,nburn
      write(*,*)txt,txtxt,labl,nburn
      labl1 = 'nburn'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,lmods
      write(*,*)txt,txtxt,labl,lmods
      labl1 = 'lmods'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,imods
      write(*,*)txt,txtxt,labl,imods
      labl1 = 'imods'
      if( labl .ne. labl1 )goto 2000

      if( imods .le. 0 .or. imods .gt. mods )then
         write(*,*)'imods ',imods
         goto 2000
      endif
c..   read cycle for model identities
      do i = 1, imods
         read (2,*)txt,labl,cmod(i)
         write(*,*)txt,txtxt,labl,cmod(i)
         labl1 = 'cmod'
         if( labl .ne. labl1 )then
            write(*,'(/4a10/)')'Expects',labl1,'got',labl
            goto 2000
         else
            inquire(file=cmod(i),exist=tobe)
            if( .not. tobe )then
               write(*,'(/5x,a8,a25)')cmod(i),
     1              ' is not in this directory,'
               write(*,'(a30/)')'Check genplot.in for typos.'
               stop'genplot: cmod'
            endif
         endif
      enddo

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

c..   dummy read
      read (2,*)text
      write(*,*)text

      close(2)

c..   BUILD NETWORK
      call build(0)

c..   size parameter for finite difference between envelopes
      epsilon = 0.005d0
      write(*,*)'setting epsilon ',epsilon

c..   construct low T opacity table lookups
      call kappinit

c..   loop to read models and extract xx, yy
      do i = 1, imods

c..   read of model files
         inquire(file=cmod(i),exist=tobe)
         if( tobe )then
c..   dummy for character checking
c..   default is flush left so first nonblank in (1:1)
            cmodd = cmod(i)
            if( cmodd(1:1) .eq. 'H' )then
c..   binary files
               write(*,*)'binary',cmodd(1:1)
cccccccccccccccccccccccccc
               open(8,file=cmod(i),form='unformatted')
               nc = 1

               call readb(note,nc)
               write(*,*)note

            else
c..   formatted files
               write(*,*)'formatted'
cccccccccccccccccccccccccccccccc
               open(8,file=cmod(i),form='formatted')
               rewind 8
               nc = 1

               read(8,'(a72)')note
               write(*,*)note
               rewind 8

               if( note(7:7) .eq. '7' .or. note(7:7) .eq. '8' )then
c..   newer formatted read
                  write(*,*)'CALLING READF'
                  call readf(note,nc)
               else
c..   older incomplete formatted read
                  write(*,*)'CALLING NREADF'
                  call nreadf(note,nc)
c..   get mixing length alphaml from header note in model file
                  write(*,*)note
                  open(22,file='scratch')
                  rewind(22)
                  write(22,'(a5)')note(35:39)
                  rewind(22)
                  read(22,'(g5.2)')alphaml
                  write(*,'(a30,5x,1pe12.3)')'alphaml is set to ',
     1                 alphaml
                  close(22)

                  if( alphaml .eq. 1.03d0 )then
                     write(*,*)'p? sequences'
                     uuml = 0.7d0
                  else
                     uuml = 1.0d0
                  endif
               endif

            endif
            close(8)

            zpop = zpop0

            write(*,'(3(a15,1pe12.3))')'alphaml',alphaml,'uuml',uuml,
     1           'zpop',zpop

            timm(i) = time
c..   define mass coordinate
            dmi(1)    = dmh(2)
            dmi(kk+1) = dmh(kk+1)
            do k = 2, kk
               xm(k)  = xm(k-1) + dmh(k)
               dmi(k) = 0.5d0*( dmh(k+1) + dmh(k) )
            enddo
            xm(kk+1) = xm(kk) + dmh(kk+1)
            xmsol(i) = xm(kk+1)/sol

            write(*,'(a10,a8,0pf10.5,a15,a10,i5,a8)')cmod(i),
     2           'MASS =',xmsol(i),
     1           ' solar masses,',' nnuc =',nnuc,' nuclei'
            write(*,'(a12,i5,2(a12,1pe12.3))')'modes',modes,
     1           'Menv/sol',dmh(kk+1)/sol,'epsilon',epsilon

c..   check network and model consistency (for abundance sanity)
            inquire(file='net.rc',exist=tobe)
            if( tobe )then
               open(9,file='net.rc',form='formatted',status='old')
 1001          continue
               read(9,'(3i5,a5,0pf10.4,1pe12.4)',end=1000)
     1              n,lz(n),ln(n),cnuc(n),qex(n),solarx(n)
               if( lz(n) .eq. 2 .and. ln(n) .eq. 2 )then
                  goto 1000
               endif
               goto 1001
            else
               write(*,*)'no file net.rc in current directory'
               stop' genplot2 error 1'
            endif
 1000       continue
            close(9)

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
                  sum = sum + x(n,k)
               enddo
               sume      = 0.0d0
               do n = 1, nnuc
                  sume  = sume + x(n,k) * dble( lz(n) )/xa(n)
               enddo
               ye(k)   = x(ndim,k)
               yion(k) = sume
               ytot    = ye(k) + yion(k)
               if( abs( sum-1.0d0 ) .gt. abs(diffold) )then
                  kworst = k
                  diffold = sum  - 1.0d0
                  difzold = sume - ye(k)
               endif

            enddo

            if( abs(diffold) .gt. 1.0d-6 .or.
     1           abs(difzold) .gt. 1.0d-6 )then
               write(*,'(a10,i5,a10,1pe12.3,a10,1pe12.3)')
     1              'kworst',kworst,'diff X',diffold,'diff Z',difzold
            endif


c..   flag which eos to use for state and fitenv
            nopeos(i) = nopaleos
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc

c..   fill kk eos arrays and define some abundance indices
            call state(2,kk,nc)

            xmet = x(ndim-2,kk) + x(ldeut,kk)
            ymet = x(ndim-1,kk) + x(lhe3,kk)
            zmet = 1.0d0 - xmet - ymet
            write(*,'(3(a5,0pf10.6),a10,1pe12.3)')
     1           'X', xmet, 'Y', ymet, 'Z', zmet,
     2           'z/zsun',
     3           zmet/0.0191d0

c..   define a few parameters
            ncymax  = 300
            ncytest = ncymax - 10
            it      = 0

            call burn(nc,0,0,2,kk)

c..gravity and area
            do k = 2, kk
               g(k) = grav*xm(k)/r(nc,k)**2
               a(k) = pi4*r(nc,k)**2
            enddo


            xmass = xm(kk)+dmh(kk+1)

c......................................................
c..   add envelope

            if( modes .eq. 2 )then

               write(*,*)'modes = ',modes,' with ',kk,' zones'

c..	l=1 to initialize fitenv for each model
               l = 1

               call fitenv(nc)
               
c..set convection variables

               call cinit(2,kk,nc)


               do j = 2,nvmax(3)-1
                  if( cvarx .eq. 'k' )then
                     xxa(j) = kk + nvmax(3) - j +1
                  elseif( cvarx .eq. 'xm' )then
                     xxa(j) = ( xm(kk+1) + vm(j,3) )/sol
                  elseif( cvarx .eq. 'r' )then
                     xxa(j) = vr(j,3)
                  endif
                  if( cvary .eq. 'm/r4' )then
                     yya(j) = (vp(j,3)-vp(j-1,3))/(vm(j,3)-vm(j-1,3))/(
     2                    0.5d0*(grav*(xm(kk+1)+vm(j,3))/
     3                    (pi4*vr(j,3)**4)
     3                    +grav*(xm(kk+1)+vm(j-1,3))/
     3                    (pi4*vr(j-1,3)**4) ) )
     4                    +1.0d0
                  endif

               enddo

c..............................
               write(*,'(a20,8a12)')'GRID of ENVELOPES','del L: 1',
     1              'del R: 2','ref.: 3','-del L: 4','-del R: 5',
     2              'value'
               write(*,'(a20,1p8e12.4)')'ref L',(vl(1,j),j=1,5),
     1              tl(2,kk)
               write(*,'(a20,1p8e12.4)')'ref R',(vr(1,j),j=1,5),
     1              r(2,kk+1)
               write(*,'(a20,1p8e12.4)')'ref Te',(vtem(1,j),j=1,5),
     1              t(2,kk+2)
               write(*,'(a20,1p8e12.4)')'ref P(kk)',(vp(nvmax(j),j),
     1              j=1,5),p(2,kk+1)
               write(*,'(a20,1p8e12.4)')'ref T(kk)',(vtem(nvmax(j),j),
     1              j=1,5),t(2,kk+1)
               write(*,'(a20,1p8e12.4)')'ref r(kk)',(vr(nvmax(j),j),
     1              j=1,5),r(2,kk)
               write(*,'(a20,1p8e12.4)')'ref rad-ad',
     1              (vnrad(nvmax(j),j)-vnad(nvmax(j),j),j=1,5),
     2              dnrad(kk)-dnad(kk)
               write(*,'(a20,1p8e12.4)')'ref cvel(kk)',
     1              (vvel(nvmax(j),j),j=1,5),h(kk),hp(kk)

               write(*,'(a20,1p8e12.4)')'ref dlnP(kk)',(
     1              grav*xm(kk)/(pi4*vr(nvmax(j),j)**4)/vp(nvmax(j),j)
     2              *dmh(kk)+ dlog( vp(nvmax(j),j) ),j=1,5), 
     3              dlog( vp(nvmax(3),3) ),
     4              dlog( p(nc,kk) ) 
               write(*,'(a20,1p8e12.4)')'ropal(jmax)',(
     1              vrho(nvmax(j),j)/(vtem(nvmax(j),j)*1.0d-6)**3,
     1              j=1,5),1.0d0/(v(nc,kk)*(t(nc,kk)*1.0d-6)**3)
c.........................................

c..   pad model with envelope values, ignore photosphere kk+2 and
c..   dummy envelope zone kk+1

c..   envelope mass xmeint
               xmeint = xm(kk+1)-xm(kk)
               write(*,*)kk,nvmax(3)
               
               dmh(kk+1) =  vm(nvmax(3)-1,3) - vm(nvmax(3),3)
               dmi(kk) = 0.5d0*( dmh(kk) + dmh(kk+1))
               write(*,'(8a13)')'xmkk','dmhkk','dmhkk+1',
     1              'dmikk','avedmhkk'
               write(*,'(1p8e13.5)')xm(kk), dmh(kk),dmh(kk+1),dmi(kk),
     1              0.5d0*( dmh(kk)+dmh(kk+1))          

               do k = kk+1, kk + nvmax(3) -1
c..   j counts from nvmax(3) (k = kk+1) down to 3
                  j = nvmax(3) -k +1 + kk
c..   define mass coordinate as interpolation variable
c..   interpolate to midpoint in mesh (between envelope boundaries vm)
                  if( j .ge. nvmax(3) )then
                     do n = 1, 4
c..   ji = nvmax(3) down to nvmax(3)-3 (toward surface at nvmax=1)
c..   this is near join boundary
                        ji(n) = nvmax(3)-n+1
                        xint(n) = vm(ji(n),3)
                     enddo
                     xxint   = 0.5d0*(xint(1)+xint(2))
                  elseif( j .lt. 3 )then
c..   ji = 4 down to 1
c..   this is near photosphere
                     do n = 1, 4
                        ji(n) = 5-n
                        xint(n) = vm(ji(n),3)
                     enddo 
                     xxint   = 0.5d0*(xint(3)+xint(4))
                  else
c..   ji = j down to j-3
c..   this is in between
                     do n = 1, 4
                        ji(n) = j-n+2
                        xint(n) = vm(ji(n),3)
                     enddo
                     xxint   = 0.5d0*(xint(2)+xint(3))
                  endif

c..sintr3 is cubic interpolation
c..sintr  is quadratic
                  do n = 1,4
                     fint(n) = vtem(ji(n),3)
                  enddo
                  call sintr3(xint,xxint,fint,t(nc,k))
                  do n = 1,4
                     fint(n) = 1.0d0/vrho(ji(n),3)
                  enddo
                  call sintr3(xint,xxint,fint,v(nc,k))

c                  do n = 1,4
c                     fint(n) = vp(ji(n),3)
c                  enddo
c..use to test accuracy of interpolation
c                  call sintr(xint,xxint,fint,p(nc,k))
c                  ff = p(nc,k)

c                  call sintr3(xint,xxint,fint,p(nc,k))
                  
c                  write(*,'(i5,1p8e12.4)')k,ff,p(nc,k),p(nc,k)/ff-1.0d0,
c        1   dlog(p(nc,k)),xxint
                  
c                  do n = 1,4
c                     fint(n) = ventr(ji(n),3)
c                  enddo
c                  call sintr3(xint,xxint,fint,entropy(k))
c                  do n = 1,4
c                     fint(n) = vsound(ji(n),3)
c                  enddo
c                  call sintr3(xint,xxint,fint,sound(k))
c                  do n = 1,4
c                     fint(n) = ve(ji(n),3)
c                  enddo
c                  call sintr3(xint,xxint,fint,e(nc,k))
                  
c                  do n = 1,4
c                     fint(n) = vak(ji(n),3)
c                  enddo
c                  call sintr3(xint,xxint,fint,ak(k))

c..   zone centers versus boundaries
c..   envelope index j-1 is boundary zone k
c..   interpolated zone center k-1/2 is between j and j-1
                  r(nc,k)  = vr(j-1,3)
                  h(k)     = vvel(j-1,3)
                  dnab(k)  = vnab(j-1,3)
                  dnrad(k) = vnrad(j-1,3)
                  dnad(k)  = vnad(j-1,3)
c..doux defined to be zero in envelope
                  scr(k)   = vnad(j-1,3)-vnab(j-1,3)
c..   constant velocity in envelope
                  u(nc,k)  = u(nc,kk)

                  if( dnrad(k) .gt. dnad(k) )then
                     ic(k) = 1
                  else
                     ic(k) = 0
                  endif
c..   update to new mass values in envelope
                  dmh(k) = vm(j-1,3) - vm(j,3)
                  dmi(k-1) = 0.5d0*(dmh(k)+dmh(k-1))
                  xm(k)  = xm(k-1) + dmh(k)
                  a(k)   = pi4*r(nc,k)**2
                  g(k)   = grav*xm(k)/r(nc,k)**2

                  tl(nc,k) = vl(j,3)
c..   zone centers versus boundaries
c..   envelope assumed homogeneous in composition
c..   so boundaries = centers
                  do n = 1, ndim
                     x(n,k) = x(n,kk)
                  enddo

               enddo
                            
c..run state over envelop interpolated values to fill in all
c..variables from state.f (p,e,ak,entropy,sound,...)
               call state(kk+1,kk+nvmax(3)-1,nc)

  
c..   define value at join from envelope integration
               j = nvmax(3) +1
               scr(kk) =  vnad(j-1,3)-vnab(j-1,3)
c..   define values in envelope
               do k = kk,  kk+nvmax(3)-2
                  if( g(k) .gt. 0.0d0 )then
                     hmp(k) = (p(nc,k)*v(nc,k) + p(nc,k+1)*v(nc,k+1))
     1                    *0.5d0/g(k)
                     nsqr(k) = g(k)*( scr(k) )/hmp(k)
                  else
                     hmp(k) = 0
                     nsqr(k) = 0
                  endif
               enddo
c..   redefine mass variables
               do k = kk, kk+nvmax(3)-2
                  dmi(k) = 0.5d0*( dmh(k+1) + dmh(k) )
               enddo

               write(*,'(8a12)')'k','dP','Gmdm/4pi r4','error'
               write(*,'(i12,1p8e12.3)')kk+1,p(nc,kk+1)-p(nc,kk+2),
     1              grav*xm(kk+1)*dmi(kk+1)/(pi4*r(1,kk+1)**4),
     2              (p(nc,kk+1)-p(nc,kk+2))/
     3              (grav*xm(kk+1)*dmi(kk+1)/(pi4*r(1,kk+1)**4))
     4              -1.0d0

               write(*,'(a6,i6,1p8e12.3)')'kk',kk,p(nc,kk)-p(nc,kk+1),
     1              grav*xm(kk)*dmi(kk)/(pi4*r(1,kk)**4),
     2              (p(nc,kk)-p(nc,kk+1))/
     3              (grav*xm(kk)*dmi(kk)/(pi4*r(1,kk)**4))
     4              -1.0d0

               write(*,'(i12,1p8e12.3)')kk-1,p(nc,kk-1)-p(nc,kk),
     1              grav*xm(kk-1)*dmh(kk-1)/(pi4*r(1,kk-1)**4),
     2              (p(nc,kk-1)-p(nc,kk))/
     3              (grav*xm(kk-1)*dmi(kk-1)/(pi4*r(1,kk-1)**4))
     4              -1.0d0

               write(*,'(8a12)')'rkk','vr(max,3)','err'
               write(*,'(1p8e12.3)')r(1,kk),vr(nvmax(3),3),
     1              (r(1,kk)-vr(nvmax(3),3))/r(1,kk)
               write(*,'(8a12)')'Lkk','L-vl(max,3)','L/Lsol'
               write(*,'(1p8e12.3)')tl(1,kk),tl(1,kk)-vl(nvmax(3),3),
     1              tl(1,kk)/sollum

               kkold(i) = kk
               if( modes .eq. 2 )then
c..   loops runt to kk+1
                  kkeff(i) = kk + nvmax(3) -2
               else
                  kkeff(i) = kk + nvmax(3) -1
               endif

c..   kk now runs over envelope values too
               kk = kkeff(i)


            else

c..   no envelope
c..   save this value for each of i models
               kkeff(i) = kk
               write(*,*)'modes = ',modex,' with ',kk,' zones'

            endif
c..............................................................


c..   x-coordinate
            call getvec(xx,xmin,xmax,cvarx,cxlabel)

c..   y-coordinate
            call getvec(yy,ymin,ymax,cvary,cylabel)

            do k = 1, kk
               ax(i,k) = xx(k)
               ay(i,k) = yy(k)
            enddo

            axmin(i) = xmin
            axmax(i) = xmax
            aymin(i) = ymin
            aymax(i) = ymax

         else
            write(*,*)'genplot2: no file ',cmod(i)
            stop'genplot2: no model file error'
         endif
      enddo

c..   find limits over all models
      do i = 1, imods-1
         xmin = amin1( xmin, axmin(i) )
         xmax = amax1( xmax, axmax(i) )
         ymin = amin1( ymin, aymin(i) )
         ymax = amax1( ymax, aymax(i) )
      enddo
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

      IF ( pgbeg(0,cdevice,1,1) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
      call PGASK (.FALSE.)
c..   roman font=2
      call PGSCF(1)
c..   font scaling (size)
      call PGSCH(1.0)
c..   new page
      call PGPAGE

c..   window
      call PGSWIN(xmin0,xmax0,ymin0,ymax0)
c..   set color
      call PGSCI(1)
c..   set line width (1-201)
      call PGSLW(5)
c..   tickmarks on x(bottom=B,top=C) and y(left)
      call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
c..   Left label
      call PGMTXT('L',2.0,0.5,0.5,cylabel)
c..   bottom label
      call PGMTXT('B',2.0,0.5,0.5,cxlabel)

      do i = 1, imods

         kk = kkeff(i)

         do k = 1, kk
            xx(k) = ax(i,k)
            yy(k) = ay(i,k)
         enddo

c..   set color
         call PGSCI(i+1)
c..   draw curve
         if( linesty .eq. 0 )then
            call PGLINE(kk,xx,yy)
         elseif( linesty .eq. 1 )then
            call PGPT(kk,xx,yy,20)
         elseif( linesty .eq. 2 )then
c..large square dots = 16
            call PGLINE(kk,xx,yy)
            call PGPT(kk,xx,yy,16)
         else
c..small dots = 20
            call PGLINE(kk,xx,yy)
            call PGPT(kk,xx,yy,20)
         endif

c..   indicate join
c..   kkold(i) is the join index (kk)
c..   note shift -1 because xx,yy arrays start at 1 not 2
         call PGSCI(1)
         call PGPT(1,xx(kkold(i)-1),yy(kkold(i)-1),-4)

c..   reset color
         call PGSCI(i+1)
c..   label to accompany model designation
         if( lmods .eq. 0 )then
c..   mass/Msol
            tmass   = xmsol(i)
         elseif( lmods .eq. 1 )then
c..   time before last model (imods) in seconds
            tmass = timm(imods) - timm(i)
         else
c..   time(sec)
            tmass = timm(i)
         endif

         call ftoc(tmass,axistag)
         ctlabel = axistag//' '//cmod(i)
         xlabel = xmin
         ylabel = ymax - float(i-1)*(ymax-ymin)*0.05

c..garbage to foil rampaging segfault
c         write(*,*)xlabel
c         write(*,*)ylabel
c         write(*,*)ctlabel
c         write(*,*)axistag
c         write(*,*)cmod(i)
c         write(*,*)i
ccccccccccccccccccccccccccccccccccc

         call PGTEXT(xlabel,ylabel,ctlabel)

      enddo

      if( device .eq. "/xwin" )then
         pause
      endif

      close(10)

      call PGEND

      stop'successful termination'

 2000 continue
      write(*,*)'genplot: error in  input file: genplot.in'
      write(*,'(4a10)')'labl',labl,'labl1',labl1
      stop'genplot error'

      end


