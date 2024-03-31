      subroutine dthoul(dthc)
c..7-30-06
c..   diffusion using Thoul subroutine: DIFFUSION(M,A,Z,X,CL,AP,AT,AX)
c--------------------------------------------------------------------

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'cconst'
      include 'cburn'
      include 'comdif'

c..   diffusion variables
c..   mdif = number of diffusion groups
      integer*4 mdif
      parameter( mdif = 6 )
c..   coefficients for diffusion in tridiagonal solver
      real*8 tdifa(ndim,kdm),tdifb(ndim,kdm),tdifc(ndim,kdm),
     1     tdifd(ndim,kdm)
      real*8 dfdyp(mdif,kdm),dfdym(mdif,kdm),dxdif(mdif,kdm)
      real*8 vdiff(mdif,kdm)
      integer*4 kcfl,icfl
      real*8 dtcfl,dtcflold

      real*8 difa(mdif),difz(mdif),difx(mdif),difcl(mdif,mdif)
      real*8 difap(mdif),difat(mdif),difax(mdif,mdif)
      real*8 dlnpdr, dlntdr, dlnCdr(mdif)
      real*8 AXsum,funits,thalf(kdm),rhohalf(kdm),xhhalf,xhehalf,xyehalf
      real*8 dxddt(ndim,kdm),fd(ndim,kdm),dxddtmx,deltar(kdm)
      real*8 xo16half,xc12half,xn14half
c..   tridiagonal
      real*8 ttta(kdm),tttb(kdm),tttc(kdm),tttd(kdm),tttu(kdm)
c..   variables to test mixing
      real*8 delxtry

      real*8 sum, sum0, dthc

      integer*4 i,k,j,n, nc,kworst,kxloop

c..   flag for each nucleus giving its diffusion group
      integer*4 jd(ndim),nbugdiff,kdxddt,kdift

c..   intermediate variables for coulomb logarithms:
      REAL*8 ZXA,AC,NI,CZ,XIJ,NE,AO,LAMBDAD,LAMBDA,difC(mdif)

      real*8 hignore

      data nbugdiff/0/,kdift/2/,hignore/1.0d-8/
c--------------------------------------------------------------

c..initialize restrictive dtdiff for diffusion to no restriction at all
c..timstp.f will ignore this unless it is reset to a lower value below
      dtdiff = dthc * 3.0d0

      if( mixmode .eq. 0 .or. mixmode .eq. 2 )then
c..   diffusion is on
c         write(*,*)"diffusion is on"
c..   set up mdif abundance groups for diffusion
         do i = 1, ndim
            jd(i) = 0
         enddo
c..   hydrogen
         jd(nnuc-1) = 1
c..   helium
         jd(nnuc) = 2
c..   carbon
         jd(lc12) = 3
c..   nitrogen
         jd(ln14) = 4
c..   oxygen
         jd(lo16) = 5
c..   electrons = ndif
         jd(nnuc+1) = mdif

         do k = 1, kk
            do j = 1, nnuc
               fd(j,k) = 0.0d0
            enddo
         enddo
         do k = 2, kk

            if( x(nnuc-1,k) .gt. hignore .and. 
     1           x(nnuc-1,k+1) .gt. hignore )then
c..   diffusive flux calculation; boundary centering for fluxes
               if( k .lt. kk )then
                  thalf(k)   = 0.5d0*( t(1,k) + t(1,k+1) )
                  rhohalf(k) = 0.5d0*(1.0d0/v(1,k) + 1.0d0/v(1,k+1) )
                  xhhalf     = 0.5d0*( x(nnuc-1,k) + x(nnuc-1,k+1) )
     1                 /xa(nnuc-1)
                  xhehalf    = 0.5d0*( x(nnuc  ,k) + x(nnuc  ,k+1) )
     1                 /xa(nnuc)
                  xyehalf    = 0.5d0*( x(nnuc+1,k) + x(nnuc+1,k+1) )
c..
                  xo16half   = 0.5d0*( x(lo16,k) + x(lo16,k+1) )
     1                 /xa(lo16)
                  xc12half   = 0.5d0*( x(lc12,k) + x(lc12,k+1) )
     1                 /xa(lc12)
                  xn14half   = 0.5d0*( x(ln14,k) + x(ln14,k+1) )
     1                 /xa(ln14)
               else
c..   values at envelope join
                  thalf(k)   = t(1,k)
                  rhohalf(k) = 1.0d0/v(1,k)
c..   envelope values
                  xhhalf     =  x(nnuc-1,k+1) /xa(nnuc-1)
                  xhehalf    =  x(nnuc  ,k+1) /xa(nnuc)
                  xyehalf    =  x(nnuc+1,k+1) 
c..
                  xo16half   =  x(lo16,k+1)   /xa(lo16)
                  xc12half   =  x(lc12,k+1)   /xa(lc12)
                  xn14half   =  x(ln14,k+1)   /xa(ln14)
               endif

c..   group charges
               difz(1)    =  1.0d0
               difz(2)    =  2.0d0
               difz(3)    =  6.0d0
               difz(4)    =  7.0d0
               difz(5)    =  8.0d0
               difz(mdif) = -1.0d0

c..   atomic weights Thoul style (for ions, not atomic masses)
c..   electrons
               difa(mdif) = egrestm*avagadro
c..   hydrogen group minus electron mass
               difa(1) = xa(nnuc-1) -difz(1)*difa(mdif)
c..   helium group minus electron mass
               difa(2) = xa(nnuc)   -difz(2)*difa(mdif)
c..   carbon group minus electron mass
               difa(3) = xa(lc12)   -difz(3)*difa(mdif)
c..   nitrogen group minus electron mass
               difa(4) = xa(ln14)   -difz(4)*difa(mdif)
c..   oxygen group minus electron mass
               difa(5) = xa(lo16)   -difz(5)*difa(mdif)

c..   mass fractions Thoul style (for ions, not atomic masses)
c..   electrons
               difx(mdif) = egrestm*avagadro*xyehalf
c..   hydrogen group
               difx(1)    = difa(1)*xhhalf
c..   helium group
               difx(2)    = difa(2)*xhehalf
c..   carbon group
               difx(3)    = difa(3)*xc12half
c..   nitrogen group
               difx(4)    = difa(4)*xn14half
c..   oxygen group
               difx(5)    = difa(5)*xo16half

c..   Ye * me/amu
c..   adjust for consistency in Thoul convention
               difx(mdif) = 0.0d0
               do j = 1, mdif-1
                  difx(mdif) = difx(mdif)+ difx(j)/difa(j)*difz(j)
               enddo
               difx(mdif) = difx(mdif)*difa(mdif)

               sum = 0.0d0
               do j = 1, mdif
                  sum = sum + difx(j)
               enddo

c..   uncomment for testing
c      if( nbugdiff .ne. 0 )then
c         write(*,*)'group normalization error ',sum-1.0d0
c         write(*,'(a10,1p8e12.3)')'difx',difx
c      endif
cccccccccccccccccccccccccccccccccccc

c..   make sure groups are renormalized!
               do i = 1,mdif
                  difx(i) = difx(i)/sum
               enddo

c..   uncomment for testing
c      sum = -1.0d0
c      do j = 1, mdif
c         sum = sum + difx(j)
c      enddo
c      if( nbugdiff .ne. 0 )then
c         write(*,*)'renormalization error ',sum
c         write(*,'(a10,1p8e12.3)')'difx',difx
c      endif
cccccccccccccccccccccccccccccccccccccccc

c..   from Thoul's README: compute coulomb logarithms
c..   calculate concentrations from mass fractions:
               ZXA = 0.0d0
               DO I = 1, Mdif-1
                  ZXA = ZXA + difZ(I)*difX(I)/difA(I)
               ENDDO
               DO I = 1, Mdif-1
                  difC(I) = difX(I)/(difA(I)*ZXA)
               ENDDO
               difC(Mdif) = 1.0d0
c..   calculate density of electrons (NE) from mass density (RHO):
               AC = 0.0d0
               DO I = 1, Mdif
                  AC = AC + difA(I)*difC(I)
               ENDDO	
               NE = rhohalf(k)/( 1.6726d-24*AC )
c     calculate interionic distance (AO): 
               NI = 0.0d0
               DO I = 1, Mdif-1
                  NI = NI + difC(I)*NE
               ENDDO
               AO = (0.23873d0/NI)**(1.0d0/3.0d0)	
c     calculate Debye length (LAMBDAD):
               CZ = 0.0d0
               DO I = 1, Mdif
                  CZ = CZ + difC(I)*difZ(I)**2
               ENDDO
               LAMBDAD = 6.9010d0*dSQRT( thalf(k)/(NE*CZ) )
c     calculate LAMBDA to use in Coulomb logarithm:
               LAMBDA = dMAX1(LAMBDAD,AO)
c     calculate Coulomb logarithms:
               DO I = 1, Mdif
                  DO J = 1, Mdif
                     XIJ = 2.3939d3*thalf(k)
     1                    *LAMBDA/ABS(difZ(I)*difZ(J))
                     difCL(I,J) = 0.81245d0*
     1                    dLOG( 1.0d0 + 0.18769d0*XIJ**1.2d0 )
                  ENDDO
               ENDDO
c..   Thoul diffroutine.f
c..   input
c     mdif = number of groups
c     atomic weight     = difa
c     charge            = difz
c     mass fraction     = difx    (sums to 1 over all groups) 
c     coulomb logrithms = difcl
c..   output
c     coefficient vectors ap, at, and coefficient matrix ax
c     (simplifies for mdif=3)

               call  DIFFUSION(mdif,difa,difz,difx,difcl,difap,difat,
     1              difax)

               if( nbugdiff .ne. 0 .and. k .eq. kdift )then
                  write(*,'(a15,i5,3(a12,1pe12.3))')'THOUL: zone ',k,
     1                 'T(K)',t(2,k),'rho(g/cc)',1.0d0/v(2,k),
     2                 'X(he,kk+1)',x(nnuc,k)*xa(nnuc)

                  write(*,'(a5,10a12)')'g','difap','difat',
     1                 'difax(g,1)','difax(g,2)','difax(g,3)',
     1                 'difax(g,4)','difax(g,5)','difax(g,6)'

                  do i = 1, mdif
                     write(*,'(i5,1p10e12.3)')i,difap(i),difat(i),
     1                    (difax(i,j),j=1,mdif)
                  enddo

                  write(*,'(a5,10a12)')'g',
     2                 'difcl(g,1)', 'difcl(g,2)', 'difcl(g,3)',
     3                 'difcl(g,4)', 'difcl(g,5)', 'difcl(g,6)'
                  do i = 1, mdif
                     write(*,'(i5,1p10e12.3)')i,
     2                    (difcl(i,j),j=1,mdif)
                  enddo


c..   check convervation consistency of thoul coefficients
c..   these should sum to zero (to roundoff accuracy)
                  sum = 0
                  do i = 1, mdif
                     sum = sum + difap(i)*difx(i)
                  enddo
                  write(*,'(a35,1pe12.3)')'residual error in AP',sum
                  sum = 0
                  do i = 1, mdif
                     sum = sum + difat(i)*difx(i)
                  enddo
                  write(*,'(a35,1pe12.3)')'residual error in AT',sum
                  
                  sum = 0
                  do i = 1, mdif
                     do j = 1, mdif
                        sum = sum + difax(i,j)*difx(i)
                     enddo
                  enddo
                  write(*,'(a35,1pe12.3)')'residual error in AX',sum
               endif


c..   k.ge.2 .and. k.le.kk
c..   differenced centered on the boundary k
               if( k .ge. 2 .and. k .lt. kk )then
c..   inner boundary k=1 is set to zero 
                  deltar(k) = 0.5d0*(r(1,k+1)-r(1,k-1))
                  dlnpdr = 2.0d0/(p(1,k)+p(1,k+1))*(p(1,k+1)-p(1,k))
     1                 /deltar(k)
                  
                  dlntdr = 2.0d0/(t(1,k)+t(1,k+1))*(t(1,k+1)-t(1,k)) 
     1                 /deltar(k)
c..   hydrogen
c..   note: xa's cancel out here and below
                  dlncdr(1) = 2.0d0/( x(nnuc-1,k+1)/x(nnuc+1,k+1) +
     1                 x(nnuc-1,k)/x(nnuc+1,k) )
     2                 *( x(nnuc-1,k+1)/x(nnuc+1,k+1) -
     3                 x(nnuc-1,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   helium
                  dlncdr(2) = 2.0d0/( x(nnuc,k+1)/x(nnuc+1,k+1) +
     1                 x(nnuc,k)/x(nnuc+1,k) )
     2                 *( x(nnuc,k+1)/x(nnuc+1,k+1) -
     3                 x(nnuc,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   c12
                  dlncdr(3) = 2.0d0/( x(lc12,k+1)/x(nnuc+1,k+1) +
     1                 x(lc12,k)/x(nnuc+1,k) )
     2                 *( x(lc12,k+1)/x(nnuc+1,k+1) -
     3                 x(lc12,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   n14
                  dlncdr(4) = 2.0d0/( x(ln14,k+1)/x(nnuc+1,k+1) +
     1                 x(ln14,k)/x(nnuc+1,k) )
     2                 *( x(ln14,k+1)/x(nnuc+1,k+1) -
     3                 x(ln14,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   o16
                  dlncdr(5) = 2.0d0/( x(lo16,k+1)/x(nnuc+1,k+1) +
     1                 x(lo16,k)/x(nnuc+1,k) )
     2                 *( x(lo16,k+1)/x(nnuc+1,k+1) -
     3                 x(lo16,k)/x(nnuc+1,k) )
     4                 /deltar(k)

               elseif( k .eq. kk )then
c..   difference to avoid "average envelope values"
c      deltar(k) = dmi(k)/( a(k)*rhohalf(k) )
      deltar(k) = r(1,k) - r(1,k-1)
c     dlnpdr = 2.0d0/(p(1,k)+p(1,k-1))*(p(1,k)-p(1,k-1))
c     1              /deltar(k)
c     
c     dlntdr = 2.0d0/(t(1,k)+t(1,k-1))*(t(1,k)-t(1,k-1)) 
c     1              /deltar(k)

c                  deltar(k) = deltar(k-1)

                  dlnpdr = 2.0d0/(p(1,k)+p(1,k-1))*(p(1,k)-p(1,k-1))
     1                 /deltar(k)
                  
                  dlntdr = 2.0d0/(t(1,k)+t(1,k-1))*(t(1,k)-t(1,k-1)) 
     1                 /deltar(k)

c..   hydrogen
c..   note: xa's cancel out here and below
                  dlncdr(1) = 2.0d0/( x(nnuc-1,k+1)/x(nnuc+1,k+1) +
     1                 x(nnuc-1,k)/x(nnuc+1,k) )
     2                 *( x(nnuc-1,k+1)/x(nnuc+1,k+1) -
     3                 x(nnuc-1,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   helium
                  dlncdr(2) = 2.0d0/( x(nnuc,k+1)/x(nnuc+1,k+1) +
     1                 x(nnuc,k)/x(nnuc+1,k) )
     2                 *( x(nnuc,k+1)/x(nnuc+1,k+1) -
     3                 x(nnuc,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   c12
                  dlncdr(3) = 2.0d0/( x(lc12,k+1)/x(nnuc+1,k+1) +
     1                 x(lc12,k)/x(nnuc+1,k) )
     2                 *( x(lc12,k+1)/x(nnuc+1,k+1) -
     3                 x(lc12,k)/x(nnuc+1,k) )
     4                 /deltar(k)

c..   n14
                  dlncdr(4) = 2.0d0/( x(ln14,k+1)/x(nnuc+1,k+1) +
     1                 x(ln14,k)/x(nnuc+1,k) )
     2                 *( x(ln14,k+1)/x(nnuc+1,k+1) -
     3                 x(ln14,k)/x(nnuc+1,k) )
     4                 /deltar(k)
c..   o16
                  dlncdr(5) = 2.0d0/( x(lo16,k+1)/x(nnuc+1,k+1) +
     1                 x(lo16,k)/x(nnuc+1,k) )
     2                 *( x(lo16,k+1)/x(nnuc+1,k+1) -
     3                 x(lo16,k)/x(nnuc+1,k) )
     4                 /deltar(k)

               else
                  write(*,*)'THOUL: diffusion k error:',k,kk
                  stop'thoul diff k error'
               endif

c..   transform gradients to Thoul units for comutation of velocities
               dlnpdr    = dlnpdr    * solrad
               dlntdr    = dlntdr    * solrad
               dlncdr(1) = dlncdr(1) * solrad
               dlncdr(2) = dlncdr(2) * solrad
               dlncdr(3) = dlncdr(3) * solrad
               dlncdr(4) = dlncdr(4) * solrad
               dlncdr(5) = dlncdr(5) * solrad

c..   Thoul's diffusion should insure that this is zero
               dlncdr(mdif) = 0.0d0

c..   change velocities to cgs units
c..   add T**2.5/rho factor
c..   scale radius, time, temperature and density
               funits =  solrad/(secpy*6.0d13) * 
     1              (thalf(k)/1.0d7 )**2.5d0 *( 1.0d2 / rhohalf(k) )

c..   construct diffusion velocity vdiff(i,k) for group i and
c..   zone boundary k, in cgs units
               do i = 1, mdif
c..   sum terms for all composition gradients
                  axsum = 0.0d0
                  do j = 1, mdif
                     axsum = axsum + difax(i,j)*dlncdr(j)
                  enddo

                  vdiff(i,k) = ( difAP(I)*dlnpdr + difAT(I)*dlntdr 
     1                 + AXsum ) * funits

               enddo


               if( nbugdiff .ne. 0 .and. k .eq. kdift )then
                  write(*,'(/a30,8f12.1)')'group charge',
     1                 (difz(i),i=1,mdif)
                  write(*,'(a30,1p8e12.3)')'pressure flux',
     1                 (difap(i)*dlnpdr,i=1,mdif)
                  write(*,'(a30,1p8e12.3)')'temperature flux',
     1                 (difat(i)*dlntdr,i=1,mdif)
                  do j = 1, mdif
                     write(*,'(a25,i5,1p8e12.3)')'composition flux',j,
     1                    (difax(i,j)*dlncdr(j),i=1,mdif)
                  enddo
                  write(*,'(a30,1p8e12.3)')'sums to vdiff-->',
     1                 (vdiff(i,k)/funits,i=1,mdif)

                  write(*,'(/a30,1p8e12.3)')'composition gradients',
     1                 (dlncdr(i),i=1,mdif)
                  write(*,'(a30,1p8e12.3)')'O16 composition coefs',
     1                 (difax(5,i),i=1,mdif)

               endif

               if( nbugdiff .ne. 0 .and. k .eq. kdift )then
                  sum = 0.0d0
                  do j=1,mdif
                     do i = 1, mdif
                        sum = sum + difax(i,j)*dlncdr(j)*difx(i)
                     enddo
                  enddo

                  write(*,'(a35,1p8e12.3)') ' residual error in AX*X',
     1                 sum
               endif

cccccccccccccccccccccccccccccccccccccc
c..   fthoul=0.8 is scaled to effective rate 
c..   from Delahaye and Pinsonneault 2005
ccccccccccccccccccccccccccccccccccccccccccc


c..   calculate diffusion fluxes
               do i = 1, ndim
                  if( jd(i) .ne. 0 )then
                     if( jd(i) .ne. mdif )then
c..   fd is the mass flux for nucleus i across zone boundary k
c..   nucleus i belongs to group jd(i)
c..   difx(i) is the mass fraction of group jd(i) at that boundary
                        fd(i,k) = vdiff(jd(i),k) * rhohalf(k)
     1                       * difx(jd(i))
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   fthoul factor multiplies all fluxes
     1                       *fthoul
                     else
c..   electrons go with nuclei (careful, this is subtle because
c..   of Thoul's choice of atomic mass convention)
                        fd(i,k) = 0.0d0
                     endif

                     if( nbugdiff .ne. 0 .and. k .eq. kdift )then   
                        if( i .eq. lc12 )write(*,'(10x,3a5,a5,9a14)')
     1                       'k','i','jd','cnuc','vdiff','fd','difx',
     2                       'dt*v/dr'

                        write(*,'(a10,3i5,a5,1p9e14.5)')'fluxes',
     1                       k,i,jd(i),cnuc(i),
     1                       vdiff(jd(i),k), fd(i,k),difx(jd(i)),
     5                       vdiff(jd(i),k)*dth(2)/(r(1,k+1)-r(1,k))
                     endif
                  endif
               enddo

c..   derivatives of flux with respect to Y(hydrogen)
c..   convert to cgs
               dfdyp(1,k) = difa(1)*difax(1,1)/deltar(k)
c     1           *( 1.0d0 - difx(1) )/(1.0d0 + difx(1))
     2              *funits*rhohalf(k)*solrad
ccccccccccccccccccccccccccccccccccccccccccccccccccc
     2              *fthoul

               dfdym(1,k) = -dfdyp(1,k)
c..   use diagonal terms for other groups
               do j = 2, mdif-1
                  dfdyp(j,k) = difa(j)*difax(j,j)/deltar(k)
     1                 *funits*rhohalf(k)*solrad
ccccccccccccccccccccccccccccccccccccccccccccccccccc
     2                 *fthoul

                  dfdym(j,k) = -dfdyp(j,k)
               enddo


               dlnpdr = 2.0d0/(p(1,k)+p(1,k-1))*(p(1,k)-p(1,k-1))
     1              /deltar(k)
               
               dlntdr = 2.0d0/(t(1,k)+t(1,k-1))*(t(1,k)-t(1,k-1)) 
     1              /deltar(k)

               if( nbugdiff .ne. 0 .and. k .eq. kdift )then
                  write(*,'(/a25,8a12)')'Test Thoul solution',
     1                 'vdiff sum','vdiff 1st'
                  sum = 0.0d0
                  do j = 1, mdif
                     sum = sum + vdiff(j,k)*difx(j)/difa(j)*difz(j)
                  enddo
                  write(*,'(a25,1p10e12.3)')'charge current',
     1                 sum,(vdiff(j,k)*difx(j)/difa(j)*difz(j), 
     2                 j=1,mdif)

                  sum = 0.0d0
                  do j = 1, mdif
                     sum = sum + vdiff(j,k)*difx(j)
                  enddo
                  write(*,'(a25,1p8e12.3)')'mass current',
     1                 sum,(vdiff(j,k)*difx(j),j=1,mdif) 

                  write(*,'(a25,8a12)')'group variables','error','sum'
                  sum = 0.0d0
                  do j = 1, mdif
                     sum = sum + difx(j)
                  enddo
                  write(*,'(a25,1p8e12.3)')'mass',
     1                 sum -1.0d0,sum
                  sum = 0.0d0
                  do j = 1, mdif
                     sum = sum +  difx(j)/difa(j)*difz(j)
                  enddo
                  write(*,'(a25,1p8e12.3)')'charge',
     1                 sum, -difx(mdif)/difa(mdif)*difz(mdif)

                  write(*,'(a25,a12,10i12)')'check flux sums','sum',
     1                 (j,j=1,mdif)
                  sum = 0.0d0
                  do j = 1,mdif
                     sum = sum + vdiff(j,k)*difx(j)
                  enddo
                  write(*,'(a25,1p8e12.3)')'new fluxes: vdiff',
     1                 sum,(vdiff(j,k)*difx(j),j=1,mdif)

                  sum = -1.0d0 
                  do i = 1, nnuc
                     if( jd(i) .ne. 0 )then
                        sum = sum + difx(jd(i))/difa(jd(i))*xa(i)
                     endif
                  enddo
                  write(*,'(a25,1p8e12.3)')'new group mass error',
     1                 sum

                  write(*,'(a25,a12,a46)')'This column ------->',  
     1                 '---------',
     1                 '  <-- should be near zero for all Thoul tests'
               endif

            else

c..   ignore diffusive settling for vanishing hydrogen
               dxddt(n,k) = 0.0d0
               fd(n,k)    = 0.0d0

            endif

         enddo
c..   boundary conditions on diffusion flux
         do i = 1, nnuc
            fd(i,kk+1) = 0.0d0
            fd(i,1)    = 0.0d0
         enddo

c..   tridiagonal coef.s for thoul diffusion (only first order)
         do j = 1, mdif-1
            dfdyp(j,kk+1) = 0.0d0
            dfdym(j,1)    = 0.0d0
            do k = 2, kk+1
               tdifa(j,k) =  dth(2)*a(k  )/dmh(k)*dfdyp(j,k  )
               tdifc(j,k) = -dth(2)*a(k-1)/dmh(k)*dfdym(j,k-1)
               tdifb(j,k) = 1.0d0 -dth(2)*a(k)/dmh(k)*dfdyp(j,k)
     1              +dth(2)*a(k-1)/dmh(k)*dfdym(j,k-1)
               tdifd(j,k) = 0.0d0
            enddo
         enddo


c..   find cfl condition
         dtcflold = 0.0d0
         dtcfl    = dtcflold
         kcfl     = 0
         icfl     = 0
c..   uses group which is fastest
         do k = 2, kk
            do i = 1,mdif-1
               dtcfl = 2.0d0*abs( vdiff(i,k)*dthc/(r(1,k+1)-r(1,k-1)) )
               if( dtcfl .gt. dtcflold )then
                  dtcflold = dtcfl
                  kcfl = k
                  icfl = i
               endif
            enddo
         enddo



c..   simple diffusion critical coef. is 0.25
c..   0.5 seems to work for the Sun
         if( dtcflold .gt. 0.1d0 )then
c..this restriction will cause timstp.f to estimate a
c..smaller timestep for the next cycle
c..unchanged if dtcfold = dtcfl(max) = 0.5
            dtdiff = dth(2) * 0.1d0/dtcflold
            dthc   = dtdiff
            write(*,'(a30,i5,2(a6,1pe12.3))')
     1           'DTHOUL: Reducing dtc in k = ',kcfl,
     1           'to',dthc,'from',dth(2)
         endif



c..   choose reference nucleus n = nnuc-1
         kxloop = kk

c..   logic assumes n < nnuc+1
         do n = 1, nnuc
            do k = 2, kk+1
c.....................diffusion   uses Thoul H flux only................
c..   explicit Thoul nucleon number time derivative * dth
c..   using mass fluxes forces mass conservation
               if( jd(n) .ge. 1 .and. jd(n) .lt. mdif )then
                  dxddt(n,k) = ( -a(k  )*fd(n,k) + a(k-1)*fd(n,k-1) )
     1                 /( dmh(k) * xa(n) )
               else
                  dxddt(n,k) = 0.0d0
               endif
            enddo
         enddo

c..   solve tridiagonal by Thomas method to get revised values
         do j = 1, mdif-1
            do k = 2, kk+1
               if( j .eq. 1 )then
                  tttd(k-1) = dxddt(nnuc-1,k)*dthc
               elseif( j .eq. 2 )then
                  tttd(k-1) = dxddt(nnuc,k)*dthc
               elseif( j .eq. 3 )then
                  tttd(k-1) = dxddt(lc12,k)*dthc
               elseif( j .eq. 4 )then
                  tttd(k-1) = dxddt(ln14,k)*dthc
               elseif( j .eq. 5 )then
                  tttd(k-1) = dxddt(lo16,k)*dthc
               endif
               ttta(k-1) = tdifa(j,k)
               tttc(k-1) = tdifc(j,k)
               tttb(k-1) = tdifb(j,k)
               tttu(k-1) = 0.0d0
            enddo

            call triad(ttta,tttb,tttc,tttd,tttu,kk)
           
            do k = 2, kk+1
c..   save implicit guesses
               dxdif(j,k) = tttu(k-1)
            enddo

         enddo

         do k = 1, kk
            do n = 1, nnuc
               x(n,k+1) = scrx(n,k+1)
            enddo
c..   explicit overwrites for diffusive changes
            
c..   hydrogen
c..   mass fractions
            x(nnuc-1,k+1) = scrx(nnuc-1,k+1) + dxdif(1,k+1)*xa(nnuc-1)
c..   c12
            x(lc12,k+1)   = scrx(lc12,k+1) + dxdif(3,k+1)*xa(lc12)
c..   n14
            x(ln14,k+1)   = scrx(ln14,k+1) + dxdif(4,k+1)*xa(ln14)
c..   o16
            x(lo16,k+1)   = scrx(lo16,k+1) + dxdif(5,k+1)*xa(lo16)
c..   helium4
            x(nnuc,k+1) = scrx(nnuc,k+1) - (
     1           + dxdif(1,k+1)*xa(nnuc-1)
     2           + dxdif(3,k+1)*xa(lc12)
     3           + dxdif(4,k+1)*xa(ln14)
     4           + dxdif(5,k+1)*xa(lo16)
     5           )

c..   reset Ye and test normalization
            x(nnuc+1,k+1) = 0.0d0
            sum = -1.0d0
            do n = 1, nnuc
               x(nnuc+1,k+1) = x(nnuc+1,k+1) 
     1              + dble( lz(n) )*x(n,k+1)/xa(n)
               sum = sum + x(n,k+1)
            enddo
         enddo

c..   sum to check total conservation of nucleus n, including envelope
c      write(*,*)
c     1        'THOUL',
c     1        ' doing conservation check of each nucleus over all k'
         do n = 1, nnuc
            sum  = 0.0d0
            sum0 = 0.0d0
            do k = 2, kk+1
               sum  = sum  +    x(n,k)*dmh(k)/xm(kk+1)
               sum0 = sum0 + scrx(n,k)*dmh(k)/xm(kk+1)
            enddo
            if( abs(sum-sum0) .gt. 1.0d-8 )then
               write(*,'(a20,a5,3(a10,1pe14.7))')'total check',cnuc(n),
     1              'new sum',sum,'old sum',sum0,
     2              'diff',sum-sum0
            endif
         enddo

c..   find fastest changing zone for each nucleus
         do n = 1, nnuc
            sum = 0.0d0
            kworst = 0
            do k = 2, kk+1
               delxtry = x(n,k) - scrx(n,k)
               if( abs(delxtry) .gt. abs(sum) )then
                  kworst = k
                  sum = delxtry
               endif
            enddo

            if( kworst .gt. 0 .and. abs(sum) .ge. 1.0d-3 )then
               write(*,'(a30,i5,a5,i5,a5,1pe12.3)')'dthoul: k fastest',
     1              kworst,'kk',kk,cnuc(n),sum
            endif
         enddo

c      write(*,*)'THOUL: doing xchecks for kk ',kk
         call xcheck(kk,x)
         call xcheck(kk,scrx)
c      write(*,*)'THOUL: xcheck over 2 to kk done, doing 2 to kk+1'
         call xcheck(kk+1,x)
         call xcheck(kk+1,scrx)

         if(  nbugdiff .ne. 0 )then
            write(*,'(a5,5a12,a45)')'k','err scrx','err x',
     1           'He','H','metals',
     1           'xchecks after thoul, before conv triad'
            do k = kk-2, kk+1
               sum  = -1.0d0
               sum0 = -1.0d0
               do n = 1, nnuc
                  sum0 = sum0 + scrx(n,k)
                  sum  = sum  +    x(n,k)
               enddo
               write(*,'(i5,1p8e12.3)')k,sum0,sum,x(nnuc,k),
     1              x(nnuc-1,k),
     2              1.0d0 -x(nnuc,k) -x(nnuc-1,k)
            enddo
         endif

c..   fluxes now determined, calculate changes in abundance due
c..   to convection and wave-induced mixing
         do k = 2, kk
            do n = 1,nnuc
               if( abs( xa(n)*dxddt(n,k) )*dthc  .gt. 2.0d-2 )then
c..   dangerously large diffusion for this timestep
                  write(*,'(a30,a5,i5)')'dangerous dth for diffusion',
     1                 cnuc(n),k
      dthc = 0.5d0*dth(2)
cccccccccccccccccccccccccccccccccccccccccccccc
                  write(*,'(6(a10,1pe12.3))')'dYddt',dxddt(n,k),
     1                 'dxddt*dt',dxddt(n,k)*dth(2)*xa(n),
     2                 'dth',dth(2),
     3                 'dthc',dthc,
     4                 'xold',xold(n,k),
     5                 'x',x(n,k)
               endif
            enddo
         enddo

         kdxddt = 0
         dxddtmx = 0.0d0
         do k = 2, kk+1
            if( abs( dxddt(nnuc-1,k) ) .gt. dxddtmx )then
               kdxddt = k
               dxddtmx = dxddt(nnuc-1,k)
            endif
         enddo

         if( kdxddt .gt. 0 .and. abs(dxddtmx) .gt. 1.0d-8 )then
            write(*,'(a25,2(a8,i5),5(a14,1pe12.3))')
     1           'THOUL: max. diff. change: ','kdxddt',kdxddt,
     2           'kk',kk,'dxddtmx',dxddtmx,
     3           'dxddtmx*dt',dxddtmx*dthc,
     4           'x(nnuc-1,k)',x(nnuc-1,kdxddt),
     5           'x(nnuc-1,k+1)',x(nnuc-1,kdxddt+1)
            write(*,*)'outer zone + envelope:'
            write(*,'(a5,8a12)')'k','dxddt H','dxddt H dt','vdiff',
     1           'fd k','fd k-1','dmi','xa'
            k = kk
            write(*,'(i5,1p8e12.3)')k,dxddt(nnuc-1,k),
     1           dxddt(nnuc-1,k)*dthc,vdiff(1,k),
     1           fd(nnuc-1,k),fd(nnuc-1,k-1),dmi(k),xa(nnuc-1)
            k = kk+1
            write(*,'(i5,1p8e12.3)')k,dxddt(nnuc-1,k),
     1           dxddt(nnuc-1,k)*dthc,vdiff(1,k),
     1           fd(nnuc-1,k),fd(nnuc-1,k-1),dmh(k),xa(nnuc-1)

         endif

c..   update
         nc   = 2

         call state(2,kk,nc)

         call cinit(2,kk,2)

      else
c..   no diffusion
         do n = 1, nnuc+1
            do k = 1, kk+1
               dxddt(n,k) = 0.0d0
               fd(n,k)    = 0.0d0
            enddo
         enddo
         nc   = 2

         call state(2,kk,nc)

         call cinit(2,kk,2)
      endif

      return
      end
