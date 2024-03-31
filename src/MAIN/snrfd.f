
      subroutine snrfd(v1,d,t, p,pt,vpv, e,et,vev, s,si,k)
      implicit none
c..   last revised 10-29-03

c..   fermi-dirac integrals n = 1/2, 3/2 for nonrelativistic
c..   electron eos

      include 'cconst'

      real*8    v1,d,t,p,pt,vpv,e,et,vev,s,si
      integer*4 k

      integer*4 nerfd,nopt
      parameter(nerfd = 29)
      real*8  bsi(nerfd),bj2(nerfd),bj3(nerfd),bphi(nerfd)
      common/nrfddata/bsi,bj2,bj3,bphi

      real*8 root2,b, b3h, xj2, xj3,phi,b5h,fact,betaf,rhoc

      parameter( root2 = 1.414213562d0, betaf = egrestm*crad2/boltz, 
     1     rhoc = pi4*(egrestm*crad/planck)**3/avagadro*2.0d0 )
c----------------------------------------------------------------
c..   si is ue/kT, s  is S/avagadro*k

      b   = betaf / t
      b3h = b*dsqrt( b )
c..   estimate fermi integral j2 for no positrons
      xj2 = b3h*d / rhoc /root2

      if( xj2 .lt. bj2( 1)        )then
         call loi(xj2,xj3,si,phi)
         nopt = 0
      elseif( xj2 .gt. bj2(nerfd) )then
         call hii(xj2,xj3,si,phi)
         nopt = 2
      else
         call midi(xj2,xj3,si,phi)
         nopt = 1
      endif


      b5h  = b*b3h
      fact = xj3*pcons/b5h *2.0d0*root2

      p    = fact/3.0d0 
      pt   = 0.5d0*p*(5.0d0 - 3.0d0*phi)/t 
      vpv  = - p * phi 
c..   energy per unit volume
      e    = 1.5d0 * p 
      et   = 1.5d0 * pt 
c..   derivative per unit mass
      vev  = 1.5d0 * vpv + e 
c..   entropy per unit Ye
      s    = 5.0d0*xj3/(3.0d0*xj2) - si 
c..   convert to state forms
      vpv  = vpv / v1
      e    = e   * v1
      et   = et  * v1

      s    =  s * v1 * d

      return
      end


      subroutine hii(xj2,xj3,si,phi)
      implicit none
c..   last revised 10-29-03

c..   degeneracy limit
c..   nonrelativistic fermi dirac functions j1/2, j3/2

      include 'cconst'

      real*8 xj2,xj3,si,phi
      real*8 pi2, pio8, u, a, si2, tsi2, bx, cx
      integer*4 n

      parameter( pi2 = pi**2, pio8 =  pi2/8.0d0 )

c..   input: xj2
c..   output: si = (mu -mc**2)/kT, xj3, phi=(3)j2**2/j1*j3
c----------------------------------------------------------------

      u = 0.0d0
      a = 1.5d0 * xj2
      u = a**(2.0d0/3.0d0)
      do n = 1, 5
         u = ( a /(1.0d0 + pio8/(u*u) ) )**(2.0d0/3.0d0)
      enddo

      si   = u
      si2  = si*si
      tsi2 = 3.0d0*si2
c..bx and cx have corrected coeficients
      bx   = 1.0d0 + 0.625d0*pi2/si2
      cx   = 1.0d0 + 0.125d0*pi2/si2
      xj3  = 0.6d0*xj2*si*bx/cx

c..   exact.coef.
      phi  =  ( 1.0d0 - pi2/tsi2 )/0.6d0

      return
      end


      subroutine midi(xj2,xj3,si,phi)
      implicit none

c..intermediate values of fermi dirac functions j1/2, j3/2
c..from table

      real*8 xj2,xj3,si,phi

      real*8 a, b, c, dj2p, dj2m, dj2, del, dj3p, dj3m
      real*8 dasip, dasim

      integer*4 i, j, k, l
      
      integer*4 nerfd
      parameter(nerfd = 29)
      real*8  bsi(nerfd),bj2(nerfd),bj3(nerfd),bphi(nerfd)
      common/nrfddata/bsi,bj2,bj3,bphi

c..input xj2, asi, bj2, bj3, bphi
c..output si, xj3, phi

      j = 1
      do k = 2, 29
        j = k
        if( xj2 .lt. bj2(k) )go to 101
      enddo
  101 continue

c..linear interpol. for (3)j2**2/(j1*j3)
      l   = j - 1
      a   = bphi(j) - bphi(l)
      b   = bj2(j)  - bj2(l)
      c   = xj2     - bj2(l)
      phi = a*c/b   + bphi(l)

      if( j .gt. 28 )j = 28
      l = j - 1
      i = j + 1
c..quadratic interpolation (17.2.81) for j3 and si
      dj2p = bj2(i) - bj2(j)
      dj2m = bj2(j) - bj2(l)
      dj2  = dj2p + dj2m
      del  = xj2  - bj2(j)
c..j3  interpolation
      dj3p  = (bj3(i) - bj3(j))/dj2p
      dj3m  = (bj3(j) - bj3(l))/dj2m
      a     = bj3(j)
      b     = (dj3m*dj2p + dj3p*dj2m)/dj2
      c     = (dj3p - dj3m)/dj2
      xj3   = a + (b + c*del)*del

      dasip = (bsi(i) - bsi(j))/dj2p
      dasim = (bsi(j) - bsi(l))/dj2m
      a     = bsi(j)
      b     = (dasim*dj2p + dasip*dj2m)/dj2
      c     = (dasip - dasim)/dj2
      si    = a + (b + c*del)*del
c...........................................
      return
      end


      subroutine loi(xj2,xj3,si,phi)
      implicit none
c..nondegenerate, nonrelativistic fermi dirac j1/2, j3/2

      real*8 xj2,xj3,si,phi

      real*8 root2, pi, rootpi
      parameter( root2 = 1.414213562d0, pi = 3.14159265359d0 )
      parameter( rootpi = 1.772453851d0 )

      real*8 esi, a, b
c---------------------------------------------------------
c..input: xj2 (J1/2)
c..output: si, xj3 (J3/2), phi (3[J1/2]**2/Jm1/2*J3/2)

      esi = 2.0d0 * xj2/( rootpi - xj2/root2 )
      si  = dlog( esi ) 
      a   = 1.0d0 - esi/root2/4.0d0
      b   = 1.0d0 - esi/root2/2.0d0
      xj3 = 1.5d0 * xj2 * a / b
      phi = 1.003205d0

c..values of phi in bphi(1)

      return
      end









