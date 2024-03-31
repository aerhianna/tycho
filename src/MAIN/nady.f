      subroutine nady(ye,rho,t,p,pt,vpv,e,et,vev,s,si,k)
      implicit none

c..   last modified wda 10-11-01

c..Nadyozhin solution for extremely relativistic e+e- pair gas,
c..both degenerate and nondegenerate
c..use modern constants
      include 'cconst'

      real*8    pairc,tmcc,rhoc,pi2,tiny
      parameter( tiny = 1.0d-10, pi2 = pi*pi )
      parameter( pairc = avagadro * evrestm * ergspermev,
     1     tmcc = pairc/rgas, 
     2     rhoc  = 2.0d0*pi4*(egrestm * crad/planck)**3/avagadro )

      integer*4 k,i
      real*8    ye,rho,t,p,pt,vpv,e,et,vev,s,si
      real*8    d,b,b2,b3,b4,w3,w4,w30,w0,piob2,w,w2,f,df,dw,v
      real*8    dmudkt,tdsdkt,dwdrho,dsdv,perd,vperdv,eerd,veerdv
c----------------------------------------------------------------

c..si is ue/kT
c..s  is S / avagadro*k
      d   = rho * ye
      b   = tmcc / t
      b2  = b*b
      b3  = b*b2
c..estimate fermi integral j2 for no positrons
      w30  = 3.0d0 * d / rhoc
      w0   = w30**(1.0d0/3.0d0)
c..iterate for chemical potential/mc**2
      piob2 = pi2/b2
      w     = w0
      do i = 1, 30
        w2    = w * w
        f     = w30 - w*( w2  + piob2 )
        df    =    -3.0d0*w2 - piob2
        dw    = -f/df
        if( dmin1( abs( f ), abs( dw) ) .lt. tiny )goto 100
        w = w + dw
      enddo
      write(*,*)'nonconvergence in nady.f, k',k
      stop'nady'
 100  continue

c..1.803**2/4 = 0.812702; 1.803 is J2(0)
      b4 = b*b3
      w2 = w*w
      w3 = w*w2
      w4 = w2*w2
      p  = ( 0.25d0*w4
     1   + (0.5d0*pi2/b2 )*w2
     2   + 11.3644d0/b4         )*pcons/3.0d0

      v  = 1.0d0/rho
c..7*pi**4/15 = 45.4576d0
      s     = rhoc /(3.0d0 * rho )*(pi2*w2/b + 4.54576d1/b3)
      e     = rgas*t*s - p/rho + w*pairc*ye
      dmudkt = - 2.0d0*pi2*w/b/( 3.0d0*w2 + pi2/b2 )
      pt    = rho*rgas*( s + ye*dmudkt )
      tdsdkt = rhoc/(3.0d0*rho)*( pi2*w2/b
     1      + 3.0d0*4.54576d1/b3 + dmudkt* 2.0d0*pi2*w/b2 )
      et    =  rgas * tdsdkt
      vpv   = -(rho*ye)**2/rhoc*3.0d0*pairc/(3.0d0*w2 + pi2/b2)
      dwdrho = 3.0d0*ye/rhoc/(3.0d0*w2 + pi2/b2)
      dsdv  = - rhoc/3.0d0*2.0d0*pi2*w/b *rho * dwdrho
     1      + s*rho
      vev = -p*v + rgas*t*v*dsdv

c..extreme relativistic, T=0 values
      perd   =  ( 0.25d0*w0**4 )*pcons/3.0d0
      vperdv = - 4.0d0/3.0d0*perd
      eerd   = - perd/rho + w0*pairc*ye
      veerdv = -perd*v

c..return only thermal part, p(T) - p(T=0) for m=0 gas
c..because chandra.f gives the T=0 part already
      p   = p   - perd
      vpv = vpv - vperdv
      e   = e   - eerd
      vev = vev - veerdv

      return
      end



