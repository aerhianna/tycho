
      subroutine ioniz(tt,d,yht,yat,xz,ynuc,yenuc,yh,ya1,ya2,yei,
     1     eion,yze,aodebye,uion1,uion2,sion,pion,yh2,
     2     weak,strong,zeta,k,n,nn)

c..   8-7-06 wda
      implicit none

      include 'cconst'

c..   ionization of h and he, he+ with coulomb (debye and ion sphere) 
c..   corrections and excluded volume

      real*8 tt,d,yht,yat,xz,yh,ya0,ya1,ya2,yei,eion,yze
      real*8 uion1,uion2,sion,pion,weak,strong
      real*8 ynuc,yenuc
      real*8 t,x1,x2,x3,fact
      real*8 fh,fhe1,fhe2, yz1, ya,t4,yzt0,fz, ye2,ye3
      real*8 fh0,fhe10,fhe20,phih20
      real*8 yzt,yzi,yeiold, top,bottom
      real*8 qh,qa1,qa2
      real*8 fion
      real*8 dtop,dbottom,dye2,dye3,delyei,fun,dfun
      real*8 t6min, t6max
      real*8 aodebye,parth,eparth,parthe,fakion

      real*8 fw,fs,fc,fw2,fs2,pch,sch,ech,pca1,eca1,sca1,pca2,eca2,sca2
      real*8 pw,ps,pc,sw,ss,sc,ew,es,ec,uw,us,uc,dd
      real*8 expon

      real*8 ccev,dlog2,trsat,tvsat,ktev,theta,mamu,mamu3h
      real*8 trot,tvib,zrot,urot,zvib,uvib
      real*8 phih2,yhn,yh2
      real*8 alpha,beta,yhg1,yhg2, qh2, etahi,etaha,etahm
c,etahe
      real*8 etahe0,etahe1,etahe2
      real*8 sh,she,shi,sha,shm,srot,svib
      real*8 ph,phe,ehi,eha,ehm,ehe
      real*8 a2ob,dadye,dbdye
      real*8 zmel,dzmeldt,smel,umel

      real*8 zeta,zetar
      real*8 wfact,xfact,sfact

      integer k,n, nn
c      integer*4 id

      real*8 yech,xch,efermi,wweak,aaodebye
c..ionization potentials (eV) for HI, HeI, HeII
c.. dissociation potential for H2 molecule
      real*8 iph1,iphe1,iphe2,dph2
      data iph1/13.59844d0/,iphe1/24.58741d0/,iphe2/54.41778d0/
      data dph2/4.4781d0/

c.. excludeded volume per mole
      real*8 volh0,volhm,volhe1,volsum

      data volh0/3.92d0/,volhm/7.85d0/,volhe1/4.97d0/

      real*8 uxh0,sxh0,pxh0
      real*8 uxhe0,sxhe0,pxhe0
      real*8 uxhe1,sxhe1,pxhe1
      real*8 uxhm,sxhm,pxhm

      real*8 mfact, lamb

      real*8 ratio,ratio2

c..   flag to turn on/off coulomb corrections in the ionization iteration
c..   ifcoul=0 is off
      integer*4 ifcoul
      data ifcoul /1/, wfact/1.0d0/,sfact/1.0d0/
      data xfact/8.0d0/,mfact/1.0d0/

c..temperature range for explicit ionization calculation
c..avoid helium burning by using 4.0d8 K in t6max (wda 9/20/04)
      data t6min/1.0d-3/, t6max/4.0d2/

c..   fakion=( 2.0d0 *pi *boltz*egrestm/ planck**2 )**1.5d0/avagadro
c..   * 1d9 for t6 variable
c..   no electron spin factor of 2 is folded in
c..   using cconst values
      parameter( fakion = 4.00971158d0 )

      data mamu/0.504d0/
      data trot/170.8d0/, tvib/6.100d3/
      real*8 oldfun,dyeiold

c-------------------------------------------------------------------
c..   input: tt, d, yht, yat, xz, ynuc, yenuc, aodebye, uion1, uion2, k, n
c..   output: yh, ya1, ya2, yei, eion, yze, yh2, eion, pion, sion

c     yht = total mole fraction of h
c     yat =                        he
c     ynuc= total mole number of nuclei
c     yenuc= Ye for electrons associated with ynuc
c     yh  = fraction of h ionized
c     yhn = fraction of h in neutral atoms
c     yh2 = fraction of h in molecules/2
c     ya0 = fraction of he unionized
c     ya1 = fraction of he at least singly ionized
c     ya2 = fraction of he at least doubly ionized
c     y   = fraction of electrons free
c     uion1 = coulomb chem pot (eV) for singly ionized
c     uion2 = coulomb chem pot (eV) for doublely ionized
c     pion = pressure from H, He, with coulomb effects
c     sion = entropy
c     eion = energy/unit mass
c
c     yei = mole fraction of free electrons
c     tt in units deg K
c     t  in units 1.0d6 deg K
c     d  in cgs units
c-------------------------------------------------------------------

c..   ccev converts eV to degrees Kelvin, consistent (1e-16) with eos
      ccev   = cergs/rgas*1.0d-6
c..   fully ionized guess for weak = I0 * aodebye; do not revise as
c..   it is needed in screening.f
c..   approximate ye if fully ionized 
c..   (accurate to mass deficit for finite metal abundance)
      yech = 0.5d0*( 1.0d0 + yht )
      xch  = ( d * yech / 0.9735d6)**(1.0d0/3.0d0)
c..eV energy units
      efermi = 0.511d6*( sqrt(1.0d0 + xch**2) -1.0d0 )
      ktev  = tt/ccev
      aaodebye = sqrt( 30.507d0 * d /ktev * zeta )
      wweak = iph1 * aaodebye * wfact 
c..   hydrogen molecule
      dlog2  = dlog(2.0d0)
      qh2    = dph2*cergs*1.0d-6
c..this is the factor in the chemical potential
c..u/kT = ln( Y * theta / (g * A**1.5 ) )
c..where the mass is in amu (and is used for nuclei and atomic chemical
c..potentials
      theta  = 3.20338727d3 * d /tt**1.5d0
c..   it is the inverse of fakion it tt-->T6

c..   saturation temperatures for H2 molecule partition functions
c..   choose saturation at T(fermi) for hydrogen (13.6ev/k=1.577e5 K)
      trsat  = 1.0d6
      tvsat  = trsat
      mamu3h = mamu**1.5d0
c..   avagadros number * ionization potential in erg/mole (previous values)
      qh  = iph1 *cergs*1.0d-6
      qa1 = iphe1*cergs*1.0d-6
      qa2 = iphe2*cergs*1.0d-6
      t      = tt * 1.0d-6

      if( t .gt. t6max  )then

c-----fully ionized--------------------------------------------
         yhn = 0.0d0
         yh  = 1.0d0
         ya0 = 0.0d0
         ya1 = 0.0d0
         ya2 = 1.0d0
         ye3 = yht
         ye2 = yat * 2.0d0
         yze = xz*yenuc
         yei = ye3 + ye2 + yze
         xch    = ( d * yei / 0.9735d6 )**(1.0d0/3.0d0)
c..eV energy units
         efermi = 0.511d6*( sqrt(1.0d0 + xch**2) -1.0d0 )
         ratio  = dexp( -(iph1 -efermi)/ktev )   
         ratio2 = dexp( -(iphe2-efermi)/ktev )
         lamb = 1.5d0 * ktev / iph1 * d * volh0 * mfact
         yhn  = (iph1 /efermi)**1.5d0 *(ratio - 1.0d0)/(ratio + 1.0d0) 
         yh2  = 0.0d0
         yh   = 1.0d0 - yhn - yh2*2.0d0

         ya1  = (iphe2/efermi)**1.5d0 *(ratio2 - 1.0d0)/(ratio2 + 1.0d0) 
         ya0  = 0.0d0
         ya2  = 1.0d0 -ya1 -ya0 
c..first order correction
         yei  = yht*yh + yat*(ya2*2.0d0 + ya1)
         xch    = ( d * yei / 0.9735d6 )**(1.0d0/3.0d0)
         strong = iph1 * 102.01d0 * xch * sfact
         eparth = 0.0d0
         nn   = 1

c         theta = 3.20338727d3 * d /tt**1.5d0
c..   free energies for weak and strong screening
         if( yht .gt. 1.0d-15 )then
c..   ionized hydrogen
            ye3 = yht*yh
            fs  = -strong * cergs *1.0d-6 *yht * yh
            fw  = - wweak * cergs *1.0d-6 *yht * yh * 2.0d0
            fw2 = fw**2
            fs2 = fs**2
            dd  = sqrt( fw2 + fs2 )
            fc  = - fw*fs/dd
            pw  = 0.5d0 * d *fw
            ps  = d * fs / 3.0d0
            sw  = 0.5d0 * fw /( rgas * tt )
            ss  = 0.0d0
            ew  = fw + sw * rgas * tt
            es  = fs
            pc  = (-fs/dd)**3 * pw + (-fw/dd)**3 * ps
            sc  = (-fs/dd)**3 * sw + (-fw/dd)**3 * ss
            ec  = (-fs/dd)**3 * ew + (-fw/dd)**3 * es
            pch    = pc
            sch    = sc
            ech    = ec

c..   pressure of H
            ph    = rgas*tt*d*yht*yh + pch
c..   hydrogen u/kT
            etahi = dlog( yht*yh *theta )
c..   extra entropy NOT in chemical potential
            sh   =  yht*yh *( 2.5d0 - etahi ) + sch
c..   ionization energy of metals ignored in eos
c..   eion goes to 0.0 when plasma is neutral atoms
            ehi  =  yht*yh * ( qh + 1.5d0*rgas*tt ) + ech
         else
            ph  = 0.0d0
            sh  = 0.0d0
            ehi = 0.0d0
            eha = 0.0d0
            ehm = 0.0d0
            ye3 = 0.0d0
         endif
c..   free energies for weak and strong screening
         if( yat .gt. 1.0d-15 )then
            pca1 = 0
            sca1 = 0
            eca1 = 0
c..   doublely ionized helium
            fw  = - wweak  * cergs *1.0d-6 *yat * ya2 * 6.0d0
            fs  = - strong * cergs *1.0d-6 *yat * ya2
     1           * 2.0d0**(5.0d0/3.0d0)
            fw2 = fw**2
            fs2 = fs**2
            dd  = sqrt( fw2 + fs2 )
            fc  = - fw*fs/dd
            pw  = 0.5d0 * fw * d
            ps  = d * fs /3.0d0
            sw  = 0.5d0 * fw /( rgas * tt ) 
            ss  = 0.0d0
            ew  = fw + sw * rgas * tt
            es  = fs
            pc  = (-fs/dd)**3 * pw + (-fw/dd)**3 * ps
            sc  = (-fs/dd)**3 * sw + (-fw/dd)**3 * ss
            ec  = (-fs/dd)**3 * ew + (-fw/dd)**3 * es

            pca2 = pc
            sca2 = sc
            eca2 = ec
c..   pressure of He 
            phe   = rgas*tt*d*yat + pca1 + pca2
c..   helium u/kT
            if( ya0 .gt. 0.0d0 )then
               etahe0 = dlog( yat*ya0    *theta )
               etahe2 = dlog( yat*ya2    *theta )
               etahe1 = dlog( yat*ya1/2.0d0   *theta )
c..   helium entropy for all ionization stages
               she   = yat*( ya1*( 2.5d0 - etahe1)
     1              +        ya2*( 2.5d0 - etahe2) 
     2              +        ya0*( 2.5d0 - etahe0 ) ) + sca1 + sca2
            else
c..protect against 0*ln(0) when neutral He --> 0 at high T
               etahe0 = 0.0d0
               etahe2 = dlog( yat*ya2    *theta )
               etahe1 = dlog( yat*ya1/2.0d0   *theta )
c..   helium entropy for all ionization stages
               she   = yat*( ya1*( 2.5d0 - etahe1)
     1              +        ya2*( 2.5d0 - etahe2) 
     2              +        ya0*( 2.5d0 - etahe0 ) ) + sca1 + sca2
            endif
c..   heII partition function constant
            ehe  = yat*( ya1*qa1 + ya2*(qa1+qa2)) + yat*1.5d0*rgas*tt
     1           + eca1 + eca2
         else
c..   no He
            phe = 0.0d0
            she = 0.0d0
            ehe = 0.0d0
            ye2 = 0.0d0
         endif
c..   free electrons
         yei = ye3 + ye2 + yze
c..   pressure of H and He
         pion = ph + phe 
c..   entropy of H and He
         sion = sh + she 
c..   energy/unit mass of H and He
         eion = ehi + ehe
         return
         
      elseif( t .lt. t6min )then

c-----un-ionized except for some metals, H is molecular-----------------
         yh  = 0.0d0
         ya1 = 0.0d0
         ya2 = 0.0d0
         yze = 0.0d0
         ye3 = 0.0d0
         ye2 = 0.0d0
         yze = 1.0d-6*yenuc
         yei = ye3 + ye2 + yze
         eparth = 0.0d0
         nn  = 1
         yhn = 0.0d0
         yh2 = 0.5d0

c..   free energies for weak and strong screening
         if( yht .gt. 1.0d-15 )then
c..no screening, neutral gas
            pch = 0.0d0
            sch = 0.0d0
            ech = 0.0d0
c..   pressure of H
            ph    = rgas*tt*d* yht*yh2

c..   hydrogen u/kT
            etahm = dlog( yht*yh2/(zrot*zvib*zmel)*mamu3h*theta )
c..   extra entropy NOT in chemical potential
            srot  = yht*yh2* trsat/(trsat + tt)
            svib  = yht*yh2* tvib/tt/
     1           (exp(tvib/tt/tvsat*(tvsat+tt))-1.0d0)
c..   rough molecular electron partition function
            zmel   = 1.0d0 + 0.2d0*exp(-3.0d3/tt) + 0.4d0*exp(-1.0d4/tt)
     1           + 0.1d0*exp(-1.3d4/tt)
            dzmeldt = ( 0.2d0*3.0d3*exp(-3.0d3/tt) 
     1           + 0.4d0*1.0d4*exp(-1.0d4/tt)
     2           + 0.1d0*1.3d4*exp(-1.3d4/tt)  )/( tt*zmel )
            smel  = yht*yh2* dzmeldt 
            shm   =  yht*yh2*( 2.5d0 - etahm ) + srot + svib + smel
c     1           -yht*yh2*dlnchi
            sh    =  shm + srot + svib + smel
c..   ionization energy of metals ignored in eos
c..   eion goes to 0.0 when plasma is neutral atoms
            ehm  =  yht*( yh2*(-qh2                 + 1.5d0*rgas*tt  
     1           + rgas*tt*( srot + svib + smel ) ) )
         else
            ph = 0.0d0
            sh = 0.0d0
            ehi = 0.0d0
            eha = 0.0d0
            ehm = 0.0d0
         endif

         return

      else

c-----general partially ionized, partially dissociated case----------

c..first guess at composition for excluded volume (all H0 and He0)
         volsum = ( yht*volh0 + yat*volhe1 )*xfact
         uxh0  = d * volsum * xfact 
         if( uxh0 .gt. 100.0d0 )then
            uxh0 = 100.0d0
         endif
         uxhe0 = uxh0
         uxhe1 = uxh0
c..   approximate ye as if fully ionized 
c..   (accurate to mass deficit for finite metal abundance)
         yech   = 0.5d0*( 1.0d0 + yht )
         xch    = ( d * yech / 0.9735d6 )**(1.0d0/3.0d0)
c..   strong shielding is slightly enhanced to make the weak-strong
c..   transition more accurate
         strong = iph1 * 102.01d0 * xch * sfact
ccccccccccccccccccccccccccccccccccccccccccccccc

c..   eV energy units
         efermi = 0.511d6*( sqrt(1.0d0 + xch**2) -1.0d0 )
         ktev   = tt/ccev
         aaodebye = sqrt( 30.507d0 * d /ktev * zeta )
         wweak  = iph1 * aaodebye * wfact
c..   partition functions saturate at temperatures trsat and tvsat
         zrot =  tt/trot *trsat/(trsat+tt) 
         urot = -ktev*dlog(zrot)
         zvib = 1.0d0/( 1.0d0
     1        -dexp( -tvib/tt*(tvsat+tt)/tvsat ) )
         uvib = -ktev*dlog(zvib)
c..   rough molecular electron partition function
         zmel    = 1.0d0 + 0.2d0*exp(-3.0d3/tt) + 0.4d0*exp(-1.0d4/tt)
     1        + 0.1d0*exp(-1.3d4/tt)
         dzmeldt = ( 0.2d0*3.0d3*exp(-3.0d3/tt) 
     1        + 0.4d0*1.0d4*exp(-1.0d4/tt)
     2        + 0.1d0*1.3d4*exp(-1.3d4/tt)  )/( tt*zmel )
         umel    = -ktev*dlog(zmel)
c..   net effect of H molecule and atom excluded volume not 0 
         uxhm = d*(volhm-2.0d0*volh0)*xfact
         phih20   = 4.0d0*mamu3h
     1        /theta*exp(-(dph2 -urot -uvib -umel )/ktev 
     2        + uxhm )

c..avoid overflow for large phi20 (in a2ob variable in iteration below)
         phih20= dmin1( phih20,1.0d50 )

         fact  = fakion * t**1.5d0/d 
c..ionization potentials for H,He+,He from Cox/Allen, p. 36
         x1 = iph1 /ktev
         x2 = iphe1/ktev
         x3 = iphe2/ktev
c..   deal with over and underflow
         x1 = dmin1( x1, 500.0d0 )
         x1 = dmax1( x1,-500.0d0 )
         x2 = dmin1( x2, 500.0d0 )
         x2 = dmax1( x2,-500.0d0 )
         x3 = dmin1( x3, 500.0d0 )
         x3 = dmax1( x3,-500.0d0 )
c..   nuclear factors cancel so only electrons must be considered
c..   HI ground state is 2
         parth = 2.0d0
c..   eparth = d ln parth / d ln T
         eparth = 0.0d0
c..   heII partition function 
         parthe = 1.0d0
c..   classical factor (no composition dependence)
c..   g*g/g factors = 1, 4, 1
         fh0    = fact*dexp(-x1) 
         fhe10  = fact*dexp(-x2) * 4.0d0
         fhe20  = fact*dexp(-x3) 
c..   set composition dependent part (weak) using full ionization approximation

         if( ifcoul .eq. 1 )then
c..   first guess to define chemical potentials
            yh  = 1.0d0
            ya1 = 0.5d0
            ya2 = 0.5d0
c..   for H
            fw  = - wweak  * cergs *1.0d-6 *yht * yh * 2.0d0
            fs  = - strong * cergs *1.0d-6 *yht * yh
            fw2 = fw**2
            fs2 = fs**2
            dd  = sqrt( fw2 + fs2 )
            fc  = - fw*fs/dd
            uw  = - 2.0d0 * wweak * (1.0d0 + yht*yh*2.0d0/(2.0d0*zeta))
            us  = - strong 
            uc  = (-fs/dd)**3 * uw + (-fw/dd)**3 * us
            expon = -uc/ktev + uxh0
            expon = dmin1(300.0d0,expon)
            fh   = fh0 * exp( expon )
c..   He+
            fw  = - wweak  * cergs *1.0d-6 *yat * ya1 * 2.0d0
            fs  = - strong * cergs *1.0d-6 *yat * ya1
            fw2 = fw**2
            fs2 = fs**2
            dd  = sqrt( fw2 + fs2 )
            fc  = - fw*fs/dd
            uw  = - 2.0d0 * wweak * (1.0d0 + yat*ya1*2.0d0/(2.0d0*zeta))
            us  = - strong 
            uc  = (-fs/dd)**3 * uw + (-fw/dd)**3 * us
            expon = -uc/ktev + uxhe0 - uxhe1
            expon = dmin1(300.0d0,expon)
            fhe1 = fhe10 * exp( expon )
c..   He++
            fw  = - wweak  * cergs *1.0d-6 *yat * ya2 * 6.0d0
            fs  = - strong * cergs *1.0d-6 *yat * ya2 
     1           * 2.0d0**(5.0d0/3.0d0)
            fw2 = fw**2
            fs2 = fs**2
            dd  = sqrt( fw2 + fs2 )
            fc  = - fw*fs/dd
            uw  = - 6.0d0 * wweak * (1.0d0 + yat*ya2*6.0d0/(2.0d0*zeta))
            us  = - strong * 2.0d0**(5.0d0/3.0d0)
            uc  = (-fs/dd)**3 * uw + (-fw/dd)**3 * us
            expon = -uc/ktev + uxhe1
            expon = dmin1(300.0d0,expon)
            fhe2 = fhe20* exp( expon )
         elseif( ifcoul .eq. 2 )then
c..   no coulomb but excluded volume is added
            fh     = fh0  *dexp( uxh0 )
            fhe1   = fhe10*dexp( uxhe0 - uxhe1)
            fhe2   = fhe20*dexp( uxhe1 )
         else
            fh     = fh0
            fhe1   = fhe10
            fhe2   = fhe20
         endif

         lamb = 1.5d0 * ktev / iph1 * d * volh0 * mfact
c..   initialize for iteration
         fh     = fh0  
         fhe1   = fhe10
         fhe2   = fhe20
         phih2  = phih20
c..   estimate free electrons prior to iteration
c..   pick dominant term for first guess for hydrogen
c..   suppose molecular hydrogen is negligible
         if( fh .lt. 1.0d12 )then
c..   y(H+) is Y(H+)/YH
            if( yht .gt. 1.0d-15 )then
c..   finite amount of hydrogen
c..   y(H+) is Y(H+)/YH
               yhg1 = ( sqrt( 4.0d0*fh*yht + fh**2 ) - fh )*0.5d0/yht
c..   suppose ionized hydrogen negligible
c..   y(H0) is Y(H0)/YH
               if( phih2 .lt. 1.0d3 )then
                  yhg2 = phih2*0.25d0*( sqrt( 1.0d0 + 8.0d0*yht/phih2 )
     1                 -1.0d0 )/yht
               else
                  yhg2 =  1.0d0 - 2.0d0*yht/phih2 
               endif
            else
               yhg1 = 1.0d0
               yhg2 = 1.0d0
            endif
         else
c..   guess it is ionized
            yhg1 = 1.0d0
            yhg2 = 0.0d0
            fh   = 1.0d12
         endif
         if( 1.0d0-yhg1 .gt. yhg2 )then
c..   use 2 , H+ is small
            yhn = yhg2
            yh  = yhg1
            yh2 = 0.5d0*(1.0d0 - yhg2/(1.0d0+yhg1) )
         else
c..   use 1, H+ is large
            yh  = yhg1
            yh2 = 0.0d0
            yhn = 1.0d0 - yhg1
         endif
c..   helium contribution
c..   try hydrogen estimate plus full ionization of He as first guess
         yei    = yht*yh + yat*2.0d0
         bottom = yei**2 + yei*fhe1 + fhe1*fhe2
c..yhe+
         ya1 = yei*fhe1/bottom
c..yhe++
         ya2 = fhe1*fhe2/bottom
c..he charge Ye(He)
         ya  = ya1 + 2.0d0*ya2
c..neutral helium
         ya0 = yei**2/bottom
         yei = yht*yh + yat*( ya1 + 2.0d0*ya2 )
c..   metal electrons
         if( xz .gt. 0.0d0 )then
            t4 = 100.0d0*t
c--------------------------------------------------------------------
            yzt0    =     (  8.3d-5/(1.0d0 + (0.4d0/t4)**15 )
     1           + 2.3d-6/(1.0d0 + (0.25d0/t4)**15 )
     2           + 1.99d-2/(2.0d0 + (3.0d0/t4)**6 )  )/0.02d0
c--------------------------------------------------------------------
c..   hydrogenic approximation
            fz    = fakion*t**1.5d0/d*dexp(-x1)
            yz1   = fion(fz,ynuc)
c..   normalize for correct high T, rho value
            yzt0  = dmin1(yzt0,yenuc)
c..   mean A is 1/ynuc
            yzt   = (yz1*yenuc + yzt0 )
            yzi   = dmin1(yzt,yenuc)
         else
            yzi = 0.0d0
         endif

c..   includes metal electrons
         yze    = xz * yzi
         yei    = yh*yht + yat*ya1 + 2.0d0*yat*ya2  + yze
         if( yht + yat .gt. 0.0d0 )then
            yeiold = yei
            fun    = 0.0d0
            delyei = 0.1d0
            yeiold = yei

c..   loop for iterative refinement
            do nn = 1, 200
c..   excluded volume for atoms and molecules; take constant 
c..   for better convergence
c               volsum = volh0 
               volsum = volh0
c..   chemical potential corrections for excluded volume / kT
               uxh0  = d * volsum * xfact 
               if( uxh0 .gt. 100.0d0 )then
                  uxh0 = 100.0d0
               endif
               uxhe0 = uxh0
               uxhe1 = uxh0
c..   add coulomb correction to chemical potential 
c..   (modifies effective masses)
               zetar = 2.0d0*yht*yh +2.0d0*yat*ya1
     1              + 6.0d0*yat*ya2
               aaodebye = sqrt( 30.507d0 * d /ktev * zetar )
               wweak    = iph1 * aaodebye *wfact
c..   delta mu/kT due to excluded volume

c..set composition dependent part (weak) 
               if( ifcoul .eq. 1 )then
c..   h+
                  uw  =  -2.0d0 * wweak 
     1                 *(1.0d0 + yht*yh*2.0d0/(2.0d0*zetar))
                  us  =  -strong * 4.0d0/3.0d0

c.. yh dependence cancels in fs/dd, fw/dd ratios, but not Zi factors
                  fw  = - wweak  * cergs *1.0d-6 *yht * yh * 2.0d0
                  fs  = - strong * cergs *1.0d-6 *yht * yh
                  fw2 = fw**2
                  fs2 = fs**2
                  dd  = sqrt( fw2 + fs2 )
                  if( dd .gt. 0.0d0 )then
                     uc  = (-fs/dd)**3 * uw + (-fw/dd)**3 * us
                  else
                     uc = uw
                  endif
                  expon =  - uc/ktev + uxh0 
                  expon = dmin1(300.0d0,expon)
                  fh   = fh0 * exp( expon )
c..   he+
                  uw  =  -2.0d0 * wweak 
     1                 * (1.0d0 + yat*ya1*2.0d0/(2.0d0*zetar))
                  us  =  -strong  * 4.0d0/3.0d0
c.. ya1 dependence cancels in fs/dd, fw/dd ratios
                  fw  = - wweak  * cergs *1.0d-6 *yat * ya1 * 2.0d0
                  fs  = - strong * cergs *1.0d-6 *yat * ya1
                  fw2 = fw**2
                  fs2 = fs**2
                  dd  = sqrt( fw2 + fs2 )
                  uc  = (-fs/dd)**3 * uw + (-fw/dd)**3 * us
                  expon = - uc/ktev + uxhe0 - uxhe1 
                  expon = dmin1(300.0d0,expon)
                  fhe1 = fhe10 * exp( expon )
c..   he++
                  uw  =  -6.0d0 * wweak
     1                 *(1.0d0 + yat*ya2*6.0d0/(2.0d0*zetar))
                  us  =  -strong * 2.0d0**(5.0d0/3.0d0) * 4.0d0/3.0d0
c.. ya2 dependence cancels in fs/dd, fw/dd ratios
                  fw  = - wweak  * cergs *1.0d-6 *yat * ya2 * 6.0d0
                  fs  = - strong * cergs *1.0d-6 *yat * ya2 
     1                 * 2.0d0**(5.0d0/3.0d0 )
                  fw2 = fw**2
                  fs2 = fs**2
                  dd  = sqrt( fw2 + fs2 )
                  uc  = (-fs/dd)**3 * uw + (-fw/dd)**3 * us
                  expon =  - uc/ktev + uxhe1 
                  expon = dmin1(300.0d0,expon)
                  fhe2 = fhe20 * exp( expon )
               elseif( ifcoul .eq. 2 )then
c..no coulomb but excluded volume is included
                  fh     = fh0  *dexp( uxh0 )
                  fhe1   = fhe10*dexp( uxhe0 - uxhe1)
                  fhe2   = fhe20*dexp( uxhe1 )
               else
c..no coulomb or excluded volume
                  fh     = fh0
                  fhe1   = fhe10
                  fhe2   = fhe20
               endif

c..   helium contribution
               top    = fhe1*(yei + 2.0d0*fhe2)
               bottom = yei**2 + yei*fhe1 + fhe1*fhe2

               if( ifcoul .eq. 1 )then
                  dtop    = fhe1
                  dbottom = 2.0d0*yei + fhe1
               elseif( ifcoul .eq. 2 )then
                  dtop    = fhe1
                  dbottom = 2.0d0*yei + fhe1

               else
                  dtop    = fhe1
                  dbottom = 2.0d0*yei + fhe1
               endif

               ye2    = yat * top/bottom
               dye2   = yat*(dtop/top -dbottom/bottom)*top/bottom

c..   hydrogen contribution
               if( yht .gt. 1.0d-15 )then
                  if( fh .lt. 1.0d90 )then
c..........................   1.0d30 gives a glitch (wda 9-27-05)
c..   normal hydrogen partitial ionization
                     alpha = (1.0d0 + fh/yei)*phih2*(fh/yei)*0.25d0
                     beta  = yht             *phih2*(fh/yei)**2*0.5d0
                     a2ob  = alpha**2/beta

                     if( ifcoul .eq. 0 )then
c..   no coulomb or excluded volume
                        dbdye = -2.0d0*beta/yei
                        dadye = -0.25d0*phih2*(fh/yei)**2
     1                       *(2.0d0/yei + 1.0d0/fh)
                     elseif( ifcoul .eq. 1 )then
c..   coulomb and excluded volume
                        dbdye = -2.0d0*beta/yei
                        dadye = -0.25d0*phih2*(fh/yei)**2
     1                       *(2.0d0/yei + 1.0d0/fh
     2                       )
                     else
c..no coulomb but excluded volume ifcoul=2
                        dbdye = -2.0d0*beta/yei
                        dadye = -0.25d0*phih2*(fh/yei)**2
     1                       *(2.0d0/yei + 1.0d0/fh )
                     endif
c..treat limiting cases sanely
                     if( a2ob .lt. 1.0d-5 )then
c..   second order expansion of sqrt term (low  ionization)
                        ye3  = sqrt(beta)*( 1.0d0-sqrt(a2ob)+0.5d0*a2ob)
                        dye3 = -sqrt(beta)/yei
                     elseif( a2ob .gt. 1.0d4 )then
c..   third order expansion of sqrt term (high ionization)
c..   trouble with 1.0d5 (slow convergence occasionally)
                        ye3  = 0.5d0*beta/alpha - 0.125d0*alpha/a2ob**2
     1                       + 0.125d0*beta/alpha/a2ob**3
c..   second order derivative in ye3
                        dye3 = 0.5d0*(dbdye - beta/alpha*dadye)/alpha
     1                       -0.125d0*dadye/a2ob**2
                     else
                        ye3  = sqrt( alpha**2 + beta ) - alpha
                        dye3 = (alpha*dadye +0.5d0*dbdye)
     1                       /sqrt(alpha**2+beta) - dadye
                     endif
c..   derivative assumed insensitive
                     yh = ye3/yht
                  else
c..   extreme ionization, so alpha and beta formulation becomes singular
c..   and we use this limiting case instead
                     ye3  = yht*( 1.0d0 - yei/fh )
c..   first order derivative in ye3
                     dye3 = -yht/fh
                     a2ob = 1.0d30
                  endif
               else
c..   neutral limit
                  ye3   = 0.0d0
                  yh    = 1.0d0
                  dye3  = 0.0d0
                  alpha = 1.0d0
                  beta  = 0.0d0
                  a2ob  = 1.0d0
                  dadye = 0.0d0
                  dbdye = 0.0d0
               endif
c..   add H, He, and metal ions to charge conservation condition
               oldfun  = fun
               dyeiold = delyei
               fun    = yei   - ( ye3  + ye2 + yze )
               dfun = 1.0d0 - ( dye3 + dye2      ) 
               delyei = - fun / dfun 

               if( abs(delyei) .gt. 5.0d-1 .and. nn .gt. 1 )then
                  write(*,'(a15,i4,1p13e11.3)')'large delyei',
     1                 nn,tt,yei,yht+2.0d0*yat,delyei,fun,dfun,
     2                 yhg1,yhg2,fh,fhe1,fhe2,phih2,yei+delyei
               endif

               if( tt .le. 8.0d0 )then
                  if( nn .le. 1 )then
                     write(*,'(a20,1p12e12.3)')'guesses',yhg1,yhg2,
     1                    phih2,yh2,yhn,yh
                     write(*,'(a5,13a11)')'nn','tt','yei','delyei',
     2                    'fun','dfun','dye3','ye3','a2ob','alpha',
     3                    'beta','d*yhn*vh0'
                  endif 

                  write(*,'(i5,1p14e11.3)')nn,tt,yei,delyei, 
     1                 fun,dfun,dye3,ye3,a2ob,alpha,beta,
     2                 d*volsum

               endif

               if( yei + delyei .gt. yht + yat*2.0d0 + yze +1.0d-15
     1              )then
                  if( yei .gt.  yht + yat*2.0d0 + yze )then
                     delyei = yht + yat*2.0d0 + yze -yei
                     stop'b'
                  else
c..   trial overshoots fully ionized limit; adjust
                     delyei = 0.5d0*( yht + yat*2.0d0 + yze - yei )
                  endif
               endif

c..   yhn and yh2 are fractions of H, as is yh (yh+yhn+yh2=1)
               if( yht .gt. 1.0d-15 )then
c..   this flip of /yht to *yht is correct
c..   ye3,yei are mole fractions, yhn,yh2 are fractions of H nuclei
c..   so there is a factor of yht floating around
                  yhn = yei*ye3/fh /yht
                  yh2 = yhn**2/phih2 *yht
               else
                  yhn = 0.0d0
                  yh2 = 0.0d0
               endif
               if( yat .gt. 1.0d-15 )then
                  ya0 = yei**2/bottom
                  ya1 = yei  * fhe1 / bottom
                  ya2 = fhe2 * fhe1 / bottom
               else
                  ya0 = 0.0d0
                  ya1 = 0.0d0
                  ya2 = 1.0d0
               endif
                  
               if( nn .ge. 190  )then
c..   convergence trouble; write diagnostic data
                  if( nn .eq. 190 )then
                     write(*,'(2a5,13a11)')'k','nn','yei','yht',
     1                    'fun','dfun','dye3','a2ob','ye3','b/2a',
     1                    '1st order','2nd order','delyei','wweak'
                  endif
                  write(*,'(3i5,1p16e11.3)')k,nn,n,yei,yht,fun,dfun,
     1                 dye3,
     1                 a2ob,ye3,beta/alpha*0.5d0,ye3-beta/alpha*0.5d0,
     2                 sqrt(alpha**2+beta)-alpha-beta/alpha*0.5d0
     3                 +0.125d0*alpha/a2ob**2,delyei,wweak,
     4                 alpha,beta
               endif

c..check for exit criteria
               if( abs(fun) .lt. 1.0d-10 .and. nn .ge. 2 )goto 100

c..   update and continue
               yeiold = yei
               yei    = yeiold + delyei

            enddo

c..   error exit
            write(*,'(a30)')"ERROR: nonconvergence in ioniz"
            write(*,'(3a6,8a12)')"k","n","nn","yei","yeiold",
     1           "tt", "d", "yht", "yat", "yz"
            write(*,'(3i6,1p8e12.3)')
     1           k,n,nn,yei,yeiold,tt,d,yht,yat,yze
            write(*,'(8a12)')'ye3','ye2','yemax','fun','dfun',
     1           'dye3','dye2'
            write(*,'(1p8e12.3)')ye3,ye2,yht+2.0d0*yat+0.5d0*xz,
     1           fun,dfun,dye3,dye2

            write(*,'(8a12)')'fhe1','fhe2','top','bottom',
     1           'dtop','dbottom'
            write(*,'(1p8e12.3)')fhe1,fhe2,top,bottom,dtop,dbottom
            write(*,'(8a12)')'dd','uc','uw','us','wweak','iph1',
     1           'aodebye','wfact'
            write(*,'(1p8e12.3)')dd,uc,uw,us,wweak,iph1,aaodebye,wfact

            write(*,'(8a12)')'efermi','ktev','zetar'
            write(*,'(1p8e12.3)')efermi,ktev,zetar

            write(*,'(8a12)')'yh','ya1','ya2','fhe10','fhe20'
            write(*,'(1p8e12.3)')yh,ya1,ya2,fhe10,fhe20

            write(*,'(2(a20,1pe12.3))')'entering T6 ',t,
     1           'yei-yeiold ',yei-yeiold
            write(*,'(a20,i5)')'ifcoul',ifcoul
            stop'ioniz.f aaa'

 100        continue

c..   set with iterated values
c..   yhn and yh2 are fractions of H, as is yh
            if( yht .gt. 1.0d-15 )then
               yhn       = yei*ye3/fh /yht
               yh2       = yhn**2/phih2 *yht
            else
               yhn = 0.0d0
               yh2 = 0.0d0
            endif
c..   unionized helium
            ya0    = yei  * yei  / bottom
c..   singly ionized helium
            ya1    = yei  * fhe1 / bottom
c..   doubly ionized helium
            ya2    = fhe2 * fhe1 / bottom
c..   renormalize for safety
            ya0 = ya0/(ya0 + ya1 + ya2)
            ya1 = ya1/(ya0 + ya1 + ya2)
            ya2 = ya2/(ya0 + ya1 + ya2)
            ya  = ya1 + ya2

         else
c..no H or He
            pion = 0.0d0
            eion = 0.0d0
            sion = 0.0d0
            ya0  = 0.0d0
            ya1  = 0.0d0
            ya2  = 1.0d0
            yh   = 1.0d0
            yh2  = 0.0d0
            yhn  = 0.0d0
            yei  = yht + yat + yze
         endif
      endif

c..   free energies for weak and strong screening

      if( yht .gt. 1.0d-15 )then

c..   fw is the part of the helmholtz free energy sum 
c..   for hydrogen coulomb interactions 
c..   the fw's sum to the helmholtz free energy
c..   Z**2 + Z = 2 here

c..   for H
         if( ifcoul .eq. 1 )then
            fw  = - wweak  * cergs *1.0d-6 *yht * yh * 2.0d0
            fs  = - strong * cergs *1.0d-6 *yht * yh
            fw2 = fw**2
            fs2 = fs**2
            dd  = sqrt( fw2 + fs2 )
            fc  = - fw*fs/dd
            pw  = 0.5d0 * d *fw
            ps  = d * fs / 3.0d0

            sw  = 0.5d0 * fw /( rgas * tt )
            ss  = 0.0d0
            ew  = fw + sw * rgas * tt
            es  = fs
            pc  = (-fs/dd)**3 * pw + (-fw/dd)**3 * ps
            sc  = (-fs/dd)**3 * sw + (-fw/dd)**3 * ss
            ec  = (-fs/dd)**3 * ew + (-fw/dd)**3 * es
            
         else
            pc = 0
            sc = 0
            ec = 0
         endif
         pch = pc
         sch = sc
         ech = ec

c..   excluded volume H0 atom
         pxh0 = rgas*tt*d *yht*yhn *d*volsum 
         sxh0 = -yht*yhn *d*volsum 
c..   excluded volume H2 molecule
         pxhm = rgas*tt*d *yht*yh2 *d*volsum
         sxhm = -yht*yh2 *d*volsum 
c..   pressure of H
         ph   = rgas*tt*d* yht*(yh+yh2+yhn) 
     1        + pch
     2        + pxh0 + pxhm
c..   theta has no statisical factors so they must be added
c..   hydrogen u/kT
c..   ionized hydrogen (g=1)
         etahi = dlog( yht*yh        *theta )
c..   atomic hydrogen (g=2=parth)
         etaha = dlog( yht*yhn/parth *theta )
c..   molecular hydrogen (add partition functions)
         etahm = dlog( yht*yh2/(zrot*zvib*zmel)*mamu3h*theta )
c..   extra entropy NOT in chemical potential
         srot  = yht*yh2* trsat/(trsat + tt)
         svib  = yht*yh2* tvib/tt/(exp(tvib/tt/tvsat*(tvsat+tt))-1.0d0)
         smel  = yht*yh2* dzmeldt 

         shi   =  yht*yh *( 2.5d0 - etahi ) + sch
         sha   =  yht*yhn*( 2.5d0 - etaha ) + yht*yhn*eparth 
         shm   =  yht*yh2*( 2.5d0 - etahm ) + srot + svib + smel
         sh    =  shi + sha + shm + sxh0 + sxhm
c..   ionization energy of metals ignored in eos
c..   eion goes to 0.0 when plasma is neutral atoms
c..   and is negative when H molecules form
         ehi  =  yht*( yh *( qh                  + 1.5d0*rgas*tt ) )
     1        + ech
         eha  =  yht*( yhn*( eparth*rgas*tt*yhn  + 1.5d0*rgas*tt ) )
         ehm  =  yht*( yh2*(-qh2                 + 1.5d0*rgas*tt  
     1        + rgas*tt*( srot + svib + smel ) ) )
c..   net energy from excluded volume is zero
      else
         ph  = 0.0d0
         sh  = 0.0d0
         ehi = 0.0d0
         eha = 0.0d0
         ehm = 0.0d0
      endif

c..   free energies for weak and strong screening
      if( yat .gt. 1.0d-15 )then
c..   unionized helium
c..   excluded volume He0
            pxhe0 = rgas*tt*d *yat*ya0 *d*volsum
            sxhe0 = -yat*ya0           *d*volsum
c..   singly ionized helium
         if( yat*ya1 .gt. 1.0d-15 )then
c..   Z**2 + Z = 2 here
c..   He+
            if( ifcoul .eq. 1 )then
               fw  = - wweak  * cergs *1.0d-6 *yat * ya1 * 2.0d0
               fs  = - strong * cergs *1.0d-6 *yat * ya1
               fw2 = fw**2
               fs2 = fs**2
               dd  = sqrt( fw2 + fs2 )
               fc  = - fw*fs/dd
               pw  = 0.5d0 * fw * d
               ps  = d * fs /3.0d0
               sw  = 0.5d0 * fw /( rgas * tt ) 
               ss  = 0.0d0
               ew  = fw + sw * rgas * tt
               es  = fs
               pc  = (-fs/dd)**3 * pw + (-fw/dd)**3 * ps
               sc  = (-fs/dd)**3 * sw + (-fw/dd)**3 * ss
               ec  = (-fs/dd)**3 * ew + (-fw/dd)**3 * es
            else
               pc = 0
               sc = 0
               ec = 0
            endif

c..excluded volume He+
            pxhe1 = rgas*tt*d *yat*ya1 *d*volsum
            sxhe1 = -yat*ya1           *d*volsum
            pca1 = pc
            sca1 = sc
            eca1 = ec
         else
c..   use benign values if small abundance (avoid overflow)
            fw    = 0.0d0
            pxhe1 = 0
            sxhe1 = 0
            pca1  = 0.0d0
            eca1  = 0.0d0
            sca1  = 0.0d0
         endif
c..   doublely ionized helium
         if( yat*ya2 .gt. 1.0d-15 )then
c..   Z**2 + Z = 6 here
c..   He++
            if( ifcoul .eq. 1 )then
               fw  = - wweak  * cergs *1.0d-6 *yat * ya2 * 6.0d0
               fs  = - strong * cergs *1.0d-6 *yat * ya2
     1              * 2.0d0**(5.0d0/3.0d0)
               fw2 = fw**2
               fs2 = fs**2
               dd  = sqrt( fw2 + fs2 )
               fc  = - fw*fs/dd
               pw  = 0.5d0 * fw * d
               ps  = d * fs /3.0d0
               sw  = 0.5d0 * fw /( rgas * tt ) 
               ss  = 0.0d0
               ew  = fw + sw * rgas * tt
               es  = fs
               pc  = (-fs/dd)**3 * pw + (-fw/dd)**3 * ps
               sc  = (-fs/dd)**3 * sw + (-fw/dd)**3 * ss
               ec  = (-fs/dd)**3 * ew + (-fw/dd)**3 * es
            else
               pc = 0
               sc = 0
               ec = 0
            endif

            pca2 = pc
            sca2 = sc
            eca2 = ec
         else
c..benign values if small abundance
            fw    = 0.0d0
            pxhe1 = 0.0d0
            sxhe1 = 0.0d0
            pca2  = 0.0d0
            eca2  = 0.0d0
            sca2  = 0.0d0
         endif


c..   pressure of He 
         phe   = rgas*tt*d* yat + pca1 + pca2 + pxhe0 + pxhe1
c..   helium u/kT
         etahe0 = dlog( yat*ya0   *theta )
         etahe2 = dlog( yat*ya2   *theta )
         etahe1 = dlog( yat*ya1/2.0d0  *theta )
c..   extra entropy NOT in chemical potential
c..   helium entropy for all ionization stages
         she   = yat*( ya1*( 2.5d0 - etahe1)
     1        +        ya2*( 2.5d0 - etahe2) 
     2        +        ya0*( 2.5d0 - etahe0 ) ) + sca1 + sca2
     3        + sxhe0 + sxhe1
c..   heII partition function constant
         ehe  = yat*( ya1*qa1 + ya2*(qa1+qa2)) + yat*1.5d0*rgas*tt
     1        + eca1 + eca2
c..   net energy term from excluded volume is zero
      else
         phe = 0.0d0
         she = 0.0d0
         ehe = 0.0d0
      endif
c..pressure of H and He
      pion = ph + phe 
c..entropy of H and He
      sion = sh + she 
c..energy/unit mass of H and He
      eion = ehi + eha + ehm + ehe

      return
      end


