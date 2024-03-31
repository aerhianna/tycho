
      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      implicit none
      save

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax

c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar

c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),zbarxx,fxtytot1,abar,zbar

      zbarxx  = 0.0d0
      fxtytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       fxtytot1    = fxtytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/fxtytot1
      zbar   = zbarxx * abar
      return
      end

      subroutine eosfxt(temp, den, zone, nc, csound)
      implicit none
      save
      include 'vector_eos.dek'

      include 'dimenfile'
      include 'timmes.inc'
      include 'comod'
      include 'ceoset'

c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

c..declare
      double precision pi,amu,kerg,clight,avo,qe,planck,ssol,asol
      parameter       (pi      = 3.1415926535897932384d0,
     1                  amu    = 1.6605402d-24,
     2                  kerg   = 1.380658d-16,
     3                  clight = 2.99792458d10, 
     4                  avo    = 6.0221367d23,
     5                  qe     = 4.8032068d-10,  
     6                  planck      = 6.6260755d-27,
     7                  ssol   = 5.67051d-5,
     8                  asol   = 4.0d0 * ssol / clight)

      integer          i,j,zone,nc
      double precision fxtx,deni,tempi,xni,
     1                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     4                 kt,ktinv,
     5                 xnem,pele,eele,sele,                 
     7                 etaele,
     8                 detadt,detadd,xnefer,dxnedt,dxnedd,s,
     9                 temp,den,fxtytot1,fxtye,csound,
     &                 sioncon,forth,forpi,kergavo,ikavo,asoli3,light2

      parameter        (sioncon = (2.0d0*pi*amu*kerg)/(planck*planck),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)

c..for the abar derivatives
      double precision dpepda,deepda,dsepda,                 
     1                 detada,dxneda


c..for the zbar derivatives
      double precision dpepdz,deepdz,dsepdz,
     1                 detadz,dxnedz


c..for the tables, in general
      integer          imax,jmax
      parameter        (imax = 211, jmax = 71)
      double precision d(imax),fxtt(jmax)

c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)

c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),
     2                 dpdfdt(imax,jmax)

c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),
     2                 efdt(imax,jmax)

c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),
     2                 xfdt(imax,jmax)

c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_tt,df_dt
      double precision fxtdth,dt2,fxtdti,dt2i,fxtdd,dd2,ddi,dd2i,xt,
     1                 xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md

c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)


c..for the coulomb corrections
      double precision a1,b1,c1,d1,e1,a2,b2,c2,
     1                 third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third = 1.0d0/3.0d0,
     9                  esqu  = qe * qe)
      logical tobe
c..for initialization
      integer          ifirst
      data             ifirst/0/ 

c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt

c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0d0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)

c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt

c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1

c..open the table
      inquire(file='helm_table.dat',exist=tobe)
      if( tobe )then
         open(unit=2,file='helm_table.dat',status='old')
      else
         write(*,*)'timmes: no file helm_table.dat in directory'
         stop'timmes: no .in file error'
      endif

c..read the helmholtz free energy table
       tlo   = 4.0d0
       thi   = 11.0d0
       tstp  = (thi - tlo)/dble(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -10.0d0
       dhi   = 11.0d0
       dstp  = (dhi - dlo)/dble(imax-1)
       dstpi = 1.0d0/dstp
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        fxtt(j) = 10.0d0**tsav
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**dsav
         read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo

c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        fxtdth      = fxtt(j+1) - fxtt(j)
        dt2         = fxtdth * fxtdth
        fxtdti      = 1.0d0/fxtdth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = fxtdth
        dt2_sav(j)  = dt2
        dti_sav(j)  = fxtdti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        fxtdd       = d(i+1) - d(i)
        dd2         = fxtdd * fxtdd
        ddi         = 1.0d0/fxtdd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = fxtdd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo

       close(unit=2)
       write(6,*)
       write(6,*) 'finished reading HELMHOLTZ eos table'
       write(6,04) 'imax=',imax,' jmax=',jmax
       write(6,03) 'temp(1)   =',fxtt(1),' temp(jmax)   =',fxtt(jmax)
       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
       write(6,*)

      end if

c..start of vectorization loop, normal execution starts here
      eosfail = .false.
c      do j=jlo_eos,jhi_eos
cc      write(*,*)temp,den
       if (temp .le. 0.0) stop 'temp less than 0 in helmeos'
       if (den  .le. 0.0) stop 'den less than 0 in helmeos'

c       temp  = temp_row(j)
c       den   = den_row(j)
c       abar  = abar_row(j)
c       zbar  = zbar_row(j)
      
       fxtytot1 = 1.0d0/abar
       fxtye    = fxtytot1 * zbar

c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt

c..electron-positron section:
c..assume complete ionization 
       xni     = avo * den * fxtytot1
       xnem    = xni * zbar

c..enter the table with fxtye*den
       din = fxtye*den

c..bomb proof the input
       if (temp .gt. fxtt(jmax)) then
        write(6,01) 'temp=',temp,' fxtt(jmax)=',fxtt(jmax)
        write(6,*) 'temp too hot, off grid'       
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (temp .lt. fxtt(1)) then
        write(6,01) 'temp=',temp,' fxtt(1)=',fxtt(1)
        write(6,*) 'temp too cold, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .gt. d(imax)) then
        write(6,01) 'den*fxtye=',din,' d(imax)=',d(imax)
        write(6,*) 'fxtye*den too big, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .lt. d(1)) then
        write(6,01) 'fxtye*den=',din,' d(1)=',d(1)
        write(6,*) 'fxtye*den too small, off grid'
        write(6,*) 'setting eosfail to true and returning'
c        call flush(6)
        eosfail = .true.
        return
       end if

c..hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))


c..access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)


c..various differences
       xt  = max( (temp - fxtt(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

c..the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)

       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)

c..derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)

       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)

       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)

       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)

c..second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)

c       ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c       ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c       ddsi2d =   ddpsi2(xd)

c       ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c       ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c       ddsi2md =  ddpsi2(mxd)


c..the free energy
       free  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density
       df_d  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..derivative with respect to temperature
       df_t = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density**2
c       df_dd = h5(iat,jat,
c     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2         ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

c..derivative with respect to temperature**2
       df_tt = h5(iat,jat,
     1       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to temperature and density
       df_dt = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)

       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)

       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)

       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)


c..derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)


c..look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

c..pressure derivative with density
       dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(fxtye * dpepdd,0.0d0)

c..look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)

c..electron chemical potential etaele
       etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
       fxtx       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = fxtye * fxtx

c..derivative with respect to temperature
       detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      detada = -fxtx * din * fxtytot1
      detadz =  fxtx * den * fxtytot1

c..look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

c..electron + positron number densities
      xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
      fxtx        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
      fxtx = max(fxtx,0.0d0)
      dxnedd   = fxtye * fxtx

c..derivative with respect to temperature
      dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      dxneda = -fxtx * din * fxtytot1
      dxnedz =  fxtx  * den * fxtytot1
       
c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of this one
       fxtx       = din * din
       pele    = fxtx * df_d
       dpepdt  = fxtx * df_dt
c       dpepdd  = fxtye * (fxtx * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/fxtye - 2.0d0 * din * df_d
       dpepda  = -fxtytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*fxtytot1*(2.0d0 * din * df_d  +  s)

       fxtx       = fxtye * fxtye
       sele    = -df_t * fxtye
       dsepdt  = -df_tt * fxtye
       dsepdd  = -df_dt * fxtx
       dsepda  = fxtytot1 * (fxtye * df_dt * din - sele)
       dsepdz  = -fxtytot1 * (fxtye * df_dt * den  + df_t)


       eele    = fxtye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = fxtx * df_d + temp * dsepdd
       deepda  = -fxtye * fxtytot1 * (free +  df_d * din) + 
     1     temp * dsepda
       deepdz  = fxtytot1* (free + fxtye * df_d * den) + temp * dsepdz


c..store this row
        pe(zone)      = pele
        pet(zone)     = dpepdt
        pev(zone)     = -dpepdd / v(nc,zone) * den
        ee(zone)      = eele
        eet(zone)     = deepdt
        eev(zone)     = -deepdd / v(nc,zone) * den
        sel(zone)     = sele/kergavo

        sele_row(zone)   = sele 
        spos_row(zone)   = 0.0d0
        dsept_row(zone)  = dsepdt 
        dsepd_row(zone)  = dsepdd 
        dsepa_row(zone)  = dsepda        
        dsepz_row(zone)  = dsepdz

        xnem_row(zone)   = xnem
        xne_row(zone)    = xnefer
        dxnet_row(zone)  = dxnedt
        dxned_row(zone)  = dxnedd
        dxnea_row(zone)  = dxneda
        dxnez_row(zone)  = dxnedz
        xnp_row(zone)    = 0.0d0

        etaele_row(zone) = etaele
        detat_row(zone)  = detadt
        detad_row(zone)  = detadd
        detaa_row(zone)  = detada
        detaz_row(zone)  = detadz
        etapos_row(zone) = 0.0d0


c..end of vectorization loop
c      enddo
      return
      end


