c<<<<<<< .mine
c
c     
      
c=======

     
      
c>>>>>>> .r98
      subroutine state(kin,ktop,nc)
c<<<<<<< .mine

c     
c     equation of state
c     
c     uses HELMHOLTZ table if T > 10**7.7 K
c     includes ideal ion gas, radiation pressure, and noninteracting
c     ferm-dirac gas of electrons and positrons (Timmes version)
c     in the equation of state
c     
c     coulomb: weak debye+salpeter strong screening 
c     via free energy interpolation
c     
c     uses OPAL EOS if T < 10**8 K
c     uses mass fraction as composition variable in x array
c     
c     interpolates OPALEOS and HELMHOLTZ if 10**7.7 K < T < 10**8 K
c     
c     OPAL opacity, both type 1 and type 2
c     schatzmann electron conduction scaled to Hubbard-Lampe
c     uses type 1 opacities (z interpolation) if nopac = 1
c     uses type 2 opacities (alpha enhanced) if nopac = 0
c     uses thomson if nopac = -1


c     ----------------
c     REVISION HISTORY:
c     -----------------
c     4-17-08 wda, interpolate HELMHOLTZ and OPALEOS
c     4-01-08 wda, add L to call sequence for opal
c
c     
c     TODO: THIS ROUTINE NEEDS SERIOUS SPRING CLEANING
c     (C.MEAKIN, 9-29-2008)
c     


c=======

c     
c     equation of state
c     
c     uses HELMHOLTZ table if T > 10**7.7 K
c     includes ideal ion gas, radiation pressure, and noninteracting
c     ferm-dirac gas of electrons and positrons (Timmes version)
c     in the equation of state
c     
c     coulomb: weak debye+salpeter strong screening 
c     via free energy interpolation
c     
c     uses OPAL EOS if T < 10**8 K
c     uses mass fraction as composition variable in x array
c     
c     interpolates OPALEOS and HELMHOLTZ if 10**7.7 K < T < 10**8 K
c     
c     OPAL opacity, both type 1 and type 2
c     schatzmann electron conduction scaled to Hubbard-Lampe
c     uses type 1 opacities (z interpolation) if nopac = 1
c     uses type 2 opacities (alpha enhanced) if nopac = 0
c     uses thomson if nopac = -1


c     ----------------
c     REVISION HISTORY:
c     -----------------
c     4-17-08 wda, interpolate HELMHOLTZ and OPALEOS
c     4-01-08 wda, add index L to call sequence for opal
c
c     
c     TODO: THIS ROUTINE NEEDS SERIOUS SPRING CLEANING
c     (C.MEAKIN, 9-29-2008)
c     


c>>>>>>> .r98
      implicit none
      save
      
      include 'dimenfile'
      include 'comod'
      include 'compu'
      include 'ceoset'
      include 'cconst'
      include 'cburn'
      include 'timmes.inc'
      include 'cgen'
      include 'cnabla'


      integer*4 kin,ktop,nc
c..   flag for opal eos initialization
      integer*4 ieos2
      character*10 cdummy
      character*11 ch11

c..   opal number fraction, mass fraction, solar equivalent mass fraction
      real*8 xnum(28),xfra(28),xelem(28)
      character*2 celem(28)
      integer*4 inum

      real*8 chemcon,chemfak,chema(ndim)

c...  used in eos subroutine
      integer*4 ikd(kdm),ikt(kdm), intab(kdm)
      common/ineos/ikd,ikt,intab

c..   local variables
      integer*4 k,kd,kkt,n,j
      
      real*8    ffkd,ffkt
      real*8    t1,v1,yht,xht,yat,ynuc,yenuc,facti,cdel,t2,v2
      real*8    yei,yei1,yei2,ya1,ya2
      real*8    eion,eion1,eion2
      real*8    yh,ddt,ddv,eiont,eionv,yeit,yeiv
      real*8    t3,t4
      real*8    uion1,uion2,sion,pion
      real*8    sion1,sion2,pion1,pion2
      real*8    ioh
c..   array for OPALEOS data before interpolation
      real*8    ppp(kdm),pppt(kdm),pppv(kdm),eee(kdm),
     1     eeet(kdm),eeev(kdm),sss(kdm)

c..   opacity scaling for solar model
      real*8    factak,akalpha,aklam,akrzero,aktzero,akgam,akscale

c..   local variables for H molecule
c..   yh2 = fraction of H which is molecular (0<yh2<0.5)
      real*8    yh2

c..   local variables for coulomb free energy
      real*8    cfw(kdm),cfs(kdm),cf(kdm),aodeb(kdm),pfrm(kdm)
      real*8    aodebye,xbhc,weak,strong

c..   interpolation variables for smooth join at table edges
      real*8    pegas,pegst,pegsv,eegas,eegst,eegsv,sen,si,fact1,fact2
      real*8    pp,ppt,vpv,eel,vev,dtdv,cc,ett
      real*8    thelm,topal

      real*8    d1,d2, xz,yze,xxc,xxo,t6,ropal
      real*8    zeta,xpf,sumex
      real*8    akc,akct,akcv,ako,akot,akov,x2,sumx,sumz,varz
      real*8    sumc,sumn,sumo,sumfe
      real*8    tydxc,tydxo
      real*8    thomson
      
      real*8    fivethrds(ndim),la(ndim),lz2(ndim),xms(ndim),yms(ndim)
      real*8    csound

      integer*4 nerfd
      parameter( nerfd = 29 )
      real*8    bsi(nerfd),bj2(nerfd),bj3(nerfd),bphi(nerfd)
      common/nrfddata/bsi,bj2,bj3,bphi

      integer*4 nii

      logical   tobe

      data lc12/-1/,lo16/-1/
      data ieos/-1/,ieos2/-1/
      data ioh/13.598d0/
c<<<<<<< .mine



c=======

c------------------------------
c     write(*,*)'ENTERING STATE ',kin,ktop,nc
c------------------------------

c>>>>>>> .r98
c------------------------------------------------------------------
c     
c..   regions for different eos procedures
c
c  1e-15   1e-10 g/cc             1e8 g/cc    1e+11 g/cc   rho*Ye
c     
c-------------------+-----------------------------+---------T=1e11 K
c          |   region 2                           |
c          |                   HELMHOLTZ table    |
c--+-------+--------------------------------------+---------T=1e8 K = thelm
c  |       |  region 1   INTERPOLATE |            .
c--+-------+--------------------------------------+-------T=1e7.7 K = topal
c  |                              |              
c  |          region 0          |                           T=1e4 K
c  |          OPALEOS table   |
c  |                        |
c--+--------------------------------------------------------T=2e3 K
c 
c  1e-15   1e-10 g/cc              1e8 g/cc    1e+11 g/cc   rho*Ye
c
c     ^
c     | log T
c     |            
c     |
c     +----log rhoz*Ye -->

c..   Timmes: 
c     log10 T(low)        = 4
c     log10 T(high)       = 11
c     71 values of T(K)
c     log10 rhoz*Ye(low)  = -10
c     log10 rhoz*Ye(high) =  11
c     211 values of rhoz*Ye(g/cc)




c..Use OPALEOS where available, HELMHOLTZ for advanced stages (log T> 8)
c..   OPALEOS (not square in T, rhoz)
c     log10 T(low)        = 3.3
c     log10 T(high)       = 8
c
c     log10 rhoz(low)     -15      8 > log T > 3.3
c     log10 rhoz(high)      7      8 > log T > 7.352
c------------------------------------------------------------------
c..   upper and lower limits in T for Helmholtz-OPALEOS interpolation
      thelm = 10.0d0**7.9d0
      topal = 10.0d0**7.7d0


c..   time this operation and accumulate it in runstate

      call seconds(runstate0)


c...........................................
c..   CHECK INPUT:
c..   catch undefines and impossible values
c...........................................
      do k = kin, ktop
         if( v(nc,k) .lt. 1.0d-20 .or. v(nc,k) .gt. 1.0d25 )then
            write(*,'(/a30,3i5,1p8e12.3)')
     1           'STATE: bad input: k, kk, nc, V',k,kk,nc,v(nc,k)
            if( k .le. kk )then
               write(*,*)'called on Henyey grid'
            else
               write(*,*)'called by Envelope: fitenv branch'
               write(*,*)kin,ktop,k,nc,v(nc,k)
            endif
            stop'STATE input V'
         endif
      enddo
      do k = kin, ktop
         if( t(nc,k) .lt. 1.0d0 .or. t(nc,k) .gt. 1.0d15 )then
            write(*,'(/a30,3i5,1p8e12.3)')
     1           'STATE: bad input: k, kk, nc, T',k,kk,nc,t(nc,k)
            if( k .le. kk )then
               write(*,*)'called on Henyey grid'
            else
               write(*,*)'called by Envelope: fitenv branch'
            endif
            stop'STATE input T'
         endif
      enddo



c.........................................................................
c..   BEGIN INITIALIZATION:
c.........................................................................
      if( ieos .ne. 0 )then
         print*,'MSG(state): INIT EOS, OPACITY, CONDUCTIVITY'

         ieos  = 0
c..   flag to initialize OPALEOS also
         ieos2 = 1
c<<<<<<< .mine
         

c.........................................
c     COMMENTED: NOT CURRENTLY BEING USED 
c     USING HELMHOLTZ
c..   
c..   read data for ndfd (nondegenerate fermi-dirac) eos
c     open(14,file='eos1.d',status='old')
c     rewind 14
c     read(14,'(f5.1,f12.5,f12.5,f12.5)') 
c     1        ( bsi(j),bj2(j),bj3(j),bphi(j), j=1, nerfd)
c     close(14)                            
c     write(*,'(/a50)')'initializing state: file eos1.d is read'
c.........................................


c......................................
c..   SETUP nucleon number if needed
c......................................
c=======
         
c>>>>>>> .r98

c..   SETUP nucleon number if needed
c......................................
         if( xa(1) .eq. 0.0d0 )then
c..   xa array is not initialized to nucleon number
            write(*,*)'xa not defined, read net.rc first?'
            stop'state: xa'
         endif
         if( qex(netsize-1) .eq. 0.0d0 )then
            write(*,*)'qex not defined'
            stop'state: qex'
         endif


c.....................................
c..   chemical potentials for nuclei
c.....................................
         chemcon = dlog(
     .        avagadro*(planck**2*avagadro/(2.0d0*pi*boltz))**1.5d0 )
c..   set up ln(g * A**3/2) for nuclei chemical potentials
         do n = 1, netsize
            chema(n) = 1.5d0 * dlog( xa(n) )
         enddo
c..   does not include nuclear spin and partition function except for H
         chema(netsize-1) = dlog( 2.0d0 )
c..   and electrons
         chema(netsize+1) = dlog( 2.0d0*(avagadro*egrestm)**1.5d0 )
c..   set up Z**5/3 factors for strong screening eos used below
         do j = 1, netsize
            fivethrds(j) = dble( lz(j) )**(5.0d0/3.0d0)            
c..   la and lz2 are set up as real*8 for timmes azbar.f
c..   which wants the atomic mass and charge 
c..   need masses (xa) for accurate fractions by mass
            la(j)  = xa(j)
            lz2(j) = lz(j)
         enddo




         if( nopac .eq. -1 )then
            write(*,*)'MSG(state): forcing thomson opacity ',nopac
         else

c......................................................................
c..   SEARCH FOR c12, o16 INDICES FOR OPAL IF NEEDED (nopac != -1)
c..   SEARCH FOR d,li7,be9,b10,b11 for loburn and rezone
c..   (should this code be put in build.f ?)
c..
            do j = 1, netsize
               if( lz(j) .eq. 6 .and. ln(j) .eq. 6 )then
                  lc12 = j
               endif
               if( lz(j) .eq. 7 .and. ln(j) .eq. 7 )then
                  ln14 = j
               endif
               if( lz(j) .eq. 8 .and. ln(j) .eq. 8 )then
                  lo16 = j
               endif
               if( lz(j) .eq. 1 .and. ln(j) .eq. 1 )then
                  ldeut = j
               endif
               if( lz(j) .eq. 2 .and. ln(j) .eq. 1 )then
                  lhe3 = j
               endif
               if( lz(j) .eq. 3 .and. ln(j) .eq. 4 )then
                  lli7 = j
               endif
               if( lz(j) .eq. 4 .and. ln(j) .eq. 5 )then
                  lbe9 = j
               endif
               if( lz(j) .eq. 5 .and. ln(j) .eq. 5 )then
                  lb10 = j
               endif
               if( lz(j) .eq. 5 .and. ln(j) .eq. 6 )then
                  lb11 = j
               endif
               if( lz(j) .eq. 13 .and. ln(j) .eq. 13 )then
                  lal26 = j
               endif
            enddo
            if( lc12 .le. 0 .or. lo16 .le. 0 )then
               write(*,*)'ERR(state): no c12/o16 for opal:',lc12,lo16
               stop'ERR(state): no c12/o16 for opal'
            endif
            if( ldeut .le. 0 .or. lli7 .le. 0 .or. lbe9 .le. 0 .or.
     1           lb10 .le. 0 .or. lb11 .le. 0 )then
               write(*,*)'WARN(state): no dlibeb for loburn ',ldeut,
     1              lli7,lbe9,lb10,lb11
c..   no dlibeb burning (loburn not called)
               dlibeb = 1
            else
c..   allow dlibeb burning (default)
               dlibeb = 0
            endif
c..

            
c......................................................................
c     SETUP OPAL OPACITY TABLES 
c     type 2 (alpha enriched)
c     FILE(S): codataa
c......................................................................
            if( nopac .eq. 0 )then
c..   check metallicity of opacity tables
               inquire(file='codataa',exist=tobe) 
               if( .not. tobe )then
                  write(*,*)'No file codataa found'
                  stop'state.f codataa'
               endif
               open(20,file='codataa',status='old')
               do j = 1, 62
                  if( j .ge. 32 .and. j .le. 59 )then
                     inum = j-31
c..   inum is charge number Z
c..   celem is symbol for element
c..   xnum is fraction of metals by number
c..   xfra is fraction of metals by mass
                     read(20,'(3x,a2,a11,27x,0pf8.6,8x,0pf8.6)'),
     .                    celem(inum),ch11,xnum(inum),xfra(inum)
                  else
                     read(20,'(a10)')cdummy
                  endif
               enddo
               read(20,'(13x,a12,26x,0pf6.5)')copal,zmetal
               print*,'MSG(state): Z IN OPAL TYPE 2 TABLES'
               write(*,'(a21,a14,a20,0pf6.5)')'OPAL opacity file 2:',
     .              copal,'metallicity =',zmetal
               close(20)
               

c.........................................................................
c     SETUP OPAL OPACITY TABLES
c     type1 (interpolae in metallicity)
c     FILE(S): hzdata
c.........................................................................
            elseif( nopac .eq. 1 )then
c..   type 1 opacity tables (interpolate in metallicity)
               inquire(file='hzdata',exist=tobe) 
               if( .not. tobe )then
                  write(*,*)'No file hzdata found'
                  stop'state.f hz'
               endif
               write(*,*)'MSG(state): opening type 1 opac file: hzdata'
               open(20,file='hzdata',status='old')
               
               do j = 1, 62
                  if( j .ge. 32 .and. j .le. 59 )then
                     inum = j-31
c..   inum is charge number Z
c..   celem is symbol for element
c..   xnum is fraction of metals by number
c..   xfra is fraction of metals by mass
                     read(20,'(3x,a2,a11,27x,0pf8.6,8x,0pf8.6)'),
     .                    celem(inum),ch11,xnum(inum),xfra(inum)
                  else
                     read(20,'(a10)')cdummy
                  endif
               enddo
               

c.......................................
c..   SETUP zmetal FOR OPAL1:
c.............................
c..   (ADJUST FOR MIXMODE)
c..   * HOW GENERAL IS THIS MIXMODE STUFF?
c..   * IS THIS THE BEST PLACE FOR THIS?
c........................................
               if( mixmode .eq. 2 .or. mixmode .eq. 0 )then
c..   diffusion, use original metallicity
                  zmetal = zpop0
                  write(*,*)'MSG(state): metal zpop0 is ',zmetal
               else
c..   use metallicity of outer zone (kk) in mesh
                  zmetal = 0.0d0
                  do j = 1, nnuc
                     if( lz(j) .gt. 2 )then
                        zmetal = zmetal + x(j,kk)
                     endif
                  enddo
                  write(*,*)'MSG(state): metal in zone kk is ',zmetal
                  write(*,*)mixmode
               endif

c....................................................
c..   GET ELEMENTAL abundances for OPAL1 TABLE DATA
c....................................................
               read(20,'(13x,a12)')copal
               print*,'MSG(state):'
               write(*,'(a14,a14,a35)')'OPAL Type 1:',copal,
     1              'interpolate H, He, metallicity'
               write(*,'(a5,a5,3a10)')'Z','el','f(num)','f(mass)',
     1              'X(mass)'
               do j = 6, 28
                  write(*,'(i5,a5,0pf10.6,0p2f10.6)')j,celem(j),xnum(j),
     .                 xfra(j),xfra(j)*zmetal
               enddo
               if( xfra(8) .ne. 0.0d0 .and. xnum(8) .ne. 0.0d0 )then
                  write(*,'(2(a20,1pe12.3))')
     .                 'X(Fe)/X(O) ',xfra(26)/xfra(8),
     .                 'N(Fe)/N(O) ',xnum(26)/xnum(8)
               else
                  stop'error in OPAL data file'
               endif
               close(20)
            endif




c.........................................................
c..   COMPARE OPAL 1 & 2 ABUNDANCES WITH SOLAR
c..   construct solar elemental abundance from solarx data 
c.........................................................
c..   (mass fractions)
            do inum = 1, 28
               xelem(inum) = 0.0d0
c..   solar mass fractions for elements (normalized to 1)
               do j = 1, nnuc
                  if( inum .eq. lz(j) )then
                     xelem(inum) = xelem(inum) + solarx(j)
                  endif
               enddo
            enddo


c....................................
c..   sum of solar mass fractions 
c..   for metals (H and He omitted)
c....................................
            varz = 0
            do j = 3, 28
               varz = varz + xelem(j)
            enddo
            print*,'MSG(state): sum of mass fracs for metals'
            write(*,'(a20,3(a10,1pe12.3))')'metallicity:',
     .           'solarx',varz
            write(*,'(2(a5,a5,2a12,a15))')
     .           'Z','El','Xopal','Xsol','opal/sol',
     .           'Z','El','Xopal','Xsol','opal/sol'

            do inum = 3, 28, 2
               write(*,'(2(i5,a5,1p2e12.3,0pf15.10))')inum,celem(inum),
     .              xfra(inum)*zmetal,xelem(inum),
     .              xfra(inum)*zmetal/xelem(inum),
     .              inum+1,celem(inum+1),
     .              xfra(inum+1)*zmetal,xelem(inum+1),
     .              xfra(inum+1)*zmetal/xelem(inum+1)
            enddo

c..   save opal values for use in defining excess C and O (type 2)
            opalc =  xfra(6)*zmetal
            opaln =  xfra(7)*zmetal
            opalo =  xfra(8)*zmetal
            print*,'MSG(state): opal values of C*Z, N*Z and O*Z:'
            write(*,'(a10,4(a12,1pe12.4))')'OPAL','zmetal',zmetal,
     .           'opalC',opalc,'opalN',opaln,'opalO',opalo


c..   check nucleon conservation
            sumx  = 0.0d0
c..   mass fraction of nuclei above helium (lithium and up)
            sumz  = 0.0d0
            sumc  = 0.0d0
            sumn  = 0.0d0
            sumo  = 0.0d0
            sumfe = 0.0d0
            do j = 1, netsize
               sumx = sumx + x(j,kk)
               if( lz(j) .gt. 2 )then
                  sumz = sumz + x(j,kk)
               endif
               if( lz(j) .eq. 6 )then
                  sumc = sumc + x(j,kk)
               endif
               if( lz(j) .eq. 7 )then
                  sumn = sumn + x(j,kk)
               endif
               if( lz(j) .eq. 8 )then
                  sumo = sumo + x(j,kk)
               endif
               if( lz(j) .eq. 26)then
                  sumfe = sumfe + x(j,kk)
               endif
            enddo
            tydxc = sumc - sumfe*xelem(6)/xelem(26)
            tydxo = sumo - sumfe*xelem(8)/xelem(26)
            sumz  = sumz - (tydxc + tydxo)
            sumc  = sumc - tydxc
            sumo  = sumo - tydxo


c     print*,'MSG(state): mass fracs above He (Li and up):'
            write(*,'(a10,6(a12,1pe12.4))')'TYCHO(kk)',
     .           'zsum(ty)',sumz,
     .           'tychoC',sumc,
     .           'tychoN',sumn,
     .           'tychoO',sumo

            write(*,'(a10,4(a12,1pe12.4))')'TY-OPAL',
     .           'z',sumz-zmetal,
     .           'C',sumc-opalc,
     .           'N',sumn-opaln,
     .           'O',sumo-opalo

c            print*,'MSG(state): tycho sumc, sumo, sumfe =',
c     .           sumc,sumo,sumfe
            print*,'MSG(state): These C,O excess based on Fe sol ratio'
            print*,'MSG(state): tycho dXC = ', tydxc
            print*,'MSG(state): tycho dXO = ', tydxo

            if( abs(1.0d0-sumx) .gt. 1.0d-2 )then
               write(*,*)'ERR(state): init error: sumx-1 ',sumx-1.0d0,kk
               stop'ERR(state): sumx'
            endif



c..........................................................
c..   CHECK CONSISTENCY BETWEEN TYCHO Z and OPAL TABLES Z
c..   NOTE: metallicity .lt. 1.0d-8 is essentially zero
c..........................................................
            varz = (zmetal-sumz) /
     .           (zmetal+sumz+1.0d-8) * 2.0d0

            if( nopac .eq. 0 )then
               
               if( mixmode .eq. 1 .or. mixmode .eq. 3 )then
c.....................
c..   OPAL TYPE 2: no gravitational settling ("Thoul diffusion")
c.....................
                  if( abs( varz ) .gt. 1.0d-1  )then
                     print*,'ERR(state): Opacity/Tycho Z Discrepancy'
                     write(*,'(/a26,1pe10.2,a14,1pe10.2,a8,1pe10.2,
     .                    a4,i5)')
     .                    'STATE: ERROR: OPAL z =',zmetal,
     .                    ' .ne. BURN z =', sumz,
     .                    ' error =',sumz-zmetal,' kk=',kk
                     write(*,*)'Fractional error in metallicity is ',
     .                    varz
                     write(*,*)'Such a large error may give garbage'
                     write(*,*)'either:' 
                     write(*,*)
     .                    '  1. choose different set of codata? files'
                     write(*,*)
     .                    '  see makecodata script',
     .                    ' in ./opal2.dir/type2tab'
                     write(*,*)
     .                    '  2. redifine imodel abundances using genex'
                     write(*,*)'see ./src/genex.f and ./genex.in'
                     write(*,*)
     .                    'if you insist, change the ',
     .                    'abs(varz) criterion'
                     write(*,*)'in ./src/state.f, near line 306'
                     write(*,*)'NOTE: mixmode 2,4 have ',
     .                    'different criteria'
                     stop'ERR(state): OPAL/burn.f abundance mismatch'
                  else
c<<<<<<< .mine
                     print*,'MSG(state): sumz(opal)-zmetal(zone kk)=',
     .                    sumz-zmetal
c                     write(*,'(a48,1pe12.3,a20)')
c     .                    'State: OPAL z agrees with actual z to',
c     .                    sumz-zmetal,' in outer zone kk'
c=======
                     print*,'MSG(state): sumz(opal)-zmetal(zone kk)=',
     .                    sumz-zmetal
c>>>>>>> .r98
                  endif

               else
c..........................
c..   OPAL TYPE 2: with gravitational settling so metallicity changes
c..........................
                  write(*,'(a48,1pe12.3,a20)')
     .                 'STATE: OPAL z agrees with actual z to',
     .                 sumz-zmetal,' in outer zone kk'
                  write(*,'(/a24,i3)')'Mixing option mixmode =',
     .                 mixmode
                  if( nopac .eq. 0 )then
                     write(*,*)
     .                    'Gravitational settling (Thoul diffusion)'
                     write(*,*)'Enhanced oxygen will be ',sumz-zmetal,
     .                    ' above ',opalo
                     if(  varz .gt. 1.0d-1  )then
                        write(*,*)
     .                       'Metallicity is too low for accurate ',
     .                       'opacity with type 2 tables'
                        write(*,*)'Fractional error is ',varz
                        write(*,*)'Use type 1 tables (nopac=1), or'
                        write(*,*)'choose an OPAL CO table with',
     .                       ' metallicity of ',sumz,' or less'
                        write(*,'(a50,a15)')
     .                       'and higher metallicity will be ',
     .                       'approximated by extra oxygen'
                        write(*,*)'see state.f, about line 443'
                        stop'Diffusive mixing/OPAL table error'
                     endif
                  endif
               endif

c..   mass fraction of nuclei above helium (lithium and up) 
c..   at center (k=2)
               sumz = 0.0d0
               do j = 1, netsize
                  sumx = sumx + x(j,2)
                  if( lz(j) .gt. 2 )then
                     sumz = sumz + x(j,2)
                  endif
               enddo
c<<<<<<< .mine
               print*,'MSG(state): mfracs above He (Li and up)'
               print*,'MSG(state): NO correction for excess over sol:'
c=======
               print*,'MSG(state): mfracs above He (Li and up)'
               print*,
     .              'MSG(state): including Z contrib and enrichment'
c>>>>>>> .r98
               write(*,'(a10,4(a12,1pe12.4))')'TYCHO(2)',
     .              'zsum(ty)',sumz,
     .              'tychoC',x(lc12,2),
     .              'tychoN',x(ln14,2),
     .              'tychoO',x(lo16,2)
            endif
         endif
      endif


c<<<<<<< .mine




c......................................................................
c=======


c......................................................................
c>>>>>>> .r98
c..   load temperature, density, Ye arrays for eos, test for bad input
c......................................................................
      do k = kin, ktop
         tem(k) = t(nc,k)
         
c..   CALCULATE Ye
c..   (sumex = excess mass due to nuclear binding, zero is c12)
         sumex          = 0.0d0
         x(netsize+1,k) = 0.0d0
         do j = 1, netsize
            x(netsize+1,k) = x(netsize+1,k) + x(j,k)/xa(j)*dble( lz(j) )
            sumex = sumex + x(j,k)/xa(j)* qex(j)/931.487d0            
         enddo
         ye(k) = x(netsize+1,k)
         
c..   see Arnett, 1996, p. 8, mass fraction representation
c..   for network:
c     rhoz(k) = (1.0d0-sumex)/v(nc,k)
c..   use mass fractions for eos
         rhom(k) = 1.0d0/v(nc,k)
         rhoz(k) = rhom(k)
c..   compute indices in tablular EOS
c..   Timmes uses mass fractions
         denye(k) = ye(k)*rhom(k)
      enddo

c..   test for unphysical input data
      do k = kin, ktop
         if( denye(k) .le. 0.0d0 
     .        .or. ye(k) .lt. 0.0d0 .or. ye(k) .gt. 1.0d0 )then
            write(*,*)'state test: bad input Ye'
            write(*,'(5a5,5a12)')'model','it','k','kin','ktop',
     .           'T','denye','ye','rhoz','rhom'
            write(*,'(5i5,1p5e12.3)')model,it,k,kin,ktop,
     .           tem(k),denye(k),ye(k),rhoz(k),rhom(k)
            stop'state input Ye'
         endif
      enddo

c<<<<<<< .mine

c...................................................................
c=======


c...................................................................
c>>>>>>> .r98
c..   OPACITY AND CONDUCTION
c...................................................................

c..................
c     THOMSON
c..................
      if( nopac .eq. -1 )then
         do k = kin, ktop
            thomson = 0.399d0
            ako     = ye(k) * thomson
            akot    = 0.0d0
            akov    = 0.0d0
            
c....................................................
c..   electron conduction, schatzman, w.dwarfs, p83
c....................................................
            yht = x(netsize-1,k)/xa(netsize-1)
            yat = x(netsize,k)  /xa(netsize)
            xz  = 1.0d0 - yht*xa(netsize-1) - yat*xa(netsize)
            x2      = ( 1.015d-6 * denye(k) )**0.666667d0
            akc     = 1.77d-6*( 1.0d0 + x2 )*( tem(k) / rhom(k) )**2
     .           * ( 0.4d0*yat*4.0d0 + yht*0.1d0 + 0.8d0*xz )
c..   last factor adjusted from hubbard-lampe (1969, apjs 18,297)
c..   tables for H, He, C
c..   save for dout3.f output
            akap(k) = akc
            akct    = 2.0d0 * akc / tem(k)
            akcv    = rhom(k)*(2.0d0*akc -
     .           1.77d-6*(tem(k)/rhoz(k))**2 / 1.5d0 * x2 )
c..   combine in parallel (add by inverses)
            ak(k)   = 1.0d0/( 1.0d0/ako + 1.0d0/akc )
            akt(k)  = ak(k)**2*(akct/akc**2 + akot/ako**2 )
            akv(k)  = ak(k)**2*(akcv/akc**2 + akov/ako**2 )
         enddo

c<<<<<<< .mine
c.................................................
      else                      !OPAL
c.................................................
c=======
c>>>>>>> .r98

c<<<<<<< .mine
c=======
c.................................................
c      else                      !OPAL
c.................................................

c>>>>>>> .r98
         do k = kin, ktop
c..   opac composition variables yht,xz,xxc,xxo
            yht = x(netsize-1,k)/xa(netsize-1)
            yat = x(netsize,k)  /xa(netsize)
            xht = x(netsize-1,k)
c..   opal table metallicity
c..   xz    = zmetal
c..   actual zone metallicity. requires deuterium, he3 in network
            xz    = 1.0d0 - x(netsize-1,k) - x(netsize,k)
     .           -x(ldeut,k) - x(lhe3,k)
            xz    = dmax1( 0.0d0, xz )
            t6    = tem(k)*1.0d-6
            ropal = rhom(k)/t6**3
            

c....................................................
            if( nopac .eq. 0 )then !...OPAL TYPE 2
c....................................................
c..   extra c and o for type 2 tables
               xxc   = x(lc12,k) + x(ln14,k)*0.5d0
     .              - opalc -0.5d0*opaln
               xxo   = x(lo16,k) + x(ln14,k)*0.5d0
c<<<<<<< .mine
     2              - opalo - 0.5d0*opaln
               if( xz .lt. zmetal )then
c=======
c     .              - opalo - 0.5d0*opaln
c              if( xz .lt. zmetal )then
c>>>>>>> .r98
c..   tabular value is lowest metallicity available with type 2 tables
                  xxc = 0.0d0
                  xxo = 0.0d0
c..   THIS IS AN ERROR, MAYBE SMALL
               else
                  xxc = dmax1( 0.0d0, xxc )
c..   extra metals are taken to be "oxygen like"
                  xxo = xz - zmetal - xxc
               endif

c..   feed in table metallicity zmetal in type 2 opal tables, 
c..   extra metals in xxc and xxo
               call opal(t6,ropal,zmetal,xht,xxc,xxo,ako,akot,akov,
     .              k,kk,nopac,l)


c........................................................     
            elseif( nopac .eq. 1 )then !...OPAL TYPE 1
c........................................................
c..   type 1 tables interpolate for any z
               xxc = 0.0d0
               xxo = 0.0d0
c..   use actual metallicity of zone k (OPAL type 1 tables)
               call opal(t6,ropal,xz,xht,xxc,xxo,ako,akot,akov,
     1              k,kk,nopac,l)


               if( nkscale .eq. 1 )then                  
c..   scale opacity in interior
                  ako  = ako *1.1d0
                  akot = akot*1.1d0
                  akov = akov*1.1d0
                  
               elseif( nkscale .eq. 2 )then
c..   bahcall basu pinsonneault serenelli
                  akalpha = 0.21d0
                  aktzero = 2.18d6
                  akgam   = 0.1d0*aktzero
                  if( k .lt. kk )then
c..   scale opacity only in radiative region
                     factak = 1.0d0 +
     1                    akalpha*akgam**2/
     2                    ( (t(1,k) - aktzero)**2 + akgam**2 )
                  else
                     factak = 1.0d0
                  endif
c..   scale opacity in interior
                  akscale = ako * factak
                  akot = akot*factak
     1                 + ako*( -2.0d0*(t(1,k)-aktzero) )* (factak-1.0d0)
     2                 /( (t(1,k) - aktzero)**2 + akgam**2 )**2
                  akov = akov*factak
c..   now update ako
                  ako  = akscale
                  

c......................................................
               elseif( nkscale .eq. 3 )then
c..   g-mode waves
c..   scale with mixing length estimate
                  akalpha = 0.05d0*alphaml
c.................................
                  akrzero = 0.713d0*solrad
                  aklam   = 0.2d0 * akrzero
                  if( r(1,k) .gt. akrzero )then
                     factak = 1.0d0 + akalpha
                  else
                     factak = 1.0d0 
     1                    + akalpha * exp( (r(1,k)-akrzero)/aklam )
                  endif
                  ako  = ako *factak
                  akot = akot*factak
                  akov = akov*factak
               endif
            else
               write(*,*)nopac
               stop'nopac error in state.f'
            endif


c....................................................
c..   ELECTRON CONDUCTION: schatzman, w.dwarfs, p83
c....................................................
            x2      = ( 1.015d-6 * denye(k)  )**0.666667d0
            akc     = 1.77d-6*( 1.0d0 + x2 )*( tem(k) / rhom(k) )**2
     1           * ( 0.4d0*yat*4.0d0 + yht*0.1d0 + 0.8d0*xz )
c..   last factor adjusted from hubbart-lampe (1969, apjs 18,297)
c..   tables for H, He, C
c..   save for dout3.f output
            akap(k) = akc
            akct    = 2.0d0 * akc / tem(k)
            akcv    = 2.0d0 * akc * rhom(k) 
     1           *(  1.0d0 - x2/( 3.0d0* (1.0d0+x2) )   )
c..   combine in parallel (add by inverses)
            ak(k)   = 1.0d0/( 1.0d0/ako + 1.0d0/akc )
            akt(k)  = ak(k)**2*(akct/akc**2 + akot/ako**2 )
            akv(k)  = ak(k)**2*(akcv/akc**2 + akov/ako**2 )
         enddo
      endif





c--------------------------------------------------------------------
c..   EQUATIONS OF STATE
c..   use opaleos at log T<8, else wda+helmholtz eos
c---------------------------------------------------------------------

c..   nuclei (except hydrogen and helium) and radiation component
      do k = kin, ktop
c..   sum number of nuclei, except for atomic and molecular hydrogen
c..   and helium, which are treated separately
         yion(k) = 0.0d0
         do n = 1, netsize-2
            yion(k) = yion(k) + x(n,k)/xa(n)
         enddo
         t1    = tem(k)
         t2    = t1*t1
         t3    = t2*t1
         t4    = t2*t2
         v1    = v(nc,k)
         pnuc(k)  = rgas*t1/v1 *yion(k)
         pnuct(k) = rgas/v1    *yion(k)
         pnucv(k) = -pnuc(k)/v1
         enuc(k)  = rgas*t1*1.5d0 *yion(k)
         enuct(k) = rgas*   1.5d0 *yion(k)
         enucv(k) = 0.0d0
         if( norad .eq. 0 )then
c..   radiation pressure (default)
            prad(k)  = arad*t4/3.0d0
            pradt(k) = 4.0d0*arad*t3/3.0d0
            pradv(k) = 0.0d0
            erad(k)  = 3.0d0*prad(k)*v1
            eradt(k) = 3.0d0*pradt(k)*v1
            eradv(k) = 3.0d0*prad(k)
            srad(k)  = pradt(k)*v1/rgas
         else
c..   set radiation pressure to zero for comparisons
            prad(k)  = 0.0d0
            pradt(k) = 0.0d0
            pradv(k) = 0.0d0
            erad(k)  = 0.0d0
            eradt(k) = 0.0d0
            eradv(k) = 0.0d0
            srad(k)  = 0.0d0
         endif

c..   nuclear nse and coulomb terms are not included here, 
c..   only ideal gas of nuclei
c..   sum nuclear entropies
c..   snuc is the entropy of the nuclear component for all nuclei,
c..   ionized or not, excepting hydrogen
         snuc(k) = 0.0d0
         chemfak =  dlog( rhom(k)/tem(k)**1.5d0 )
         do n = 1, netsize-2
c..   sum all nuclei in atoms or ions 
c..   but atomic and molecular hydrogen and helium
c..   which are treated separately
            if( x(n,k) .gt. 1.0d-16 .and. n .ne. netsize-1 )then
c..   pay 9-4-06               .......
c..   etanuc is chemical potential/kT
c..   gspin is 2J+1 for nuclear ground stateg
               etanuc(n,k) = dlog( x(n,k)/xa(n)/gspin(n) ) +chemfak 
     1              -chema(n) +chemcon 
               snuc(k) = x(n,k)/xa(n)*(2.5d0 - etanuc(n,k)) + snuc(k)
            endif
         enddo
c..   add hydrogen and helium nuclei in ions or atoms
c..   electrons added below in ionization
c     do n = netsize-1,netsize
c     if( x(n,k) .gt. 0.0d0 )then
c..   etanuc is chemical potential/kT
c     etanuc(n,k) = dlog( x(n,k)/gspin(n) ) +chemfak 
c     1              -chema(n) +chemcon 
c     snuc(k) = x(n,k)*(2.5d0 - etanuc(n,k)) + snuc(k)
c     endif
c     enddo
      enddo





c......................................................................
c..   LOAD TEMP, DENS, Ye ARRAYS FOR EOS
c...  IDENTIFY EOS REGIONS
c..   (uses denye(k), as defined last above)
c......................................................................
      do k = kin, ktop

c..   0.1 is log10 increment in density per table element;
c..   -10.0 is log10 of first d value in table
         
         ffkd   = ( dlog10( denye(k) ) + 10.0d0 )*10.0d0 + 1.0d0
         ikd(k) = int( ffkd )
c..   0.1 is log10 increment in temperature per table element;
c..   4.0 is log10 of first t value in table
         ffkt   = ( dlog10( tem(k) ) - 4.0d0 )*10.0d0 + 1.0d0
         ikt(k) = int( ffkt )
         kd     = ikd(k)
         kkt    = ikt(k)
                  
         if( kkt .lt. 38 )then
c..   low temperature, pure OPALEOS (38--> log T = 7.7000)
            intab(k) = 0
            
         elseif( kkt .ge. 40 )then
c..   high temperature, pure helmholtz (40--> log T = 7.9000)
c..   high temperature, pure helmholtz (41--> log T = 8.0000)
            intab(k) = 2
            
         else
c..   interpolate
            intab(k) = 1
         endif
         
      enddo
c.........................................
c..   END DATA LOADING / REGION IDENTIFY
c.........................................






c......................................
c..   IONIZATION AND COULOMB EFFECTS
c......................................
      do k = kin, ktop
         t1 = tem(k)
         v1 = v(nc,k)
c..   composition variables (netsize = nnuc .le. ndim-1)
         yht = x(netsize-1,k)/xa(netsize-1)
         yat = x(netsize,k)/xa(netsize)
         xz  = 1.0d0 - x(netsize-1,k) - x(netsize,k)
c..   sum number of nuclei and associated electrons/nucleus
c..   protons and alphas excepted
         ynuc  = 0.0d0
         yenuc = 0.0d0
         do n = 1, netsize-2
            yenuc = yenuc + x(n,k)/xa(n) * dble( lz(n) )
            ynuc  = ynuc  + x(n,k)/xa(n)
         enddo
c..   this defines yenuc for ioniz.f
c..   yenuc is the Ye of the nuclei, ignoring H1 and He4
c..   ynuc is the mole fraction of those nuclei
c..   zbar is xz*yenuc/ynuc so
c..   abar is xz/ynuc
         if( xz .gt. 0.0d0 )then
            yenuc = yenuc/xz
         else
            yenuc = 0.0d0
         endif
c..   use debye theory for ionization corrections
c..   add coulomb effects to ionic eos
c..   construct average Z**2 + Z for debye approximation
c..   overestimates zeta for partial ionization
         zeta   = 0.0d0
         do n = 1, netsize
            zeta = zeta + x(n,k)/xa(n)*dble( lz(n)*(1+lz(n)) )
         enddo
c..   save for weak screening in rate.f
c     wszeta(k) = zeta

c..   debye (weak screening) approximation to coulomb effects
c..   debye sphere model (unbound nondegenerate electrons)
c..   aodebye is bohr radius/debye length
         aodebye = sqrt( 3.5402d5 * rhom(k) / tem(k) * zeta )
         weak    = ioh * aodebye
         xpf     = ( denye(k) / 0.9735d6 )**0.3333333333d0
c..   strong screening approximation to coulomb effects
         xbhc   = 0.0d0
         do n = 1, netsize
c..   fivethrds is Z**(5/3) for strong screening eos
            xbhc = xbhc + x(n,k)/xa(n)*fivethrds(n)
         enddo
c..   dimensionless fermi momentum xpf
         xpf    = ( denye(k) / 0.9735d6 )**0.3333333333d0
c..   102.01 * xpf is (9/5)a0*me*c/h(4pi/3)**(2/3)*xpf = kappas*a0
         strong = 13.6d0 * 102.01d0 * xpf
c..   force weak value if fw*fs/(fw+fs)

c     strong = ioh * ( 1.0d0 + 102.01d0 * xpf ) 
c     ecouls = -strong * xbhc * avagadro * ergspermev *1.0d-6

c..   define the chemical potential correction due to screening
c..   singly ionized chemical potential correction in eV
c..   note factor of Z**2 + Z = 2 for H+ He+, 6 for He++, 0 for neutrals
         uion1 = - weak*2.0d0
c..   doubly ionized
         uion2 = - weak*6.0d0

         if( nocoul .eq. 1 )then
c..   no coulomb correction in ionization
            aodebye = 0
         else
c..   avoid debye expansion nonconvergence (a0/R(D) << 1)

         endif
         aodeb(k) = aodebye
         pfrm(k)  = xpf
c..   trouble if cdel .gt. 1d-4 for convergence in convection
         cdel  = 1.0d-4
         facti = 1.0d0 + cdel
         t2    = t1*facti
         v2    = v1*facti
         if( noion .ne. 0 )then
c..   complete ionization
            yei   = ye(k)
            yei1  = yei
            yei2  = yei
            ya1   = 0
            ya2   = 1
            yh    = 1
            eion  = 0
            eion1 = eion
            eion2 = eion
         else

c..   default ionization-eos calculation
            do n = 1,3
               if( n .eq. 1 )then
                  d1 = 1.0d0/v1
                  call ioniz(t2,d1,yht,yat,xz,ynuc,yenuc,yh,ya1,ya2,
     .                 yei1,eion1,yze,aodebye,uion1,uion2,sion1,
     .                 pion1,yh2,weak,strong,zeta,k,n,nii)

               elseif( n .eq. 2 )then
                  d2 = 1.0d0/v2
                  call ioniz(t1,d2,yht,yat,xz,ynuc,yenuc,yh,ya1,ya2,
     .                 yei2,eion2,yze,aodebye,uion1,uion2,sion2,
     .                 pion2,yh2,weak,strong,zeta,k,n,nii)

               elseif( n .eq. 3 )then
                  d1 = 1.0d0/v1
                  call ioniz(t1,d1,yht,yat,xz,ynuc,yenuc,yh,ya1,ya2,
     .                 yei,eion,yze,aodebye,uion1,uion2,sion,
     .                 pion,yh2,weak,strong,zeta,k,n,nii)

                  if( pion .lt. 0.0d0 )then
                     write(*,*)'STATE: Pion < 0'
                     write(*,'(a5,12a12)')'k','T','rho','Pion',
     .                    'Pnuc+rad','yei','yh','ya1','ya2','yze'
                     write(*,'(i5,1p12e12.3)')k,
     .                    tem(k),rhom(k),pion,pnuc(k)+prad(k),
     .                    yei,yh,ya1,ya2,yze
                     stop'ERR(state): pion < 0.d0'
                  endif
               else
                  write(*,*)'ERR(state): n= ',n
                  stop'ERR(state): err with n'
               endif
            enddo
         endif

         niit(k)   = nii
         ddt       = t2-t1
         ddv       = v2-v1
         eiont     = (eion1 - eion)/ddt
         eionv     = (eion2 - eion)/ddv
         yeit      = (yei1  - yei )/ddt
         yeiv      = (yei2  - yei )/ddv
         yef(k)    = yei
         yeft(k)   = yeit
         yefv(k)   = yeiv
         eions(k)  = eion
         eionst(k) = eiont
         eionsv(k) = eionv

         edis(k)   = eion
         edist(k)  = (eion1 - eion)/ddt
         edisv(k)  = (eion2 - eion)/ddv
         pdis(k)   = pion
         pdist(k)  = (pion1 - pion)/ddt
         pdisv(k)  = (pion2 - pion)/ddv
         sdis(k)   = sion
         
         chemfak   = dlog( rhom(k)/tem(k)**1.5d0 )
c..   fraction of hydrogen which is ionized
         yhi(k)    = yh
c..   fraction of hydrogen which is molecular (0<yh2<0.5)
         yhm(k)    = yh2
c..   fraction of helium which is singly ionized
         yhe1(k)   = ya1
c..   fraction of helium which is doubly ionized
         yhe2(k)   = ya2
c..   electron entropy, corrected for partial ionization
c     if( yef(k) .gt. 1.0d-16 )then
         etanuc(netsize+1,k) = dlog(yef(k)) +chemfak 
     .        -chema(netsize+1) +chemcon
         sel(k)    = yef(k)*(2.5d0 - etanuc(netsize+1,k))
c..   adjust and save for electron eos below
         denye(k) = rhom(k)*yef(k)
      enddo
c.........................................
c..   END IONIZATION AND COULOMB EFFECTS
c.........................................




c.............................................................
c      print*,'DBG(state): OPALEOS tables'
c.............................................................
c..   REGION 0 and 1                       OPALEOS table   
c..   2 is pure HELMHOLTZ
c.............................................................
c..   
      do k = kin, ktop
         if( intab(k) .ne. 2 )then
c..   0 and 1
            d1  = rhom(k)
            t1  = tem(k)
            if( t1 .lt. 2.0d3 )then
c<<<<<<< .mine
               write(*,'(2(a12,i5),2(a12,1pe12.3))')'k',k,
     .              'ktop',ktop,'T',t1,'rho',d1
               write(*,'(1p8e12.4)')dmh(k),dmh(k+1),xm(k),r(nc,k),
     .              r(nc,k+1)
               write(*,'(2(a12,i5))')'L',l,'model',model
               stop'state T: opaleos'
c=======
               write(*,*)'STATE: LOW T ERROR'
               write(*,'(2(a12,i5),2(a10,1pe12.4))')'k',k,
     .              'ktop',ktop,'T',t1,'rho',d1
               write(*,'(5(a12,1pe12.3))')'dmh(k)',dmh(k),
     1          'dmh(k+1)',dmh(k+1),'xm(k)',xm(k),'r(nc,k)',r(nc,k),
     2          'r(nc,k+1)', r(nc,k+1),'tl(nc,k)',
     3          tl(nc,k),'tl(nc,k+1)',tl(nc,k+1)
               write(*,'(2(a12,i5))')'step',l,'model',model

               write(*,'(1p8e12.3)')sigma*t(nc,k)**4*pi4*r(nc,k)**2,
     1  (tl(nc,k)/(sigma*pi4*r(nc,k)**2))**0.25d0

               stop'STATE: T is below opaleos table'
c>>>>>>> .r98
            endif
            if( d1 .lt. 1.0d-15 )then
c<<<<<<< .mine
               write(*,'(a12,i5,2(a12,1pe12.3))')'k',k,'T',t1,
     .              'rho',d1
               stop'state d: opaleos'
c=======
               write(*,*)'STATE: LOW DENSITY ERROR'
               write(*,'(a12,i5,2(a12,1pe12.4))')'k',k,'T',t1,
     .              'rho',d1
               stop'STATE:  density is below opaleos table'
c>>>>>>> .r98
            endif

c..   hydrogen mass fraction
            xht = x(netsize-1,k)
            
c..   metallicity
c..   chooses EOS table using zpop0 which is passed as xz
c..   xz = zpop0
c..   allows change in metalicity with evolution (diffusion)
            xz = zmetal


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   forces average metallicity in EOS to equal that in opal
c..   EOS tables use single metallicity!!!!!
c..   it does interpolate in H,He so most of settling is accounted for
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..ieos2 is flag for opal eos to be read (if ieos2 .ne. 0)
c..set to 1 in initialization above, set to 0 in eosxtrin.f
            
            call eosxtrin(t1,d1,xht,xz,ppp(k),pppt(k),pppv(k),
     .           eee(k),eeet(k),eeev(k),sss(k),norad,k,ieos2)
            
            
            if( intab(k) .eq. 0 )then
c..   use OPAL estimate, 
c..   else wda estimate from above not overwritten
c..   mole number of free electrons = PV/RT - Yions            
               yef(k) =  ppp(k)/(rgas*t1*d1)
     .              -(xht/xa(nnuc-1)+(1.0d0-xht-xz)/xa(nnuc))
c..   keep yef(k) sane in electron degenerate conditions
               yef(k) = dmax1(yef(k),0.0d0)
               yef(k) = dmin1(yef(k),ye(k))
            endif
         endif
      enddo

c..   REGION 0
      do k = kin, ktop
c..   no radiation yet; added below
c..   put p,e etc in arrays; adjust arrays in interpolation region 1
         if( intab(k) .eq. 0 )then
            p(nc,k)    = ppp(k)
            pt(k)      = pppt(k)
            pv(k)      = pppv(k)
            e(nc,k)    = eee(k)
            et(k)      = eeet(k)
            ev(k)      = eeev(k)
            entropy(k) = sss(k)
         endif
      enddo




c----------------------------------------------------------------------
c..   REGION 1 and 2             HELMHOLTZ table
c..   0 is pure OPALEOS
c----------------------------------------------------------------------
      do k = kin, ktop
         if( intab(k) .ne. 0)then
            
            do j = 1, netsize
               xms(j) = x(j,k)
            enddo
            call azbar(xms,la,lz2,netsize,yms,abar,zbar)
c..   only ionized electrons
            zbar = zbar * yef(k)/ye(k)
            call eosfxt(tem(k),rhom(k),k,nc,csound)
c..   P,E, etc are in pe(k),ee(k),etc arrays
         endif
      enddo


c..   REGION 1
c..   interpolate HELMHOLTZ AND OPALEOS
      do k = kin, ktop
         if( intab(k) .eq. 1 )then
c..   topal and thelm are set at beginning of this subroutine
c..   interpolaton factors in temperature
            fact2 = (tem(k)-topal)
     .           /(thelm-topal)
c            fact2 = (dlog(tem(k))-dlog(topal))/(dlog(thelm)-dlog(topal))
            fact1 = 1.0d0 - fact2
            
            p(nc,k)    = (pnuc(k) + pe(k)  + pdis(k) )*fact2
     .           + ppp(k)*fact1
            pt(k)      = (pnuct(k)+ pet(k) + pdist(k))*fact2
     .           + pppt(k)*fact1
            pv(k)      = (pnucv(k)+ pev(k) + pdisv(k))*fact2
     .           + pppv(k)*fact1
            e(nc,k)    = (enuc(k) + ee(k)  + edis(k) )*fact2
     .           + eee(k)*fact1
            et(k)      = (enuct(k)+ eet(k) + edist(k))*fact2
     .           + eeet(k)*fact1
            ev(k)      = (enucv(k)+ eev(k) + edisv(k))*fact2
     .           + eeev(k)*fact1
            entropy(k) = (snuc(k) + sel(k) + sdis(k) )*fact2
     .           + sss(k)*fact1
         endif
      enddo




c..............................................................
c..   REGION 2 pure HELMHOLTZ
c..   collect components for wda-helm EOS, excepting radiation
c..............................................................
      do k = kin, ktop
         if( intab(k) .eq. 2 )then
            p(nc,k)    = pnuc(k) + pe(k)  + pdis(k)
            pt(k)      = pnuct(k)+ pet(k) + pdist(k)
            pv(k)      = pnucv(k)+ pev(k) + pdisv(k)
            e(nc,k)    = enuc(k) + ee(k)  + edis(k)
            et(k)      = enuct(k)+ eet(k) + edist(k)
            ev(k)      = enucv(k)+ eev(k) + edisv(k)
            entropy(k) = snuc(k) + sel(k) + sdis(k)
         endif
      enddo

c.........................................
c..   ALL REGIONS
c..   add radiation pressure and energy
c.........................................
      do k = kin, ktop
         p(nc,k)    = p(nc,k)    + prad(k)  
         pt(k)      = pt(k)      + pradt(k) 
         pv(k)      = pv(k)      + pradv(k) 
         e(nc,k)    = e(nc,k)    + erad(k)  
         et(k)      = et(k)      + eradt(k) 
         ev(k)      = ev(k)      + eradv(k) 
         entropy(k) = entropy(k) + srad(k)  
      enddo




c..............................................................
c..   CALCULATE SOUND SPEED, sqrt(cc), and SANITY CHECK: c*c>0
c..............................................................
      do k = kin, ktop
         dtdv = - (p(nc,k) + ev(k))/et(k)
         gam  = - v(nc,k)/p(nc,k)*( pv(k) + pt(k)*dtdv )
         gamma1(k) = gam
     
    
c..   special relativistically correct
         cc   = gam * crad2 /
     .        ( 1.0d0 + (e(nc,k) + crad2)/(p(nc,k) * v(nc,k)) )
c..   
         if( cc .gt. 0.0d0 )then
            sound(k) =  dsqrt( cc )
         else
c..   trouble in EOS
            write(*,*)'STATE: c*c<0 error in eos at k =',k
            write(*,'(5(a11,1pe11.3))')'T',t(nc,k),'V',v(nc,k),
     .           'Tem',tem(k),'rhoz',rhoz(k)

            do n = 1, nnuc+1
               write(*,'(5(a5,1pe11.3))')cnuc(n),x(n,k)
            enddo

            write(*,'(3(a11,i11),2(a11,1pe11.3))')"intab(k)",intab(k),
     .           "ikt(k)", ikt(k), "ikd(k)", ikd(k)

            write(*,*)'gam is  - v/p*( pv + pt*dtdv )'
            write(*,'(4(a11,1pe11.3))')"gam",gam,"v/p*pv",
     .           v(nc,k)/p(nc,k)*pv(k),
     .           "pv+pt*dtdv",pv(k)+pt(k)*dtdv,"dtdv",dtdv

            write(*,'(6(a11,1pe11.3))')"p",p(nc,k),"pt",pt(k),
     .           "pv",pv(k),'entropy',entropy(k)
            write(*,'(6(a11,1pe11.3))')"e",e(nc,k),"et",et(k),
     .           "ev",ev(k),"cfw",cfw(k),"cfs",cfs(k),"cf",cf(k),
     .           "pcs",cfs(k)/(3.0d0*v(1,k)),"rgas*t*v",
     .           rgas*t(1,k)*v(1,k),'yef/ye',yef(k)/ye(k),
     .           "pfrm",pfrm(k),"pdeg",6.0d22*1.6d0*pfrm(k)**5,
     .           "ef/kt",pfrm(k)**2/2.0d0*(11.605d9*0.511d0/t(1,k))
            write(*,'(6(a11,1pe11.3))')"pnuc",pnuc(k),"prad",prad(k),
     .           "pe",pe(k),"pcoul",pcoul(k),"pdis",pdis(k)
            write(*,'(6(a11,1pe11.3))')"enuc",enuc(k),"erad",erad(k),
     .           "ee",ee(k),"ecoul",ecoul(k),"eions",eions(k)
            write(*,'(6(a11,1pe11.3))')"snuc",snuc(k),"srad",srad(k),
     .           "sel",sel(k),"scou",scou(k),"sdis",sdis(k)
            write(*,'(6(a11,1pe11.3))')"pnuct",pnuct(k),
     .           "pradt",pradt(k),"pet",pet(k),
     .           "pe*yet/ye", pe(k)*yeft(k)/yef(k),
     .           "pct",pct(k),"pdist",pdist(k)
            write(*,'(6(a11,1pe11.3))')"pnucv",pnucv(k),
     .           "pradv",pradv(k),"pev",pev(k),
     .           "pe*yev/ye", pe(k)*yefv(k)/yef(k),
     .           "pcv",pcv(k),"pdisv",pdisv(k)
            write(*,'(6(a11,1pe11.3))')"enuct",enuct(k),
     .           "eradt",eradt(k),"eet",eet(k),
     .           "edist",edist(k)
            write(*,'(6(a11,1pe11.3))')"enucv",enucv(k),
     .           "eradv",eradv(k),"eev",eev(k),
     .           "edisv",edisv(k)
            write(*,'(6(a11,1pe11.3))')"pdis",pdis(k),"pdist",pdist(k),
     .           "pdisv",pdisv(k),"Sdis",sdis(k)
            write(*,'(6(a11,1pe11.3))')"edis",edis(k),"edist",edist(k),
     .           "edisv",edisv(k)
            write(*,'(6(a11,1pe11.3))')"ye",ye(k),"yef",yef(k),
     .           "yhf",yhf(k),"(1-yhf)/2",0.5d0*(1.0d0-yhf(k)),
     .           "yeft/yef",yeft(k)/yef(k),"yefv/yef",yefv(k)/yef(k)
            write(*,'(13a11)')'pt','pv','et','ev',
     .           'tpt/p','vpv/p','tet/e','vev/e'
            write(*,'(1p13e11.3)')pt(k),pv(k),et(k),ev(k),
     .           t(nc,k)*pt(k)/p(nc,k), v(nc,k)*pv(k)/p(nc,k),
     .           t(nc,k)*et(k)/e(nc,k), v(nc,k)*ev(k)/e(nc,k)
            write(*,'(a5,8a12)')'k','T','rhom','cfw','cfs','cf',
     .           'yef','ye','R*S-E/T'
            write(*,'(i5,1p8e12.3)')k,tem(k),rhom(k),cfw(k),cfs(k),
     .           cf(k),yef(k),ye(k),rgas*entropy(k)-e(nc,k)/tem(k)
            write(*,'(a10,1p8e12.3)')'Ropal',rhom(k)/(tem(k)*1.0d-6)**3
            write(*,'(2(a10,1pe12.3))')'sound',sound(k),'sound**2',cc
            stop'state: sound speed negative'
         endif
      enddo



c..   timing state, sum of elapsed time in runstate
      call seconds(runstate1)
      runstate = runstate + runstate1-runstate0



      
c     SUCCESS
      return
c     
      end
      
