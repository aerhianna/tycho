c     ----------------------------------------------------
c     TYCHO - A Stellar Evolution Program
c     ----------------------------------------------------
c     Primary website:
c     http://chandra.as.arizona.edu/~dave/tycho-intro.html
c     ----------------------------------------------------
      
      program tycho
      
      
c     MAJOR REVISIONS
c     revised 12-11-06
c     1-d hydrodynamic stellar evolutionary code
c     general reaction network
c     helmholtz equation of state
c     OPAL  and Anderson-Ferguson opacities
c     includes hydrogen molecule and pressure ionization/dissociation
c     and envelope-interior join correction
c     Thomas-Kippenhahn rotation (centrifugal force)
c     uses mass fractions for x(n,k)
c     Richardson convection criterion
c     Turbulent entrainment
      
      
      
c..   TYCHO VERSION LOG
c..   6.9  on website for 535 spring 2006
c..   6.91 dthoul.f tests for no Hydrogen
c..        gtintp.f has mass and pressure error tests, radius stepsize
c..        cgtintp  has larger dimension for try(i,k) zoning tests
c..        solven.f only write negflag for repeated problem
c..        lintrx.f uses linear interpolation
c..   7.0  uses both binary and formatted models, 
c..   complete and backward compatible
c..   

      
c..   ADDITIONAL MAINTENANCE LOG
c..   15 September 2008: Code cleaned up and more
c..   orderly comment format started (C.A.Meakin)
c..   

      implicit none

      include 'dimenfile'
      include 'comod'
      include 'compu'
      include 'ctmstp'
      include 'cqdv.h'
      include 'cruntm.h'
      include 'cbug'
      include 'czone'
      include 'cnabla'
      include 'cphot'
      include 'cgtr.h'
      include 'cconst'
      include 'conline'
      include 'caeps'
      include 'csurface'
      include 'cenv'
      include 'cburn'
      include 'cgen'
      include 'cpgplot'
      include 'cenrchk'
      include 'crate'
      include 'ceoset'

      integer*4 idebug, nupdat, no, ncycle, ktype, kkman, ks
      integer*4 mix, kxmx, nc,  k, j, i, n


      real*8 gamma,del
      real*8 epssum, epsg, fak, runc, runs, stepsec, dum
      real*8 tlum0, tlum1, tlum2, tlfak, tluma
      real*8 radius0,radius1,radius2
      real*8 tim2,period,dthc, area,flux,aog,aogmax
      
c..   hermite variables
      real*8 f0her,f1her,d0her,d1her,dmher,fhher,dther
      
      real*8 sum
      integer*4 kaogmax,naogmax
      logical tobe
      data ctstamp/"TYCHO 8.00"/


c---------------------------------------------------------------------

      write(*,*)ctstamp

c---------------------------------------------------------------------
c     
c     --------------
c     time centering
c     --------------
c     
c     for the variable u
c     dimensioning subscript 2 permits a time subscript to take on the
c     values 1 or 2..1 implies time(n-1/2), 2 implies time(n+1/2)
c     for all other time dependent variables
c     values 1 or 2..1 implies time(n), 2 implies time(n+1)
c     only the boundary variable u(2,k) is taken at n-1/2,n+1/2,n+3/2
c     other time dependent variables(boundary or zone) are at n-1,n,n+1 
c     
c     dth(1)=time(n)-time(n-1)
c     dth(2)=time(n+1)-time(n)
c     dti(1)=time(n+1/2)-time(n-1/2)
c     dti(2)=time(n+3/2)-time(n+1/2)
c     time(n-1)....
c     
c     values 1 or 2.. 1 implies time(n-1/2), 2 implies time(n+1/2) for u 
c     values 1 or 2.. 1 implies time(n), 2 implies time(n+1) for all els
c     
c     .............time(n+3/2)...velocity
c     .
c     dti(2).            time(n+1)...thermo variables..
c     .
c     ^        ......u......time(n+1/2)..velocity..........dth(2)
c     |        .
c     |   dti(1).          time(n).....thermo vars.(r,p,t,v,d,du,s,)
c     |         .
c     |         ....u......time(n-1/2)..velocity............dth(1)
c     
c     
c     mode = 0, 1, 2, 3
c     u(n,k) is at t(n), not t(n-1/2)
c     
c     
c     ---------------
c     space centering
c     ---------------
c     kk is the no. of boundaries
c     kk-1 is the no. of interior zones
c     there are two exterior zones
c     hence, kk+1 is the total no. of zones
c     
c     boundary,zone orientation:
c     zone 1    bound 1   zone 2    bound 2--------bound kk  zone
c     1/2        1        3/2        2             kk      kk+1/2
c     .          1         .         2-------------kk        .
c     1          1         2         2             kk      kk+1
c     1                    2                 kk
c     1                    2                 kk
c     
c     
c     [NEW DEPICTION OF ZONE INDEXING: Meakin 2008]
c     
c     .      |-------|-------|  ...  |-------|======|======|
c     . k=   1,      2,      3, ..., kk-1,  kk,     kk+1,   kk+2   (boundary)
c     . k=(1)    2,      3,     ...,    kk,    kk+1,   kk+2        (interior) 
c     
c     [NOTE: Interior zone k=1 seems to be a placeholder ]
c     k is index which keeps track of zone and boundary computation
c     a boundary variable at physical boundary k has subscipt k
c     a zone variable at physical zone (k-1/2) has subscript k
c     
c     * see comments in dimension area describing which variables are
c     computed on boundaries, and which are computed at zone points
c     
c-------------------------------------------------------------------


c..   initialize arrays to zero
      call zrow
      
      idebug = 0
      nbug   = 0
      it     = 0
      nupdat = 0
      l      = 0
      no     = 0

      work   = 0.0d0
      elum   = 0.0d0
      esou   = 0.0d0
      enucl  = 0.0d0

      runburn0  = 0
      runburn1  = 0
      runburn   = 0
      runstate0 = 0
      runstate1 = 0
      runstate  = 0
      runenvel0 = 0
      runenvel1 = 0
      runenvel  = 0

      tim2   = 0.0d0
      temin  = 0.0d0
      period = 0.0d0
      ncycle = 0
      ktype  = 13+2
      kkman  = 1
      ks     = 2
      mix    = 0
      kxmx   = 0

c..   negative value is a flag for lack of redefinition later
      zpop = -1.0
      
      
c     -----------------
c..   FILE UNITS for IO
c     -----------------
c     5 stdin (don't use))
c     6 stdout (don't use!)
c
c     LOCAL ONLY:
c..   2  params.d
c..   2  timmes helm_table.dat
c..   2  tycho params.d
c..   2  readco
c..   2  readop
c..   10 abinit (ssystem.dat)
c..   10 gen scratch
c..   10 pgplot scratch.pg
c..   12 masses.f
c..   12 partfun netwinv
c..   12 tycho soltest.file
c..   12 getfkt isotope.lib
c..   14 deo eos.d, eos1.d, ceosfile
c..   14 state eos1
c..   19 weakread.f netweak
c..   20 state codataa
c..   30 abinit (net.rc)
c..   30 kapp ku2data
c..   31 kappinit alxtab
c..   30 new1kapp,new2kapp alxdata
c..   33 nse  data.a
c..   33 sqse data.a
c..   55  getfkt netsu

c..   GLOBAL (do not open or close these):
c..   3  tycho out.file
c..   8  gen   filout
c..   8  dout5 filout
c..   11 gen   filecv
c..   11 cvplt filecv
c..   77  gen   filehr
      


c     ---------------------------------------
c     INPUT PARAMETERS TO CONTROL COMPUTATION
c     ---------------------------------------
c     
      inquire(file='params.d',exist=tobe)
      if( tobe )then
         open(2,file='params.d',form='formatted',status='old')
      else
         write(*,*)'tycho: no file params.d in this directory'
         stop'tycho'
      endif
      


c..   output unit for output to logfile for current run 
c..   not permanent, will be overwritten by next run
      open(3,file='logfile',form='formatted',status='unknown')



c     -------------------------------------------
c..   BUILD NUCLEAR REATION NETWORK, 
c     reads: isotope.lib, netsu, netwinv, netweak
c     -------------------------------------------
      call build(idebug)

      

c     -------------------------------
c     READ IN INITIAL MODEL
c     READS: params.d, imodel, 
c     SUB abinit IS CALLED AND READS:
c     net.rc, ssystem.dat
c     -------------------------------
      call gen

      
c     -------------------------------------
c..   construct low T opacity table lookups
c     reads alxtab
c     -------------------------------------
      call kappinit
      
      xlol  = 0
      telog = 0


      
c..   set up secondary arrays for n=1 level in time
      nc   = 1
      

c..   equation of state, reads eos1.d, helm_table.dat
c..   and for OPAL type 2 (nopac=0): 
c..   codataa, codatab, codatac, codatad, codatae
c..   for OPAL type 1 (nopac=1):
c..   hzdata
c..   for OPALEOS (in eosxtrin.f) EOSdata.148 or EOSdata.188
      


      
      
c..   THIS NEEDS TO BE UPDATED TO INTERPOLATION IN z (metallicity)
      
      call state(2,kk,nc)
      
c..   save abundances from start of time step; this is required to
c..   restart aborted attempt at convergence
      do k = 1,kk+1
         do j = 1, ndim
            xold(j,k) = x(j,k)
         enddo
      enddo


c..   einitialize burn and update energy generation only 
c..   (nc,1,...) for value 1
      call burn(nc,1,idebug,2,kk)
      

c..   compute general relativistic correction factors
      call gtrs(p,q,e,v,r,u,xm,dmi,dmh,kk,nc,newt,mode)
      
c..   initialize convection variables
c      call cinit(2,kk,nc)
c	   do this in hydro.f


c..   initialize energy checks
      call enrchk(no)


c..   insure that timestep conditions arre satisfied on first step
c..   for mode = 0 uses lesser of courant and dth1 from params.d
      call timstp(dthc,period,ktype,kkman,1,kxmx)
      
c..   always reset past timestep for initial step
c..   to guard against large correction in timstp.f
      dth(1) = dth(2)
      
      boflag = 0  !for AMANDA BlackmanOwen Routine
      zamstime = 0.00
c     ---------------------
c..   end of initialization
c     ---------------------
      




c-----------------------------------------------|
c     BEGIN TIME EVOLUTION LOOP                 |
c     (logic returns here after each time step) |
c-----------------------------------------------|
 505  continue
      model = model + 1

c..   l is index for time steps attempted
      l = l+1
      

c..   save abundances from start of time step; this is required to
c..   restart aborted attempt at convergence
      do k = 1,kk+1
         do j = 1, ndim
            xold(j,k) = x(j,k)
         enddo
      enddo


c..   calculate next time step for momentum equation
      dti(1) = 0.5d0*(dth(1) + dth(2))


c..   inner boundary condition
      call bound

c     ------------------
c..   AMANDA 10/03/16-5/26/17
c..   Blackman & Owen magnetic field mass loss
c     ------------------
c      call blackowen(1)


c..   outer boundary condition: mass accretion and loss
      call advec(nc,kk)


c..   update to n=2; with updated abundances and convective velocities
      call hydro(ks,no,mix,kxmx,dthc)

	  call blackowen(1)

      if( no .ne. 0 )then
c..   count nonconvergences
         nupdat = nupdat + 1
      endif


c..   number of network cycles (steps) computed
      ncycle = ncycle + it


      if( no .eq. 0 .and. mset .eq. 0 )then
c..   update time and angular velocity on successful convergence
         time = time + dth(2)
         tim2 = tim2 + dth(2)
c..   angular velocity
         do k = 2, kk
            omeg(k) = ajay(k)/r(2,k)**2
         enddo
         if( r(2,kk+1) .gt. 0.0d0  )then
            omeg(kk+1) = ajay(kk+1)/r(2,kk+1)**2
         endif
         omeg(1) = omeg(2)
      endif




c     ------------------
c..   SOME IO & GRAPHICS
c     ------------------
c..   0) analysis of new time step for on-line output
      call online(tim2,gamma,del,no,kkman,ktype,ks)
      
c..   1) pgplot monitors evolution
      call pgplot(no,ktype,kkman)
      
c..   2) hrplt gives location in HR diagram
      call hrplt(no,ks)
      
c..   3) cvplt give plot of evolution of convection zones
      call cvplt(no)

      
      
c..   estimate best new time step
      call timstp(dthc,period,ktype,kkman,2,kxmx)




c     ----------------------------
c     TEST FOR STOPPING CONDITIONS
c     ----------------------------
c     1) Time step size
      if( dth(2) .lt. 1.0d-20 )then
         ll = l
         write(*,*)'TYCHO: time step too small: ',dth(2)
      endif
      
c     2) stop when core begins to collapse
      if( 1.0d0/v(1,2) .gt. 2.0d10 )then
         ll = l
         write(*,*)'TYCHO: central density limit exceeded'
      endif
      
      if( mode .eq. 1 )then
         
c     3) stop when gamma - 4/3 negative (and at high density)
         if( del .lt. 0.0d0 .and. 1.0d0/v(1,2) .gt. 1.0d8)then
            ll = l
            write(*,*)'mode 1: gamma less than 4/3 limit exceeded'
         endif
         
c     4) stop when surface temperature drops
         if( t(1,kk) .lt. 1.0d3 )then
            ll = l
            write(*,*)'mode 1: surface tempterature too low'
         endif
      endif
      
c     5) test for supereddington luminosity on grid
      if( mode .ne. 0 )then
         aogmax  = 0.0d0
         kaogmax = 1
         naogmax = 0
c     join boundary radiative flux must be defined differently
         do k = 2, kk-1
            area = pi4*r(2,k)**2
c     kappa omitted in flux definition because it cancels in aog [a/g]
            flux = arad*crad/3.0d0 * (t(2,k)**4-t(2,k+1)**4)
     .           * area/dmi(k)
            aog = area*flux/xm(k)/(pi4*grav*crad) 
c     count zones in excess of Eddington luminosity
            if( aog .gt. 1.0d0 .and. ic(k) .eq. 0 )then
c     radiative zones only (they should become convective)
               naogmax = naogmax + 1
            endif
c     find worst zone
            if( aog .gt. aogmax )then
               aogmax  = aog
               kaogmax = k
            endif
         enddo
         if( aogmax .gt. 0.9d0 )then
            write(*,'(a25,4(a10,i5),a10,1pe12.3)')
     .           'Tycho: SuperEdd > 0.9:','kk',kk,'worst k',kaogmax,
     .           'ic(max)',ic(kaogmax),'num zns',naogmax,
     .           'a/g max',aogmax
            write(*,'(6(a10,1pe12.3))')'nab',dnab(kaogmax),
     .           'nad',dnad(kaogmax),'doux',doux(kaogmax),
     .           'n-a-x',dnab(kaogmax)-dnad(kaogmax)-doux(kaogmax),
     .           'nrad',dnrad(kaogmax),
     .           'alt ar/g',2.0d0*arad/3.0d0*(
     .           t(2,kaogmax)**4/p(2,kaogmax)+
     .           t(2,kaogmax+1)**4/p(2,kaogmax+1) )*dnab(kaogmax)
         endif
         if( naogmax .ge. kk/10 )then
c     about 10 percent of zones
c            ll = l
            write(*,*)'Eddington limit exceeded at step ',ll,
     .           ' model ',model,' in naogmax zones'
         endif
      endif

      
c..   6) only nup update defaults allowed
c..   nup = 1 continues until a nonconvergence (nupdat = 1)
      if( nupdat .ge. nup )then
         ll = l
      endif

      
c..   7) check machine time used
      call seconds(runt)
      if( runtmx .gt. 0.0d0 ) then
         if( runt .ge. runtmx - 30.0d0 )  ll = l
      endif
      
      
c..   8) check existence of lockfile= "starlock"
      call haltwda(l,ll)
      

c..   9) stop at star time
      if( stime .ne. 0.0d0 )then
         if( stime-time .le. dth(2) .and. 
     .        stime-time .gt. 0.0d0 )then
c..   next step will overshoot target unless dth is adjusted
c..   accurate to 0.0001
c..   avoid negative dth(2)
            dth(2) = 1.0001d0*(stime-time)
            dth(1) = dth(2)
         endif
      endif
      if( stime .ne. 0.0d0 .and. time .ge. stime )then
         ll = l
         
         write(*,'(a25,a6,1pe12.4,a6,a10,1pe12.4,a6)')
     .        '--- startime exceeded -->','age',time/secpy,' years ',
     .        'target',stime/secpy,' years '
c     write(*,'(a25,i5,a12,0pf10.3)')'mixmode',mixmode,
c     1        'alpha(ml)',alphaml
         write(*,*)'Basu & Antia(1997) Rcz/Rsolar = (7.13+-0.001)e-01'
         write(3,'(a25,a6,1pe12.4,a6,a10,1pe12.4)')
     .        '--- startime exceeded ---','age',time/secpy,' years ',
     .        'target',stime/secpy
         
c..   solar test
         if( nsoltest .gt. 0 )then
            open(12,file='sol.test.file')
            write(12,'(a25,a6,1pe12.4,a6,a10,1pe12.4)')
     .           '--- startime exceeded ---','age',time/secpy,' years ',
     .           'target',stime/secpy
            write(12,'(a15,i5,a15,1pe12.3)')'mixmode',mixmode,
     .           'alpha(ml)',alphaml
         endif
c..   find convective depth
         do i = 2, kk
            k = kk+1-i
            if( ic(k) .eq. 0 )then
               goto 1234
            endif
         enddo
 1234    continue
         if( k .lt. kk-1 )then
            write(*,'(4(a10,1pe12.3))')
     .           'R(sun) ', r(2,kk+1)/solrad,
     .           'R(conv)',r(2,k)/solrad,
     .           'L(sun)', tl(2,kk)/sollum,
     .           'M(sun)', xm(kk+1)/sol
            if( nsoltest .gt. 0 )then
               write(12,'(4(a10,1pe12.3))')
     .              'R(sun) ', r(2,kk+1)/solrad,
     .              'R(conv)',r(2,k)/solrad,
     .              'L(sun)', tl(2,kk)/sollum,
     .              'M(sun)', xm(kk+1)/sol
            endif
         else
c..   search for convective zone in envelope
            do i = 1, jmaxz
               k = jmaxz +1 - i
               if( znrad(k) -znad(k) .gt. 0.0d0 )goto 2345
            enddo
 2345       continue
            write(*,'(4(a10,1pe12.3))')
     .           'R(sun) ', r(2,kk+1)/solrad,
     .           'R(conv)',zr(k+1)/solrad,
     .           'L(sun)', tl(2,kk)/sollum,
     .           'M(sun)', xm(kk+1)/sol
            if( nsoltest .gt. 0 )then
               write(12,'(4(a10,1pe12.3))')
     .              'R(sun) ', r(2,kk+1)/solrad,
     .              'R(conv)',zr(k+1)/solrad,
     .              'L(sun)', tl(2,kk)/sollum,
     .              'M(sun)', xm(kk+1)/sol
            endif
         endif
         write(*,'(4(a10,1pe12.3))')
     .        'T(c)  ', t(2,2),
     .        'rho(c)', 1.0d0/v(2,2),
     .        'H1(c) ', x(nnuc-1,2),
     .        'He4(c)', x(nnuc,2),
     .        'C12(c) ', x(lc12,2),
     .        'N14(c)', x(ln14,2),
     .        'O16(c)', x(lo16,2),
     .        'P(c)  ', p(2,2)
         
         if( nsoltest .gt. 0 )then
            write(12,'(4(a10,1pe12.3))')
     .           'T(c)  ', t(2,2),
     .           'rho(c)', 1.0d0/v(2,2),
     .           'H1(c) ', x(nnuc-1,2),
     .           'He4(c)', x(nnuc,2),
     .           'C12(c) ', x(lc12,2),
     .           'N14(c)', x(ln14,2),
     .           'O16(c)', x(lo16,2),
     .           'P(c)  ', p(2,2)
            close(12)
         endif
      endif
      
      
c     10) Stop when fuel in core is exhausted
      if( ixstop .gt. 0 .and. ixstop .le. ndim )then
         if( x(ixstop,2) .lt. xstop )then
            ll = l
            write(*,*)'--- species ',cnuc(ixstop),' .lt. ',xstop,
     .           ' at core ---'
            write(3,*)'--- species ',cnuc(ixstop),' .lt. ',xstop,
     .           ' at core ---'
         endif
      endif
      
      
c     11) Stop when star swings unrealistically far to red
      if( mode .eq. 1 )then
         if( telog .le. 3.0d0 .and. telog .ne. 0.0d0 )then
            if( mloss .eq. 0 )then
c     ignore if in mass accretion mode, L may be negative
               ll = l
            endif
            write(*,*)'mode 1: log Teff < limit, = ',telog
         endif
      endif

      
c     12) Stop at ZAMS
      if( mode .eq. 1 .and. izams .eq. 1)then
c     nuclear luminosity
         epssum = 0
         do k = 2,kk
            epssum = epssum + s(5,k)*dmh(k)
         enddo
c     luminosity from gravity (work done by gravitational force in
c     hydrostatic balance)
         epsg = 0.0d0
         do k = 2, kk
            fak = e(2,k)-e(1,k) + 0.5d0*(p(2,k)+p(1,k))*(v(2,k)-v(1,k)) 
            epsg = epsg + dmh(k)*fak/dth(2)
         enddo
         write(*,'(a15,5x,3(a10,1pe12.4))')'ZAMS search:',
     1        'nuclear L',epssum,'surface L',tl(1,kk),'gravity L',epsg
         if( epssum .gt. 1.0d2*abs(epsg) )then
            write(*,*)'ZAMS found!'
            ll = l
         endif
      endif
      



c     -------
c     MORE IO
c     -------
      
c..   1) the variable l3 is control index for dout5. 
c     Subroutine dout5 dumps a model every l3 timesteps
      if( l3 .gt. 0 .and. l .gt. 0 )then
         if( mod(l,l3) .eq. 0  .or. l .ge. ll )then
            call dout5(2)
         endif
      endif
      
c     2) call dout3(...): currently a stub
      if( l .ge. ll .or. mod(l,l2) .eq. 0 )then
         nc = 2
         call dout3(nc,period)
      endif
      
c     3) Always enter edit, it writes summary data to out.file 
c..   every L2 time steps
      call edit(gamma,del,tim2,no)
      


c     ----------------------------------
c     TIMING: cpu runtime per iter cycle
c     ----------------------------------
      if( ncycle .gt. 0 )then
         if( l .eq. ll .or. mod(l,l2) .eq. 0 )then
            runc = runt / dble( ncycle )
            runs = runt / dble( l )
            if( runs .gt. 0.0d0 )then
               stepsec = 1.0d0/runs
            else
               stepsec = 0
            endif
            write(3,16) runc, runs,nupdat,stepsec,runt,runt/60.0d0
            write(6,16) runc, runs,nupdat,stepsec,runt,runt/60.0d0
 16         format( 1pe10.2,' sec/iter', 1pe10.2,
     1           ' sec/step', i6, ' defaults',
     2           1pe10.2, ' steps/s',1pe10.2,' sec total', 
     3           1pe10.2, ' min total')

            write(*,'(3(a15,1pe12.3))')'time(burn)',runburn,
     1           'time(state)',runstate,'time(envel)',runenvel
            write(*,*)'note: time(envel) and time(state) overlap'
         endif
      endif
      


      

c     ---------------------------------
c     EXIT FOR COMPLETION OF TIMESTEPS:
c     SHUTDOWN SEQUENCE
c     ---------------------------------
      if( l .ge. ll )then
         
         if( nupdat .ge. nup )then
            write(*,*)'TYCHO: nupdat= ',nupdat,' .ge. nup= ',nup
         else
            write(*,*)'TYCHO: normal shutdown sequence '
         endif

         
c     1) calculate abundance anomolies
         call anomaly

         
         write(*,'(a12,2(1pe12.3,a8))')'Stellar age',time,'seconds',
     1        time/secpy,'years'
         write(*,*)note
         if( mode .eq. 0 )then
            write(*,'(a20,i5)')'hydrodynamic mode',mode
         else
            write(*,'(a20,i5)')'quasistatic mode',mode
         endif
         if( nburn .ne. -1 )then
           if( jnb .ne. 0 )then
              write(*,'(a20,i5)')'Using JNB rates',jnb
           else
              write(*,'(a25)')'Using Rauscher/FKT rates'
           endif
         else
              write(*,'(a20)')'No burning'
         endif
         if( nopaleos .ne. 0 )then
            write(*,'(a20,i5)')'Using OPALEOS',nopaleos
         else
            write(*,'(a20,i5)')'Using HELMHOLTZ EOS',nopaleos
         endif
         if( nopac .eq. 0 )then
            write(*,'(a20,i5,a15)')'OPAL opacity type 2',nopac,copal
         elseif( nopac .eq. 1)then
            write(*,'(a20,i5,a15)')'OPAL opacity type 1',nopac,copal
         endif
         if( nkscale .ne. 0 )then
            write(*,'(a20,i5)')'Opacity scaling',nkscale
         else
             write(*,'(a20,i5)')'No Opacity scaling',nkscale
         endif
         if( fthoul .gt. 0.0d0 )then
            write(*,'(a20,1pe12.3)')'Diffusion scaling',fthoul
         endif

         write(*,'(a20,i5,a15,1pe12.3)')'mixmode',mixmode
         write(*,'(a13,1pe12.3,a15,1pe12.3)')
     1        'alpha(ml)',alphaml,'uuml',uuml

         write(*,'(2(a20,1pe12.3))')'zmetal',zmetal,'zpop(tab)',zpop
         
         write(*,'(a20,3(a5,1pe12.3))')'surface','H',x(nnuc-1,kk+1),
     1        'He',x(nnuc,kk+1),'z',1.0d0-x(nnuc,kk+1)-x(nnuc-1,kk+1)
         
         write(*,'(a20,i12)')'mvmax',nvmax(3)


         write(*,'(a20,1pe12.4)')'Mtot/Msun',xm(kk+1)/sol
         
         write(*,'(a20,1pe12.4)')'T(kk)',t(nc,kk)

         if( modes .eq. 2 )then
            write(*,'(a20,1pe12.4)')'Menv/M',dmh(kk+1)/xm(kk+1)
           write(*,'(a20,1pe12.4)')'Radius',r(nc,kk+1)
         else
           write(*,'(a20,1pe12.4)')'Radius',r(nc,kk)
         endif

c     2) Close file 10 (Q:What is this? A:various scratch files)
         close(10)
         
c     3) Remove pgplot scratch file
         call system('rm scratch.pg')
         

c     4) Close file 99: opacity scratch file
         close(99)
c     call system('rm newkapp.dummy')
         
         if( ipause .ge. 0 )then
c            pause
         endif
         
c     5) CLOSE PGPLOT GRAPHICS IF OPEN
         if( igraf .eq. 0 )then
c..   scratch file for character conversion
c     close(10)
c     call system('rm scratch.pg')
c..   keeps last model on graphics screen; asks for response
            call pgslct(idpg1)
            call pgclos
            call pgslct(idpg2)
            call pgclos
            call pgslct(idpg3)
            call pgclos
         endif

         
c     6) STOP PROGRAM! :-)
         stop'standard tycho exit'
      endif
      


      
c     --------------------
c     CHECK ERROR FLAG: no
c     --------------------
c     1) NON-CONVERGENCE
      
      if( no .ne. 0 )then
         
c     update for no success....
         write(*,10)l,model,no
         write(3,10)l,model,no
 10      format('=========== no update at l=',i6,', model=',i6,
     1        ', no=',i3)
c..   reset count of successful steps
         model = model - 1

c     reduce time step dth(2) by factor dum, update dth(1)
         dth(1) = dth(2)
         dum    = 0.5d0
c     idt read in by subroutine gen
         if( idt .ne. 0 ) dum = 0.9d0
         dth(2) = dth(2)*dum
c     zero derivatives after nonconvergence in hydro
c     leaves old (n=1) values unchanged
c     resets abundances xold-->x
c     does not update time, energy check, or call timstp or rezoning
         
         do k = 1, kk+1
            do n = 1, ndim
               x( n,k)  = xold(n,k)
               xd(n,k)  = 0.0d0
            enddo
         enddo
         
c..   zeros convective luminosity if time step is small
         if( dth(1) .lt. 1.0d0 )then
            do k = 1, kk+1
               dtl(k)  = 0
c     tl(2,k) = 0
               h(k)    = 0
               b(k)    = 0
            enddo
         endif
         
      else



c     2) CONVERGENCE
c..   reset new values to start next step: (2,k)-->(1,k), etc.
         
         do k = 1, kk+1
            h(k) = hp(k)
            if( h(k) .lt. 0.0d0 ) h(k) = 0.0d0
         enddo
         
c..   damp large mismatch at join (see fitenv.f too)
         if( modes .eq. 2 )then
            if( l .le. 1 )then
c..   tl prior to beginning of dth(1), tl(0,kk) initialized with no change
               tlum0   = tl(1,kk)
               radius0 = r(1,kk)
            else
c..   rotate to save previous luminosities
               tlum0 = tlum1
               radius0 = radius1
            endif
            tlum1   = tl(1,kk)
            tlum2   = tl(2,kk)
            radius1 = r(1,kk)
            radius2 = r(2,kk)
            
c     if( abs( tlum2-tlum1 ) .gt. 0.5d0*abs(tlum1) )then
c..   very large change so limit its size:        
c..   for large time steps relative to helmholtz-kelvin,
c..   merge toward latest value
c     tlfak = 0.05d0
c     tluma = tlfak * 0.5d0 *( tlum2 + tl(2,kk-1))
c     1              + (1.0d0-tlfak)*tlum1
c..   relax to thermal balance at join
c     write(*,'(a40)')
c     1              'tycho: L change HEAVILY damped'
c     write(*,'(3(a20,1pe12.3))')'tlum0',tlum0,
c     1              'tlum1',tlum1,'tlum2',tlum2
c     write(*,'(3(a20,1pe12.3))')'(L2-L1)/L1',
c     1              ( tlum2-tlum1 )/abs(tlum1),
c     2              '(L2-2L1+L0)/L1',
c     3              (tlum2-2.0d0*tlum1+tlum0)/abs(tlum1),
c     4              'L(est)',tluma
c     tl(2,kk) = tluma
c     endif
            
            if( tlum2 .gt. vl(nvmax(1),1) .or.
     1           tlum2 .lt. vl(nvmax(4),4) )then
               write(*,'(a30, 4(a7,1pe11.3) )')'tycho: L out of box: ', 
     2              'Lmax',vl(nvmax(1),1),'L',tlum2,
     3              'Lmin',vl(nvmax(4),4),'dlnL',tlum2/tlum1-1.0d0
            endif
            
            if( radius2 .gt. vr(nvmax(2),2) .or.
     1           radius2 .lt. vr(nvmax(5),5) )then     
               write(*,'(a30, 4(a7,1pe11.3) )')'tycho: R out of box: ', 
     2              'Rmax',vr(nvmax(2),2),'R',radius2,
     3              'Rmin',vr(nvmax(5),5),'dlnR',radius2/radius1-1.0d0
            endif
            
         endif
         
         do k = 1, kk+1
            r(1,k)  = r(2,k)
            tl(1,k) = tl(2,k)
            u(1,k)  = u(2,k)
            v(1,k)  = v(2,k)
            t(1,k)  = t(2,k)
            dv(1,k) = dv(2,k)
            dt(1,k) = dt(2,k)
            q(1,k)  = q(2,k)
            p(1,k)  = p(2,k)
            e(1,k)  = e(2,k)

            if( mset .eq. 0 )then
c..   update with new changes means accept x(n,k) which are updated
               do n = 1, ndim
c..   no negative x(n,k)
                  if( x(n,k) .le. 0.0d0 ) x(n,k) = 0.0d0
               enddo
            else
c..   freeze composition
               do n = 1, ndim
                  x(n,k) = xold(n,k)
               enddo
            endif
            
         enddo
         
         p(2,kk+1) = p(1,kk+1)
         t(2,kk+1) = t(1,kk+1)

c..   quick mixing for 87A limiting case
c     if( mode .eq. 0 )then
c     if( l .le. 1 )write(*,*)'QUIKMIX ENABLED'
c     call quikmix(kk,p,q,v,x,dmh)
c     endif


c..   map envelope onto grid........................................
         if( mapenv .eq. 1 )then
            write(*,'(/a10,a40)')'TYCHO:',
     1           'Mapping envelope onto Henyey grid'
            write(*,*)'Generate new stencil of envelope models'

c..   reset envelope integration so central value 3 exactly matches
c..   this is cleaner than interpolating existing values

            call fitenv(1)

            do k = 2, nvmax(3)-1
               write(*,'(i5,1p12e12.4)')k,vr(k,3),
     1              xm(kk+1)+vm(k,3),vp(k,3),
     1              grav*(xm(kk+1)+vm(k,3))/(pi4*vr(k,3)**4)/
     2              (vp(k+1,3)-vp(k-1,3))*(vm(k-1,3)-vm(k+1,3))
     3              -1.0d0
            enddo

            modes = 0
c            mode  = 0
c            write(*,*)'Reset hydro flag to mode = ',mode
            write(*,*)'Reset envelope flag to modes = ',modes
            write(*,*)'initial kk = ', kk
            write(*,*)'desired kk = ', kk +nvmax(3) -2
c..   save total mass to construct envelope mass coordinate
            xm(kk+nvmax(3)-1) = xm(kk+1)
            write(*,*) 'nvmax(3) ', nvmax(3),' xm(kk+1) ',xm(kk+1),
     1           ' xm(kk+1)/sol ',xm(kk+1)/sol
c..   test dimensions
            if( kk + nvmax(3 ) .ge. ktot )then
               write(*,*)'tycho: need more zones for envelope merge: '
               write(*,'(2(a20,i5))')'kk_nvmax',kk + nvmax(3 ),
     1              ' >= ktot ', ktot
               stop'tycho merge'
            endif
c..   avoid pkk+1 which is not correctly defined for this test
            write(*,'(a5,8a12)')'k','r','xm','pk','du/g'
            do k = kk-3,kk-1
               write(*,'(i5,1p8e12.3)')k,r(1,k),xm(k),p(1,k),
     1              a(k)*(p(1,k)-p(1,k+1))/dmi(k)/g(k) - 1.0d0
            enddo
c..   0.5*dmi becauseof extrapolation to zone boundary not zone center
            write(*,'(a25,1pe12.4)')'interior P bnd at kk',
     2           p(1,kk)-grav*xm(kk)/(pi4*r(1,kk)**4)*0.5d0*dmi(kk)
            write(*,'(a25,1pe12.4)')'exterior P bnd at kk',
     2           vp(nvmax(3),3)
            write(*,'(a25,1pe12.4)')'fractional error',
     2           (p(1,kk)-grav*xm(kk)/(pi4*r(1,kk)**4)*0.5d0*dmi(kk))
     2           /vp(nvmax(3),3)-1.0d0

c..   map outer zone abundances onto envelope (kk+1)
            do k = kk+1, kk + nvmax(3) - 1
               do n = 1, ndim
                  x(n,k) = x(n,kk)
               enddo
            enddo

c..add envelope to henyey grid
            do j = 1, nvmax(3)-1
               k = kk+j
               i = nvmax(3)-j
               xm(k)   = xm(kk+nvmax(3)-1) + vm(i,3)
               r(1,k)  = vr(i,3)
               a(k)    = pi4*r(1,k)**2
               g(k)    = grav*xm(k)/r(1,k)**2
               tl(1,k) = vl(i,3)
               u(1,k)  = u(1,kk)
c..   initial guess for temperature; to be iterated
               t(1,k)  = 0.5d0*(vtem(i,3)+vtem(i+1,3))

c..   hermite interpolation for pressure
               f0her = vp(i,3)
               f1her = vp(i+1,3)
               d0her = grav*( xm(kk+nvmax(3)-1) + vm(i,3) )
     1              /(pi4*vr(i,3)**4)
               d1her = grav*( xm(kk+nvmax(3)-1) + vm(i+1,3) )
     1              /(pi4*vr(i+1,3)**4)
               dmher = vm(i,3) - vm(i+1,3)

               call hermite( f0her,f1her,d0her,d1her,dmher,fhher)


               dmh(k)  = vm(i,3)-vm(i+1,3)
c..   specific volume from mass conservation
               v(1,k) = (pi43*(vr(i,3)**3 - vr(i+1,3)**3))/dmh(k)
c..   iterate T for desired pressure
               do n = 1, 30

                  call state(k,k,1)

                  dther  = (fhher-p(1,k))/pt(k)
c..slow convergence
                  if( n .gt. 10 )then
                     dther = dther*0.5d0
                  endif
                  t(1,k) = t(1,k) + dther

                  if( n .ge. 25 )then
                     write(*,'(2i5,1p8e13.5)')k,n,p(1,k),fhher,
     1                    fhher/p(1,k)-1.0d0,t(1,k),dther/t(1,k),
     2                    pt(k),yef(k)
                  endif
                  if( abs(dther/t(1,k) ) .lt. 1.0d-5 .and.
     1                 abs( fhher/p(1,k) -1.0d0) .lt. 1.0d-5 )goto 100
               enddo

               write(*,*)'tycho.f error in T iteration in envelope map'
               write(*,'(5(a12,i5))')'k',k,'n',n
               write(*,'(3(a12,1pe12.3))')'P(1,k)',p(1,k),
     1              'p target',fhher,'fract P',p(1,k)/fhher-1.0d0,
     2              'PT',pt(k)
               write(*,'(3(a12,1pe12.3))')'T(1,k)',t(1,k),
     1              'dT',dther,'fract T',dther/t(1,k)
                  
               stop'TYCHO: envelope mapping'
               
 100           continue
               
               h(k)    = vvel(i,3)
c..   not centered, will be recalculated from t and v below
c     p(1,k)  = vp(i,3)
               dnab(k) = vnab(i,3)
               dnrad(k) = vnrad(i,3)
               dnad(k)  = vnad(i,3)
               doux(k) = 0.0d0
               sound(k) = vsound(i,3)
               
c..   fill in abundances; envelope = outer henyey zone
c     do n = 1, ndim
c     x(n,k) = x(n,kk)
c     enddo
            enddo
            kk = kk + nvmax(3)-2
            write(*,*)'new kk = ', kk
            write(*,*)'nc = ',nc
            
c..   reset it (iteration index) and xold(n,kk) for tests in burn
            it = 0
            do n = 1, ndim
               xold(n,kk) = x(n,kk)
            enddo
            
c..   fill n=2 values over whole new grid
            do k = 1, kk+1
               r(2,k)  = r(1,k)
               tl(2,k) = tl(1,k)
               u(2,k)  = u(1,k)
               dr(k)   = 0.0d0
               dtl(k)  = 0.0d0
               t(2,k)  = t(1,k)
               v(2,k)  = v(1,k)
               dt(2,k) = 0.0d0
               dv(2,k) = 0.0d0
               dt(1,k) = 0.0d0
               dv(1,k) = 0.0d0
            enddo
            
c..   mass coordinate
            dmi(1)    = dmh(2)
            dmi(kk+1) = dmh(kk+1)
            do k = 2, kk
               xm(k)  = xm(k-1) + dmh(k)
               dmi(k) = 0.5d0*( dmh(k+1) + dmh(k) )
            enddo
c     if( modes .eq. 2 )then
c..   for consistency with hstat and fitenv
            dmi(kk) = dmh(kk)
c     endif
            xm(kk+1) = xm(kk) + dmh(kk+1)
            
c     call state(2,kk+1,1)
            call state(2,kk+1,2)
            
            call gtrs(p,q,e,v,r,u,xm,dmi,dmh,kk,nc,newt,mode)
            
c     do k = 2, kk
c     g(k) = grav*xm(k)/r(nc,k)**2
c     a(k) = pi4*r(nc,k)**2
c     enddo
c     if( r(nc,1) .gt. 0.0d0 )then
c     g(1) = grav*xm(1)/r(nc,1)**2
c     a(1) = pi4*r(nc,1)**2
c     else
c     g(1) = 0.0d0
c     a(1) = 0.0d0
c     endif
            
            do k = 1, kk+1
               if( k .ge. 2 .and. k .le. kk )then
                  du(k) = -g(k) +a(k)*(p(nc,k)-p(nc,k+1))/dmi(k)
               else
                  du(k) = 0.0d0
               endif
               if( g(k) .gt. 0.0d0 )then
                  write(*,'(i5,1p8e12.4)')k,du(k),g(k),du(k)/g(k)
               endif
            enddo
            
            call cinit(2,kk,nc)
            
c..   avoid convective constraint
            dthc = 0.0d0
c..   estimate best new time step, allowing courant time
            
            call timstp(dthc,period,ktype,kkman,2,kxmx)
            
            dth(1) = dth(2)
            dti(1) = dth(2)
            dti(2) = dth(2)
            time = 0.0d0
            ll = l+1
            
c..   keep velocities consistent with convection flags
            do k = 2, kk
               if( ic(k) .eq. 0 )then
                  h(k) = 0.0d0
               endif
            enddo
            
            write(*,'(2(a10,i5),a10,1pe12.3)')
     1           'mode =',mode,'modes =',modes,
     2           'dth =',dth(2)
            write(*,'(a30//)')'TYCHO: envelope now on grid'
            
c     go to 505
            call dout5(2)
            stop'envelope mapping'
            
         elseif( mapenv .eq. 2 )then
            write(*,*)'converting to radiative outer boundary'
            do k = 2,5
               write(*,'(2i5,1p12e11.3)')k,ic(k),xm(k),dmh(k),r(1,k),
     1              u(1,k),tl(1,k),t(1,k),v(1,k),p(1,k),g(k),
     2              a(k)*(p(1,k)-p(1,k+1))/dmi(k)-g(k),
     3              arad*t(1,k)**4/(3.0d0*p(1,k))
            enddo
            write(*,*)kk
            do k = kk-4,kk+1
               write(*,'(2i5,1p12e11.3)')k,ic(k),xm(k),dmh(k),r(1,k),
     1              u(1,k),tl(1,k),t(1,k),v(1,k),p(1,k),g(k),
     2              a(k)*(p(1,k)-p(1,k+1))/dmi(k)-g(k),
     3              arad*t(1,k)**4/(3.0d0*p(1,k))
            enddo
            stop'ty2'
cccccccccccccccccc
            
            go to 505
            
         endif




         
c     --------------
c     REZONE SECTION
c     --------------

         if( ktot .ne. 0 .and. mod(l,l4) .eq. 0 )then
            if( l .ne. 0 )then
c..   
               nc  = 1
               jj = kk
c..   kk may be changed upon return
               
c     debug: will print a dump file ??dmp?? 
c     where the leading ?? is the prefix and
c     the trailing ascends from 00 to the number of such dumps 
c     call dout5(1)
c     debug
               
c..   puts rezoned model in arrays (1,k), nc=1
c..   arrays (2,k) are not used until overwritten for next step
               
               call rezone(nc)
               
c..   interpolate for boundary conditions for envelope
c..   update kk
               if( njj .ne. jj .and. njj .ne. 0 )then
c..   rezoning has occurred, reset kk
                  kk = njj
               endif
               
c..   uniform rezoning, force no new zoning
               if( ktot .lt. 0 .and. jj .ne. kk ) ktot = 0
               
               s(1,kk+1) = 0.0d0
               s(2,kk+1) = 0.0d0
               s(2,kk  ) = 0.0d0
               
c..   reset it and xold(n,kk) for tests in burn
               it = 0
               do n = 1, ndim
                  xold(n,kk) = x(n,kk)
               enddo
               
               call state(2,kk,nc)
               
c..   no composition change due to burning 
c..   only nuclear and neutrino heating and cooling 
c..   burn(...,1,...)
c     call burn(nc,1,idebug,2,kk)
               
               call gtrs(p,q,e,v,r,u,xm,dmi,dmh,kk,nc,newt,mode)
               
               do k = 2, kk
                  g(k) = grav*xm(k)/r(nc,k)**2
                  a(k) = pi4*r(nc,k)**2
               enddo
               if( r(nc,1) .gt. 0.0d0 )then
                  g(1) = grav*xm(1)/r(nc,1)**2
                  a(1) = pi4*r(nc,1)**2
               else
                  g(1) = 0.0d0
                  a(1) = 0.0d0
               endif
               
               call cinit(2,kk,nc)
c..   
               if( ismoo .lt. 0  )then
                  write(*,'(a5,10a11)')'k','xm','dmi','dmh','g',
     1                 'du','A/dmi','pk','pk+1','rk'
                  dmh(kk+1) = dmh(kk)
                  dmi(kk+1) = dmh(kk+1)
                  do k = 1, kk+1
                     if( r(1,k) .gt. 0.0d0 )then
                        a(k)  = pi4*r(1,k)**2
                        g(k)  = grav*xm(k)/r(1,k)**2
                        du(k) = (p(1,k)-p(1,k+1))*a(k)/dmi(k) - g(k)
                     else
                        a(k)  = 0.0d0
                        g(k)  = 0.0d0
                        du(k) = 0.0d0
                     endif
                     
                     write(*,'(i5,1p10e11.3)')k,xm(k),dmi(k),dmh(k),
     1                    g(k),du(k),a(k)/dmi(k),p(1,k),p(1,k+1),r(1,k)
                  enddo
                  write(*,'(a5,10a11)')'k','xm','dmi','dmh','g','du',
     1                 'A/dmi','pk','pk+1','rk'
                  write(*,'(a50)')
     1                 'equal mass rezone: stop and dump new model'
                  
                  write(*,'(3(a10,i7))')'ismoo',ismoo,'kk',kk,'L',l
                  ll = l  +1
c..   change surface boundary condition to mesh with new zoning
                  modes = 0
                  kk = kk-1
               endif   
               
               
c     debug
c     call dout5(1)
c     debug
               
            endif
         endif
      endif
      
      sum = 0.0d0
      do k = 2, kk
         sum = sum + dmh(k)*ss(k)
      enddo
      
      
c     --------------------------------------------------------
c     new values defined or original values used with reduced 
c     time step go back to solution loop to construct state at 
c     next time step
c     
c     
      go to 505                                     
c     ------------------- END EVOLUTION LOOP -----------------

      
      
      
      
c     
      end






