      subroutine gtintp

c     this subroutine contains the logic for addition and deletion
c     of zones, and defines new mass coordinate arrays
c     qqn at boundaries, qq2n at centers ("n" for new, "o" for old)
c     jj is old kk
c     jjn is new kk
c     netadd is net number of added zones (negative for net deletion)

c..   rewritten 5/28/01, modified 1/4/03, 9/4/04, 9-28-05, 2-16-06,
c..   3-29-06
c..   trymin(i): delete a boundary if ALL functions try(i,j) are less
c..   trymax(i): add a boundary in a zone center if ANY function 
c..   exceeds the corresponding limit
c..   shell flame must have no larger just ahead


      implicit none

      include 'dimenfile'
      include 'comod'
      include 'compu'
      include 'czone'
      include 'cburn'
      include 'cnabla'
      include 'cenv'
      include 'cgtintp'
      include 'cconst' 
c..   ntests: number of rezoning criteria

c..   at boundaries:
c..   trial 1: d( ln rho )/dk
c..   trial 2: d( ln T   )/dk
c..   trial 3: zone ratio (or increase rate)
c..   trial 4: he4 gradient

c..   at zone centers:
c..   trial 5: flame zone
c..   trial 6: zone size (global function)
c..   trial 7: V linearity
c..   trial 8: P linearity
c..   trial 9: step size in radius

c..   added constraint on zones ahead of active flame
      integer*4 ks5max, kshell,ks5lo,ks5hi
      real*8    dms5kmax,s5kmax

c..   old and new mass coordinates:
c     qqo(kdm),qqn(kdm),qq2o(kdm),qq2n(kdm),dqqn(kdm)

      real*8    dqqs(kdm),dqqsn(kdm)
      integer*4 ispan

c..   nadd is number of zones added
c..   ndel is number of zones deleted
c     integer*4 jj,njj,nadd
      integer*4 i,j,ndel,nadd, n

c..   flags for beginning and ending of convective regions
c..   moved to cnabla
c      integer*4 kbeg(kdm),kend(kdm),nrczones

c..   nrczones is number of separate convective regions
c..   kcbeg and kcend are zones indices of a given convective region
c..   icold is updating index for convective region search
      integer*4 kcbeg,kcend,icold
c..   xmconv: working variable for convective region mass (not saved)
      real*8 xmconv

c..   jover is number of zones over that allowed by dimension
c..   jdel  is update index used by logic to keep jover in bounds
c..   jn is working index for new mass array qqn
      integer*4 jover,jdel,jn
c..   junder is number of zones allowed for deletes
      integer*4 junder,jddl
      
c..   xlnuc is summation variable for nuclear luminosity
c..   xlnuc is used for flame rezone
      real*8    xlnuc

c..   fact, rat, fak are working variables
      real*8    fact,rat,fak

c..   facte defined by gen, denotes fraction of nuclear luminosity
c..   allowed in a zone

c..   dlnt is d ln T per boundary
      real*8    dlnt

c..   working variables for abundance gradient criterion to save zone
c..   boundaries from being split
      real*8    delylim, delymin, delyfak

c..   xmchand is mass of core in which zoning is refined
c..   fmcore is desired maximum fractional mass
c..   of a zone in this region
c..   dlnp    : no added zones due to density jump if 
c     d ln P/dk is less 
c..   dmdlnp is desired range from log gradient of pressure
c..   tnucool is boundary temperature separating regimes
c..   of cooling by neutrinos (.gt.) and photons (.le.)
      real*8    xmchand, fmcore, dmdlnp, tnucool

      parameter( xmchand = 1.5d0*sol, fmcore = 0.01d0,
     1     tnucool = 1.0d9 )
c..   avoid flame division for Ne O Si burning

c..   pseudo-diffusion coefficient for zone mass adjustment
      real*8 dcoe(kdm),gradx(kdm)

c..   free particles per nucleon
      real*8 ytot(kdm),dlydm
c..   array which returns qqo index given the qqn index
      integer*4 ntoo(kdm),jp1,jp0

      character*3 cout(2*ntests+2),cblank

      data delylim / 1.0d-3 /, delymin/ 1.0d-4 /
      data ctry/'rho', 'Tem','shl', 'Yto', 'fla', 'siz', 'vli',
     1     'pli','dlr'/
      data cblank/'   '/
      data drmax/1.0d-2/
c----------------------------------------------------------------
c..   input:
c     qqo,jj, all comod variables
c..   output:
c     qqn,qq2o,qq2n,njj,nadd
c----------------------------------------------------------------

      if( mode .ne. 1 .and. mode .ne. 2 .and. mode .ne. 3 )then
         write(*,*)'Attempting to use hydrostatic rezone in mode',mode
         write(*,*)'change ktot or mode in input file'
c..   or modify to add hydrodynamic rezoning here
         stop'gitintp'
      endif

c..   set up ytotal array on qqo 
      do i = 2, jj
         ytot(i) = 0.0d0
         do n = 1, nnuc
            ytot(i) = ytot(i) + x(n,i)*(1.0d0 + dble(lz(n)))
         enddo
      enddo

c..   interpolate for old half interval mass coordinates
      do j = 1,jj-1
         qq2o(j+1) = 0.5d0 * (qqo(j) + qqo(j+1) )
      enddo
      qq2o(1) = -qq2o(2)

c..   constraints on log gradients
      dlnt = dlnv/3.0d0
c..   defines edges of acceptable boundary
      trymin(1) = 0.3d0
      trymax(1) = 1.0d0
c..   allow density inversion at join by trymax(1) .gt. 1
      trymin(2) = 0.3d0
      trymax(2) = 1.0d0
c..   thin shell (default is del)
      trymin(3) = 0.0d0
      trymax(3) = 1.0d0

c..   prevent merger if abundance gradient is greater than this
      trymin(4) = 0.03d0
c..   divide zones if greater
c..   out if reach so do not ask for rezone due to abundance gradient
      trymax(4) = 1.0d10

c..   boundaries
c..   flame
      trymin(5) = 0.1d0
      trymax(5) = 1.0d0
c..   size (smooth increase in mass from analytic function)
cccccccccccccxmchand not used
      trymin(6) = 0.4d0
      trymax(6) = 1.0d0
c..   v linearity
      trymin(7) = 0.02d0 * vline
      trymax(7) = 0.06d0 * vline
c..   p linearity
      trymin(8) = 0.005d0 * pline
      trymax(8) = 0.06d0  * pline

c..   zone radius difference (fractional); drmax from params.d
      trymin(9) = 0.1d0*drmax
      trymax(9) =       drmax

c.....find radiative-convective boundaries........................
      nrczones = 0
      xmconv = 0
      kcbeg = 0
      kcend = 0
      icold = ic(1)
      ic(1) = 0
      do j = 2, jj
         if( ic(j-1) .ne. 1 .and. ic(j) .eq. 1 )then
c..   begin cz
            nrczones = nrczones + 1
            kcbeg = j-1
            kbeg(nrczones) = kcbeg
c      write(*,*)'beg ',xmconv,j-1,kcbeg,ic(j-1),ic(j)
         endif
         if( ic(j) .eq. 1 )then
c..   sum mass in current cz
            xmconv = xmconv + dmh(j)
         endif
         if( (ic(j-1) .eq. 1 .and. ic(j) .ne. 1) )then
c     1        .or. (j .eq. kk) )then
c..   past zone convective but current zone is not
c..   ending cz
            xmconv = xmconv + dmh(j)
            kcend = j
            kend(nrczones) = kcend
c      write(*,*)'end ',xmconv,j,kcend,ic(j-1),ic(j)
         endif
      enddo
c      write(*,*)'gtintp kcbeg,kcend ',kcbeg,kcend,ic(kk-1),ic(kk)
      ic(1) = icold
      if( nrczones .gt. 0 )then
         if( kcbeg .ge. kcend )then
c..   unfinished outer convection zone -- extends beyond join
c..   just extend to join in logic here (kk; kk+1 is whole envelope)
            xmconv = xmconv + dmh(kk)
            kcend = kk
            kend(nrczones) = kcend
         endif
      endif

      if( nrczones .gt. 1 )then
         write(*,*)'GTINTP: ',nrczones,
     1        ' convection zones on grid'

         do j = 1, nrczones
            write(*,'(i5,2(i5,1pe12.4))')j,kbeg(j),xm(kbeg(j))/sol,
     1           kend(j),xm(kend(j))/sol
         enddo
      endif

c.....end of find radiative-convective boundaries..................

c..   flame zone constraint: calculate total nuclear luminosity
c..   for photon cooled regimes
      xlnuc = 0.0d0
      do j = 2, jj
         if( t(1,j) .le. tnucool )then
            xlnuc = xlnuc + ss(j) * dmh(j)
         endif
      enddo
c..   maximun flame
      ks5max = 0
      s5kmax = 0.0d0
      do j = 2, jj
         if( ss(j) .gt. s5kmax )then
            s5kmax = ss(j)
            ks5max = j
         endif
      enddo

c..   small or no flame if in helmholtz-kelvin contraction
      if( t(1,2) .lt. 1.0d7 )then
         xlnuc = dmax1( xlnuc, tl(1,jj) )
      endif

c..   define test functions "try"
      do j = 2, jj
c..   log gradients in density and temperature
c..   use rho*Ye to avoid excess rezoning at H/He interface
c..   scale as dv/v + dt/t = dp/p
         if( j .lt. jj )then
            try(1,j) = abs( dlog10( v(1,j+1)/v(1,j) ) /dlnv )
            try(2,j) = abs( dlog10( t(1,j)/t(1,j+1) ) /dlnt )
         else
c..   avoid crossing join
            try(1,j) = 0.0d0
            try(2,j) = 0.0d0
         endif

c..   zone ratio
c     try(3,j) = dmh(j+1)/dmh(j) 

c..   abundance jumps, uses Ytot = sum Yi for all nuclei and electrons
         if( j .eq. 2 .or. j .eq. jj )then
            try(4,j) = 0.0d0
         else
            try(4,j) = abs( ytot(j+1) - ytot(j) )
         endif

c..   flame zone
         if( facte*xlnuc .gt. 0.0d0 .and. t(1,j) .le. tnucool )then
            try(5,j) = dmh(j)*s(5,j)/(facte * xlnuc)
         else
            try(5,j) = 0.0d0
         endif

c..   use smooth interpolation in mass
c..   smooth increase in zone mass from origin
         rat  = dmh(2)*1.1d0**(j-2)
c..   smooth transition to global constraint
c..   rat is the tentative zone mass extrapolated out from origin
c..   fact is the final interpolated (tentative) zone mass

         fact = 1.0d0/(
     1        1.0d0/rat
     2        + 1.0d0/( dmmax*xm(jj+1) ) )
ccccccccccccc xmchand not used here anymore
         try(6,j) = dmh(j)/fact *0.75d0

c         write(*,'(i5,1p8e12.3)')j,try(6,j),fact,rat,dmh(j),dmh(kk),
c     1        xm(kk)/dmh(j)
cccccccccccccccc

c..   v linearity
         if( j .eq. 2 .or. j .eq. jj )then
c..   leave these zones alone
            try(7,j) = trymin(7) * 1.01d0
         else
c..   same definition as getvec.f (for "r3")
            try(7,j) = pi43*abs( 
     1           r(1,j+1)**3/(dmh(j+1)*v(1,j+1))
     2           - r(1,j  )**3/(dmh(j+1)*v(1,j+1)) 
     3           - r(1,j  )**3/(dmh(j  )*v(1,j  )) 
     4           + r(1,j-1)**3/(dmh(j  )*v(1,j  )) 
     5           )
         endif

c..   p linearity
         if( j .eq. 2 .or. j .eq. jj )then
c..   leave these zones alone
            try(8,j) = trymin(8) * 1.01d0
         else
            try(8,j) = abs( grav/pi4/p(1,j)*( xm(j)*dmi(j)/r(1,j)**4
     1           -xm(j-1)*dmi(j-1)/r(1,j-1)**4 ) )
         endif

c         write(*,'(i5,1p8e12.3)')j,try(8,j)-trymin(8),
c     1        try(8,j)-trymax(8),
c     1        try(8,j),xm(j-1)/sol,xm(j)/sol

c..   radius width (d ln r)
         if( j .eq. 2 .or. j .eq. jj )then
c..   leave these zones alone
            try(9,j) = trymin(9) * 1.01d0
         else
            try(9,j) = (r(1,j)-r(1,j-1))/r(1,jj)
         endif

      enddo

c      stop'aaa gtintp'
cccccccccccccccccc

c..   set up logic flags from possible rezones
c..   initialize
      do j = 1, kdm
         iadd(j) = 0
         idel(j) = 0
      enddo
      do j = 1, kdm
         do i = 1,ntests
            iad(i,j) = 0
            ide(i,j) = 0
         enddo
      enddo

c..   force single zone for kforce .ne. 0
      if( kforce .eq. 0 )then
c..   multiple zone (usual case)
c..   adds
         do j = 2, jj
            do i = 1, ntests
               if( try(i,j) .gt. trymax(i) )then
                  iad(i,j) = 1
               endif
            enddo
         enddo
c..   dels
         do j = 2, jj
            do i = 1, ntests
               if( try(i,j) .le. trymin(i) )then
                  ide(i,j) = 1
               endif
            enddo

c            ispan = 0
c            do i = 1, ntests
c               ispan = ispan + ide(i,j)
c            enddo           
c            write(*,'(i5,9i3,i5)')j,(ide(i,j),i=1,ntests),ispan
c            if( ispan .eq. ntests )then
c               write(*,'(2i5,1p8e12.3)')j,ispan,xm(j)/sol,dmh(j)/sol
c            endif
ccccccccccccccccccccccccccccccccc
         enddo


c..   look for any condition for adds
         do j = 2, jj
            if( iadd(j) .eq. 0 )then
c..   adjust according to tests
               iadd(j) = iad(1,j)
               do i = 2, ntests
                  iadd(j) = iadd(j) + iad(i,j)
               enddo
               if( iadd(j) .ge. 1 )then
                  iadd(j) = 1
c..   split zones on both sides of a boundary for conditions 1 through 4
                  if( iad(1,j) .eq. 1 .or. iad(2,j) .eq. 1 .or.
     1                 iad(3,j) .eq. 1 .or. iad(4,j) .eq. 1 )then
                     iadd(j+1) = 1
c..   avoid extreme zone mass differences: split only the larger zone
                     if( dmh(j) .lt. 0.3d0*dmh(j+1) )then
                        iadd(j)   = 0
                        iadd(j+1) = 1
                     elseif( dmh(j+1) .lt. 0.3d0*dmh(j) )then
                        iadd(j)   = 1
                        iadd(j+1) = 0
                     endif
                  endif
               else
                  iadd(j) = 0
               endif
            else
c..   already demanded so set to unity
               iadd(j) = 1
            endif
         enddo

c..   do not rezone next to convective boundary for a big convective zone
c..   useful for dredge up
c         if( nrczones .gt. 0 )then
c..   there are convective zones
c            do j = 1, nrczones
c..   use only convective regions with many zones
c               if( kend(j) - kbeg(j) .gt. 10 )then
c                  do i = 2, jj
c                     if( iadd(i) .ne. 0 )then
c                        if( i - kbeg(j) .ge. -1 .and. 
c     1                       i - kbeg(j) .le. 2 )then
c                           write(*,'(a30,5(a10,i5))')
c     1                          'GTINTP: no add near conv bnd ',
c     2                          'k',i,'nzone',j,'begin',kbeg(j),
c     3                          'kk',jj
c..   too near a convective boundary
c                           iadd(i) = 0
c                        endif
c                        if( i - kend(j) .ge. -1 .and. 
c     1                       i - kend(j) .le. 2 )then
c                           write(*,'(a30,5(a10,i5))')
c     1                          'GTINTP: no add near conv bnd ',
c     2                          'k',i,'nzone',j,
c     3                          'end',kend(j),'kk',jj
c..   too near a convective boundary
c                           iadd(i) = 0
c                        endif
c                     endif
c                  enddo
c               endif
c            enddo
c         endif

c..   additional criteria for thin shell sources
c..   do not rezone most active zones unless in core
c..   avoid origin (100 seems ok)
         kshell = 100
cccccccccccccccccccccccccccccc use mass estimate (mchanra, m core expans.)
c..   if most vigorous flame is away from origin
         if( ks5max .ge. kshell )then
            do j = kshell, jj
               if( s(5,j) .gt. 0.3d0*s5kmax )then
c..   no deletes or adds near flame
                  idel(j) = 0
                  iadd(j) = 0
               endif
            enddo
         endif

c..   keep zones ahead of shell smaller than zone with peak flame
c..   shell source
         do j = kshell, jj
c            if( j .gt. ks5max .and. j .le. ks5max+20 )then
            if( j .gt. ks5max .and. j .le. ks5max+30 )then

c..   add well ahead of flame (0.75 looks ok)
               if( dmh(j) .gt. dmh(ks5max) .and. 
     1              ss(j) .lt. 0.75d0*s5kmax )then
                  iadd(j)  = 1
                  iad(3,j) = 1
               else
c..no delete before flame
                  idel(j)  = 0
                  ide(3,j) = 0
               endif
            endif
         enddo

c..   restrict rezoning in composition gradients
c..   takes priority over previous criteria
         do j = 2, jj
            if( abs(ytot(j+1)-ytot(j)) .gt. 1.0d-3 .and. 
     1           iadd(j) .ne. 0 )then
               if( xm(ks5max) .lt. 2.0d0*xmchand )then
c..   shell burning in less massive stars
                  iadd(j) = 0
               endif
            endif
         enddo
c..   no delete of zones with significant composition gradients
c..   uses all nuclei
         do j = 2, jj
c..   compute abs sum of gradients 
            delyfak = 0.0d0
            do n = 1, nnuc
               if( x(n,j) .gt. delymin )then
                  delyfak = delyfak + abs( xa(n)*(x(n,j+1)-x(n,j) ) )
               endif
            enddo
c..   delete only if abs gradient sum is small
            if( delyfak .gt. delylim )then
               idel(j) = 0
            endif
         enddo

c..   do not rezone near semiconvective zones
c..   relies on cinit.f definition of ic(j) 
c..   which is Richardson with convective shear
c..   add at zone center
         do j = 2, jj
            if( iadd(j) .eq. 1 )then
               if( ic(j-1) .eq. 2 .or. ic(j) .eq. 2 )then
                  iadd(j) = 0
               endif
            endif
         enddo
c..   del at zone boundary
         do j = 2, jj-1
            if( idel(j) .eq. 1 )then
               if( ic(j-1) .eq. 2 .or. ic(j) .eq. 2 .or.
     1              ic(j+1) .eq. 2 )then
                  idel(j) = 0
               endif
            endif
         enddo

c..   shell flame thinning
c..   overrides previous criteria
         if( t(2,ks5max) .lt. 0.5d9 )then
c..   avoid delete in flame zone
            ks5lo = max0(2,ks5max-30)
            ks5hi = min0(jj,ks5max+30)
            do j = ks5lo,ks5hi
               idel(j) = 0
            enddo
c..   only for intermediate mass evolution (H He burning)
            if( tl(2,ks5max) .gt. 0.0d0 )then
               if( dmh(ks5max)*ss(ks5max)/tl(2,ks5max) .gt. 0.08d0 )then
                  do j = ks5max,ks5max+30
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc0.1
c..   add if zones are not yet dominant
                     if( ic(j) .eq. 0 .and. ss(j)/s5kmax .lt. 0.2d0 
     1                    .and. dmh(j) .ge. 0.25d0*dmh(ks5max) )then
                        write(*,'(a25,3i5,1p8e12.3)')
     1                       'forcing add ahead of shell',j,ic(j),
     1                       ks5max,ss(j),ss(j)/s5kmax,
     1                       dmh(j)/sol,0.5d0*dmh(ks5max)/sol
                        iadd(j) = 1
                        goto 120
                     endif
                  enddo
c..   only one added for each rezone to reduce impact
 120              continue
               endif
            endif
         endif

c..   look for all conditions for dels
         do j = 2, jj
            idel(j) = ide(1,j)
            do i = 2, ntests
               idel(j) = idel(j) + ide(i,j)
            enddo
            if( idel(j) .eq. ntests )then
               idel(j) = 1
            else
               idel(j) = 0
            endif
         enddo
c..   no delete if it worsens the mass ratio
         do j = 2, jj
c..   centered on boundary j
            if( idel(j) .eq. 1 )then
               if( dmh(j)+dmh(j+1) .gt. 1.5d0*(dmh(j-1)+dmh(j+2)) )then
                  idel(j) = 0
                  ide(3,j) = 0
               endif
            endif
         enddo
c..   absolute limit of minimum zone size: add no more zones
         do j = 2, jj
            if( iadd(j) .eq. 1 )then
               if( dmh(j) .lt. 1.0d-7*xm(jj) )then
                  iadd(j) = 0
               endif
            endif
         enddo

         if( modes .eq. 2 )then
c..   envelope fitting
c..   use next to last zone in envelope integration to avoid
c..   truncation to fit in mass of envelope
c..   does not depend upon envelope solution
            dmdlnp = 0.01d0*p(2,kk)* pi4*r(2,kk)**4/(grav*xm(kk))
         else
            dmdlnp = 0.1d0*( xm(jj+1) - xm(jj) )
            dmdlnp = dmin1( dmdlnp, dmh(jj-1) )
         endif

c         write(*,'(a20,i5,1p8e12.3)')'gtintp',jj,dmh(jj),dmh(jj-1),
c     1        dmdlnp,xm(jj+1)-xm(jj)
c         stop'gtintp bbb'
cccccccccc

c..   outer zone on grid
         if( dmh(jj) .gt. dmdlnp )then
c..   outer zone too big, split it
            iadd(jj)   = 1
         elseif( dmh(jj) .lt. 0.4d0*dmdlnp )then
c..   outer zone too small, merge it
            idel(jj-1) = 1
         else
            idel(jj-1) = 0
         endif

c..   never delete boundary at join
         idel(jj) = 0

c..   limit on inner zone mass as fraction of total
         if( qqo(2)-qqo(1) .gt. 1.0d-6*qqo(jj) )then
            iadd(2) = 1
         elseif( qqo(2)-qqo(1) .lt. 1.0d-7*qqo(jj) )then
            idel(2) = 1
         else
c..   force no change
            idel(2) = 0
            iadd(2) = 0
         endif

c..   adds have priority; no dels within two zones of an add
         do j = 2, jj
            if( iadd(j) .eq. 1 )then
               if( j .eq. 2 )then
                  idel(j-1) = 0
                  idel(j) = 0
                  idel(j+1) = 0
               elseif( j .eq. jj )then
                  idel(j-2) = 0
                  idel(j-1) = 0
                  idel(j)   = 0
               else
c..   centered on k-1/2="K" for zone center
                  idel(j-2) = 0
                  idel(j-1) = 0
                  idel(j)   = 0
                  idel(j+1) = 0
               endif
            endif
         enddo

c..   no neighboring boundary deletes
         do j = 2, jj
            if( idel(j) .eq. 1 )then
               idel(j+1) = 0
            endif
         enddo
c..   neighboring boundary adds allowed

      elseif( kforce .gt. 0 .and. kforce .lt. jj )then
c..   force single add
         iadd(kforce) = 1

      elseif( kforce .lt. 0 .and. -kforce .lt. jj )then
c..   force single delete
         idel(kforce) = 1

      endif


c..   summarize changes
      nadd = 0
      ndel = 0
      do j = 2, jj
         if( iadd(j) .eq. 1 )then
            nadd = nadd + 1
         elseif( idel(j) .eq. 1 )then
            ndel = ndel + 1
         endif
      enddo

c..   restrict to maximum number of zones
c..   restrict number of added zones per rezone
c..   3 extra zones to pad array for envelope values
      jover =  jj + nadd - ndel - ktot + 3
      if( nadd .gt. jj/10 )then
c..   at most 5 percent new zones
         jover = nadd - jj/10
         write(*,*)'adjusting nadd to smaller value ',nadd,jover
      endif
      jdel = 0
      if( jover .gt. 0 )then
         write(*,*)jover,'  excess zones requested in gtintp'
         write(*,'(5(a7,i5))')'jj',jj,'add',nadd,'del',ndel,
     1        'ktot',ktot,'jover',jover
c..   use outside-in loop to negate outer request first
         do i = 2, jj-1
            j = jj + 1 - i
            if( iadd(j) .ne. 0 )then
               iadd(j) = 0
               jdel = jdel + 1
            endif
            if( jdel .ge. jover )goto 100
         enddo
 100     continue
      endif

c..   restrict number of deleted zones per rezone
      junder = ndel - jj/10
c      write(*,*)nadd,ndel,jj,junder
      jddl = 0
      if( junder .gt. 0 )then
         write(*,*)junder,' excess dels requested in gtintp'
         write(*,'(5(a7,i5))')'jj',jj,'add',nadd,'del',ndel,
     1        'ktot',ktot,'junder',junder
c..   use outside-in loop to erase outer request first
         do i = 2, jj-1
            j = jj + 1 - i
            if( idel(j) .ne. 0 )then
               idel(j) = 0
               jddl    = jddl + 1
            endif
            if( jddl .ge. junder )goto 101
         enddo
 101     continue
      endif


c..   net change in zone number
      netadd = nadd - ndel - jdel +jddl


c..   initialize new mass coordinates
      do j = 1, kdm
         qqn(j)   = 0.0d0
         qq2n(j)  = 0.0d0
         dqqsn(j) = 0.0d0
      enddo

c..   define new array with add/deletes
      qqn(1) = qqo(1)
      jn     = 2
      do j = 2, jj
         if( iadd(j) .eq. 1 )then
c..   add new boundary at midpoint of zone
c..   avoid extrapolation over join boundary
c..   1--->2
            qqn(jn)   = 0.5d0*( qqo(j) + qqo(j-1) )
            qqn(jn+1) = qqo(j)
            jn = jn + 2
         elseif( idel(j) .eq. 1 )then
            if( iad(5,j) .eq. 1 .or. iad(6,j) .eq. 1 )then
c..   merge boundaries (omit zone)
c..   3--->2
               qqn(jn) = qqo(j+1)
               qqn(jn-1) = 0.5d0*( qqo(j+1)+qqo(j-1) )
            else
c..   delete boundary
c..   2--->1
               qqn(jn) = qqo(j+1)      
            endif
         else
c..   no new zone
            qqn(jn) = qqo(j)
            jn = jn + 1
         endif
      enddo

      do j = 1, jj
         do i = 1, 2*ntests+2
            cout(i) = cblank
         enddo
         if( iadd(j) .eq. 1 .or. idel(j) .eq. 1 )then
            if( iadd(j) .eq. 1 )cout(1) = 'add'
            if( idel(j) .eq. 1 )cout(2) = 'del'
            do i = 1, ntests
               if( iad(i,j) .eq. 1 )then
                  cout(i+2) = ctry(i)
               endif
            enddo

            do i = 1, ntests
               if( ide(i,j) .eq. 1 )then
                  cout(i+2+ntests) = ctry(i)
               endif
            enddo
            if( kforce .gt. 0 )then
               cout(3) = 'f+'
            elseif( kforce .lt. 0 )then
               cout(3) = 'f-'
            endif
c..   trigger matrix for adds and dels
         endif
      enddo

c..   new "kk" index
      njj = jn - 1
      if( njj .ne. jj )then
c         write(*,'(65x,9a4,5x,9a4)')ctry,ctry
c         write(*,'(a6,5a11,a8,a36,a5,a36)')'old k','mass','dmass',
c     1        'h','ss','dYtot',
c     1        'action',
c     1        'addition criteria triggered','_____',
c     2        'deletion criteria triggered'
      endif


      qqn(njj+1) = qqo(jj+1)

c--------------------------------------------------------------
c..   smooth zone mass irregularities
c..   ntoo array contains index of qqo given index in qqn
c..   returns 0 if there is no such array (qqo is split by new bnd.)
      do i = 1, kdm
         ntoo(i) = 0
      enddo

      i = 1
      do j = 2, jj
         i = i + 1 + iadd(j) - idel(j)
         ntoo(i) = j

      enddo

c..   set up ytotal array on qqo
      do i = 2, jj
         ytot(i) = 0.0d0
         do n = 1, nnuc
            ytot(i) = ytot(i) + x(n,i)*(1.0d0 + dble(lz(n)))
         enddo
      enddo

c..   set up zone mass array
      do j = 2, njj
         dqqn(j) = qqn(j) - qqn(j-1)
      enddo
      dqqn(1) = dqqn(2)
c..   avoid origin and join zones and neighbors in qqn array
c..   relax halfway to uniform zoning
      fact = 0.25d0
      do j = 3, njj-2 
         fak = (qqn(j+1)-2.0d0*qqn(j)+qqn(j-1))/(qqn(j+1)-qqn(j-1))
         jp0 = ntoo(j)
         jp1 = ntoo(j+1)
c..   fix zeros from zone splits
         if( jp0 .eq. 0 )then
            jp0 = ntoo(j-1)

         elseif( jp1 .eq. 0 )then
            jp1 = ntoo(j+2)

         endif

         dlydm = (ytot(jp1)-ytot(jp0))/(ytot(jp1)+ytot(jp0))
     1        *2.0d0*qqo(jj)/(qqo(jp1)-qqo(jp0))

         qqn(j) = qqn(j)*(1.0d0 - 2.0d0*fact) +
     1        fact*(qqn(j+1) + qqn(j-1))

      enddo

c      write(*,*)'NEW ARRAY QQN IS SMOOTHED IN GTINTP.F'

c-----------------------------------------------------------------

c..   adjust next to join zone, over-riding previous conditions
c..   this is second-order accurate
      qqn(njj-1) = 0.5d0*( qqn(njj)+qqn(njj-2) )

c..   grid adjustments are done, get ready to interpolate
c..   new half interval mass coordinate
      do j = 1, njj-1
         qq2n(j+1) = 0.5d0 * (qqn(j) + qqn(j+1) )
      enddo
      qq2n(1) = -qq2n(2)

      do j = 2, njj
         dqqn(j) = qqn(j) - qqn(j-1)
      enddo
      dqqn(1) = dqqn(2)

c..only force zoning once
      kforce = 0

c--------------------------------------------------------------
c      write(*,*)'LEAVING GTINTP ',nadd,ndel,netadd,jdel
c--------------------------------------------------------------
      return
      end


