

      subroutine opal(t6,ropal,z,xh,xxc,xxo,ak,akt,akv,k,kk,nopac,l)
      implicit none

c..   newopal.f  last modified 5-8-06  darnett
c..   added L (step number) to allow sparse output when debugging
c..   in call sequence above and in call sequence for opacgn93

c..   using payoung's newopal.f, 
c..   added effective metalicity zalex for alexander opacities
c..   revised boundary and out-of-table procedures
c..   revised logic for Klein-Nishina limit
c..   added nopac = 1 option to use type 1 OPAL tables with opacgn93.f

c..   wrapper for OPAL 1996 and Alexander opacities
c..   input t6, r, z, xh, xxc, xxo
c..   output ak,atk,akv
c..   ropal = rho/t6**3

      real*8 t6,ropal,z,xh,xxc,xxo,ak,akt,akv
      integer*4 k,kk,l

      real*8 rlg,tlg
      real*8 tlow1,thigh1,tlow2,thigh2, rlghi,rlglo, nolghi, nolglo
      real*8 tlglo1,tlghi1,tlglo2,tlghi2, ratrlg,akcoef
      real*8 zalex

      real*8 dopact2,dopact1,xtlg,onemf

c..   tlow1  temperature below which atmospheric opacities are used
c..   now uses alexander for atmospheric opacities
c..   thigh1 temperature above which atmospheric opacities are not used
c..   morph in between
c..   tlow2  temperature above which we morph to thomson/klein-nishina
c..   thigh2 temperature above which thomson/klein-nishina is used

      parameter( tlow1 = 0.007d0, thigh1 = 0.009d0,
     1     tlow2 = 90.0d0, thigh2 = 200.0d0 )

      parameter( tlglo1 = 3.845d0, tlghi1 = 3.954d0,
     1     tlglo2 = 7.954d0, tlghi2 = 8.301d0 )

c..   limits of log10 Ropal in OPAL96 tables
      parameter( rlghi = 1.0d0, rlglo = -8.0d0 )

      real*8 opact2,dopacr2,opact1,dopacr1
      real*8 theta,dlgt,fact,akr

      integer*4 kopl,nopac

      real*8 opact,dopact,dopacr,dopactd
      common/e/ opact,dopact,dopacr,dopactd

      save


c...............................................................
c..   maximum Ropal in OPAL 1996 tables is log10 Ropal = +1.0 
c..   minimum Ropal in OPAL 1996 tables is log10 Ropal = -8.0
c...............................................................


      if( ropal .gt. 0.0d0 )then
         rlg = dlog10( ropal )
      else
         write(*,*)'OPAL.f: input error, Ropal ',Ropal
         stop'opal 1'
      endif


      if( t6 .gt. 0.0d0 )then
         tlg = dlog10( t6 ) + 6.0d0
      else
         write(*,*)'OPAL.f: input error, T6 ',t6
         stop'opal 2'
      endif


      kopl = 0

c<<<<<<< .mine

      nolghi = 10.0d0**rlghi
      nolglo = 10.0d0**rlglo
c=======

c...............................................................
c>>>>>>> .r98
c..   OPAL tables have alpha element enhancements xxc and xxo
c..   use effective metalicities for alexander tables
c..   zalex = z + xxc + xxo
c..   forces consistency with OPAL type 2 tables; no z change
c...............................................................
      zalex = z
         
      
c..   limit size of ropal to keep sane at advanced stages
c..   extrapolate low ropal by shifted lower temperature
      if( rlg .lt. rlglo  )then
c..   high entropy, low Ropal

         
         if( tlg .lt. tlglo1 )then
c..   low T, use atmospheric opacities (alexander)
c..   variable kopl signals which logic branch was used (debugging)
            
            kopl = 21
            call kapp(tlg,rlg,ak,akr,akt,zalex,xh)
            
         elseif( tlg .ge. tlglo1 .and. tlg .le. tlghi1 )then
c..   lower T interpolation between OPAL and alexander
            kopl = 22
c<<<<<<< .mine
            nolglo = 10.0d0**rlglo
c=======
            
c>>>>>>> .r98
c..   evaluate at constant T (high and low) for given Ropal
            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, thigh1, nolglo )
            else
               call opacgn93( z, xh, thigh1, nolglo,k,l )
            endif
            
            opact2  = opact 
            dopacr2 = dopacr
            dopact2 = dopact
            
c..   evaluate at constant T (high and low) for given Ropal
            call kapp(tlglo1,rlg,ak,akr,akt,zalex,xh)
            
            opact1  = ak
            dopacr1 = akr
            dopact1 = akt
            

c..   interpolation variable in temperature
            xtlg  = (tlg - tlglo1)/(tlghi1 - tlglo1)

c..   linear in log T, at fixed ropal
            fact = 1.0d0 - xtlg
            ak   =  opact2 *xtlg +  opact1*fact
            akt  =  dopact2*xtlg + dopact1*fact 
            akr  =  dopacr2*xtlg + dopacr1*fact


         elseif( tlg .ge. tlglo2 .and. tlg .le. tlghi2 )then
c..   higher T interpolation, between klein-nishina and OPAL
            kopl = 24

            theta   = 2.2d-3*t6
            opact2  = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
     .           /(1.0d0 + theta) )
            dopacr2 = 0.0d0

            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, tlow2, nolglo )
            else
               call opacgn93( z, xh, tlow2, nolglo,k,l )
            endif
            opact1  = opact
            dopacr1 = dopacr
            
            dlgt   = log10(thigh2/tlow2)
            fact   = log10(t6/tlow2)/dlgt
            ak  =  opact1*(1.0d0 - fact) + opact2*fact
            akt = (opact2-opact1)/dlgt
            akr =  dopacr1*(1.0d0 - fact) + dopacr2*fact


         elseif( tlg .gt. tlghi2 )then
c..   high T, klein-nishina
            kopl = 25
            
            theta = 2.2d-3*t6
            ak = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
     .           /(1.0d0 + theta) )
            akt = - theta/(1.0d0 + theta)
            akr = 0.0d0         
            
         else
c..   inside OPAL table boundaries in T6, but Ropal is too small
            kopl = 23
            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, t6, nolglo )
            else
               call opacgn93( z, xh, t6, nolglo,k,l )
            endif
c..opact is a function of T only, as ropal is fixed
            ak  = opact
            akt = dopact
            akr = dopacr

         endif
         
         if( ak .gt. 20.0d0 )then
c..   excessively large opacity
            ak = 20.0d0
            akt = 0.0d0
            akr = 0.0d0
         endif

c..   construct partials with V, not Ropal, constant
         akv =            -akr
         akt = akt - 3.0d0*akr
         
c..   reset to cgs units
         ak  = 10.0d0**ak
         akt = ak*akt/(t6*1.0d6)
         akv = ak*akv*(ropal*t6**3)
         
         return


      elseif( rlg .gt. rlghi )then
c..   low entropy, high Ropal
         
         if( tlg .lt. tlglo1 )then
c..   low T, use atmospheric opacities (alexander)
c..   variable kopl signals which logic branch was used (debugging)
            kopl = 11
            
            call kapp(tlg,rlghi,ak,akr,akt,zalex,xh)

c..   designed for OPAL96 boundary (ropcut = 10.0d0, log10 => +1)
            ratrlg = log10(ropal)-rlghi
            akcoef = 0.5d0
            if( ratrlg .gt. 0.0d0 )then
c..   ak is function of log T plus part linear in log10(ropal)
               ak  = ak  + akcoef * ratrlg
               akr =       akcoef
            else
               akr = 0.0d0
            endif
            
         elseif( tlg .ge. tlglo1 .and. tlg .le. tlghi1 )then
c..   lower T interpolation between OPAL and alexander
            kopl = 12
c<<<<<<< .mine
            nolghi = 10.0d0**rlghi
c=======
            
c>>>>>>> .r98
c..   evaluate at constant T (high and low) for given Ropal
            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, thigh1, nolghi )
            else
               call opacgn93( z, xh, thigh1, nolghi,k,l )
            endif
            opact2  = opact 
            dopacr2 = dopacr
            dopact2 = dopact

c..   evaluate at constant T (high and low) for given Ropal
            call kapp(tlglo1,rlghi,ak,akr,akt,zalex,xh)
            
            opact1  = ak
            dopacr1 = akr
            dopact1 = akt
            
c..   interpolation variable in temperature
            xtlg  = (tlg - tlglo1)/(tlghi1 - tlglo1) 
c..   linear in log T, at fixed ropal
            fact = 1.0d0 - xtlg
            ak   =  opact2 *xtlg +  opact1*fact
            akt  =  dopact2*xtlg + dopact1*fact 
            akr  =  dopacr2*xtlg + dopacr1*fact

c..   designed for OPAL96 boundary (ropcut = 10.0d0, log10 => +1)
            ratrlg = log10(ropal)-rlghi
            akcoef = 0.5d0
            if( ratrlg .gt. 0.0d0 )then
c..   join is linear in log10(ropal)
c..   extrapolates ak as ratrlg, approximately krammers
               ak  = ak  + akcoef * ratrlg
c     akt = akt, so unchanged
               akr = akcoef
c..   
            endif


            
         elseif( tlg .ge. tlglo2 .and. tlg .le. tlghi2 )then
c..   higher T interpolation, between klein-nishina and OPAL
            kopl = 14
            
            theta   = 2.2d-3*t6
            opact2  = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
     .           /(1.0d0 + theta) )
            dopacr2 = 0.0d0
            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, tlow2, nolghi )
            else
               call opacgn93( z, xh, tlow2, nolghi,k,l )
            endif
            opact1  = opact
            dopacr1 = dopacr
c..   designed for OPAL96 boundary (ropcut = 10.0d0, log10 => +1)
            ratrlg = log10(ropal)-rlghi
            akcoef = 0.5d0

            if( ratrlg .gt. 0.0d0 )then
c..   join is linear in log10(ropal)
c..   extrapolates ak as ratrlg, approximately krammers
               opact1  = opact  + akcoef * ratrlg
               dopacr1 = dopacr + akcoef
            endif

            dlgt   = log10(thigh2/tlow2)
            fact   = log10(t6/tlow2)/dlgt
            ak  =  opact1*(1.0d0 - fact) + opact2*fact
            akt = (opact2-opact1)/dlgt
            akr =  dopacr1*(1.0d0 - fact) + dopacr2*fact

         elseif( tlg .gt. tlghi2 )then
c..   high T, klein-nishina
            kopl = 15

            theta = 2.2d-3*t6
            ak = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
     .           /(1.0d0 + theta) )
            akt = - theta/(1.0d0 + theta)
            akr = 0.0d0

         else
            
c..   inside OPAL table boundaries in T6, but Ropal is too big
            kopl = 13

            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, t6, nolghi )
            else
               call opacgn93( z, xh, t6, nolghi,k,l )
            endif

c..   opact is a function of T only, as ropal is fixed
            ak  = opact
            akt = dopact
            akr = dopacr
c..   designed for OPAL96 boundary (ropcut = 10.0d0, log10 => +1)
            ratrlg = log10(ropal)-rlghi
            akcoef = 0.5d0
            
            if( ratrlg .gt. 0.0d0 )then
c..   ak is function of log T plus part linear in log10(ropal)
               ak  = ak  + akcoef * ratrlg
               akr =       akcoef
            else
               akr = 0.0d0
            endif
            
         endif

         
         if( ak .gt. 20.0d0 )then
c..   excessively large opacity
            ak = 20.0d0
            akt = 0.0d0
            akr = 0.0d0
         endif

c..   construct partials with V, not Ropal, constant
         akv =            -akr
         akt = akt - 3.0d0*akr
         
c..   reset to cgs units
         ak  = 10.0d0**ak
         akt = ak*akt/(t6*1.0d6)
         akv = ak*akv*(ropal*t6**3)

         return

      else
c..   medium entropy, Ropal in OPAL table



         if( tlg .lt. tlglo1 )then
c..   low T, use atmospheric opacities (alexander)
c..   variable kopl signals which logic branch was used (debugging)
            kopl = 1
            
            call kapp(tlg,rlg,ak,akr,akt,zalex,xh)
            
            
         elseif( tlg .ge. tlglo1 .and. tlg .le. tlghi1 )then
c..   lower T interpolation between OPAL and alexander
            kopl = 2


c..   evaluate at constant T (high and low) for given Ropal
            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, thigh1, ropal )
            else
               call opacgn93( z, xh, thigh1, ropal,k,l )
            endif
            opact2  = opact
            dopacr2 = dopacr
            dopact2 = dopact

c..   evaluate at constant T (high and low) for given Ropal
            call kapp(tlglo1,rlg,ak,akr,akt,zalex,xh)
            opact1  = ak
            dopacr1 = akr
            dopact1 = akt
c..   interpolation variable in temperature
            xtlg  = (tlg - tlglo1)/(tlghi1 - tlglo1) 
c..   linear in log T
            fact = 1.0d0 - xtlg
            onemf = 1.0d0 - fact
            ak   =  opact2 *onemf +  opact1*fact
            akt  =  dopact2*onemf + dopact1*fact 
            akr  =  dopacr2*onemf + dopacr1*fact


         elseif( tlg .ge. tlglo2 .and. tlg .le. tlghi2 )then
c..   higher T interpolation, between klein-nishina and OPAL
            kopl = 4

            theta   = 2.2d-3*t6
            opact2  = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
     1           /(1.0d0 + theta) )
            dopacr2 = 0.0d0
            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, tlow2, ropal )
            else
               call opacgn93( z, xh, tlow2, ropal,k,l )
            endif
            opact1  = opact
            dopacr1 = dopacr

            dlgt   = log10(thigh2/tlow2)
            fact   = log10(t6/tlow2)/dlgt
            ak  =  opact1*(1.0d0 - fact) + opact2*fact
            akt = (opact2-opact1)/dlgt
            akr =  dopacr1*(1.0d0 - fact) + dopacr2*fact

         elseif( tlg .gt. tlghi2 )then
c..   high T, klein-nishina
            kopl = 5

            theta = 2.2d-3*t6
            ak = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
     1           /(1.0d0 + theta) )
            akt = - theta/(1.0d0 + theta)
            akr = 0.0d0                     

         else
c..   inside OPAL table boundaries in T6, Ropal
            kopl = 3

c            if( l .ge. 39 .and. k .eq. 155 )then
c               write(*,'(a20,3i5,1p8e12.3)')'newopal 2',l,k,kopl,
c     1              z,xh,t6,ropal
c            endif
cccccccccccccccccccccccccccc


            if( nopac .eq. 0 )then
               call OPAC( z, xh, xxc, xxo, t6, ropal )
            else
               call opacgn93( z, xh, t6, ropal,k,l )
            endif

            ak  = opact
            akt = dopact
            akr = dopacr
         endif


         if( ak .gt. 20.0d0 )then
c..   excessively large opacity
            ak = 20.0d0
            akt = 0.0d0
            akr = 0.0d0
         endif

c..   construct partials with V, not Ropal, constant
         akv =            -akr
         akt = akt - 3.0d0*akr

c..   reset to cgs units
         ak  = 10.0d0**ak
         akt = ak*akt/(t6*1.0d6)
         akv = ak*akv*(ropal*t6**3)

         return
      endif


c..avoided for correct logic
      write(*,*)'newopal.f logic is in error!'
      stop'opal 0'


      end

