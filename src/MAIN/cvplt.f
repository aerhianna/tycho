      subroutine cvplt(no)

c..   plots evolution of convective regions
c..   saves new convective segments in file called cv.'prefix'
c..   12-28-02 (does not consider semiconvective zones convective)

      implicit none

      include 'dimenfile'
      include 'conline'
      include 'cpgplot'
      include 'cgen'
      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'cenv'
      include 'cconst'

      integer*4 j, lw, no, iloop, ncv, nbs, lastmodel, istat
      integer*4 pgopen

      real*4 xx
      real*4 xmin,xmax,ymin,ymax
      real*4 fsize,winsize

      integer*4 k

      character*40 toplab

c..   need some dummy values (time,dth,model,it) to avoid 
c..   disasterous overwrite! 
      real*8    timed,dthd
      integer*4 inmodel, itd, kdum

      real*8    fact
c..ncv is number of convective zones
c..last value of ncr is envelope

c..   cvrb ConVective Radius Begins
c..   cvre ConVective Radius Ends
c..   cvmb ConVective Mass Begins
c..   cvme ConVective Mass Ends
      real*4    cvrb(kdm),cvre(kdm),cvmb(kdm),cvme(kdm)
      real*4    burnr(2),burnm(2)

c..   burning variables
      integer*4 neps, ne
      parameter( neps = 2 )
      real*8    epmax(neps)
      real*8 tepshi(neps), tepslo(neps)
      integer*4 necol(neps)

      real*4    epxm(neps),eplgr(neps)
      integer*4 kepmax(neps)

      data lw/3/,fsize/1.5/,winsize/4.0/
      data xmin/5.2/,xmax/3.3/,ymin/-1.0/,ymax/7.0/
      data tepslo/ 0.0d0, 5.0d7/,
     1     tepshi/ 5.0d7, 1.0d10/
c-------------------------------------------------------
      lastmodel = 0
      if( l .le. 0 )then
         write(*,*)'cvplt: l error ',l
         stop'cvplt 0'
      elseif( l .eq. 1 )then

c..   get initial values from cv.?? file and reposition
c..   this gets starting values 
         read(11,72,end=160,iostat=istat)inmodel,nbs,itd,ncv,timed,dthd
         if( istat .ne. 0 )then
            write(*,*)'read error in cvplt.f'
            write(*,*)'istat ',istat
            write(*,*)'lastmodel ',lastmodel
            stop'cvplt.f, first'
         else
            lastmodel = inmodel
         endif
         rewind 11
         goto 170
 160     continue
c..   end of file encountered on first entry
c..   make current model the first one as no hr.?? exists
         timed   = time
         inmodel = model
 170     continue

         if( igraf .eq. 0 )then
c..   default igraf=0 gives online graphics
            if( ixpg3 .eq. 0 )then
c..   time axis
               fact = time - timed + dble( ll )*dth(2)
               xmin = timed                  - 0.05d0*fact
c..   extrapolate to end of this run assuming constant time steps
               xmax = time + dble(ll)*dth(2) + 0.05d0*fact
            elseif( ixpg3 .eq. 1)then
c..   model number (time step number) axis
               fact = model - inmodel + ll
               xmin = inmodel     - int(0.05d0*fact)
               xmax = model + ll  + int(0.05d0*fact)
            else
               write(*,*)'cvplt: error ixpg3 ',ixpg3
            endif


            if( iypg3 .eq. 0 )then
c..   mass coordinate in solar units on y-axis
               ymax = xm(kk+1)/sol + (xm(kk+1)-xm(1))/sol *0.05d0
               ymin = xm(1)   /sol - (xm(kk+1)-xm(1))/sol *0.05d0
            elseif( ixpg3 .eq. 1)then
c..   log10 radius coordinate on y-axis
               ymin = 8.0
               ymax = dlog10( r(1,kk) ) 
               ymax = ymax + (ymax-ymin)*0.05
            else
               write(*,*)'cvplt: error iypg3 ',iypg3
            endif

c..   override limits for plot, from gen.f, which reads params.d
            if( pg3xmin .ne. pg3xmax )then
               xmin = pg3xmin
               xmax = pg3xmax
            endif
            if( pg3ymin .ne. pg3ymax )then
               ymin = pg3ymin
               ymax = pg3ymax
            endif

            idpg3 = pgopen(device3)
            if( idpg3 .LE. 0 )stop' pgopen error in cvplt.f'

c..   no query for device
            call PGASK (.FALSE.)

c..   new page
            call PGPAGE
c..   scale window size (adjust pgscal for your preference;
c..   see data statement in pgplot.f)
            call PGPAP(winsize*pgscal,1.0)

c..   roman font = 2, sans serif = 1
            call PGSCF(1)
c..   font scaling (size)  (see data statement in pgplot.f)
            call PGSCH(fsize)
c..   window
            call PGSWIN(xmin,xmax,ymin,ymax)
c..   reset color
            call PGSCI(1)
c..   line width
            call PGSLW(lw)
c..   tickmarks on x(bottom=B,top=C) and y(left)
            call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
c..   Left label, dispacement outside in character heights,
c..   location in fractions of edge, justification(0.5=centered)
c..   and character string
c..   label axes
            if( iypg3 .eq. 0 )then
               call PGMTXT('L',2.0,0.5,0.5,'M/M(sol)')
            else
               call PGMTXT('L',2.0,0.5,0.5,'log R(cm)')
            endif
            if( ixpg3 .eq. 0 )then
               call PGMTXT('B',2.0,0.5,0.5,'time(sec)')
            else
               call PGMTXT('B',2.0,0.5,0.5,'step number')
            endif
c..   header label over graph
            toplab = 'Convection '//prefix(1)//prefix(2)
            call PGMTXT('T',2.0,0.5,0.5,toplab)

c..   set color for old data (15 = grey)
            call pgsci(15)

         endif

c..   data handling for file cv.??   ...............................
c..   position to desired timestep to avoid overwrite
         read(11,72,end=150,iostat=istat) inmodel,nbs,itd,ncv,timed,dthd
         write(3,*)'reading cv file at l =',l,' model =',model
         write(*,*)'reading cv file at l =',l,' model =',model, 
     1        ' inmodel =', inmodel
         if(model .gt. 1)then
c         if(inmodel .ge. model-1)then
c         if( iloop .eq. 1 )then
         rewind 11
c         endif
         if(inmodel .eq. model)then
         goto 150
         endif
         iloop = 1
         lastmodel = 0
 140     continue
c..   dummy reads to reposition at new "end" of file
         read(11,72,end=150,iostat=istat) inmodel,nbs,itd,ncv,timed,dthd
         if( istat .ne. 0 )then
            write(*,*)'second error read in cvplt.f'
            write(*,*)'istat ',istat
            write(*,*)'lastmodel ',lastmodel
            stop'second cvplt.f'
          endif
c         else
c            lastmodel = inmodel
c         endif
         if( nbs .eq. 1 .or. nbs .eq. 2 )then
c..   the sequence does write burning shell coordinates
            do j = 1, nbs
               read(11,71,end=150)kdum,burnr(j),burnm(j)
            enddo
         endif

         if( ncv .ge. 1 )then
c..   the model does have some convection zones
            do j = 1, ncv
               read(11,71,end=150)kdum,cvrb(j),cvre(j),cvmb(j),cvme(j)
               if( igraf .eq. 0 )then
                  if( j .lt. ncv )then
c..medium gray
                     call pgsci(15)
                  else
c..dark gray for envelope
                     call pgsci(14)
                  endif
c..   plot previously saved points
                  if( ixpg3 .eq. 0 )then
                     xx = timed
                  else
                     xx = inmodel
                  endif
                  if( iypg3 .eq. 0 )then
                     call pgmove(xx,cvmb(j))
                     call pgdraw(xx,cvme(j))
                  else
                     call pgmove(xx,cvrb(j))
                     call pgdraw(xx,cvre(j))
                  endif
               endif
            enddo
            if( nbs .eq. 1 .or. nbs .eq. 2 )then
c..   the sequence does write burning shell coordinates
c..   add (by overwriting) points (yellow=7,orange=8) for burning zone
               do j = 1, nbs
                  call pgsci(7+j-1)
                  if(  iypg3 .eq. 0 )then
                     call pgpt1(xx,burnm(j),1 )
                  else
                     call pgpt1(xx,burnr(j),1 )
                  endif
               enddo
c..color for old convective zones (15 = grey)
               call pgsci(15)
            endif

         endif
         
c..   overwrite later models
         if( inmodel .ge. model-1 )then
            if( iloop .eq. 1 )then
c..   first model in hr is larger that read by gen
c..   do not proceed further so subsequent write will overwrite 
c..   later entries in cv.?? file
               rewind 11
               write(3,*)' initializing cv'
            endif
            write(*,*)'cvplt: overwriting models beyond ',model-1
            goto 150
c            iloop = 1
         endif
         iloop = iloop + 1
         goto 140

 150     continue
         else
            write(*,*)'model = ',model,' so cv.?? read is skipped'
c..   insure starting at next model
            rewind 11
         endif

      endif

c..   add entry to cv.?? file for every converged time step
cccccccccccccccccccccccccccccccccccccccccc
c..   rewrite this (put before writes), modify format for burn zones 
c..   to find most active flame zone
c..   use maximum precision for elapsed time
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do ne = 1, neps
            epmax(ne) = 0.0d0
            kepmax(ne) = 0
            do k = 2, kk
               if( t(1,k) .lt. tepshi(ne) .and. 
     1              t(1,k) .ge. tepslo(ne) )then
                  if( s(5,k) .gt. epmax(ne) )then
                     kepmax(ne) = k
                     epmax(ne) = s(5,k)
                  endif
               endif
            enddo
            epxm(ne) = 0.0
            eplgr(ne) = 0.0d0
            if( iypg3 .eq. 0 )then
               epxm(ne) = xm( kepmax(ne) )/sol
            else
               eplgr(ne) = dlog10( r(1, kepmax(ne)) )
            endif
         enddo
         nbs = 0
         if( kepmax(1) .ge. 2 .and. kepmax(1) .le. kk )then
            burnm(1) = xm(kepmax(1))/sol
            burnr(1) = r(1,kepmax(1))/solrad
            nbs = nbs + 1
            if( kepmax(2) .ge. 2 .and. kepmax(2) .le. kk )then
               burnm(2) = xm(kepmax(2))/sol
               burnr(2) = r(1,kepmax(2))/solrad
               nbs = nbs + 1
            endif
         endif

c..now do convection
      ncv   = 0
      k     = 0
c..   boundary at origin for this logic
      ic(1) = ic(2)
      if( no .eq. 0 )then
 200     continue

         if( k .ge. kk-1 )then
c..   exit, do not draw envelope
            goto 201
         endif
         k = k+1

c        if( ic(k) .ne. 0 )then
         if( ic(k) .eq. 1 )then
c         if( h(k) .gt. 1.0d5 )then
ccccccccccccccc
c..   convective
c..   begin convective zone
            ncv  = ncv + 1
            if( r(1,k) .gt. 0.0d0 )then
               cvrb(ncv) =  dlog10( r(1,k) )
            else
               cvrb(ncv) = 8.0d0
            endif
            cvmb(ncv) =  xm(k)/sol

c..   look for end of convective zone
            do j = k+1, kk

               if( ic(j) .eq. 0 .or. ic(j) .eq. 2 )then
c..nonconvective or semiconvective
c..   outside convective zone, take next inner boundary
                  cvre(ncv) = dlog10( r(1,j-1) )
                  cvme(ncv) = xm(j-1)/sol
c..   update for zones already tested
                  k    = j - 1
c..   escape to look for beginning of next convective zone
c                 write(*,'(4i5,1p8e12.3)')ncv,k,ic(k),ic(j),
c    1   cvmb(ncv),cvme(ncv),h(k),h(j)
c                 write(*,'(3i5,1p8e12.3)')kbeg(1),kend(1),nrczones,
c    1  xm(kbeg(1))/sol,xm(kend(1))/sol
c             do k = j-20, j+1
c               write(*,'(2i5,1p8e12.3)')k,ic(k),dnabv(k),dnad(k),
c    1    doux(k),dnabv(k)-dnad(k),dnabv(k)-dnad(k)-doux(k),
c    2    h(k),x(nnuc,k),dmom(k)
c             enddo
c        stop'cvplt'
cccccccc

                  goto 202

               endif
            enddo
c..   at join boundary, still convective, escape search entirely
            cvre(ncv) = dlog10( r(1,kk) )
            cvme(ncv) = xm(kk)/sol
            goto 201
         endif
 202     continue

         if( k .lt. kk )goto 200
 201     continue
c..   henyey grid done, now add envelope
         if( modes .eq. 2 )then
c..   this treats whole envelope as convective
            ncv = ncv + 1
            cvrb(ncv) = dlog10( vr( nvmax(3),3) )
c..   set beginning slighlty outside join
            cvmb(ncv) = (xm(kk) + 0.5d0*dmh(kk))/sol
            cvre(ncv) = dlog10( vr( 1,3 ) )
            cvme(ncv) = xm(kk+1)/sol

         else
            ncv = ncv + 1
            cvrb(ncv) = dlog10( r(1,kk) )
            cvmb(ncv) = xm(kk)/sol
            cvre(ncv) = dlog10( r(1,kk) )
            cvme(ncv) = xm(kk)/sol
         endif

ccccccc

c..   always write this entry
         write(11,72)model, nbs,it, ncv, time, dth(2)

         if( nbs .gt. 0 )then
c..   burning regions have been identified
            do j = 1, nbs
               write(11,71)j,burnr(j),burnm(j)
            enddo
         endif
         if( ncv .gt. 0 )then
c..   if there are convective regions, write this entry
            do j = 1, ncv
               write(11,71)j,cvrb(j),cvre(j),cvmb(j),cvme(j)
            enddo
         endif
 70      format(3i6,1pe14.6,1pe12.4)
 71      format(i6,1p4e13.5)
 72      format(i6,i2,i4,i6,1pe14.6,1pe12.4)


         if( igraf .eq. 0 )then
c..   write new convective points on cvplt screen
            call pgslct(idpg3)
c..   set color for new data ( 2 = red)
c            call pgsci(2)
            if( ncv .gt. 0 )then
c..   some convection zones
               if( ixpg3 .eq. 0 )then
                  xx = time
               else
                  xx = real( model )
               endif
               do j = 1, ncv
c..new models: convection is red=2
                  if( j .lt. ncv )then
                     call pgsci(2)
                  else
c..envelope 12=deep red
                     call pgsci(13)
                  endif

                  if( iypg3 .eq. 0 )then
c..   mass is y coordinate
                     call pgmove(xx,cvmb(j))
                     call pgdraw(xx,cvme(j))
                  elseif( r(1,k) .gt. 0.0d0 )then
c..   log10 radius is y coordinate
                     call pgmove(xx,cvrb(j))
                     call pgdraw(xx,cvre(j))
                  endif
               enddo
            endif

c..   plot region of most energetic burning

            if( ixpg3 .eq. 0 )then
               xx = time
            else
               xx = real( model )
            endif

               if( epmax(1) .gt. epmax(2) )then
                  necol(1) = 3
                  necol(2) = 4
               else
                  necol(1) = 4
                  necol(2) = 3
               endif

            do ne = 1, neps
                  call pgsci( necol(ne) )
               if( iypg3 .eq. 0 )then
                  call pgpt1(xx,epxm(ne),1)
               else
                  call pgpt1(xx,eplgr(ne),1)
               endif
            enddo
         endif
      endif

      if( igraf .eq. 0 )then
c..   return control to main window
c..   no further ploting on this window until it is selected again
         call pgslct(idpg1)
      endif

      return
      end

