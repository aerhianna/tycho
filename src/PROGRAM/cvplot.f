      program cvplot

c..   plots evolution of convective regions
c..   read from file called cv.'prefix'
c..   8-5-02

      implicit none

      include 'dimenfile'
      include 'cpgplot'
      include 'comod'

c..ncv = number of convection regions
c..nbs = number of burning shells
      integer*4 j, lw, loop,iloop, ncv, nbs

      real*4 xx
      real*4 xmin,xmax,ymin,ymax
      real*4 xmin0,xmax0,ymin0,ymax0
      real*4 rmin,rmax,xmmin,xmmax
      real*4 fsize,winsize

      character*80 toplab

c..   need some dummy values (time,dth,model,it)
c..   to limit x axis
      real*8    timed,dthd
      integer*4 inmodel, itd, kdum
c..   final values
      real*8    timef
      integer*4 modelf
c..   initial value
      real*8    timez
c..   log values
      real*8    tzlog, tflog

c..   working scalar
      real*8    fact

c..   y-axis variables
c..   cvrb ConVective Radius Begins
c..   cvre ConVective Radius Ends
c..   cvmb ConVective Mass Begins
c..   cvme ConVective Mass Ends
      real*4    cvrb(kdm),cvre(kdm),cvmb(kdm),cvme(kdm)
c..   burnr radius of burning shell
c..   burnm mass of burning shell
      real*4    burnr(2),burnm(2)

      integer*4 linesty
      integer*4 pgbeg

      character*5  cvfile
      character*5  device
      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt

      logical tobe

      data lw/3/,fsize/1.5/,winsize/4.0/
      data xmin/5.2/,xmax/3.3/,ymin/-1.0/,ymax/7.0/
      data rmin/7.0/,rmax/7.0/,xmmin/0.0/,xmmax/0.0/
      data tflog/0.0d0/,tzlog/0.0d0/
      data txtxt/'  '/
c-------------------------------------------------------
c..   read graphics options

      inquire(file='cvplot.in',exist=tobe)
      if( tobe )then
         open(2,file='cvplot.in')
      else
         write(*,*)'cvplot: no file cvplot.in in directory'
         stop'cvplot: no cvplot.in input file error'
      endif

c..   dummy reads for labelling header lines

      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text

      read (2,*)txt,labl,device
      write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,cvfile
      write(*,*)txt,txtxt,labl,cvfile
      labl1 = 'cvfile'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ixpg3
      write(*,*)txt,txtxt,labl,ixpg3
      labl1 = 'ixpg3'
      if( labl .ne. labl1 )goto 2000
      if( ixpg3 .lt. 0 .or. ixpg3 .gt. 3 )then
cccccccccccccccccccccc
         write(*,*)'ixpg3 error: ',ixpg3
         goto 2000
      endif

      read (2,*)txt,labl,iypg3
      write(*,*)txt,txtxt,labl,iypg3
      labl1 = 'iypg3'
      if( labl .ne. labl1 )goto 2000
      if( iypg3 .lt. 0 .or. iypg3 .gt. 1 )then
         write(*,*)'iypg3 error: ',iypg3
         goto 2000
      endif

      read (2,*)txt,labl,xmin0
      write(*,*)txt,txtxt,labl,xmin0
      labl1 = 'xmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmax0
      write(*,*)txt,txtxt,labl,xmax0
      labl1 = 'xmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymin0
      write(*,*)txt,txtxt,labl,ymin0
      labl1 = 'ymin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymax0
      write(*,*)txt,txtxt,labl,ymax0
      labl1 = 'ymax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,linesty
      write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

c..   dummy reads for labelling trailer lines
      read (2,*)text
      write(*,*)text

      close(2)

c..   hr.? data files
      inquire(file=cvfile,exist=tobe)
      if( tobe )then
         open(11,file=cvfile)
      else
         write(*,*)'WARNING: ',' no file ',
     1        cvfile,' in directory'
         write(*,*)'Please generate this file or adjust cvplot.in'
      endif

c..   get initial values from cv.?? file and reposition
c..   this gets starting values 

c..   inmodel = model number of entry
c..   itd     = number of iterations needed 
c..   ncv     = number of convection zones
c..   timed   = time elapsed in seconds
c..   dthd    = time step in second

      read(11,72,end=160)inmodel,nbs,itd,ncv,timed,dthd
      if( nbs .lt. 0 .or. nbs .gt. 2 )then
         write(*,*)
     1        'CVPLOT: error in nbs, number of burning shells'
         write(*,72)inmodel,nbs,itd,ncv,timed,dthd
         stop'cvplot.in error 1'
      endif

c..save initial time
      timez = timed
      rewind 11
      goto 170

 160  continue
c..   end of file encountered on first entry
      write(*,*)'end of file at beginning of ',cvfile
      stop'cvplot: eof'
 170  continue

c..   find last values, to define x axis
c..   get minimum and maximun of y axis
      write(*,*)'reading cv file at l =',l,' model =',model
      write(*,*)'to find the limits of x and y axes'
      rewind 11
      iloop = 1
 180  continue

c..   dummy reads to reposition at new "end" of file
      read(11,72,end=190) modelf, nbs, it, ncv, timef, dth(1)

      if( nbs .lt. 0 .or. nbs .gt. 2 )then
         write(*,*)
     1        'CVPLOT: error in nbs, number of burning shells'
         write(*,72)inmodel,nbs,itd,ncv,timed,dthd
         stop'cvplot.in error'
      elseif( nbs .ne. 0 )then
c..this sequence does have burning shell data
         do j = 1, nbs
            read(11,71,end=190)kdum,burnr(j),burnm(j)
         enddo
      endif

      if( ncv .ge. 1 )then
c..   the model does have some convection zones
c..   kdum  = dummy index, equal to j, hidden from do loop
c..   cv?b, cv?e are beginning and ending of convective zone
         do j = 1, ncv
            read(11,71,end=190)kdum,cvrb(j),cvre(j),cvmb(j),cvme(j)

            rmin  = amin1( rmin,cvrb(j))
            rmax  = amax1( rmax,cvre(j))
            xmmin = amin1(xmmin,cvmb(j))
            xmmax = amax1(xmmax,cvme(j))
         enddo
      endif      
      iloop = iloop + 1
      goto 180
 190  continue

      if( ixpg3 .eq. 0 )then
c..   time axis
         fact = timef - timez 
         xmin = timez - 0.05d0*fact
         xmax = timef + 0.05d0*fact
      elseif( ixpg3 .eq. 1)then
c..   model number (time step number) axis
         fact = modelf  - inmodel
         xmin = inmodel - int(0.05d0*fact)
         xmax = modelf   + int(0.05d0*fact)
      elseif( ixpg3 .eq. 2 .or. ixpg3 .eq. 3 )then
c..   time axis using time left 
c..   this expands the axis for later stages
         fact = timef - timez

         write(*,'(a20,1p3e18.8)')'t0,tf,tdif',timez,timef,fact
cccccccccc
         if( ixpg3 .eq. 2 )then
            if( fact .gt. 0.0d0 )then
               fact = dlog10( fact )
            else
               write(*,*)' timef - timez =',fact
               write(*,*)'cannot take log10'
               goto 2000
            endif
            if( timez .le. 0.0 )then
               tzlog = dlog10( dth(1) )
            else
               tzlog = log10( timez )
            endif
            if( timef .le. 0.0 )then
               write(*,*)' timef =',timef
               write(*,*)'cannot take log10'
               goto 2000
            else
               tflog = log10( timef )
            endif
            xmin = tflog + 0.05d0*fact
            xmax = tzlog - 0.05d0*fact
            xmax = 0.0
            write(*,'(a20,1p8e12.3)')'xmin,xmax,fact',xmin,xmax,fact
ccccccccccccc
         else
            xmin = -1.05d0*(timef-timez)
            xmax =  0.05d0*(timef-timez)
         endif

      else
         write(*,*)'cvplt: error ixpg3 ',ixpg3
         goto 2000
      endif

      if( iypg3 .eq. 0 )then
c..   mass coordinate in solar units on y-axis
         ymax = xmmax
         ymin = xmmin
      elseif( iypg3 .eq. 1)then
c..   log10 radius coordinate on y-axis
         ymin = rmin
         ymax = rmax
      else
         write(*,*)'cvplt: error iypg3 ',iypg3
         goto 2000
      endif
c..   add offset for legibility
      fact = ymax - ymin
      ymin = ymin - 0.05*fact
      ymax = ymax + 0.05*fact

c..   override limits for plot, from gen.f, which reads params.d
      if( xmin0 .ne. xmax0 )then
         xmin = xmin0
         xmax = xmax0
         write(*,*)'forcing x limits ',xmin,xmax
      endif
      if( ymin0 .ne. ymax0 )then
         ymin = ymin0
         ymax = ymax0
         write(*,*)'forcing y limits ',ymin,ymax
      endif

c..initialize pgplot graphics
      if ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'

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
c..   reset color (1=white)
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
      elseif( ixpg3 .eq. 1 )then
         call PGMTXT('B',2.0,0.5,0.5,'step number')
      elseif( ixpg3 .eq. 2 )then
          call PGMTXT('B',2.0,0.5,0.5,'log(t(f)-t)')
       else
          call PGMTXT('B',2.0,0.5,0.5,'t(final)-t')
      endif
c..   header label over graph
      toplab = 'Convection '//cvfile
      call PGMTXT('T',2.0,0.5,0.5,toplab)

c..   set color for data (15 = grey, 2=red, 3 = green )
      call pgsci(2)

c..   data handling for file cv.??   ...............................
c..   position to desired timestep to avoid overwrite
      write(*,*)'reading cv file at l =',l,' model =',model
      write(*,*)'to plot plot variables'
      rewind 11
      do loop = 1,iloop-1
c..   dummy reads to reposition at new "end" of file
c         read(11,70) model, it, ncv, time, dth(1)

c..   dummy reads to reposition at new "end" of file
         read(11,72) model, nbs, it, ncv, time, dth(1)
         if( nbs .lt. 0 .or. nbs .gt. 2 )then
            write(*,*)
     1           'CVPLOT: error in nbs, number of burning shells'
            write(*,72)inmodel,nbs,itd,ncv,timed,dthd
            stop'cvplot.in error'
         elseif( nbs .ne. 0 )then
c..   this sequence does have burning shell data
            do j = 1, nbs
               read(11,71)kdum,burnr(j),burnm(j)
            enddo
         endif

         if( ncv .ge. 1 )then
c..   the model does have some convection zones
c..   kdum  = dummy index, equal to j, hidden from do loop
c..   cv?b, cv?e are beginning and ending of convective zone
            do j = 1, ncv
               read(11,71)kdum,cvrb(j),cvre(j),cvmb(j),cvme(j)
c..   plot points for convection zones
               if( ixpg3 .eq. 0 )then
                  xx = time
               elseif( ixpg3 .eq. 1 )then
                  xx = model
               elseif( ixpg3 .eq. 2 )then

                  xx = timef - time
                  if( xx .gt. 0.0 )then
                     xx = log10( xx )
                  else
                     xx = log10( dth(1) )
                  endif
               else
                  xx = -timef + time
               endif

               if( iypg3 .eq. 0 )then
                  call pgmove(xx,cvmb(j))
                  call pgdraw(xx,cvme(j))
               else
                  call pgmove(xx,cvrb(j))
                  call pgdraw(xx,cvre(j))
               endif
            enddo
         endif

c..add (by overwriting) points (yellow=7,orange=8) for burning zone
         if( nbs .ne. 0 )then
            do j = 1, nbs
            call pgsci(7+j-1)
            if(  iypg3 .eq. 0 )then
               call pgpt1(xx,burnm(j),20 )
            else
               call pgpt1(xx,burnr(j),20 )
            endif
            enddo
            call pgsci(2)
         endif

      enddo

c..   pause if xwindow is used, else go on to clean-up
      if( device .eq. '/xwin' )then
         pause
      endif

c........................................................
      write(*,*)'normal termination'

      call PGEND
      stop

 70   format(3i6,1pe14.6,1pe12.4)
 71   format(i6,1p4e13.5)
 72   format(i6,i2,i4,i6,1pe14.6,1pe12.4)

 2000 continue

      write(*,*)'cvplot: error in  input file: cvplot.in'
      stop'cvplot error'

      end

