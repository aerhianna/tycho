      subroutine hrplt(no,ks)

c..   plots evolution in HR diagram
c..   saves new HR points in file called hr.'prefix'

      implicit none

      include 'dimenfile'
      include 'conline'
      include 'cpgplot'
      include 'cgen'
      include 'comod'
      include 'cphot'
      include 'cconst'
      include 'cenv'
      include 'csurface'
      include 'cadvec'

      integer*4 lw, pgopen
      integer*4 istat,lastmodel
c..   ks = index of phospheric boundary
      integer*4 ks

      real*4 xx,yy
      real*4 xmin,xmax,ymin,ymax
      real*4 ymaxp
      real*4 fsize,winsize

      integer*4 j, n, no, iloop
      real*8    xejdummy,rhoph,speryear,bodummy

      character*40 toplab

      real*8 xej(ndim)

      real*8 Tnot,lx,omegadot,Br_star,Mdot_star,rossby
      real*8 pnot,bo_omega,bo_omega2,tauc,taug,vconv
      real*8 xbosurfomeg

c..   need some dummy values (time,dth,model,it) to avoid disasterous 
c..   overwrite! 
      real*8    timed,dthd
      integer*4 inmodel, itd

      data lw/3/,fsize/1.5/,winsize/4.0/
      data xmin/5.2/,xmax/3.3/,ymin/-1.0/,ymax/7.0/
c-------------------------------------------------------
      lastmodel = 0
      if( l .le. 0 )then
         write(*,*)'hrplt: l error ',l
      elseif( l .eq. 1 )then

c..   get initial values from hr.?? file and reposition
c..   this gets starting values 
         read(77,70,end=160,IOSTAT=istat) inmodel, itd, timed, dthd,
     1           rphot, rhoph, uphot, xlol, telog
         if( istat .ne. 0 )then
            write(*,*)'read error in hrplt.f'
            write(*,*)'istat ',istat
            write(*,*)'lastmodel ',lastmodel
            stop'hrplt.f, first'
         else
            lastmodel = inmodel
         endif
         rewind 77
         goto 170
 160     continue
c..   end of file encountered on first entry
c..   make current model the first one as no hr.?? exists
         timed   = time
         inmodel = model
 170     continue

c.................................................................
         if( igraf .eq. 0 )then
c..   activate plot limits override from params.d
            if( gtmin .ne. gtmax )then
c..   reverse for astronomers convention
               xmax = gtmin
               xmin = gtmax
               ymin = glmin
               ymax = glmax
            endif

            idpg2 = pgopen(device2)
            if( idpg2 .LE. 0 )stop' pgopen error in hrplt.f'

c..   no query for device
            call PGASK (.FALSE.)

c..   new page
            call PGPAGE
            call PGPAP(winsize*pgscal,1.0)

c..   roman font = 2
            call PGSCF(1)
c..   font scaling (size)
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
            call PGMTXT('L',2.0,0.5,0.5,'log L/sol')
            call PGMTXT('B',2.0,0.5,0.5,'log Te(K)')
            ymaxp = ymax + 0.01*(ymax-ymin)
            call pgptxt(4.6,ymaxp,0.0,0.5,'O5')
            call pgptxt(4.447,ymaxp,0.0,0.5,'B0')
            call pgptxt(4.182,ymaxp,0.0,0.5,'B5')
            call pgptxt(3.996,ymaxp,0.0,0.5,'A0')
            call pgptxt(3.869,ymaxp,0.0,0.5,'F0')
            call pgptxt(3.780,ymaxp,0.0,0.5,'G0')
            call pgptxt(3.690,ymaxp,0.0,0.5,'K0')
            call pgptxt(3.542,ymaxp,0.0,0.5,'M0')
            call pgptxt(3.4,ymaxp,0.0,0.5,'RNS')

            toplab = 'Sequence '//prefix(1)//prefix(2)

            call PGMTXT('T',2.0,0.5,0.5,toplab)
c..   set color for previous HR evolution
            call pgsci(5)
         endif
c.......................................................................

c..   hr.?? data file handling...........................
c..   position file/hr at end of file to avoid overwrite
         write(*,*)'reading hr file at l =',l,' model =',model, 
     1     ' inmodel =',inmodel
c         if(inmodel .ge. model-1)then
         if(model .gt. 1)then
c         if( iloop .eq. 1 )then
            rewind 77
         if(inmodel .eq. model) then
         goto 150
         endif
            iloop = 1
            lastmodel = 0
 140        continue
c..   dummy reads to reposition at "end" of file
            read(77,70,end=150,IOSTAT=istat) inmodel, itd, timed, dthd,
     1           rphot, rhoph, uphot, xlol, telog
c            write(*,'(3i5,1p8e12.3)') iloop,inmodel, itd, timed, dthd,
c     1           rphot, rhoph, uphot, xlol, telog

            if( istat .ne. 0 )then
               write(*,*)'read error in hrplt.f '
               write(*,*)'istat ',istat,' lastmodel ',lastmodel
               write(*,*)'inmodel ',inmodel
               stop'hrplt.f'
            endif
c..save last successfully read model number
            lastmodel = inmodel
c..   move to end of lastmodel data
            do j = 1,ndim/10+1
               read(77,72,IOSTAT=istat) xejdummy
            enddo
            read(77,73,IOSTAT=istat) bodummy
 72         format(1p10e11.3)
            if( istat .ne. 0 )then
               write(*,*)'read error in data field '
               write(*,*)'istat ',istat,' lastmodel ',lastmodel
               stop'hrplt.f'
            endif
            if( igraf .eq. 0 )then
c..   write previously saved HR points
               xx = sngl( telog )
               yy = sngl( xlol  )
               call pgpt1(xx,yy,-1)
            endif
c..   overwrite later models
            if( inmodel .ge. model-1 )then
               if( iloop .eq. 1 )then
c..   first model in hr is larger that read by gen, overwrite hr
                  rewind 77
                  write(3,*)' initializing hr'
                  write(6,*)' initializing hr'
               endif
               write(*,*)'hrplt: overwriting models beyond ',model-1
               goto 150
            endif
            iloop = iloop + 1
            goto 140

 150        continue
         else
            write(*,*)'model = ',model,' so hr.?? read is skipped'
c..   insure starting at next model
            rewind 77
         endif

      endif

c..   writes time sequence every edit except startup

      if( no .eq. 0 )then
c..   hr.? log file for HR diagram plots
c..   only converged steps written

c..   integrated mass loss rate in solar masses
         speryear = abs( peryear*dth(2)/secpy )
c..   density at envelope boundary (photosphere)
         if( v(2,ks) .le. 0.0d0 )then
            rhoph = 0
         else
            rhoph = 1.0d0/v(2,ks)
         endif

c..   abundances in ejecta
         do n = 1, ndim
            xej(n) = x(n,kk)
         enddo

c..   model    corresponding model number
c..   it       iterations used
c..   time     seconds elapsed in whole sequence
c..   dth(2)   time step
c..   rphot    photospheric radius
c..   rhoph    photospheric density (g/cc)
c..   uphot    photospheric velocity (cm/s)
c..   xlol     log( L / Lsun )              log => base 10
c..   telog    log( Teff )                  ln  => base e
c..   xm(kk)   mass (gm) interior to outer fitting radius r(1,kk)
c..   peryear  mass loss(-) or gain(+) in solar masses per year
c..   speryear solar masses ejected in this time step
c..   omeg(kk) angular rotational velocity (1/s)
c..   xej      ejected abundances (nucleon fraction) array, size ndim

         open(59,file='blackowen_vals.txt')
         read(59,73) Tnot,lx,omegadot,Br_star,Mdot_star,rossby,
     1               pnot,bo_omega,bo_omega2,tauc,taug,vconv,
     2               xbosurfomeg
         close(59)

         write(77,70) model, it, time, dth(2),
     1        rphot, rhoph, uphot, xlol, telog, xmfinl, peryear,
     2        vinf,omeg(kk),dmh(kk+1),omeg(kk+1)
         write(77,71) xej
         write(77,73) Tnot,lx,omegadot,Br_star,Mdot_star,rossby,
     1                pnot,bo_omega,bo_omega2,tauc,taug,vconv,
     2                xbosurfomeg
      endif
 70   format(2i6,1pe14.6,1p8e12.4,1p4e11.3)
 71   format(1p10e11.3)
 73   format(1p13e11.3)

      if( igraf .eq. 0 )then

c.................................................................
c..   write new HR points on HR screen
         call pgslct(idpg2)

         xx = sngl( telog )
         yy = sngl( xlol  )
c..plots interpolation box
c         x2(1) = xx
c         x2(2) = xx + 0.25d0*dlog10(1.0d0 + epsill) 
c     1        - 0.5d0*dlog10( 1.0d0 + epsirr )
c         y2(1) = yy
c         y2(2) = yy
c         call pgsci(13)
c         call pgline(2,x2,y2)
c         x2(2) = xx
c         y2(2) = yy + dlog10(1.0d0 + epsill) 
c         call pgline(2,x2,y2)
ccccccccccccccccccccccccccccccccc
         call pgsci(7)
         xx = sngl( telog )
         yy = sngl( xlol  )
         call pgpt1(xx,yy,-1)

c..   return control to main window
         call pgslct(idpg1)
c.....................................................................
      endif


c      stop'end of hrplt'

      return
      end

