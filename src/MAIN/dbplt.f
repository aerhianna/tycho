

      subroutine dbplt(no)

c..   plots evolution of iterations

      implicit none

      include 'dimenfile'
      include 'conline'
      include 'cpgplot'
      include 'cgen'
      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'cconst'
      include 'cburn'
      include 'caeps'

      real*8 fact

      integer*4 j, lw, no
      integer*4 pgopen
      integer*4 nvar, itcol,nvarlim
      parameter( nvar=5 )

      real*4 xx(kdm),yy(kdm),x2(2),y2(2)
      real*4 xmin,xmax,ymin,ymax
      real*4 fsize,winsize

      real*8 tlmax,tlmin

      integer*4 k

      character*40 toplab

      character*5 chit,chkk,chmodel

      data lw/3/,fsize/1.5/,winsize/4.0/
      data xmin/-5.0/,xmax/300.0/,ymin/-1.0/,ymax/1.0/
c-------------------------------------------------------

      if( igraf .ne. 0 )return

c..
      tlmin = tl(2,2)
      tlmax = tl(2,2)
      do k = 2, kk
         tlmin = dmin1(tlmin,tl(2,k))
         tlmax = dmax1(tlmax,tl(2,k))
      enddo
c..symmetry if large negative luminosity
      tlmax = dmax1(tlmax,-tlmin)

      if( l .le. 0 )then
         write(*,*)'dbplt: l error ',l
         stop'dbplt 0'

      elseif( l .ge. 1 )then

         if( it .eq. 1 )then
c..   setup cycle, first time called
            if( ixpg4 .eq. 0 )then
               xmin = 0. - 0.05*float(kk)
               xmax = 1.05*float(kk)
            else
               write(*,*)'dbplt: error ixpg4 ',ixpg4
            endif
            if( iypg4 .ge. 0 .and. iypg4 .le. 9 )then
               ymax = 0.05
               ymin = - ymax
            else
               write(*,*)'dbplt: error iypg4 ',iypg4
            endif
c..   override limits for plot, from gen.f
            if( pg4xmin .ne. pg4xmax )then
               xmin = pg4xmin
               xmax = pg4xmax
            endif
            if( pg4ymin .ne. pg4ymax )then
               ymin = pg4ymin
               ymax = pg4ymax
            endif

            if( l .le. 1 )then
c..   open device, for first call
               idpg4 = pgopen(device4)
               if( idpg4 .LE. 0 )stop' pgopen error in dbplt.f'
c..   no query for device
               call PGASK (.FALSE.)
               call PGPAP(2.0*winsize*pgscal,0.8)
               call PGSUBP(1,4)

            else
c..   subsequent calls, select device
               call pgslct(idpg4)  
            endif

c..   new page
            call PGPAGE
c..   roman font = 2, sans serif = 1
            call PGSCF(1)
c..   font scaling (size)
            call PGSCH(2.0*fsize)

c..   loop over 4 panels, 1=top to 4=bottom

            do j = 1, 4
               call PGPANL(1,j)

c..   window limits
               if( j .eq. 1 )then
c..   top panel
c..   flag read in from params.d
                  if( iypg4 .eq. 0 )then
c..   window
                     call PGSWIN(xmin,xmax,ymin,ymax)
c..1111111111111111111111111111111111111111111111111
                  elseif( iypg4 .eq. 1 )then
c..   window for nablas
c                     call PGSWIN(xmin,xmax,-0.1,1.0)
                     call PGSWIN(xmin,xmax,-10.0,10.0)
c..keep consistent with call to pgswin below........
c..1111111111111111111111111111111111111111111111111
                  elseif( iypg4 .eq. 2 )then
c..   window for conv. velocity
                     call PGSWIN(xmin,xmax,-5.0,0.0)
c..1111111111111111111111111111111111111111111111111
                  elseif( iypg4 .eq. 3 )then
c..   window for energy generation
                     call PGSWIN(xmin,xmax,0.0,10.0)
c..1111111111111111111111111111111111111111111111111
                  elseif( iypg4 .eq. 4 )then
c..   window for H consumption
                     call PGSWIN(xmin,xmax,-5.0,0.0)
c..1111111111111111111111111111111111111111111111111
c..keep consistent with call to pgswin below........
                  elseif( iypg4 .eq. 5 )then
c..   window for He consumption
                     call PGSWIN(xmin,xmax,-2.0,0.0)
c                     call PGSWIN(xmin,xmax,-5.0,0.0)
c..1111111111111111111111111111111111111111111111111
c..keep consistent with call to pgswin below........
                  elseif( iypg4 .eq. 6 )then
                     call PGSWIN(xmin,xmax,-1.0,3.0)
c..1111111111111111111111111111111111111111111111111
c..keep consistent with call to pgswin below........
                  elseif( iypg4 .eq. 7 )then
                     call PGSWIN(xmin,xmax,-1.0,3.0)
c..1111111111111111111111111111111111111111111111111
c..keep consistent with call to pgswin below........
                  elseif( iypg4 .eq. 8 )then
                     call PGSWIN(xmin,xmax,-1.0,1.0)
c..1111111111111111111111111111111111111111111111111
c..keep consistent with call to pgswin below........
                 elseif( iypg4 .eq. 9 )then
c                    tlmin = tl(2,2)
c                    tlmax = tlmin
c                    do k = 2, kk
c                       tlmin = dmin1(tlmin,tl(2,k))
c                       tlmax = dmax1(tlmax,tl(2,k))
c                    enddo
c                    tlmax = dmax1(tlmax,-tlmin)
 
                     call PGSWIN(xmin,xmax,-1.0,1.0)
c..1111111111111111111111111111111111111111111111111
c..keep consistent with call to pgswin below........
                  endif

ccccccccccccccccccccccccccccccccc
               elseif( j .eq. 2 )then
c..   window for Luminosity
c                  call PGSWIN(xmin,xmax,-40.0,40.0)
                  call PGSWIN(xmin,xmax,ymin,ymax)
c..222222222222222222222222222222222222222222222222

               elseif( j .eq. 3 )then

c..keep consistent with call to pgswin below........
                  if( iypg4 .eq. 1 )then
c..convection windows
                     call PGSWIN(xmin,xmax,-5.0,0.0)
                  else
                     call PGSWIN(xmin,xmax,ymin,ymax)
                  endif
c..3333333333333333333333333333333333333333333333333
               else

c..keep consistent with call to pgswin below........
                  if( iypg4 .eq. 1 )then
c..convection windows
c                     call PGSWIN(xmin,xmax,-0.1,0.1)
                     call PGSWIN(xmin,xmax,-2.0,2.0)
c..this is the bottom panel
                  else
                     call PGSWIN(xmin,xmax,ymin,ymax)
                  endif
c..44444444444444444444444444444444444444444444444444
               endif

c..   reset color
               call PGSCI(1)
c..   line width
               call PGSLW(lw)
c..   tickmarks on x(bottom=B,top=C) and y(left)
               call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
c..   Left label, dispacement outside in character heights,
c..   location in fractions of edge, justification(0.5=centered)
c..   and character string
               if( j .eq. 1 )then
                  if( iypg4 .eq. 0 )then
                     call PGMTXT('L',2.0,0.5,0.5,'delta log R')
                  elseif( iypg4 .eq. 1 )then
                     call PGMTXT('L',2.0,0.5,0.5,'rad-doux-ad')
                  elseif( iypg4 .eq. 2 )then
                     call PGMTXT('L',2.0,0.5,0.5,'log10(vconv/sound)')
                  elseif( iypg4 .eq. 3 )then
                     call PGMTXT('L',2.0,0.5,0.5,'log10(eps)')
                  elseif( iypg4 .eq. 4 )then
                     call PGMTXT('L',2.0,0.5,0.5,'delta log H1')
                  elseif( iypg4 .eq. 5 )then
                     call PGMTXT('L',2.0,0.5,0.5,'delta log He4')
                  elseif( iypg4 .eq. 6 )then
                     call PGMTXT('L',2.0,0.5,0.5,'conv flag ic')
                  elseif( iypg4 .eq. 7 )then
                     call PGMTXT('L',2.0,0.5,0.5,'dmom/sol')
                  elseif( iypg4 .eq. 8 )then
                     call PGMTXT('L',2.0,0.5,0.5,'dif')
                  elseif( iypg4 .eq. 9 )then
                     call PGMTXT('L',2.0,0.5,0.5,'L')

                  endif
                  call PGMTXT('B',2.0,0.5,0.5,'zone number k')
               elseif( j .eq. 2 )then
                  call PGMTXT('L',2.0,0.5,0.5,'delta log L')
               elseif( j .eq. 3 )then
                  if( iypg4 .eq. 1 )then
                     call PGMTXT('L',2.0,0.5,0.5,'log10(vconv/sound)')
                  else
                     call PGMTXT('L',2.0,0.5,0.5,'delta log T')
                  endif
               elseif( j .eq. 4 )then
                  if( iypg4 .eq. 1 )then
                     call PGMTXT('L',2.0,0.5,0.5,'doux')
                  else
                     call PGMTXT('L',2.0,0.5,0.5,'d ln V')
                  endif
               endif               
            enddo
c..   end of panel loop

         endif
      endif

c..   all cycles
c..   write new db points on db screen

      do k = 1, kk
         xx(k) = k
      enddo
      yy(1) = 0.

c..   set color
      if( 2+it .le. 15 )then
         itcol = 2+it
      else
         itcol = it-15
      endif
      call pgsci(itcol)

c..   limiting condition
      fact    = 0.0d0
      nvarlim = 0
      do j = 1, 4
         if( abs( wtest(j) ) .gt. abs(fact) )then
            fact = wtest(j)
            nvarlim = j
         endif
      enddo
      if( nvarlim .le. 0 .or. nvarlim .gt. 4 )then
         write(*,*)nvarlim
         stop'dbplt: nvarlim'
      endif

c..   loop over 4 panels
      do j = 1, 4 

         call pgpanl(1,j)
         

         if( j .eq. 1 )then
c..top panel
            call modflg(it,chit)
            call pgtext(xmin+5*it,ymax*1.05,chit)

            if( iypg4 .eq. 0 )then
               do k = 2, kk
                  yy(k) = dr(k)/r(1,k)
               enddo
               call PGSWIN(xmin,xmax,ymin,ymax)
c..1111111111111111111111111111111111111111111111111111111111111
            elseif( iypg4 .eq. 1 )then
               do k = 2, kk
c                  yy(k) = dnrad(k)/dnad(k)-1.0d0
                  yy(k) = dnrad(k)-doux(k)-dnad(k)
               enddo                
               yy(1) = yy(2)
c...............................................................
               call PGSWIN(xmin,xmax,-10.0,10.0)
c               call PGSWIN(xmin,xmax,-1.0,1.0)
c..keep consistent with  call PGSWIN in initialization above....
c..1111111111111111111111111111111111111111111111111111111111111
            elseif( iypg4 .eq. 2 )then
c..uses hp(k) which changes in cmix.f
               do k = 2, kk
                  if( hp(k) .gt. 1.0d-8 )then
                     yy(k) = dlog10( hp(k)/sound(k) )
                  else
                     yy(k) = -8.
                  endif
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-5.0,0.0)
c..1111111111111111111111111111111111111111111111111111111111111
               enddo

            elseif( iypg4 .eq. 3 )then
               do k = 2, kk
c                  if( s(5,k) .gt. 1.0d-8 )then
                  if( ss(k) .gt. 1.0d-8 )then
                     yy(k) = dlog10( ss(k) )
                  else
                     yy(k) = -8.
                  endif
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,0.0,10.0)
c..1111111111111111111111111111111111111111111111111111111111111
               enddo
            elseif( iypg4 .eq. 4 )then
               do k = 2, kk
                  if( aex(ndim-2,k) .gt. 1.0d-8 )then
                     yy(k) = dlog10( aex(ndim-2,k) )
                  else
                     yy(k) = -8.
                  endif
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-5.0,0.0)
c..keep consistent with call to pgswin above........
c..1111111111111111111111111111111111111111111111111111111111111
               enddo
           elseif( iypg4 .eq. 5 )then
               do k = 2, kk
                  if( aex(ndim-1,k) .gt. 1.0d-8 )then
                     yy(k) = dlog10( aex(ndim-1,k) )
                  else
                     yy(k) = -8.
                  endif
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-2.0,0.0)
c                  call PGSWIN(xmin,xmax,-5.0,0.0)
c..keep consistent with call to pgswin above........
c..1111111111111111111111111111111111111111111111111111111111111
               enddo
           elseif( iypg4 .eq. 6 )then
               do k = 2, kk
                  yy(k) = ic(k)
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-1.0,3.0)
c..keep consistent with call to pgswin above........
c..1111111111111111111111111111111111111111111111111111111111111
               enddo
           elseif( iypg4 .eq. 7 )then
               do k = 2, kk
                  yy(k) = dmom(k)/sol
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-1.0,3.0)
c..keep consistent with call to pgswin above........
c..1111111111111111111111111111111111111111111111111111111111111
               enddo

           elseif( iypg4 .eq. 8 )then
               do k = 2, kk
                  yy(k) = dif(k)
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-1.0,1.0)
c..keep consistent with call to pgswin above........
c..1111111111111111111111111111111111111111111111111111111111111
               enddo

           elseif( iypg4 .eq. 9 )then
               do k = 2, kk
                  yy(k) = dtl(k)/tlmax
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-1.0,1.0)
c..keep consistent with call to pgswin above........
c..1111111111111111111111111111111111111111111111111111111111111
               enddo

            endif

c..   top label
            toplab = 'Debug '//prefix(1)//prefix(2)
c..   white
            call pgsci(1)
            call PGMTXT('T',0.5,0.5,0.5,toplab)
c..reset color to denote iteration number
            call pgsci(itcol)

         elseif( j .eq. 2)then
            do k = 2, kk
               yy(k) = dtl(k)/tl(1,k)
            enddo
            yy(1) = 0.0d0
c            call PGSWIN(xmin,xmax,-40.0,40.0)
            call PGSWIN(xmin,xmax,ymin,ymax)
c..222222222222222222222222222222222222222222222222222222222222
         elseif( j .eq. 3)then
            if( iypg4 .eq. 1 )then
c..uses hp(k) which changes in cmix.f
               do k = 2, kk
                  if( hp(k) .gt. 1.0d-8 )then
                     yy(k) = dlog10( hp(k)/sound(k) )
                  else
                     yy(k) = -8.
                  endif
                  yy(1)=yy(2)
                  call PGSWIN(xmin,xmax,-5.0,0.0)
c..333333333333333333333333333333333333333333333333333333333333
               enddo
            else
               do k = 2, kk
                  yy(k) = dt(2,k)/t(1,k)
               enddo
               call PGSWIN(xmin,xmax,ymin,ymax)
c..444444444444444444444444444444444444444444444444444444444444
            endif
         elseif( j .eq. 4)then
            if( iypg4 .eq. 1 )then
c..convective variables
               do k = 2, kk
                  yy(k) = doux(k)
               enddo                
               yy(1) = yy(2)
c               call PGSWIN(xmin,xmax,ymin,ymax)
               call PGSWIN(xmin,xmax,-2.0,2.0)
c..keep consistent with  call PGSWIN in initialization above....
            else
               do k = 2, kk
                  yy(k) = dv(2,k)/v(1,k)
               enddo
               call PGSWIN(xmin,xmax,ymin,ymax)
            endif
         endif

         if( ymax .gt. 0.0 .and. ymin .lt. 0.0 )then
c..   horizontal line to aid the eye
            call pgsci(15)
            x2(1) = xmin
            x2(2) = xmax
            y2(1) = 0.
            y2(2) = 0.
            call pgline(2,x2,y2)
            call pgsci(itcol)
         endif

         call pgtext(float(nwtest(j)),ymax*0.8,cwtest(j))

         call pgline(kk,xx,yy)

         if( j .eq. 1 )then
            if ( iypg4 .eq. 1)then
c..denote which zones are convective
               do k = 2, kk
                  if( ic(k) .eq. 1 )then
                     call pgsci(2)
                     call PGPT(1,xx(k),0.5,ichar('1'))
                     call pgsci(itcol)
                  elseif( ic(k) .eq. 2 )then
                     call pgsci(3)
                     call PGPT(1,xx(k),0.5,ichar('2'))
                     call pgsci(itcol)
                  endif
               enddo
            elseif( iypg4 .eq. 2 )then
c..   denote which zones are convective
               do k = 2, kk
                  if( ic(k) .eq. 1 )then
                     call PGPT(1,xx(k),-0.5,ichar('1'))
                  elseif( ic(k) .eq. 2 )then
                     call PGPT(1,xx(k),-0.5,ichar('2'))
                  endif
               enddo
            endif
         endif

         if( nvarlim .eq. j )then
c..   vertical line to mark worst zone and variable
            call pgsci(2)
            x2(1) = nwtest(nvarlim)
            x2(2) = x2(1)
            y2(1) = ymin
            y2(2) = ymax
            call pgline(2,x2,y2)
            call pgsci(itcol)
            call pgtext(x2(1),y2(2)*0.6,chit)
         endif

      enddo
c..   end of panel loop

c..   reset color
      call pgsci(1)

      call modflg(kk,chkk)
      call modflg(l,chmodel)

      call pgtext(xmax*0.9,ymax*1.05,chkk)
      call pgtext((xmax+xmin)*0.5,ymax*1.05,chmodel)

      return
      end

