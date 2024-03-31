c     
c     This subroutine calculates abundance anomolies:
c     (Relative to solar???)
c     
      subroutine anomaly

      implicit none
      
c..   input
c..   zpop = metallicity for OPAL
c..   x    = abundance
c..   solarx = abundance in solar system
c..   dmh = zone mass
c..   xa = atomic number A = Z+N
c..   cnuc = character symbol of nucleus
c..   maximum of interior zone loop
c..output
c..   none, only writes to *

      include 'dimenfile'
      include 'cburn'
      include 'cgen'
      include 'comod'
      include 'compu'

      real*8 xsur(ndim), xxsur(ndim)

      real*8 zfact, zsurf
      integer*4 k,i,nsur

      character*5 chsur(ndim)
c---------------------------------------------------------------

      if( zpop .gt. 0.0d0 )then
c..z=0.0189 in full anders-grevesse table
c..z=0.0148 in full lodders table
c         zfact = zpop / 0.0191d0
c         zfact = zpop / 0.0148d0
c         zfact = 0.0d0
c         do i = 1, nnuc
c            if( lz(i) .ne. 1 .and. lz(i) .ne. 2 )then
c               zfact = zfact + solarx(i)
c            endif
c         enddo
         zfact = zsol
         zfact = dmin1(zfact,1.0d0)
         if( zfact .gt. 0.0d0 )then
            zfact = zpop/zfact
         else
            zfact = 1.0d0
         endif
         
c         if(  1.0d0 - solarx(nnuc) - solarx(nnuc-1) .gt. 0.0d0 )then
c            zfact = zpop / ( 1.0d0 - solarx(nnuc) - solarx(nnuc-1) )
c         else
c            zfact = 1.0d0
c         endif

      elseif( zpop .eq. 0.0d0 )then
         zfact = 1.0d0
         write(*,*)'BIG BANG composition, using solar Fe as basis'
      else
         write(*,*)'anomaly.f: zpop in error, ',zpop
      endif

      write(*,'(a26,0pf9.4,a10,1pe11.3,a10,0pf6.4,a10,0pf9.4)')
     1     'ANOMALY: metallicity(OPAL)',zpop,', which is ',zfact,
     2     ' of solar(',zsol,'), [Fe/H]=',dlog10( zfact )

      write(*,'(a25,f10.5)')'initial Z/H factor',zpop0/zhyd0

c..   find abundance anomalies in envelope
      if( modes .eq. 2 )then
         k = kk+1
      else
         k = kk
      endif
      write(*,'(a30,a8,3(a10,a8))')'Solar Surface fraction: H1',
     1     '0.7392'     , ' He4', '0.2486','Z','0.0122',
     2     'Z/X ','0.0165'
      zsurf = 0.0d0
      do i = 1, nnuc
         xxsur(i) = x(i,k)
         if( lz(i) .gt. 2 )then
            zsurf = zsurf + x(i,k)
         endif
      enddo

      write(*,'(a30,0pf8.4,3(a10,0pf8.4))')'Model Surface fraction: H1',
     1     x(nnuc-1,k), ' He4', x(nnuc,k),
     2     'Z surf',zsurf,'Z/H surf',zsurf/(x(nnuc-1,k))

c      nsur = 0
c      do i = 1, nnuc
c         if( solarx(i) .gt. 1.0d-12 .and. xxsur(i) .gt. 0.0d0)then
c            if( abs(  dlog10( xxsur(i)/(zfact*solarx(i) )) ) .gt. 
c     1           0.05d0 )then
c               nsur = nsur + 1
c               chsur(nsur) = cnuc(i)
c               xsur(nsur)  = dlog10( xxsur(i)/(zfact*solarx(i)) ) 
c            endif
c         endif
c      enddo
c      if( nsur .gt. 0 )then
c         write(*,*)nsur,' Anomalies (log10) at photosphere, [?/Fe]'
c         write(*,'(5(a9,0pf8.3))')(chsur(i),xsur(i),i=1,nsur)
c      endif

      nsur = 0
      do i = 1, nnuc
         if( solarx(i) .gt. 1.0d-12 .and. xxin(i) .gt. 0.0d0)then
            if( abs(  dlog10( xxin(i)/(solarx(i) )) ) .gt. 
     1           0.05d0 )then
               nsur = nsur + 1
               chsur(nsur) = cnuc(i)
               xsur(nsur)  = dlog10( xxin(i)/(solarx(i)) ) 
            endif
         endif
      enddo
      if( nsur .gt. 0 )then
         write(*,*)nsur,' Initial abundances differ from solar',
     1        ' log[xxin/solarx]'
         write(*,'(5(a9,0pf8.3))')(chsur(i),xsur(i),i=1,nsur)
      endif         

      nsur = 0
      do i = 1, nnuc
         if( xxin(i) .gt. 1.0d-12 .and. xxsur(i) .gt. 0.0d0)then
            if( abs(  dlog10( xxsur(i)/(xxin(i) )) ) .gt. 
     1           0.05d0 )then
               nsur = nsur + 1
               chsur(nsur) = cnuc(i)
               xsur(nsur)  = dlog10( xxsur(i)/(xxin(i)) ) 

            endif
         endif
      enddo
      if( nsur .gt. 0 )then
         write(*,*)'Changes from initial abundances'
         write(*,*)nsur,' Anomalies (log10) at photosphere, [?/Fe]'
         write(*,'(5(a9,0pf8.3))')(chsur(i),xsur(i),i=1,nsur)
      endif

c..   find abundance anomalies for interior
c..   weight each zone by its fraction of total mass of interior
      do i = 1, nnuc
         xxsur(i) = x(i,2)*dmh(2)/(xm(kk)-xm(2))
         do k = 3,kk
            xxsur(i) = xxsur(i) + x(i,k)*dmh(k)/(xm(kk)-xm(2))
         enddo
      enddo
      nsur = 0
      do i = 1, nnuc
         if( solarx(i) .gt. 1.0d-12 .and. xxsur(i) .gt. 0.0d0 )then
            if( abs(  dlog10( xxsur(i)/xxin(i) ) ) .gt. 
     1           0.05d0 .and. xxin(i) .gt. 0.0d0 )then
               nsur = nsur + 1
               chsur(nsur) = cnuc(i)
               xsur(nsur)  = dlog10( xxsur(i)/(xxin(i)) ) 
            endif
         endif
      enddo
      if( nsur .gt. 0 )then
         write(*,*)nsur,' Anomalies (log10) in interior, X/Xinitial'
         write(*,'(5(a9,0pf8.3))')(chsur(i),xsur(i),i=1,nsur)
      endif

      return
      end

