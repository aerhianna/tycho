
      subroutine kapp(ttlg,rlg,a,dadx,dady,z,xh)

c..   last revised 5-8-06 wda
c..   payoung 5-24-02; revised darnett 5-26-02
c..   Rosseland mean opacities from Alexander
c..   input
c..   ttlg  log10 T(K)
c..   lrop  log10 Ropal
c..   where Ropal = rho/t6**3
c..   reads unit 30 for opacity table derived from Alexander
c..   z nucleon fraction of metals
c..   xh nucleon fraction of hydrogen
c..   output
c..   a  log10 kappa
c..   dadx d log kappa/d log ropal
c..   dady d log kappa/d log T(K)


      implicit none

      include 'calex'

      integer*4 mm, nn
      parameter( mm = 85, nn=19 )
      real*8 aa(nn,mm), lropal(nn), tlg(mm),
     1       aa2(nn,mm), lropal2(nn), tlg2(mm),
     1       aa3(nn,mm), aa4(nn,mm), lropal3(nn), tlg3(mm),
     1       lropal4(nn), tlg4(mm)
      real*8 ttlg,ttlg0,lrop,rlg
      real*8 a,dadx,dady, a1h,a2h, dadx1h,dadx2h,dady1h,dady2h
      real*8 a1,dadx1,dady1,a2,dadx2,dady2,a3,dadx3,dady3,
     1       a4,dadx4, dady4
      real*8 z,xh

      real*8 ztest,htest

c.. Construct lookup tables for the compositions
c.. in the Alexander opacity tables.
c      integer*4 pp, qq
c      parameter( pp = 7, qq = 14)
c      real*8 xhtabs(pp), ztabs(qq)
c      data xhtabs(1)/0.0d0/, xhtabs(2)/0.1d0/, 
c     1          xhtabs(3)/0.2d0/, xhtabs(4)/0.35d0/,
c     1          xhtabs(5)/0.5d0/, xhtabs(6)/0.7d0/,
c     1          xhtabs(7)/0.8d0/
c      data ztabs(1)/0.1/, ztabs(2)/0.07d0/, 
c     1          ztabs(3)/0.04d0/,
c     1          ztabs(4)/0.03d0/, ztabs(5)/0.02d0/,
c     1          ztabs(6)/0.01d0/, ztabs(7)/0.004d0/,
c     1          ztabs(8)/0.002d0/, ztabs(9)/0.001d0/,
c     1          ztabs(10)/0.0003d0/, ztabs(11)/0.0001d0/,
c     1          ztabs(12)/0.00003d0/, ztabs(13)/0.00001d0/,
c     1          ztabs(14)/0.0d0/

      character*77 dummy
      character*68 ctabl

c..identifies tables as they are read in; debugging

c      integer*4 hsize, zsize
c      parameter(hsize = 434, zsize = 31)

      integer*4 i,j,jt,n,zindex,hindex,skip,zskip
      integer*4 iread, hindexold, zindexold
      logical tobe
      data iread/1/
      data hindexold/0/
      data zindexold/0/

c      integer*4 k

      save
c---------------------------------------------------------------------
c..   smallest log10 Ropal is -8.0 in OPAL, -8.0 here in alxdata
c..   limit for sanity, but do not pass back overwritten value
      lrop = dmax1( -8.0d0, rlg )
c..   save for scaling T at low Ropal
      ttlg0 = ttlg

c.. Find indices of composition tables for opal.f values
c.. indices limited to edges of tables
c.. variable xh lies between index hindex and hindex+1
      call locate(xhtabs,pp,xh,hindex)
c.. variable z lies between index zindex and zindex+1
      call locate(ztabs,qq,z,zindex)
      if(hindex .gt. pp)hindex=pp
      if(zindex .gt. qq)zindex=qq

c.. Checks to see if table for composition already
c.. in memory; avoids unnecessary reading of opacity file

      if(hindex .ne. hindexold)then
         iread = 1
      elseif(zindex .ne. zindexold)then
         iread = 1
      endif

      if(ztabs(zindex) .ge. z)then
c..   zskip = zindex-1
         if(zindex .lt. qq)then
           zskip = zindex+1
         else
           zskip = zindex
        endif
      else
c..   zskip = zindex-2
         write(*,'(2(a20,1pe26.16))')'ztabs(zindex)',ztabs(zindex),
     1        'z',z
         stop'newkapp.f: zskip error after locate'
      endif

      if( iread .eq. 1 )then
         write(*,*)'entering newkapp, iread=1'

         write(*,*)'xhtabs= H1 values in Alexander tables'
         write(*,'(1p10e10.2)')(xhtabs(j),j=1,pp)
         write(*,*)'ztabs= z values in Alexander tables'
         write(*,'(1p10e10.2)')(ztabs(j),j=1,qq)

         write(*,'(a20,1pe12.3)')'z in initial model',z
         write(*,'(a20,1pe12.3,i5)')'ztabs(zindex)',ztabs(zindex),
     1        zindex
         write(*,'(a20,1pe12.3,i5)')'ztabs(zskip)',ztabs(zskip),zskip
         write(*,'(a20,1pe12.3)')'xh',xh
         write(*,'(a20,1pe12.3,i5)')'xhtabs(hindex)',xhtabs(hindex),
     1        hindex
         write(*,'(a20,1pe12.3,i5)')'xhtabs(hindex+1)',xhtabs(hindex+1),
     1        hindex+1
      endif

      if( iread .ne. 0 )then
c..   read alexander data from opal-type table
         inquire(file='alxdata',exist=tobe)
         if( .not. tobe )then
            write(*,*)'file alxdata does not exist in this directory'
            stop'kapp.f'
         else
            open(30,file='alxdata')
            rewind(30)
         endif

c..   Skips to table for appropriate composition
c..   lowest XH and highest z are the lesser indices
         skip = (hindex-1)*hsize +(zindex-1)*zsize 

c         write(*,'(4(a20,i5))')'hsize',hsize,'zsize',zsize
c         write(*,'(4(a20,i5))')'skip',skip,'hindex-1',hindex-1,
c     1        'hi-1*hsize',(hindex-1)*hsize,'zskip*zsize',zskip*zsize
         
         if(skip .ne. 0)then
           do i = 1,skip
             read(30,'(a77)')dummy
c             write(*,*)i,dummy
           enddo
         endif
         
c..   Header
c         read(30,'(a77)')dummy
         read(30,'(a68)')ctabl


c..   notify that new table is being used
         write(*,'(2x,a68,a30,3i5)')ctabl,
     1        ' lropal; reading new alxdata',hindex,zskip,skip
   
         open(99,file='newkapp.dummy')
c..insure that dummy file begins at the beginning
         rewind 99

c..   verify that this is the correct table
         write(99,*)ctabl(61:68)
         rewind 99
         read(99,*)ztest
         rewind 99
         write(99,*)ctabl(45:52)
         rewind 99
         read(99,*)htest
         rewind 99
c..   clean up temporary file
         close(99)
         call system('rm newkapp.dummy')

         if( z .gt. ztabs(zindex) .or. z .lt. ztabs(zskip) )then
            write(*,'(1p8e12.3)')ztabs(zindex),z,ztabs(zskip)
            write(*,*)'ERROR in newkapp.f, z'
         endif
         if( xh .lt. xhtabs(hindex) .or. xh .gt. xhtabs(hindex+1) )then
            write(*,'(1p8e12.3)')xhtabs(hindex),xh,xhtabs(hindex+1)
            write(*,*)'ERROR in newkapp.f, xh'
         endif

         do i = 1, 2
            read(30,'(a77)')dummy
         enddo

c..   Ropal variable
         read(30,'(6x,0p19f7.3)')lropal
c         read(30,'(a77)')dummy
c         read(30,'(a77)')dummy

c..   TLG and log kappa(ropal,TLG)
         do jt = 1,mm
            read(30,'(0pf5.3,1x,0p19f7.3)')tlg(jt),(aa(n,jt),n=1,nn)
         enddo

c.. Get table from next higher metallicity to perform z
c.. interpolation
         
c..   Header
c         read(30,'(a77)')dummy
         read(30,'(a68)')ctabl
         write(*,'(2x,a68,a20)')ctabl,' lropal3'

         do i = 1, 2
            read(30,'(a77)')dummy
         enddo

c..   Ropal variable
         read(30,'(6x,0p19f7.3)')lropal3
c         read(30,'(a77)')dummy
c         read(30,'(a77)')dummy

c..   TLG and log kappa(ropal,TLG)
         do jt = 1,mm
            read(30,'(0pf5.3,1x,0p19f7.3)')tlg3(jt),(aa3(n,jt),n=1,nn)
         enddo
c.. Get table from next higher hydrogen fraction, same z, to
c.. perform interpolation
         if(hindex .lt. pp)then
         do i = 1,hsize-2*zsize
             read(30,'(a77)')dummy
         enddo
         endif
c         read(30,'(a77)')dummy
         read(30,'(a68)')ctabl
c..   notify that new table is being used
         write(*,'(2x,a68,a30)')ctabl,' lropal2; reading new alxdata'
c..
         do i = 1, 2
            read(30,'(a77)')dummy
         enddo
c..   Ropal variable
         read(30,'(6x,0p19f7.3)')lropal2
c         read(30,'(a77)')dummy
c         read(30,'(a77)')dummy
         
c..   TLG and log kappa(ropal,TLG)
         do jt = 1,mm
            read(30,'(0pf5.3,1x,0p19f7.3)')tlg2(jt),(aa2(n,jt),n=1,nn)
         enddo
         
c.. Get table from next higher metallicity to perform z
c.. interpolation

c..   Header
c         read(30,'(a77)')dummy
         read(30,'(a68)')ctabl
         write(*,'(2x,a68,a20)')ctabl,'lropal 4'
c..
         do i = 1, 2
            read(30,'(a77)')dummy
         enddo

c..   Ropal variable
         read(30,'(6x,0p19f7.3)')lropal4
c         read(30,'(a77)')dummy
c         read(30,'(a77)')dummy
         
c..   TLG and log kappa(ropal,TLG)
         do jt = 1,mm
            read(30,'(0pf5.3,1x,0p19f7.3)')tlg4(jt),(aa4(n,jt),n=1,nn)
         enddo

         close(30)
         iread = 0
         hindexold = hindex
         zindexold = zindex
      endif

         
c      endif
c.. find location in table of Ropal and logT. Only needs to be done
c.. once since tables are all one size.

      call locate(lropal,nn,lrop,i)

      if( i .lt. 1 )then
c..lower limit of ropal (-8 to 1 for alexander tables)
c..extrapolate behavior as shifted temperature for lropa<-8
         i=1
         lrop = lropal(1)
         ttlg0 = ttlg - 0.02d0*( rlg-lropal(1) )

         call locate(tlg,mm,ttlg0,j)

c..keep point in tables
         if( j .gt. mm-1 ) j =  mm-1
         if( j .le. 0    ) j =  1
      else

         call locate(tlg,mm,ttlg,j)

c..keep point in tables
         if( j .gt. mm-1 ) j =  mm-1
         if( j .le. 0    ) j =  1
      endif

      if( i .gt. nn-1)i=nn-1

c.. Interpolate within opacity tables to get a, dadx, dady for
c.. each composition
      call bilen(lropal ,tlg ,aa ,lrop,ttlg0,a1,
     1     dadx1,dady1,nn,mm,i,j)
      call bilen(lropal3,tlg3,aa3,lrop,ttlg0,a3,
     1     dadx3,dady3,nn,mm,i,j)
      call bilen(lropal2,tlg2,aa2,lrop,ttlg0,a2,
     1     dadx2,dady2,nn,mm,i,j)
      call bilen(lropal4,tlg4,aa4,lrop,ttlg0,a4,
     1     dadx4,dady4,nn,mm,i,j)

c.. Interpolate between a, dadx, dady for the two metallicities.
c..   Only linear since there are but two points.
     
      a1h = a3 - (a3-a1)*
     1     (z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1))
c      write(*,*)(z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1)),
c     1      z, ztabs(zindex+1),ztabs(zindex),zindex
      a2h = a4 - (a4-a2)*
     1     (z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1))
      dadx1h = dadx3 - (dadx3-dadx1)*
     1     (z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1))
      dadx2h = dadx4 - (dadx4-dadx2)*
     1     (z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1))
      dady1h = dady3 - (dady3-dady1)*
     1     (z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1))
      dady2h = dady4 - (dady4-dady2)*
     1     (z-ztabs(zindex+1))/(ztabs(zindex)-ztabs(zindex+1))

c.. Interpolate between a, dadx, dady for the two hydrogen 
c.. fractions. Only linear since there are but two points.

      if(hindex .lt. pp)then
      a = a1h + (a2h-a1h)*
     1     (xh-xhtabs(hindex))/(xhtabs(hindex+1)-xhtabs(hindex))
      dadx = dadx1h + (dadx2h-dadx1h)*
     1     (xh-xhtabs(hindex))/
     1     (xhtabs(hindex+1)-xhtabs(hindex))
      dady = dady1h + (dady2h-dady1h)*
     1     (xh-xhtabs(hindex))/
     1     (xhtabs(hindex+1)-xhtabs(hindex))
  
      else
         a = a1h
      dadx = dadx1h
      dady = dady1h
      endif
      return
      end


