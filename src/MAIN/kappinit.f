c
c
c     
      
      subroutine kappinit


c     Constructs lookup tables for 
c     Alexander & Ferguson opacities
c     
c     Tests for consistency with actual 
c     alxdata in local directory xhtabs 
c     and ztabs are tablular values 
c     from file alxtab check to see if 
c     they correspond to headers in alxdata
c     

c     -----------------------------
c     REVISION HISTORY:
c     5-16-05
c     9-29-08 (C.Meakin, comments)
c     -----------------------------
      

      implicit none

      include 'calex'
  
      real*8 ztest,htest

      integer*4 i, initial, ih, iz, numerr
      integer*4 hindex, zindex, zskip, skip
      character*77 dummy
      character*68 ctabl


      logical tobe
      data initial/1/,numerr/0/
      
c     

      inquire(file='alxtab',exist=tobe)
      if( .not. tobe )then
         write(*,*)'file alxtab does not exist in this directory'
         stop'kappinit.f'
      else
         open(31,file='alxtab')
      endif

      read(31,'(2x,i4,7x)')pp
      do i = 1, pp
         read(31,'(f7.5)')xhtabs(i)
      enddo

      read(31,'(2x,i4,7x)')qq
      do i = 1, qq
         read(31,'(f7.5)')ztabs(i)
      enddo

      hsize = zsize * qq
    


      if( initial .eq. 1 )then
         
c..   test alxdata and alxtab consistency on initial call
         write(*,'(a20,2i5)')'testing kappinit',initial
c..   set flag to avoid repeating this code

         initial = 0
         write(*,'(a10,i5)')'xhtabs',pp
         write(*,'(1p10e10.2)')(xhtabs(i),i=1,pp)
         write(*,'(a10,i5)')'ztabs',qq
         write(*,'(1p10e10.2)')(ztabs(i),i=1,qq)

c..   read alexander data from opal-type table
         inquire(file='alxdata',exist=tobe)
         if( .not. tobe )then
            write(*,*)'file alxdata does not exist in this directory'
            stop'kapp.f'
         else
            write(*,*)'kappinit: opening alxdata'
            open(30,file='alxdata')
         endif
         

c..   step through H tables
         numerr = 0
         do ih = 1, pp
            if( ih .eq. 1 )then
               skip = 0
            else
               skip = hsize-1
            endif
            hindex = ih
            zindex = 1
            zskip  = zindex + 1
c..   Steps to table for appropriate composition
            if(skip .ne. 0)then
               do i = 1,skip
                  read(30,'(a77)')dummy
               enddo
            endif
c..   Header
c            read(30,'(a77)')dummy
            read(30,'(a68)')ctabl
c..   notify that a new table is being used
            write(*,'(2x,a68,a30,3i5)')ctabl,
     1           ' lropal; reading new alxdata',hindex,zskip,skip
c..   verify that this is the correct table
            open(99,file='cinit.dummy')
            write(99,'(a7)')ctabl(61:68)
            rewind 99
            read(99,'(f7.5)')ztest
            rewind 99
            write(99,'(a7)')ctabl(45:52)
            rewind 99
            read(99,'(f7.5)')htest
            rewind 99
            if( ztest .ne. ztabs(zindex)
     1           .or. htest .ne. xhtabs(hindex) )then
               numerr = numerr+1
               write(*,'(i5,1p8e12.3)')ih,ztest,ztabs(zindex), 
     1              htest,xhtabs(hindex)
               write(*,*)'ERROR in newkapp.f'
            endif
c            write(*,'(i5,4(a10,1pe10.2))')ih,'ztest',ztest,
c     1           'z error',ztest-ztabs(zindex),
c     2           'htest',htest,'xh error',htest-xhtabs(hindex)
         enddo
         if( numerr .ne. 0 )then
            write(*,*)pp,' H entries checked, ',
     1           numerr,' errors found'
            write(*,*)'alexander opacity tables:'
            write(*,*)'files alxdata and alxtab are not consistent'
            stop'kappinit: ERROR in data'
         else
            write(*,'(i5,a40)')ih-1,
     1           'H entries checked, 0 errors found'
         endif


c..   step through z tables
         rewind 30
         numerr = 0
         do iz = 1, qq
            if( iz .eq. 1 )then
               skip = 0
            else
               skip = zsize-1
            endif
            hindex = 1
            zindex = iz
            zskip  = zindex + 1
c..   Steps to table for appropriate composition
            if(skip .ne. 0)then
               do i = 1,skip
                  read(30,'(a77)')dummy
               enddo
            endif
c..   Header
c            read(30,'(a77)')dummy
            read(30,'(a68)')ctabl

c..   notify that new table is being used
c            write(*,'(2x,a60,a30,3i5)')ctabl,
c     1           ' lropal; reading new alxdata',hindex,zskip,skip
c..   verify that this is the correct table
            open(99,file='cinit.dummy')
            write(99,'(a7)')ctabl(61:68)
            rewind 99
            read(99,'(f7.5)')ztest

            rewind 99
            write(99,'(a7)')ctabl(45:52)
            rewind 99
            read(99,'(f7.5)')htest
            rewind 99
            if( ztest .ne. ztabs(zindex)
     1           .or. htest .ne. xhtabs(hindex) )then
               numerr = numerr+1
               write(*,'(i5,1p8e12.3)')ih,ztest,ztabs(zindex), 
     1              htest,xhtabs(hindex)
               write(*,*)'ERROR in newkapp.f'
            endif
c            write(*,'(i5,4(a10,1pe10.2))')iz,'ztest',ztest,
c     1           'z error',ztest-ztabs(zindex),
c     2           'htest',htest,'xh error',htest-xhtabs(hindex)
         enddo

         if( numerr .ne. 0 )then
            write(*,*)iz-1,' z entries checked, ',
     1           numerr,' errors found'
            write(*,*)'alexander opacity tables:'
            write(*,*)'files alxdata and alxtab are not consistent'
            stop'kappinit: ERROR in data'
         else
            write(*,'(i5,a40)')iz-1,
     1           'z entries checked, 0 errors found'
         endif

c..   clean up
         close(99)
         call system('rm cinit.dummy')


c..reset to avoid repeat of tests for subsequent calls
         initial = 0

      else
         write(*,*)'kappinit called more than once'
         stop'kappinit recall'
      endif
      rewind 30

      close(31)  


c     SUCCESS
      return
      
      end


