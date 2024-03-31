c**********************************************************************
      subroutine readco1

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

c..... The purpose of this subroutine is to read the data tables
      save
      parameter (ismdata=0)   ! modified
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      character*1 dumarra(250)
      common/aa1/ q(4),h(4),xxh

      common/a1/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)

 
      common/b1/ itab(mx,mz),nta(nr),x(mx,mz),y(mx,mz),
     . zz(mx,mz)
      common/e/ opact,dopact,dopacr,dopactd
      common/ee1/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/alink1/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),xzff(100,nr)  
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
      COMMON/CST1/RLS,tmax,NRL,nset
c-------------------------------------------------------------
        if (itimeco .ne. 12345678) then
        do i=1,mx
          do j=1,mz 
            do k=1,nt
              do l=1,nr
                xz(i,j,k,l)=1.e+35
              enddo
            enddo
          enddo
        enddo
        itimeco=12345678
        endif
c
      close (2)
c..... read  tables
c..   Grevesse93 tables
c      open(2, FILE='GN93hz')
c..   Asplund05 tables
      open(2, FILE='hzdata')

c     read header
      do i=1,240
         read (2,'(a)') dumarra(i)
      enddo
c
      do 3 m=1,mx
      do 2 i=1,n(m)
c
      read(2,'(f10.5)') dum
      read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4)')
     1     itab(m,i),x(m,i),y(m,i),zz(m,i)
      read(2,'(f10.5)') dum,dum,dum
      read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
      read(2,'(f10.5)') dum
        do k=1,ntm
        read(2,'(f4.2,19f7.3)') alt(k),(xzf(k,l), l=1,nrm)
        alt(k)=alt(k)-6.
        if (isett6 .ne. 1234567) then
        t6listf(k)=10.**alt(k)
        t6arr(k)=t6listf(k)
        endif 
          do ll=1,nrm   ! modified
          xzff(k,ll)=xzf(k,ll)
          enddo
        enddo
        isett6=1234567

       if (ismdata .eq. 0) then
        tmax=10.   ! modified
        nset=65
        RLS=-8.
        nsm=1 
          RLE=1.
          nrlow=1
          nrhigh=2*(RLE-RLS)+1

        call opaltab1    !modified
       endif

 1010  continue

      ll=1
      do 110 kk=1,nre
      alr(ll)=alrf(kk)
        do k=1,nt
        t6list(k)=t6listf(k+ntb-1)
        if(ismdata .eq. 0) then
c           Following skip required because, due to missing data,
c           the X=0  low T data cannot be smoothed
          if ((m  .eq. 1) .and. (k .le. 9)) then
            xz(m,i,k,ll)=xzf(k+ntb-1,kk)
          else
            xz(m,i,k,ll)=xzff(k+ntb-1,kk)
          endif
        else
         xz(m,i,k,ll)=xzf(k+ntb-1,kk)
        endif
        enddo
  110 ll=ll+1

    2 continue
    3 continue
 
      do 12 i=2,nt
   12 dfs(i)=1./(alt(i)-alt(i-1))
      do 13 i=2,nr
   13 dfsr(i)=1./(alr(i)-alr(i-1))
      do i=2,mx-1
      dfsx(i)=1./(xx(i)-xx(i-1))
      enddo
      do i=2,mz
      dfsz(i)=1./(zza(i)-zza(i-1))
      enddo
      return
      end
c
c***********************************************************************
