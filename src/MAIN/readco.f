      subroutine readco
c..... The purpose of this subroutine is to read the data tables
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      save
      parameter (ismdata=0)   ! modified
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb
     . ,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo)
     . ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx)
     . ,nc,no
c      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101),
c     . t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr)
c     . ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) 
c     .,t6listf(ntm),opk2(nt,nr),dfsx(mx)

      logical tobe

ccccccccccccccccccccccccc

c THIS SUBROUTINE USES "INT" AS A VARIABLE; IT IS AN INTRISIC FUNCTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr),
     . t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr)
     . ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) 
     .,t6listf(ntm),opk2(nt,nr),dfsx(mx), index(101)

c      common/b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs),
c     . zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
ccccccccccccccccccc
      common/b/ x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),
     1     xca(mx,ntabs),xoa(mx,ntabs),itab(mx,ntabs),nta(nrm)

      common/alink/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),coff(100,nr)
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
ccccccccccccccccc
      COMMON/CST/RLS,tmax,NRL,nset

      common/e/ opact,dopact,dopacr,dopacrd
      character*1 dumarra(250)
      common/recoin/ itimeco,mxzero

c--------------------------------------------------------
        if (itimeco .ne. 12345678) then
        do i=1,mx
          do j=1,mc
            do k=1,mo
              do l=1,nt
                do mq=1,nr
                  co(i,j,k,l,mq)=1.d+35
                enddo
              enddo
            enddo
          enddo
        enddo
        do i=1,mx
          do j=1,mc
            do l=1,nt
              do mq=1,nr
                diag(i,j,l,mq)=1.d+35
                diago(i,j,l,mq)=1.d+35
              enddo
            enddo
          enddo
        enddo
        itimeco=12345678
        endif
      do 20 j=1,nc-1
       do 21 i=1,nc
         if(xcd(j) .ge. xc(i)) then
           n(m,j)=i+1
           if(xcd(j) .lt. xc(i)+1.d-6) n(m,j)=i
         endif
   21  continue
   20 continue
      n(m,nc)=0
c
      close (2)
c..... read X=0.0 tables
      if(m .eq. 1)then
         inquire(file='codataa',exist=tobe)
         if( .not. tobe )then
            write(*,*)'No file codataa found'
            stop'readco.f'
         endif
         open(2, FILE='codataa')
      endif
c..... read X=0.03 tables
      if(m .eq. 2)then
         inquire(file='codatab',exist=tobe)
         if( .not. tobe )then
            write(*,*)'No file codatab found'
            stop'readco.f'
         endif
         open(2, FILE='codatab')
      endif
c..... read X=0.10 tables
      if(m .eq. 3)then
         inquire(file='codatac',exist=tobe)
         if( .not. tobe )then
            write(*,*)'No file codatac found'
            stop'readco.f'
         endif
         open(2, FILE='codatac')
      endif
c..... read X=0.35 tables
      if(m .eq. 4)then
         inquire(file='codatad',exist=tobe)
         if( .not. tobe )then
            write(*,*)'No file codatad found'
            stop'readco.f'
         endif
         open(2, FILE='codatad')
      endif
c.....read X=0.70 tables
      if(m .eq. 5)then
         inquire(file='codatae',exist=tobe)
         if( .not. tobe )then
            write(*,*)'No file codatae found'
            stop'readco.f'
         endif
         open(2, FILE='codatae')
      endif

c      read header
c      write(*,*)'readco ',m
ccccccccccccccccccccccccccccccc

      read(2,'(a)') (dumarra(i),i=1,240)
c
      int=0
      do 1 j=1,no-1
      do 2 i=1,n(m,j)
        int=int+1
c
        read(2,'(f10.5)') dum
        read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)')
     x  itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
        xca(m,int)=min(xca(m,int),1.d0-x(m,int)-zz(m,int)-xoa(m,int))

        read(2,'(f10.5)') dum,dum,dum
        read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
        read(2,'(f10.5)') dum
          do k=1,ntm
            read(2,'(f4.2,19f7.3)') altin,(cof(k,l), l=1,nrm)

            do ll=1,nrm   ! modified
              coff(k,ll)=cof(k,ll)
            enddo
          enddo
          if (isett6 .ne. 1234567) then
          do k=1,ntm
            t6arr(k)=t6list(k)
          enddo
          endif
          isett6=1234567

       if (ismdata .eq. 0) then
         if ((nrm .ne. nr) .or. (ntm .ne. nt)) then
           write (*,'("Not set up to smooth data with reduced ",
     x                "T-Rho grid")')
           stop
         endif
        tmax=10.d0   ! modified
        nset=67    ! 65 in earlier version
        RLS=-8.d0
        nsm=1
          RLE=1.0d0
          nrlow=1
          nrhigh=2*idint(RLE-RLS)+1

        call opaltab    !modified
       endif

      ll=1
      do 110 kk=nrb,nre
      alr(ll)=alrf(kk)
        do k=1,nt
          if (ismdata .eq. 0) then
            if ((m .eq. 1) .and. (k .le. 9)) then
              co(m,i,j,k,ll)=cof(k+ntb-1,kk)
            else
              co(m,i,j,k,ll)=coff(k+ntb-1,kk)
            endif
          else
            co(m,i,j,k,ll)=coff(k+ntb-1,kk)
          endif
        enddo
  110 ll=ll+1
    2 continue
    1 continue
      if(x(m,1) .ne. xa(m)) then
      write(*,'(" X in the codata? file does not match xa(m)")')
      stop
      endif
c
      do i=1, nc-1
       do k=1,nt
        do l=1,nr
          diag(m,i,k,l)=co(m,n(m,i),i,k,l)
        enddo
       enddo
      enddo
c
      do 6 j=1,no-1
        int=int+1

c        write(*,*)'readco ',m,int,nc,no
cccccccccccccccccccccccccccccccccccc

        read(2,'(f10.5)') dum
        read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)')
     x  itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
        read(2,'(f10.5)') dum,dum,dum
        read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
        read(2,'(f10.5)') dum
         do k=1,ntm
           read (2,'(f4.2,19f7.3)') dum, (cof(k,l), l=1,nrm)
c     set up to smooth final "diago" opacity tables
           do l=1,nrm
              coff(k,l)=cof(k,l)
           enddo
        enddo
c     smooth final "diago" opacity tables too!
       if (ismdata .eq. 0) then
        tmax=10.d0   ! modified
        nset=67 !65 in earlier version
        RLS=-8.d0
        nsm=1
          RLE=1.0d0
          nrlow=1
          nrhigh=2*idint(RLE-RLS)+1
        call opaltab    !modified
        do k=3,NTEMP-2
          do ll=nrlow,nrhigh
c           Following skip required because, due to missing data,
c           the X=0  low T data cannot be smoothed
            if ((m .eq. 1) .and. (k .le. 9)) then
              cof(k,ll)=cof(k,ll)
            else
              cof(k,ll)=coff(k,ll)
            endif
          enddo
        enddo
       ll=1
       do kk=nrb,nre
         do k=1,nt
           diago(m,j,k,ll)=cof(k+ntb-1,kk)
         enddo
       ll=ll+1
       enddo
       endif
    6 continue

      do i=2,nt
        dfs(i)=1.d0/(alt(i)-alt(i-1))
      enddo
      do i=2,nr
       dfsr(i)=1.d0/(alr(i)-alr(i-1))
      enddo
      istep=-1
      if (mx .gt. 1 ) then
        istep=1
        do i=2,mx,istep
          dfsx(i)=1.d0/(xx(i)-xx(i-1))
        enddo
      endif
      return
      end
c
