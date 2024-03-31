      subroutine getfkt(ibug)

c..reads Thielemann file format
c..identifies symbol in Z,N
c..and echo to standard i/o
c---------------------------------------------------------------
c..to do: read netwinv to get symbols, Z, N, q, w and part. func.
c..and avoid a lot of this searching
c---------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'dimenfile'

c..these values apply to fkt reaclib 1988
c      parameter(niso = 806)
c      parameter( nzmax=36,nnmax=73)

c..these values apply to latest reaclib
      parameter(niso = 7852)
      parameter(nzmax = 112, nnmax =233)
c      parameter(niso = 5410)
c      parameter( nzmax=112,nnmax=225)
      integer*4 ncharge
      dimension p(7)
      character*5 lkh
      character*5 nam(6),inam(niso),blank,rnam(6)
      character*1 vw,nr

      dimension isoz(niso), ison(niso)
      dimension izr(6),inr(6),id(6)
      dimension nnn(nzmax,nnmax)

      logical tobe

      include 'crate'
      include 'comcsolve'

c===============================================================
      data blank/'     '/
c............................................
      do j=1,nnmax
        do i=1,nzmax
          nnn(i,j) = 0
        enddo
      enddo
c............................................
c..get isotope names, Z and N
      inquire(file='isotope.lib',exist=tobe)
      if( .not. tobe )then
         write(*,*)'No file isotope.lib found'
         stop'getfkt.f'
      endif
      open(12,file='isotope.lib')
      i = 0
10    continue
      i = i + 1
      read(12,'(a5,2i4)',end=11)inam(i),isoz(i),ison(i)
c      write(*,'(a5,2i4)')inam(i),isoz(i),ison(i)
      goto 10
11    continue
      i = i-1
c      write(*,*)'isotope.lib read,',i,' isotopes'
c      write(3,*)'isotope.lib read,',i,' isotopes'
      close(12)


c..read thielemann coef.s

      open(55,file='netsu')
c..friedel-style format
300   FORMAT(i2,3x,6a5,8x,a5,a1,a1,2x,1pe12.5,i5)
310   FORMAT(4e13.6)

      kk = 0
      ktype = 0
20    continue
      READ(55,300,END=90) k,(nam(j),j=1,6),lkh,nr,vw,qqf,iffnu
      READ(55,310) (p(j),j=1,7)

      if( k .ne. 0 )then
c..flag for type of reaction (collected in "decks")
        ktype = ktype + 1
c        write(*,*)'deck is ',ktype
        do j =1,6
         nam(j) = blank
        enddo
      else
c..store values in arrays for rate subroutine
        kk = kk+1
        do j = 1,6
          rname(j,kk) = nam(j)
        enddo
        do j = 1,7
          rcoef(j,kk) = p(j)
c          if(ktype .eq. 1)then
c             write(*,*)j,rcoef(j,kk)
c          endif
        enddo
        qval(kk) = qqf
        rlkh(kk) = lkh
        rnr(kk)  = nr
        rvw(kk)  = vw
        ideck(kk) = ktype
        iffn(kk) = iffnu
c..identify all participating nuclei
        call find(nam(1),id(1),inam,blank)
        call find(nam(2),id(2),inam,blank)
        call find(nam(3),id(3),inam,blank)
        call find(nam(4),id(4),inam,blank)
        call find(nam(5),id(5),inam,blank)
        call find(nam(6),id(6),inam,blank)
      endif

c      write(*,*)'before decks '

c--deck 1--------------------------------------------
      if( ktype .eq. 1 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1)  = isoz(id(1))
           inr(1)  = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         ncharge = izr(2) - izr(1)
         nnucleon = ncharge + inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        ncharge,nnucleon,lkh,nr,vw
        endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype,kk
           stop
         endif
         if( ncharge .eq. 0 )then
           write(*,*)'charge nonconservation',ktype,kk,izr(1),izr(2)
           stop
         endif
c         if( ncharge .eq. 1 )then
c           if( lkh .ne. 'bet-' )then
c            if( lkh .ne. ' ffn' )then
c             write(*,*)'error in beta decay'
c             write(*,'(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
c     1            rnam(1),izr(1),inr(1),id(1),
c     1            rnam(2),izr(2),inr(2),id(2),
c     1            ncharge,nnucleon,lkh,nr,vw
             
c             stop
c            endif
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                if( lkh .ne. ' ffn')then
c                 write(*,*)'error in positron decay'
c                 stop
c                endif
c               endif
c             endif
c           endif
         endif
      endif

c--deck 2--------------------------------------------
      if( ktype .eq. 2 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         ncharge = izr(3) + izr(2) - izr(1)
         nnucleon = ncharge + inr(3) + inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'2error in beta decay, ncharge=',ncharge
         endif
c         if( ncharge .eq. 1 )then
c           if( lkh .ne. 'bet-' )then
c             write(*,*)'2error in beta decay'
c             stop
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                 write(*,*)'2error in positron decay'
c                 stop
c               endif
c             endif
c           endif
         endif
      endif
c
c--deck 3--------------------------------------------
      if( ktype .eq. 3 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         ncharge  = izr(4) + izr(3) + izr(2) - izr(1)
         nnucleon = ncharge + inr(4) + inr(3) + inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'3error in beta decay, ncharge=',ncharge
         endif
c        if( ncharge .eq. 1 )then
c          if( lkh .ne. 'bet-' )then
c            write(*,*)'3error in beta decay'
c            stop
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                 write(*,*)'3error in positron decay'
c                 stop
c               endif
c             endif
c           endif
         endif
      endif
c
c--deck 4--------------------------------------------
      if( ktype .eq. 4 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         ncharge = izr(3) - izr(2) - izr(1)
         nnucleon = ncharge + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'4error in beta decay, ncharge=',ncharge
         endif
c         if( ncharge .eq. 1 )then
c           if( lkh .ne. 'bet-' )then
c             write(*,*)'4error in beta decay'
c             stop
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                 write(*,*)'4error in positron decay'
c                 stop
c               endif
c             endif
c           endif
         endif
      endif
c
c--deck 5--------------------------------------------
      if( ktype .eq. 5 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         ncharge  = izr(4) + izr(3) - izr(2) - izr(1)
         nnucleon = ncharge + inr(4) + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
         write(*,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         write(3,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'5error in beta decay, ncharge=',ncharge
         endif
c         if( ncharge .eq. 1 )then
c           if( lkh .ne. 'bet-' )then
c             write(*,*)'5error in beta decay'
c             stop
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                 write(*,*)'5error in positron decay'
c                 stop
c               endif
c             endif
c           endif
         endif
      endif
c
c--deck 6--------------------------------------------
      if( ktype .eq. 6 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         ncharge  = izr(5) + izr(4) + izr(3) - izr(2) - izr(1)
         nnucleon = ncharge + inr(5) + inr(4) + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(5(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(5(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'6error in beta decay, ncharge=',ncharge
         endif
c         if( ncharge .eq. 1 )then
c           if( lkh .ne. 'bet-' )then
c             write(*,*)'6error in beta decay'
c             stop
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                 write(*,*)'6error in positron decay'
c                 stop
c               endif
c             endif
c           endif
         endif
      endif
c--deck 7--------------------------------------------
      if( ktype .eq. 7 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         if( id(6) .gt. 0 .and. id(6) .le. niso )then
           rnam(6) = inam(id(6))
           izr(6) = isoz(id(6))
           inr(6) = ison(id(6))
         else
           rnam(6) = blank
           izr(6) = 0
           inr(6) = 0
         endif
         ncharge  = izr(6) + izr(5) + izr(4)
     1            + izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(6) + inr(5) + inr(4)
     2            + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        rnam(6),izr(6),inr(6),id(6),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        rnam(6),izr(6),inr(6),id(6),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype,nnucleon,ncharge,
     1     rnam(1),rnam(2),rnam(3),rnam(4),rnam(5),rnam(6)
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'7error in beta decay, ncharge=',ncharge
         endif
c         if( ncharge .eq. 1 )then
c           if( lkh .ne. 'bet-' )then
c             write(*,*)'7error in beta decay'
c             stop
c           endif
c         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
c           if( lkh .ne. 'bet+')then
c             if( lkh .ne. ' bec')then
c               if( lkh .ne. '  ec')then
c                 write(*,*)'7error in positron decay'
c                 stop
c               endif
c             endif
c           endif
         endif
      endif

c--deck 8--------------------------------------------
      if( ktype .eq. 8 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
c Deck 8 redefined by reaclib to e1 + e2 + e3 -> e4. Implemented 7/13
c         if( id(5) .gt. 0 .and. id(5) .le. niso )then
c           rnam(5) = inam(id(5))
c           izr(5) = isoz(id(5))
c           inr(5) = ison(id(5))
c         else
c           rnam(5) = blank
c           izr(5) = 0
c           inr(5) = 0
c         endif
         ncharge  = izr(4) - izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(4) - inr(3) - inr(2) - inr(1)
         if( nbug .ne. 0 )then
          write(*,'(5(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(5(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'8error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
c             write(*,*)'8error in beta decay'
c             stop
           endif
         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
c                 write(*,*)'8error in positron decay'
c                 stop
               endif
             endif
           endif
         endif
      endif
c--deck 9--Added decks 9, 10, 11. PAY 7/13
      if( ktype .eq. 9 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         ncharge  = izr(5) + izr(4)
     1            - izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(5) + inr(4)
     2            - inr(3) - inr(2) - inr(1)
         if( nbug .ne. 0 )then
          write(*,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'9error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
c             write(*,*)'9error in beta decay', lkh
c             stop
           endif
         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
c                 write(*,*)'9error in positron decay', lkh
c                 stop
               endif
             endif
           endif
         endif
      endif
c--deck 10--------------------------------------------
      if( ktype .eq. 10 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         if( id(6) .gt. 0 .and. id(6) .le. niso )then
           rnam(6) = inam(id(6))
           izr(6) = isoz(id(6))
           inr(6) = ison(id(6))
         else
           rnam(6) = blank
           izr(6) = 0
           inr(6) = 0
         endif
         ncharge  = izr(6) + izr(5) - izr(4)
     1            - izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(6) + inr(5) - inr(4)
     2            - inr(3) - inr(2) - inr(1)
         if( nbug .ne. 0 )then
          write(*,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        rnam(6),izr(6),inr(6),id(6),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        rnam(6),izr(6),inr(6),id(6),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'10error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
c             write(*,*)'10error in beta decay', lkh
c             stop
           endif
         endif
         if( ncharge .eq. -1 )then
            ec(k)=1
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
c                 write(*,*)'10error in positron decay', lkh
c                 stop
               endif
             endif
           endif
         endif
      endif
c--deck 11--Added decks 9, 10, 11. PAY 7/13
      if( ktype .eq. 11 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif

         ncharge  = izr(5) - izr(4)
     1            - izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(5) - inr(4)
     2            - inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          write(*,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           write(*,*)'11error in beta decay, ncharge=',ncharge
         endif
      endif

      do nc = 1, 6
      if( nam(nc) .ne. blank )then

c        if(  nnn(izr(nc),inr(nc))  .eq. 0 )then
c          if( izr(nc) .ne. 0 .and. inr(nc) .ne. 0)then

        if( izr(nc) .ne. 0 .and. inr(nc) .ne. 0)then
          if(  nnn(izr(nc),inr(nc))  .eq. 0 )then
            nnn(izr(nc),inr(nc)) = id(nc)
          endif
        endif
      endif

      enddo

c-------------------------------------------------------

      goto 20

90    continue
      close(55)
      write(*,*)' getfkt has finished reading netsu',
     1      kk,' reaction entries'

c      write(3,*)' getfkt has finished reading netsu',
c     1      kk,' reaction entries'

c--------------------------------------------------------------------
c--------------------------------------------------------------------

      number = 0
      do j=1,nnmax
        do i=1,nzmax
c..do not count alpha as nucleus
          if(nnn(i,j) .ne. 0 .and. nnn(i,j) .ne. 6 )then
            number = number + 1
          endif
        enddo
      enddo
      itot = number + 3

      write(*,*)itot,' nuclei',nzmax,nnmax

      if( itot+1 .gt. ndim )then
        write(*,'(a20,2(a10,i5))')'getfkt:','itot+1',itot+1,'ndim',ndim
        write(*,*)'network too large for dimensions'
        write(*,*)'in dimenfile,'
        write(*,*)'change ndim to .ge. ',itot+1
        write(*,*)'change nreac to .ge. ',kk

        stop'getfkt 1'
      endif
      inuc = number
      ireac = kk


c..setup nz,nn,xid vectors for burn
      l = 0
      do k = 1, niso
        do j=1,nnmax
          do i=1,nzmax
c..do not count alpha as nucleus
            if( k .eq. nnn(i,j) .and. k .ne. 6 )then
              l = l+1
              nz(l) = isoz(k)
              nn(l) = ison(k)
              xid(l) = inam(k)
              if( ibug .ne. 0 )then
                write(*,'(3i5,2x,a5)')l,nz(l),nn(l),xid(l)
              endif
            endif
          enddo
        enddo
      enddo
c==============================================
c 1,2, and 6 are n,p, and alphas in isotope.lib
c==============================================
      xid(number+1) = inam(1)
      nz (number+1) = isoz(1)
      nn (number+1) = ison(1)
      xid(number+2) = inam(2)
      nz (number+2) = isoz(2)
      nn (number+2) = ison(2)
      xid(number+3) = inam(6)
      nz (number+3) = isoz(6)
      nn (number+3) = ison(6)
c electrons are in number + 4 = ndim usually
      xid(number+4) = '   Ye'

      if( ibug .ne. 0 )then
        do l = number+1, number+3
              write(*,'(3i5,2x,a5)')l,nz(l),nn(l),xid(l)
        enddo
      endif
c
c..count number of reactions of each type (deck)
      do j = 1, 11
        ndeck(j) = 0
      enddo
      do i = 1, ireac
        do j = 1,11
        if( ideck(i) .eq. j )ndeck(j) = ndeck(j) + 1
        enddo
      enddo

      write(*,*)
     1 'number of reactions ndeck(n) in each deck (n=1 to 11)'
      write(*,'(8i5)')ndeck
      if( ndeck(2) .ne. ndeck(4) )then
        write(*,'(a50)')
     1 'WARNING: 1-->2+3 does not match inverse 2+3-->1'
        write(*,'(a50)')'         because ndeck(2) .ne. ndeck(4)'
      endif
      if( ndeck(3) .ne. ndeck(8) )then
        write(*,'(a50)')
     1 'WARNING: 1-->2+3+4 does not match inverse 2+3+4-->1'
        write(*,'(a50)')'         because ndeck(3) .ne. ndeck(8)'
      endif
      if( ndeck(6) .ne. ndeck(9) )then
        write(*,'(a55)')
     1 'WARNING: 1+2-->3+4+5 does not match inverse 3+4+5-->1+2'
        write(*,'(a55)')'         because ndeck(6) .ne. ndeck(9)'
      endif
      if( ndeck(7) .ne. ndeck(10) )then
        write(*,'(a55)')
     1 'WARNING: 1+2-->3+4+5+6 does not match inverse 3+4+5+6-->1+2'
        write(*,'(a55)')'         because ndeck(7) .ne. ndeck(10)'
      endif
      if( mod( ndeck(5),2 ) .ne. 0 )then
        write(*,'(a50)')
     1 'WARNING: 1-->2+3 does not match inverse 2+3-->1'
        write(*,'(a50)')'         because ndeck(2) .ne. ndeck(4)'
      endif

      return
      end



