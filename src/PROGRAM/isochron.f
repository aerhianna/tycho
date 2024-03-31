      program isochron

      integer*4 nnuc, nstep, nmod, maxmod, intin(2), i, j, k

      parameter (nnuc = 178, maxmod = 100)
      real*8 surfin(12), nucin(nnuc), age, iso(2,maxmod), myear,
     1       agediff, agediffold
      character*5 modlist(maxmod)
      character*8 filein
      character*12 fileout

      data myear/3.15d13/

      write(*,*)'Enter isochrone age in My: '
      read(*,*)age
      write(*,*)'Enter name of list of hr files: '
      read(*,'(a8)')filein

      nmod=1

      open(8,file=filein)

 100  continue
      read(8,'(a5)',end=101)modlist(nmod)

      open(9,file=modlist(nmod))
      agediffold = 9.9d99

 102  continue
      read(9,'(2i6,1pe14.6,1p8e12.4,1p3e11.3)',end=103)intin(1),
     1    intin(2),surfin
      read(9,'(1p10e11.3)')nucin
c      write(*,*)surfin
      agediff = dabs(age*myear - surfin(1))
      if(agediff .lt. agediffold)then
         agediffold = agediff
      else
         goto 103
      endif
      goto 102

 103  continue
      close(9)
      iso(1,nmod) = surfin(6)
      iso(2,nmod) = surfin(7)
c      write(*,*)iso(1,nmod),iso(2,nmod),nmod
      nmod = nmod+1
      goto 100

 101  continue
      close(8)
      fileout = filein//'.iso'
      open(10,file=fileout)
      do i = 1, nmod-1
         write(10,'(1p2e11.4)')iso(1,i), iso(2,i)
      enddo

      close(10)

      end
