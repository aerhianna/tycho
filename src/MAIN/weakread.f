      subroutine weakread

c..   read netweak and set up arrays interpolation in t9, rho

      implicit none

      include 'dimenfile'
      include 'crate'
      include 'comcsolve'

      integer*4 rloc, i,j,k, kup, klo
c..   defined in dimenfile      parameter(tsize = 13, rhosize = 11)
      real*8 ratin(tsize,rhosize), 
     1       enuin(tsize,rhosize)
      
      character*72 null
c      character*5 ci1,ci2
c      character*4 clrat
c-----------------------------------------------------------------

      kup = k1deck(1)
      klo = k2deck(2)
      do k = kup, klo
c..loop over decks 1 and 2 which have weak rates
         if( iffn(k) .gt. 0 )then

c..   location in netweak data set
            rloc = iffn(k) 
c..   open file starts at the beginning
            open(19,file='netweak')

            do i = 1, rloc-1
               read(19,*)null
               do j = 1,32
                  read(19,*)null
               enddo
            enddo   
            read(19,*)null

c..   weak rate and neutrino energy emission arrays (in log10)
c..   netweak data is listed as
c     x(t1,d1),x(t2,d1),...,x(t13,d1)
c..   x(t1,d2),x(t2,d2),...,x(t13,d2)
c..   ...                   x(t13,d11)

            read(19,*)((ratin(i,j),enuin(i,j),i=1,tsize),
     1           j=1,rhosize)

            do i=1, tsize
               do j=1,rhosize
                  ratarray(rloc,i,j) = ratin(i,j)
                  enuarray(rloc,i,j) = enuin(i,j)
               enddo
            enddo

            close(19)

         endif            
      enddo

      return
      end

      
