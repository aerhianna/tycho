      subroutine modflg(model,flag)
      implicit none

c..input - model number

c..output - corresponding a-format representation
c..replaces leading zeros with blanks

      integer*4    model, i, k, l, irem, nscale
      character*5  flag
      character    c, cc(5), zero, blank

      data zero/'0'/, blank/' '/
c-----------------------------------------------
c  integer remainder = irem
      irem = model

      do  k = 1, 5
        l = 5 + 1 - k
        nscale = 10**(l-1)
        i = irem/nscale

            call itoa(i,c)

        cc(k) = c
        irem = irem - i*nscale

      enddo

c..replace leading zeros with blanks
      do k = 1,4
        if( cc(k) .ne. zero )then
          goto 100
        else
          cc(k) = blank
        endif
      enddo

 100  continue

      flag = cc(1)//cc(2)//cc(3)//cc(4)//cc(5)

      return
      end
