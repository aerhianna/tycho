c     
c     
      subroutine modflgo(model,flag)
      
c..   input   - model number
c..   output  - corresponding a-format representation
c..   does not replace leading zeros with blanks
      
      implicit none
      integer*4    model, i, k, l, irem, nscale
      character*5  flag
      character    c, cc(5)

      
      irem = model
      
      do  k = 1, 5
         l = 5 + 1 - k
         nscale = 10**(l-1)
         i = irem/nscale
         
         call itoa(i,c)
         
         cc(k) = c
         irem = irem - i*nscale
         
      enddo
      
      flag = cc(1)//cc(2)//cc(3)//cc(4)//cc(5)

      
c     SUCCESS
      return
c     
      end
