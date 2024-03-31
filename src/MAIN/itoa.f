c
c
c     
      subroutine itoa(i,c)
      
      implicit none

c     translates number (integer i) to character (a code)
c     0 .le. number .le. 9
      
      integer i, k
      integer icode(10)
      character c
      character ccode(10)
      
      data ccode/'1','2','3','4','5','6','7','8','9','0'/
      data icode/1,2,3,4,5,6,7,8,9,0/
c     
      c = ccode(10)
      
      do k = 1, 9
         if( i .eq. icode(k) ) c = ccode(k)
      enddo
      
c     
      return
c     
      end
