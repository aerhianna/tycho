
      subroutine sintr3(xold,xnew,f,fnew)
      implicit none

c..   for single point
c--------------------------------------------------------------
c     interpolates for  xnew   in old  xold(4) array
c     to find new       fnew from old    f(4)
c     does cubic  interpolation; 
c     will extrapolate below xold(1) or above xold(4)

      integer*4 k
      real*8    xold(4),f(4),xnew,fnew
      real*8    fact1,fact2,fact3,fact0

c--------------------------------------------------------------
      k = 1
c..   cubic interpolation
      fact0 = ( (xnew   - xold(k+1))*(xnew      - xold(k+2)) 
     1     *    (xnew   - xold(k+3))  ) /
     2     ( (xold(k)   - xold(k+1))*(xold(k)   - xold(k+2))
     3     * (xold(k)   - xold(k+3)) )

      fact1 = ( (xnew   - xold(k)  )*(xnew      - xold(k+2))
     1     *    (xnew   - xold(k+3))  ) /
     2     ( (xold(k+1) - xold(k)  )*(xold(k+1) - xold(k+2))
     3     * (xold(k+1) - xold(k+3))  )

      fact2 = ( (xnew   - xold(k+1))*(xnew      - xold(k)  )
     1     *    (xnew   - xold(k+3)) ) /
     1     ( (xold(k+2) - xold(k+1))*(xold(k+2) - xold(k)  )
     3     * (xold(k+2) - xold(k+3)) )

      fact3 = ( (xnew   - xold(k))*(xnew      - xold(k+1)  )
     1     *    (xnew   - xold(k+2)) ) /
     1     ( (xold(k+3) - xold(k))*(xold(k+3) - xold(k+1)  )
     3     * (xold(k+3) - xold(k+2)) )


         fnew = fact0*f(k) + fact1*f(k+1) + fact2*f(k+2) + fact3*f(k+3)
     
      return
      end



