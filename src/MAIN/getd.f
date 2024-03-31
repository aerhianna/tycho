
      SUBROUTINE GETD(F,N,D,FP1,FPN)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
C
C  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS
C  OF UNITY.
C
C
      DIMENSION F(N),D(N),T(85)
C
      FP1=(-11.0d0*F(1)+18.0d0*F(2)-9.0d0*F(3)+2.0d0*F(4))/6.0d0
      FPN=(11.0d0*F(N)-18.0d0*F(N-1)+9.0d0*F(N-2)-2.0d0*F(N-3))/6.0d0
C
      D(1)=-0.5d0
      T(1)=0.5d0*(-F(1)+F(2)-FP1)
C
      DO 10 J=2,N-1
         D(J)=-1.0d0/(4.0d0+D(J-1))
         T(J)=-D(J)*(F(J-1)-2.d0*F(J)+F(J+1)-T(J-1))
   10 CONTINUE
C
      D(N)=(FPN+F(N-1)-F(N)-T(N-1))/(2.0d0+D(N-1))
C
      DO 20 J=N-1,1,-1
         D(J)=D(J)*D(J+1)+T(J)
   20 CONTINUE
C
      RETURN
      END
