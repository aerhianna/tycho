      FUNCTION DILUT(U)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C          CALCULATES DILUTION FACTOR
      IF(U .LT. 1.d-3) GO TO 1
      X = 1.0d0-U*U
      DILUT = 0.5d0*(1.0d0-SQRT(X))
      RETURN
 1         DILUT = 0.25d0*U*U
      RETURN
      END


