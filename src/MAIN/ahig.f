      FUNCTION AHIG(DELTA)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C          CALCULATES BETACRIT FOR ALPHA LE 0.65
      IF (DELTA .LT. 0.1d0) GO TO 1
      AHIG = 0.7d0
      RETURN
 1         IF (DELTA .LT. 0.085d0) GO TO 2
      AHIG = -20.0d0*(DELTA-0.1d0)+0.7d0
      RETURN
 2         IF (DELTA .LT. 0.06d0) GO TO 3
      AHIG = 1.0d0
      RETURN
 3         IF (DELTA .LT. 0.05d0) GO TO 4
      AHIG = -100.0d0*(DELTA-0.06d0)+1.0d0
      RETURN
 4         AHIG = 2.0d0
      RETURN
      END

