      FUNCTION ALOW(DELTA)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C          CALCULATES BETACRIT FOR ALPHA GT 0.75
      IF (DELTA .LT. 0.1d0) GO TO 1
      ALOW = 0.25d0
      RETURN
 1         IF (DELTA .LT. 0.085d0) GO TO 2
      ALOW = -30.0d0*(DELTA-0.1d0)+0.25d0
      RETURN
 2         IF (DELTA .LT. 0.06d0) GO TO 3
      ALOW = 0.7d0
      RETURN
 3         IF (DELTA .LT. 0.05d0) GO TO 4
      ALOW = -30.0d0*(DELTA-0.06d0)+0.7d0
      RETURN
 4         IF (DELTA .LT. 0.035d0) GO TO 5
      ALOW = 1.0d0
      RETURN
 5         IF (DELTA .LT. 0.025d0) GO TO 6
      ALOW = -100.0d0*(DELTA-0.035d0)+1.0d0
      RETURN
 6         ALOW = 2.0d0
      RETURN
      END


