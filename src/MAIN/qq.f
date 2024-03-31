      FUNCTION QQ(BETA,DELTA)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ++++ CALCULATES Q AS FUNCTION OF BETA AND DELTA ++++
C ++++ Q = A(BETA)**DELTA-1. ++++++++
C + A(BETA) LINEAR BETWEEN A(2) = 22.1,A(1) = 7.5,A(0.7) = 4.0,++
C + A(0.5) = 2.5,A(0.25) = 1.18 ++++++
      IF(BETA .LT. 0.2d0) GO TO 2
      IF(BETA .LT. 2.5d0) GO TO 1
      WRITE(*,*)'from FUNCTION QQ: BETA TOO LARGE'
      STOP
 1         IF(BETA.LT.1.0d0) GO TO 11
      A = 15.0d0*(BETA-1.0d0)+7.5d0
      GO TO 4
 11       IF(BETA.LT.0.7d0) GO TO 12
      A = 11.66667d0*(BETA-0.7d0)+4.0d0
      GO TO 4
 2         WRITE(*,*)'from FUNCTION QQ: BETA TOO SMALL'
      STOP
 12       IF(BETA.LT.0.5d0) GO TO 13
      A = 7.5d0*(BETA-0.5d0)+2.5d0
      GO TO 4
 13       IF(BETA.LT.0.2d0) GO TO 2
      A = 5.28d0*(BETA-0.25d0)+1.18d0
 4         CONTINUE
      QQ = A**DELTA-1.0d0
      RETURN
      END


