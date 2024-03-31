      FUNCTION BETCR(ALPHA,DELTA)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ***** CALCULATES BETA VALUE FOR UCRIT ACC.TO PAULDRACH *
C       IT HAS BEEN MODIFIED TO TO BE SMOOTH
      IF (ALPHA.LE.0.65d0) GO TO 2
      IF (ALPHA.GT.0.75d0) GO TO 1
      BLOW = ALOW(DELTA)
      BHIG = AHIG(DELTA)
      BETCR = 10.0d0*(BHIG-BLOW)*(ALPHA-0.65d0)+BLOW
      RETURN
 1         BETCR = AHIG(DELTA)
      RETURN
 2         BETCR = ALOW(DELTA)
      RETURN
      END


