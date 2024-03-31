      FUNCTION FCCOR(U,VQ,Y,ALPHA,RSTCM,IAPPR)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C  FINITE CONE ANGLE CORRECTION FACTOR FOR GIVEN U,Y,V,RST(CM)
C
C             for IAPPR = 0  we adopt h = 0
C                       1             = calculated correctly
      IF (U .LT. 1.d-5) THEN
        FCCOR = 1.0d0
        RETURN
      ELSE
      CONTINUE
      ENDIF
      IF (IAPPR .EQ. 1) THEN
        XL = U*U*(1.0d0-RSTCM*VQ/U/Y)
      ELSE
        XL = U*U
      ENDIF
      ABXL = ABS(XL)
      IF (ABXL .LE. 1.d-03) THEN
        FCCOR = 1.0d0-ALPHA*XL/2.0d0
        fcccc = fccor
      ELSE
      AP1 = ALPHA+1.0d0
      FCCOR = (1.0d0-(1.0d0-XL)**AP1)/AP1/XL
      fcccc = fccor
      ENDIF
      RETURN
      END


