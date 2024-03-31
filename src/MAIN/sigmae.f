      FUNCTION SIGMAE(YPS,XIHE)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C +++++ SIGMAE = THOMSON ABS.COEFF. DIVIDED BY DENSITY
C +++++ XIHE = NUMBER OF ELECTRONS PROVIDED BY HE NUCLEUS
C +++++ YPS = N(HE)/N(H)
      SIGMAE = (1.0d0+XIHE*YPS)/(1.0d0+4.0d0*YPS)
      SIGMAE = SIGMAE*0.3975d0
      RETURN
      END


