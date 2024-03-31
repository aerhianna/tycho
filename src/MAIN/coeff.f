      SUBROUTINE COEFF(ALPHA,BETA,A0,A1,B1,B2)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ++++ COEFFICIENTS A0,A1,B1,B2 FOR APROX. TERM. VELOCITY ++++
      HILFEX = 1.0d0-ALPHA
      HILFEX = 1.0d0/HILFEX
      ALAAF = (1.0d0/(1.0d0+ALPHA))**HILFEX
      A0 = (1.0d0+BETA-ALAAF)/BETA
      A1 = A0-ALAAF
      B1 = ALPHA*HILFEX/(2.0d0*BETA)
      B2 = B1*(BETA+1.0d0)
      RETURN
      END










