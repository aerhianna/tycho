      FUNCTION XFCNEW(XMDCA,UCRIT,YCRIT,VC,CNW11,ALPHA,DELTA,RSTCM,
     &  IAPPR)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C               CALCULATES M_dot finite cone
C               using new formula of Kudritzki 97, July 4
      ALPRIME = ALPHA-DELTA
      EX1 = DELTA/ALPRIME
      EX2 = 1.0d0/ALPRIME
      EX3 = ALPHA/ALPRIME
      W = DILUT(UCRIT)
      FAC1 = CNW11*UCRIT*UCRIT/VC/W
      xlfac1 = log10(fac1)
      xlfac1 = xlfac1*ex1
      FAC2 = FCCOR(UCRIT,VC,YCRIT,ALPHA,RSTCM,IAPPR)
      xlfac2 = log10(fac2)
      xlfac3 = log10(xmdca)
      xlfac2 = xlfac2*ex2
      xlfac3 = xlfac3*ex3
      XFCNEW = XLFAC1 + XLFAC2 + XLFAC3
      xfcnew = 10.0d0**xfcnew
      RETURN
      END


