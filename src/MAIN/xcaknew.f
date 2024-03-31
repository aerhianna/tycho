      FUNCTION XCAKNEW(XK,ALPHA,ALEFF,STLUM,SIG,CA,CAFAC,BFAC,VPROT,
     &     VIERPI)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C     ++++ XMDCAK = MDOT AFTER CAK
C     BUT NOW INCLUDING V_SOUND/V_ESC AND V_SOUND/V_C TERMS AT
c     CRITICAL POINT
      B  = BFAC
      AS = CAFAC
C     END OF TEST
      HILFEX  = 1.0d0/ALPHA
      HILF    = 1.2762d023*XK*STLUM
      hilf    = log10(hilf)
      fac1    = HILF*HILFEX
      fac2    = ALEFF/VPROT/B
      fac2    = log10(fac2)
      HILFEX  = HILFEX*(1.0d0-ALPHA)
      HI      = SIG*(1.0d0-ALEFF)/(VIERPI*CA*AS)
      hi      = log10(hi)
      fac3    = hi*hilfex
      xlcak   = fac1+fac2+fac3
      xcaknew = 10.0d0**xlcak
      RETURN
      END


