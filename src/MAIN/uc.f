      FUNCTION UC(ALPHA,DELTA,VSQUAR)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ******** ESTIMATE OF RECIPROCAL CRITICAL POINT******
      SIGMA = ALPHA*ALPHA*(1.0d0-ALPHA)/VSQUAR
      XPHI  = PHICR(VSQUAR)
      XPHI  = XPHI*XPHI
      SIGMA = SIGMA/XPHI
      BBB   = BETCR(ALPHA,DELTA)
      Q     = QQ(BBB,DELTA)
      XC    = -2.0d0*Q/SIGMA/3.0d0
      HILF  = 1.46d0/SIGMA
      HILF  = HILF*(ALPHA+1.0d0)*(ALPHA**0.6d0)
      HILF  = HILF**0.33333333333d0
      XC    = HILF + XC
      xc    = dmax1(xc, 1.03d0)
      UC    = 1.0d0/XC
      RETURN
      END


