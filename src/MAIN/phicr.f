      FUNCTION PHICR(VSQUAR)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ****** CALCULATES PAULDRACH PHI-FUNCTION AT UCRIT****
      VSQ   = SQRT(VSQUAR)
      XL    = log10(VSQ)
      XL    = 0.36d0+XL
      XL    = 0.3d0*XL
      XL    = VSQ**XL
      PHICR = 3.0d0*XL
      RETURN
      END

