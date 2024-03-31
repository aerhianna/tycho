      FUNCTION VAUCR(VSOUND,VESC,ALPHA)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ***** COMPUTES VELOCITY AT CRITICAL POINT ****
C ***** ACCORDING TO PAULDRACH *****
      VSQ   = VSOUND*VSOUND/VESC/VESC
      PHI   = PHICR(VSQ)
      HILF  = PHI/(1.0d0-ALPHA)
      HILF  = SQRT(HILF)
      EXXX  = 2.0d0/ALPHA
      H     = (1.0d0-ALPHA)**EXXX
      H     = 1.0d0-H
      HILF  = HILF/H
      VAUCR = VSOUND*HILF
      RETURN
      END


