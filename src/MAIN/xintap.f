      FUNCTION XINTAP(ALPHA,BETA,DELTA,UCRIT)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C ++++++ APPROX. WIND INTEGRAL = EQ.(60) +++++++
      Q = QQ(BETA,DELTA)
      CALL COEFF(ALPHA,BETA,A0,A1,B1,B2)
      BET1 = 1.0d0/(1.0d0+BETA)
      BET2 = BET1*BET1
      BET3 = BET2*BET1
      BET4 = BET3*BET1
      BET5 = BET4*BET1
      UC1 = UCRIT
      UC2 = UC1*UC1
      UC3 = UC2*UC1
      UC4 = UC3*UC1
      HILF = 1.0d0/(1.0d0-ALPHA)
      X = A0*(UC1-BET1+HILF*Q*(UC3-BET3)/3.0d0)
      X = X-A1*((UC2-BET2)/2.0d0+HILF*Q*(UC4-BET4)/4.0d0)
      X = X+BET1
      X = X+B1*(BET2/2.0d0+HILF*Q*BET4/4.0d0)
      X = X-B2*(BET3/3.0d0+HILF*Q*BET5/5.0d0)
      X = X+HILF*Q*BET3/3.0d0
      G = 1.0d0/(A0-A1*UC1)
      G = G*((1.0d0/(Q*UC2+1.0d0))**HILF)
      Z = (1.0d0/G)
C ***************** NEW *************************************
      IF(Z.LT.1.0d0) GO TO 1
      Z = G
      GO TO 2
C **************** UNTIL HERE *******************************
 1         Z = 2.0d0/ALPHA*(1.0d0-Z)
      Z = G*(1.0d0+SQRT(Z))
 2         XINTAP = Z*X
      RETURN
      END


