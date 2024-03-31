      SUBROUTINE INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

C
C  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE
C  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN
C  Numerical Recipes, PP. 118 TO 120
C
C
      PARAMETER(IPR=20)
      COMMON/CF/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      DIMENSION B(16)
      LOGICAL IERR
C
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
ccccccccccccccccc
      COMMON/CST/RLS,tmax,NRL,nset

C
C  EXTREME LIMITS ALLOWED ARE:-
C     (3.800-0.0125) TO (8.000+0.0125) FOR LOG10(T)
C     (RLS-0.125) TO (RLE+0.1254) FOR LOG10(R)
C     (ALLOWING FOR SMALL EXTRAPOLATIONS BEYOND TABULAR VALUES)
C
C  FUNCTION DEFINITIONS FOR CUBIC EXPANSION
C
      FF(S,T)=    B( 1)+T*(B( 2)+T*(B( 3)+T*B( 4)))
     +   +S*(     B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     +   +S*(     B( 9)+T*(B(10)+T*(B(11)+T*B(12)))
     +   +S*(     B(13)+T*(B(14)+T*(B(15)+T*B(16))) )))
C
      FFX(S,T)=   B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     +   +S*(  2*(B( 9)+T*(B(10)+T*(B(11)+T*B(12))))
     +   +S*(  3*(B(13)+T*(B(14)+T*(B(15)+T*B(16)))) ))
C
      FFY(S,T)=   B( 2)+S*(B( 6)+S*(B(10)+S*B(14)))
     +   +T*(  2*(B( 3)+S*(B( 7)+S*(B(11)+S*B(15))))
     +   +T*(  3*(B( 4)+S*(B( 8)+S*(B(12)+S*B(16)))) ))
C
      FFXY(S,T)=  B( 6)+T*(2*B( 7)+3*T*B( 8))
     +   +S*(   2*B(10)+T*(4*B(11)+6*T*B(12))
     +   +S*(   3*B(14)+T*(6*B(15)+9*T*B(16)) ))
C
C
      IERR=.FALSE.
C
      X=20.d0*(FLT-3.800d0)+1
      FLR=FLRHO+18.d0-3.d0*FLT
      Y=2*( FLR - RLS )+1
C
      IF(X.LT.2.)THEN
         IF(X.LT.0.75d0)THEN
            IERR=.TRUE.
            i=0
         ELSE
            I=1
         ENDIF
      ELSEIF(X.GT.84)THEN
         IF(X.GT.85.25d0)THEN
            IERR=.TRUE.
            i=85
         ELSE
            I=84
         ENDIF
      ELSE
         I=int(X)
      ENDIF
      U=X-I
C
      IF(Y.LT.2.)THEN
         IF(Y.LT.0.75d0)THEN
            IERR=.TRUE.
            j=0
         ELSE
            J=1
         ENDIF
      ELSEIF(Y.GT.NRL-1)THEN
         IF(Y.GT.NRL+0.25d0)THEN
            IERR=.TRUE.
            j=nrl
         ELSE
            J=NRL-1
         ENDIF
      ELSE
         J=int(Y)
      ENDIF
      V=Y-J
C
      IF(IERR)THEN
         G=9.999d0
         DGDT=9.999d0
         DGDRHO=9.999d0
         RETURN
      ENDIF
C
C
C  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
      B(1)=F(I,J)
      B(2)=FY(I,J)
      B(3)=3*(-F(I,J)+F(I,J+1))-2*FY(I,J)-FY(I,J+1)
      B(4)=2*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1)
C
      B(5)=FX(I,J)
      B(6)=FXY(I,J)
      B(7)=3*(-FX(I,J)+FX(I,J+1))-2*FXY(I,J)-FXY(I,J+1)
      B(8)=2*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1)
C
      B(9)=3*(-F(I,J)+F(I+1,J))-2*FX(I,J)-FX(I+1,J)
      B(10)=3*(-FY(I,J)+FY(I+1,J))-2*FXY(I,J)-FXY(I+1,J)
      B(11)=9*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     + +6*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))
     + +4*FXY(I,J)
     + +3*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))
     + +2*(FXY(I,J+1)+FXY(I+1,J))
     + +FXY(I+1,J+1)
      B(12)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     + +4*(-FX(I,J)+FX(I,J+1))
     + +3*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))
     + +2*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))
     + -FXY(I+1,J)-FXY(I+1,J+1)
C
      B(13)=2*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J)
      B(14)=2*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J)
      B(15)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     + +4*(-FY(I,J)+FY(I+1,J))
     + +3*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))
     + +2*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))
     + -FXY(I+1,J+1)-FXY(I,J+1)
      B(16)=4*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     + +2*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)
     +    +FY(I,J)-FY(I+1,J)-FY(I+1,J+1)+FY(I,J+1))
     + +FXY(I,J)+FXY(I+1,J)+FXY(I+1,J+1)+FXY(I,J+1)
C
C  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),
C      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)
      G=FF(U,V)
      DGDT=20.0d0*FFX(U,V)-6.0d0*FFY(U,V)
      DGDRHO=2.0d0*FFY(U,V)
C
C
      RETURN
      END
