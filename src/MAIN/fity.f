      SUBROUTINE FITY
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
C
C  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS
C  FY AND FXY
C
C
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
ccccccccccccccccc
      COMMON/CST/RLS,tmax,NRL,nset

C
      PARAMETER(IPR=20)
      DIMENSION A(IPR),B(IPR),AD(IPR),BD(IPR)
      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
C
      DO 30 I=1,nset   ! modified
         DO 10 J=1,NRL
            A(J)=F(I,J)
            B(J)=FX(I,J)
   10    CONTINUE
C
         CALL GETD(A,NRL,AD,AP1,APN)
         CALL GETD(B,NRL,BD,BP1,BPN)
C
         FY(I,1)=AP1
         FY(I,NRL)=APN
         FXY(I,1)=BP1
         FXY(I,NRL)=BPN
         DO 20 J=2,NRL-1
            FY(I,J)= -A(J)+A(J+1)-2.0d0*AD(J)-AD(J+1)
            FXY(I,J)=-B(J)+B(J+1)-2.0d0*BD(J)-BD(J+1)
   20    CONTINUE
   30 CONTINUE
C
      RETURN
      END

