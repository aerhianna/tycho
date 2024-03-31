      SUBROUTINE FITX
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
C
C  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.
C  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.
C
C
      PARAMETER(IPR=20)
      DIMENSION A(85),D(85)
C
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
ccccccccccccccccc
      COMMON/CST/RLS,tmax,NRL,nset

      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
C
      DO 30 J=1,NRL
         DO 10 I=1,nset ! modified
            A(I)=F(I,J)
   10    CONTINUE
         CALL GETD(A,nset,D,AP1,APN)  ! modified
         FX(1,J)=AP1
         FX(nset,J)=APN   ! modified
         DO 20 I=2,nset-1  ! modified
            FX(I,J)=-A(I)+A(I+1)-2.0d0*D(I)-D(I+1)
   20    CONTINUE
   30 CONTINUE
C
      RETURN
      END

