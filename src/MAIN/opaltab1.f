      subroutine opaltab1

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

C
C  CODE FOR FITTING AND SMOOTHING OPAL DATA. ADAPTED FROM A CODE
C     WRITTEN BY MIKE SEATON(obtained june 1993)
C
C     OPAL DATA.
C     ASSUMES FIRST T6=0.006, LAST T6=10.OR 0.04). Depending on position
C     in the table. 
C     USES RECTANGULAR ARRAY FOR VARIABLES T6 AND LOG10(R)
C
C     (1) NSM=NUMBER OF PASSES THROUGH SMOOTHING FILTER.
C     USE OF NSM=1 OR 2 IS RECOMMENDED.
C     NO SMOOTHING WITH NSM=0
C     (2) RANGE FOR LOG10(R),
C     RLS=FIRST VALUE, RLE=LAST VALE
C     (RLS MUST BE FIRST VALUYE IN TABLE)
C
C  SUBROUTINE INTERP
C     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF
C               SUBROUTINE INTERP
C     WHICH GIVES LOG(ROSS) AND TWO FIRST DERIVATIVES FOR ANY
C     VALUES OF LOG(T) AND LOG(RHO). SEE BELOW FOR FURTHER
C     EXPLANATION.
C
C  OUTPUT FOR THE CASE OF NSM.GT.0.
C     INTERP IS USED TO OBTAIN SMOOTHED DATA INTERPOLATED
C     BACK TO THE ORIGINAL OPAL MESH. TWO FILES ARE WRITTEN.
C
C
C  THE SUBROUTINES SPLINE AND SPLINT ARE ADAPTED FROM THOSE GIVE BY
C  W.H. Press, S.A. Teulolsky, W.T. Vettering and B.P. Flannery,
C  "Numerical Recipes in FORTRAN", 2nd edn., 1992, C.U.P.
C  OTHER REFERENCES ARE MADE TO METHODS DESCRIBED IN THAT BOOK.
C
      PARAMETER(IP=100,IPR=20)
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      DIMENSION U(IP),ROSSL(IP,IPR),V(IP),V2(IP)
      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      CHARACTER*1 HEAD(100)
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
      COMMON/CST/RLS,tmax,NRL,nset
      common/alink/ N,NSM,nrlow,nrhigh,RLE,t6arr(100),xzff(100,nr)  
      LOGICAL IERR

C
      NRL=2*(RLE-RLS)+1
C
C     STORE LOG10(T) IN U AND LOG10(ROSS) IN ROSSL
C     CHECK FIRST VALUE OF T6
      T6=t6arr(1)
      do j=1,NRL
      ROSSL(1,j)=xzff(1,j)
      enddo

      if (abs(T6-.0056341325) .lt. 1.e-8) then
         U(1)=6.+LOG10(T6)
      ENDIF
C     SET ROSSL UP TO T6=t6arr(nset)
      I=1
    5 I=I+1
      T6=t6arr(I)
      do j=1,NRL
      ROSSL(I,j)=xzff(I,j)
      enddo
         U(I)=6+LOG10(T6)
         IF(T6.LT.tmax)GOTO 5
      N=I
      IF(N.GT.IP)THEN
         PRINT*,' REQUIRE PARAMETER IP OF AT LEAST ',N
         STOP
      ENDIF
C
C
C     DEFINE VARIABLES
C         X=20.0*(LOG10(T)-3.80)+1
C         Y=2.0*(LOG10(R)-RLS)+1
C     USE INDICES I=1 TO nset AND J=1 TO NRL
C     X AND Y ARE SUCH THAT, ON MESH-POINT (I,J), X=I AND Y=J
C     OBTAIN:-
C         F(I,J)=LOG10(ROSS)
C         FX(I,J)=dF/dX
C         FY(I,J)=dF/dY
C         FXY(I,J)=ddF/dXdY
C
C
C     FIRST GET F AND FX, INTERPOLATING FROM OPAL T6 TO
C     INTERVAL OF 0.05 IN LOG10(T).
      DO 40 J=1,NRL
C        FOR EACH LOG10(R), STORE LOG10(ROSS) IN V(I)
         DO 20 I=1,N
            V(I)=ROSSL(I,J)
   20    CONTINUE
C
C        GET FIRST DERIVATIVES AT END POINTS
C
C        GET SECOND DERIVATIVES FOR SPLINE FIT
         CALL SPLINE(U,V,N,V2)
C
C        INTERPOLATE TO LOG10(T)=FLT, FLT=3.8(0.05)8.0
         DO 30 I=1,nset ! modified
            FLT=3.75+0.05*I
            CALL SPLINT(U,V,N,V2,FLT,F(I,J),FX(I,J))
   30    CONTINUE
C
   40 CONTINUE
C
C
C  OPTION FOR SMOOTHING
      IF(NSM.GT.0)THEN
         DO 35 NS=1,NSM
            CALL SMOOTH
   35    CONTINUE
         CALL FITX
      ENDIF
C
C
C  GET FY AND FXY
      CALL FITY
C
C  THE ARRAYS F, FX, FY AND FXY ARE NOW STORED
C
C  CAN NOW DO INTERPOLATIONS USING
C       CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
C       INPUT IS FLT=LOG10(T), FLRHO=LOG10(RHO)
C       OUTPUT IS G=LOG10(ROSS)
C              DGDT=dG/d(LOG10(T))
C            DGDRHO=dG/d(LOG10(RHO))
C              IERR=.TRUE. IF INPUT FLT, FLRHO ARE OUT-OF-RANGE,
C                          ELSE IERR=.FALSE.
C
C INTERPOLATE BACK TO OPAL POINTS
      IF(NSM.GT.0)THEN
         do l=1,NRL
         xzff(1,l)=ROSSL(1,l)
         enddo

         DO 70 K=2,N
            FLT=U(K)
            DO 50 L=nrlow,nrhigh
               FLR=RLS+.5*(L-1)
               FLRHO=FLR-18.+3.*FLT
               CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
               IF(IERR)THEN
               ENDIF
               V(L)=G
   50       CONTINUE
            T6=t6arr(K)
            do l=nrlow,nrhigh
            xzff(K,l)=V(l)

            enddo

   70    CONTINUE
      ENDIF
C
C
 1000 FORMAT('  SMOOTHED OPAL DATA')
 1100 FORMAT('  OPAL DATA, (SMOOTHED-ORIGINAL)')
 2000 FORMAT(100A1)
 2222 FORMAT(F8.3,20F7.3)
 6000 FORMAT(/' FIRST T6=',1P,E10.3,', SHOULD BE 0.006')
 6003 FORMAT(/' !!! OUT-OF-RANGE !!!'/' FLT=',1P,E10.3,', FLRHO=',E10.3,
     + ', FLR=',E10.3)
C
      END
C
