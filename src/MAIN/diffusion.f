C*************************************************************
C This routine was written by Anne A. Thoul, at the Institute
C for Advanced Study, Princeton, NJ 08540.
C See Thoul et al., Ap.J. 421, p. 828 (1994)
C The subroutines LUBKSB and LUDCMP are from Numerical Recipes.
C*************************************************************
C This routine inverses the burgers equations.
C
C The system contains N equations with N unknowns. 
C The equations are: the M momentum equations, 
C                    the M energy equations, 
C                    two constraints: the current neutrality 
C                                     the zero fluid velocity.
C The unknowns are: the M diffusion velocities,
C                   the M heat fluxes,
C                   the electric field E
C                   the gravitational force g.
C
C**************************************************
      SUBROUTINE DIFFUSION(M,A,Z,X,CL,AP,AT,AX)

C The parameter M is the number of species considered.
C
C Fluid 1 is the hydrogen
C Fluid 2 is the helium
C Fluids 3 to M-1 are the heavy elements
C Fluid M is the electrons
C
C The vectors A,Z and X contain the atomic mass numbers, 
C the charges (ionization), and the mass fractions, of the elements.
C NOTE: Since M is the electron fluid, its mass and charge must be
C      A(M)=m_e/m_u
C      Z(M)=-1.
C
C The array CL contains the values of the Coulomb Logarithms.
C The vector AP, AT, and array AX contains the results for the diffusion 
C coefficients.

      IMPLICIT NONE

      INTEGER*4 M,N,I,J,L,MMAX,NMAX
      PARAMETER (MMAX=20,NMAX=42)
      INTEGER*4 INDX(NMAX)
      REAL*8 A(M),Z(M),X(M),AP(M),AT(M),AX(M,M),CL(M,M)
      REAL*8 C(MMAX),CC,AC,XX(MMAX,MMAX),Y(MMAX,MMAX),YY(MMAX,MMAX),
     $     K(MMAX,MMAX)
      REAL*8 ALPHA(NMAX),NU(NMAX),GAMMA(NMAX,NMAX),DELTA(NMAX,NMAX),
     $     GA(NMAX)
      REAL*8 TEMP,KO,D

C The vector C contains the concentrations
C CC is the total concentration: CC=sum(C_s)
C AC is proportional to the mass density: AC=sum(A_s C_s)
C The arrays XX,Y,YY and K are various parameters which appear in 
C Burgers equations.
C The vectors and arrays ALPHA, NU, GAMMA, DELTA, and GA represent
C the "right- and left-hand-sides" of Burgers equations, and later 
C the diffusion coefficients.
      
C Initialize parameters:

      KO=2.0d0
      N=2*M+2
      DO I=1,M
         C(I)=0.0d0
      ENDDO
      CC=0.0d0
      AC=0.0d0
       
C Calculate concentrations from mass fractions:

      TEMP=0.0d0
      DO I=1,M-1
         TEMP=TEMP+Z(I)*X(I)/A(I)
      ENDDO
      DO I=1,M-1
         C(I)=X(I)/A(I)/TEMP
      ENDDO
      C(M)=1.0d0

C Calculate CC and AC:
         
      DO I=1,M
         CC=CC+C(I)
         AC=AC+A(I)*C(I)
      ENDDO

C Calculate the mass fraction of electrons:

      X(M)=A(M)/AC

C Calculate the coefficients of the burgers equations

      DO I=1,M
         DO J=1,M
            XX(I,J)=A(J)/(A(I)+A(J))
            Y(I,J)=A(I)/(A(I)+A(J))
            YY(I,J)=3.0d0*Y(I,J)+1.3*XX(I,J)*A(J)/A(I)
            K(I,J)=1.0d0*CL(I,J)*
     $           SQRT(A(I)*A(J)/(A(I)+A(J)))*C(I)*C(J)*
     $           Z(I)**2*Z(J)**2
         ENDDO
      ENDDO

C Write the burgers equations and the two constraints as
C alpha_s dp + nu_s dT + sum_t(not 2 or M) gamma_st dC_t 
C                     = sum_t delta_st w_t

      DO I=1,M
         ALPHA(I)=C(I)/CC
         NU(I)=0.0d0
         DO J=1,M
            GAMMA(I,J)=0.0d0
         ENDDO
         DO J=1,M
            IF ((J.NE.2).AND.(J.NE.M)) THEN
               GAMMA(I,J)=-C(J)/CC+C(2)/CC*Z(J)*C(J)/Z(2)/C(2)
               IF (J.EQ.I) THEN
                  GAMMA(I,J)=GAMMA(I,J)+1.0d0
               ENDIF
               IF (I.EQ.2) THEN
                  GAMMA(I,J)=GAMMA(I,J)-Z(J)*C(J)/Z(2)/C(2)
               ENDIF
               GAMMA(I,J)=GAMMA(I,J)*C(I)/CC
            ENDIF
         ENDDO

         DO J=M+1,N
            GAMMA(I,J)=0.0d0
         ENDDO
      ENDDO
      
      DO I=M+1,N-2
         ALPHA(I)=0.0d0
         NU(I)=2.5d0*C(I-M)/CC
         DO J=1,N
            GAMMA(I,J)=0.0d0
         ENDDO
      ENDDO
      
      ALPHA(N-1)=0.0d0
      NU(N-1)=0.0d0
      DO J=1,N
         GAMMA(N-1,J)=0.0d0
      ENDDO
      
      ALPHA(N)=0.0d0
      NU(N)=0.0d0
      DO J=1,N
         GAMMA(N,J)=0.0d0
      ENDDO
      
      DO I=1,N
         DO J=1,N
            DELTA(I,J)=0.0d0
         ENDDO
      ENDDO
      
      DO I=1,M
         DO J=1,M
            IF (J.EQ.I) THEN
               DO L=1,M
                  IF(L.NE.I) THEN
                     DELTA(I,J)=DELTA(I,J)-K(I,L)
                  ENDIF
               ENDDO
            ELSE
               DELTA(I,J)=K(I,J)
            ENDIF
         ENDDO
         
         DO J=M+1,N-2
            IF(J-M.EQ.I) THEN
               DO L=1,M
                  IF (L.NE.I) THEN
                     DELTA(I,J)=DELTA(I,J)+0.6d0*XX(I,L)*K(I,L)
                  ENDIF
               ENDDO
            ELSE
               DELTA(I,J)=-0.6d0*Y(I,J-M)*K(I,J-M)
            ENDIF
         ENDDO
         
         DELTA(I,N-1)=C(I)*Z(I)
         
         DELTA(I,N)=-C(I)*A(I)
      ENDDO
      
      DO I=M+1,N-2
         DO J=1,M
            IF (J.EQ.I-M) THEN
               DO L=1,M
                  IF (L.NE.I-M) THEN
                     DELTA(I,J)=DELTA(I,J)+1.5d0*XX(I-M,L)*K(I-M,L)
                  ENDIF
               ENDDO
            ELSE
               DELTA(I,J)=-1.5d0*XX(I-M,J)*K(I-M,J)
            ENDIF
         ENDDO
         
         DO J=M+1,N-2
            IF (J-M.EQ.I-M) THEN
               DO L=1,M
                  IF (L.NE.I-M) THEN
                     DELTA(I,J)=DELTA(I,J)-Y(I-M,L)*K(I-M,L)*
     $                    (1.6d0*XX(I-M,L)+YY(I-M,L))
                  ENDIF
               ENDDO
               DELTA(I,J)=DELTA(I,J)-0.8d0*K(I-M,I-M)
            ELSE
               DELTA(I,J)=2.7d0*K(I-M,J-M)*XX(I-M,J-M)*Y(I-M,J-M)
            ENDIF
         ENDDO
         
         DELTA(I,N-1)=0.0d0
         
         DELTA(I,N)=0.0d0
      ENDDO
      
      DO J=1,M
         DELTA(N-1,J)=C(J)*Z(J)
      ENDDO
      DO J=M+1,N
         DELTA(N-1,J)=0.0d0
      ENDDO
      
      DO J=1,M
         DELTA(N,J)=C(J)*A(J)
      ENDDO
      DO J=M+1,N
         DELTA(N,J)=0.0d0
      ENDDO




C Inverse the system for each possible right-hand-side, i.e.,
C if alpha is the r.h.s., we obtain the coefficient A_p
C if nu    ---------------------------------------- A_T
C if gamma(i,j) ----------------------------------- A_Cj
C 
C If I=1, we obtain the hydrogen diffusion velocity
C If I=2, ------------- helium   ------------------
C If I=3,M-1, --------- heavy element -------------
C If I=M, ------------- electrons -----------------
C For I=M,2M, we get the heat fluxes
C For I=N-1, we get the electric field
C For I=N, we get the gravitational force g

      CALL LUDCMP(DELTA,N,NMAX,INDX,D)
      
      CALL LUBKSB(DELTA,N,NMAX,INDX,ALPHA)
      CALL LUBKSB(DELTA,N,NMAX,INDX,NU)
      DO J=1,N
         DO I=1,N
            GA(I)=GAMMA(I,J)
         ENDDO
         CALL LUBKSB(DELTA,N,NMAX,INDX,GA)
         DO I=1,N
            GAMMA(I,J)=GA(I)
         ENDDO
      ENDDO

C The results for the coefficients must be multiplied by p/K_0:

      DO I=1,M
         ALPHA(I)=ALPHA(I)*KO*AC*CC
         NU(I)=NU(I)*KO*AC*CC
         DO J=1,M
            GAMMA(I,J)=GAMMA(I,J)*KO*AC*CC
         ENDDO
      ENDDO

      DO I=1,M
         AP(I)=ALPHA(I)
         AT(I)=NU(I)
         DO J=1,M
            AX(I,J)=GAMMA(I,J)
         ENDDO
      ENDDO
      

      RETURN
      
      END                                                             

*********************************************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER*4 N,NP
C     ..
C     .. Array Arguments ..
      REAL*8 A(NP,NP),B(N)
      INTEGER*4 INDX(N)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      INTEGER*4 I,II,J,LL
C     ..
      II = 0
      DO 12 I = 1,N
          LL = INDX(I)
          SUM = B(LL)
          B(LL) = B(I)
          IF (II.NE.0) THEN
              DO 11 J = II,I - 1
                  SUM = SUM - A(I,J)*B(J)
   11         CONTINUE

          ELSE IF (SUM.NE.0.0d0) THEN
              II = I
          END IF

          B(I) = SUM
   12 CONTINUE
      DO 14 I = N,1,-1
          SUM = B(I)
          IF (I.LT.N) THEN
              DO 13 J = I + 1,N
                  SUM = SUM - A(I,J)*B(J)
   13         CONTINUE
          END IF

          B(I) = SUM/A(I,I)
   14 CONTINUE
      RETURN

      END

*********************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      IMPLICIT NONE

C     .. Parameters ..
      INTEGER*4 NMAX
      REAL*8 TINY
      PARAMETER (NMAX=177,TINY=1.0D-40)
C     ..
C     .. Scalar Arguments ..
      REAL*8 D
      INTEGER*4 N,NP
C     ..
C     .. Array Arguments ..
      REAL*8 A(NP,NP)
      INTEGER*4 INDX(N)
C     ..
C     .. Local Scalars ..
      REAL*8 AAMAX,DUM,SUM
      INTEGER*4 I,IMAX,J,K
C     ..
C     .. Local Arrays ..
      REAL*8 VV(NMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..

c..to avoid g77 warning (wda 9/9/04)
      imax = 0
      D = 1.0d0
      DO 12 I = 1,N
          AAMAX = 0.0d0
          DO 11 J = 1,N
              IF (ABS(A(I,J)).GT.AAMAX) AAMAX = ABS(A(I,J))
   11     CONTINUE
          IF (AAMAX.EQ.0.0d0)then
             write(*,*)'equation ',i,' of ',n
             stop 'LUDCMP: Singular matrix.'
          endif
          VV(I) = 1.0d0/AAMAX
   12 CONTINUE
      DO 19 J = 1,N
          IF (J.GT.1) THEN
              DO 14 I = 1,J - 1
                  SUM = A(I,J)
                  IF (I.GT.1) THEN
                      DO 13 K = 1,I - 1
                          SUM = SUM - A(I,K)*A(K,J)
   13                 CONTINUE
                      A(I,J) = SUM
                  END IF

   14         CONTINUE
          END IF

          AAMAX = 0.0d0
          DO 16 I = J,N
              SUM = A(I,J)
              IF (J.GT.1) THEN
                  DO 15 K = 1,J - 1
                      SUM = SUM - A(I,K)*A(K,J)
   15             CONTINUE
                  A(I,J) = SUM
              END IF

              DUM = VV(I)*ABS(SUM)
              IF (DUM.GE.AAMAX) THEN
                  IMAX = I
                  AAMAX = DUM
              END IF

   16     CONTINUE
          IF (J.NE.IMAX) THEN
              DO 17 K = 1,N
                  DUM = A(IMAX,K)
                  A(IMAX,K) = A(J,K)
                  A(J,K) = DUM
   17         CONTINUE
              D = -D
              VV(IMAX) = VV(J)
          END IF

          INDX(J) = IMAX
          IF (J.NE.N) THEN
              IF (A(J,J).EQ.0.0d0) A(J,J) = TINY
              DUM = 1.0d0/A(J,J)
              DO 18 I = J + 1,N
                  A(I,J) = A(I,J)*DUM
   18         CONTINUE
          END IF

   19 CONTINUE
      IF (A(N,N).EQ.0.0d0) A(N,N) = TINY
      RETURN

      END










