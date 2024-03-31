C**** Incorporates changes givXENX in Remark by Einarsson
C**** CACM 17(4) April 1974 p.225
C      REAL FUNCTION WEW_A ( XINPX, XENX )
      SUBROUTINE PRODLOG(XINPX,XOPX)

C      IMPLICIT NONE
C
C  ITERATIVE SOLUTION OF XINPX = W * EXP ( W ) WHERE XINPX IS GIVEN.
C  (NOVEMBER 1970)
C  (REVISED - SEPTEMBER 1971)
C  VERSION A -- CDC 6600 MACHINE ACCURACY.
C
C  INPUT PARAMETER:
C    XINPX  ARGUMXENXT OF W(XINPX)
C
C  OUTPUT PARAMETERS:
C    WEW  THE DESIRED SOLUTION.
C    XENX   THE LAST RELATIVE CORRECTION TO W(XINPX).
C
C  SET CONSTANTS...
C     .. Scalar ArgumXENXts ..
      REAL*8 XENX,XINPX
C     ..
C     .. Local Scalars ..
C**** Next 2 statement replaced in remark by Einarsson
C     REAL C1,C2,C3,C4,FLOGIN,TMPRY,TMPRY2,XOPX,Y,ZZNZ
C     INTEGER NEWE
      REAL*8 FLOGIN,TMPRY,TMPRY2,XOPX,YYY,ZZNZ
C**** End of replaced statements
C     ..
C     .. Intrinsic Functions ..
C      INTRINSIC LOG
C**** Next 10 statements deleted in remark by Einarsson
C     ..
C     .. Save statement ..
C      SAVE C1,C2,C3,C4,NEWE
C      DATA NEWE / 1 /
C      IF ( NEWE ) 10, 20, 10
C10    NEWE = 0
C      C1 = 4. / 3.
C      C2 = 7. / 3.
C      C3 = 5. / 6.
C      C4 = 2. / 3.
C**** End of replaced statements
C
C  COMPUTE INITIAL GUESS...
20    FLOGIN = LOG(XINPX)
      IF ( XINPX - 6.46d0 ) 30, 30, 40
C**** Next statement replaced in remark by Einarsson
C30    XOPX = XINPX * ( 1. + C1 * XINPX ) / ( 1. + XINPX * ( C2 + C3 * XINPX ) )
30    XOPX=XINPX*(3.0d0+4.0d0*XINPX)/
     1     (3.0d0+XINPX*(7.0d0+2.5d0*XINPX))
C**** End of replaced statements
      ZZNZ = FLOGIN - XOPX - LOG(XOPX)
      GO TO 50
40    XOPX = FLOGIN
      ZZNZ = -LOG(XOPX)
50    CONTINUE
C
C  ITERATION ONE...
      TMPRY = 1.d0 + XOPX
C**** Next statememt replaced in remark by Einarsson
C     YYY = 2. * TMPRY * ( TMPRY + C4 * ZZNZ ) - ZZNZ
      YYY=2.d0*TMPRY*(TMPRY+ZZNZ/1.5d0)-ZZNZ
C**** End of replaced statements
      XOPX=XOPX*(1.d0+ZZNZ*YYY/(TMPRY*(YYY-ZZNZ)))
C
C  ITERATION TWO...
      ZZNZ = FLOGIN - XOPX - LOG(XOPX)
      TMPRY = 1.d0 + XOPX
**** Next statement replaced in remark by Einarsson
C     TMPRY2 = TMPRY + C4 * ZZNZ
      TMPRY2 = TMPRY + ZZNZ/1.5d0
      XENX=ZZNZ*TMPRY2/(TMPRY*TMPRY2-0.5d0*ZZNZ)
      XOPX = XOPX * ( 1.d0 + XENX )
C
C  RETURN...
      RETURN
      END
