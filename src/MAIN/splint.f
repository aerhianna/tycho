      SUBROUTINE SPLINT(XA,YA,N,Y2A,X,Y,YP)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.0d0) stop 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0d0
      YP=0.05d0*  (  (-YA(KLO)+YA(KHI))/H
     +   +      ( -(3*A**2-1)*Y2A(KLO)
     +            +(3*B**2-1)*Y2A(KHI) )*H/6.0d0 )
      RETURN
      END

