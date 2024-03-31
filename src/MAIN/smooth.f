      SUBROUTINE SMOOTH
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
C
C  THIS SUBROUTINE USES A 2-DIMENSIONAL GENERALISATION OF THE SMOOTHING 
C  TECHNIQUES DESCRIBED ON PP. 644 TO 649 OF Numerical Recipes.
C
C  CONSIDER THE 25 POINTS DEFINED BY
C       I+n, n=-2,-1,0,1,2 AND J+m, m=-2,-1,0,1,2.
C  THE FUNCTION TO BE SMOOTHED IS FITTED TO A BI-CUBIC, INVOLVING
C  16 COEFFICIENTS, USING TECHNIQUES OF LEAST-SQUARES. THE SMOOTHED
C  FUNCTION (TEMPORARILY STORED IN FXY) IS GIVEN BY THE FITTED VALUE
C  AT THE POINT I AND J.
C
C  THE FITTING IS SHIFTED FOR POINTS CLOSE TO BOUNDARIES.
C
C
      PARAMETER(IPR=20)
C
c      COMMON/CST/NRL,RLS,nset,tmax  ! modified
ccccccccccccccccc
      COMMON/CST/RLS,tmax,NRL,nset

      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
C
      DIMENSION GAM(6)
      DATA GAM/+0.0073469388d0,-0.0293877551d0,-0.0416326531d0,
     +         +0.1175510204d0,+0.1665306122d0,+0.2359183673d0/
      DIMENSION BET(11)
      DATA BET/
     + -0.0048979592d0,-0.0661224490d0,-0.0293877551d0,+0.0195918367d0, 
     +  0.2644897959d0,+0.1175510204d0,-0.0783673469d0,+0.0277551020d0, 
     +  0.3746938776d0,+0.1665306122d0,-0.1110204082d0/
      DIMENSION ALP(11)
      DATA ALP/
     + -0.0844897959d0,-0.0048979592d0,+0.0073469388d0,+0.0012244898d0, 
     +  0.3379591837d0,+0.0195918367d0,-0.0293877551d0,+0.4787755102d0, 
     +  0.0277551020d0,-0.0416326531d0,-0.0069387755d0/
C
C
      DO 20 I=3,nset-2
C
         J=1
         FXY(I,J)=
     +    ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     +   +ALP(2)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J+3)+F(I+2,J+3)
     +          +F(I-1,J+4)+F(I+1,J+4) )
     +   +ALP(3)*( F(I-2,J+2)+F(I+2,J+2) )
     +   +ALP(4)*( F(I-2,J+4)+F(I+2,J+4) )
     +   +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     +   +ALP(6)*( F(I-1,J+1)+F(I+1,J+1)+F(I-1,J+3)+F(I+1,J+3) )
     +   +ALP(7)*( F(I-1,J+2)+F(I+1,J+2) )
     +   +ALP(8)*  F(I  ,J  )
     +   +ALP(9)*( F(I  ,J+1)+F(I  ,J+3) )
     +   +ALP(10)* F(I  ,J+2) +ALP(11)*F(I  ,J+4)
C
         J=2
         FXY(I,J)=
     +    BET(1)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J+3)+F(I+2,J+3) )
     +   +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     +   +BET(3)*( F(I-2,J+1)+F(I+2,J+1) )
     +   +BET(4)*( F(I-2,J+2)+F(I+2,J+2)+F(I-1,J-1)+F(I+1,J-1)
     +            +F(I-1,J+3)+F(I+1,J+3) )
     +   +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     +   +BET(6)*( F(I-1,J+1)+F(I+1,J+1) )
     +   +BET(7)*( F(I-1,J+2)+F(I+1,J+2) )
     +   +BET(8)*( F(I  ,J-1)+F(I  ,J+3) )
     +   +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J+1) +BET(11)*F(I  ,J+2)
C
         DO 10 J=3,NRL-2
            FXY(I,J)=
     +         GAM(1)*( F(I-2,J-2)+F(I-2,J+2)+F(I+2,J-2)+F(I+2,J+2) )
     +        +GAM(2)*( F(I-2,J+1)+F(I-2,J-1)+F(I-1,J-2)+F(I-1,J+2)
     +                 +F(I+1,J-2)+F(I+1,J+2)+F(I+2,J-1)+F(I+2,J+1) )
     +        +GAM(3)*( F(I-2,J  )+F(I+2,J  )+F(I  ,J-2)+F(I  ,J+2) )
     +        +GAM(4)*( F(I-1,J-1)+F(I-1,J+1)+F(I+1,J-1)+F(I+1,J+1) )
     +        +GAM(5)*( F(I-1,J  )+F(I  ,J-1)+F(I  ,J+1)+F(I+1,J  ) )
     +        +GAM(6)*  F(I  ,J  )
   10    CONTINUE
C
         J=NRL-1
         FXY(I,J)=
     +     BET(1)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J-3)+F(I+2,J-3) )
     +    +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     +    +BET(3)*( F(I-2,J-1)+F(I+2,J-1) )
     +    +BET(4)*( F(I-2,J-2)+F(I+2,J-2)+F(I-1,J+1)+F(I+1,J+1)
     +             +F(I-1,J-3)+F(I+1,J-3) )
     +    +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     +    +BET(6)*( F(I-1,J-1)+F(I+1,J-1) )
     +    +BET(7)*( F(I-1,J-2)+F(I+1,J-2) )
     +    +BET(8)*( F(I  ,J+1)+F(I  ,J-3) )
     +    +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J-1) +BET(11)*F(I  ,J-2)
C
         J=NRL
         FXY(I,J)=
     +     ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     +    +ALP(2)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J-3)+F(I+2,J-3)
     +             +F(I-1,J-4)+F(I+1,J-4) )
     +    +ALP(3)*( F(I-2,J-2)+F(I+2,J-2) )
     +    +ALP(4)*( F(I-2,J-4)+F(I+2,J-4) )
     +    +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     +    +ALP(6)*( F(I-1,J-1)+F(I+1,J-1)+F(I-1,J-3)+F(I+1,J-3) )
     +    +ALP(7)*( F(I-1,J-2)+F(I+1,J-2) )
     +    +ALP(8)*  F(I  ,J  )
     +    +ALP(9)*( F(I  ,J-1)+F(I  ,J-3) )
     +    +ALP(10)* F(I  ,J-2) +ALP(11)*F(I  ,J-4)
C
   20 CONTINUE
C
      DO 40 I=3,nset-2   ! modified
         DO 30 J=1,NRL
            F(I,J)=FXY(I,J)
   30    CONTINUE
   40 CONTINUE
C
      RETURN
      END

