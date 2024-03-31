      FUNCTION EDGAMMA(STLUM,STMASS,SIG)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

C +++++ EDDINGTON LIMIT GAMMA
C ++++ SIG COMPUTED BY FUNCTION SIGMAE
      EDGAMMA = 7.6395d-5*SIG*STLUM/STMASS
      RETURN
      END


