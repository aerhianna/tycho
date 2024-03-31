      SUBROUTINE RADWIND(XLOGL,STMASS,TEFF,YPS,XMDCM,VINF,ZMET)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-m)

c calculates mass-loss rates and terminal velocity as function
c of stellar parameters

c SECOND version with simplified treatment of helium ionization
c(not important)
c and first empirical fits to formce multipliers
c AND METALLICITY DEPENDENCE

c note: no attempt so far for smooth transitions in
c       wind properties when going from O- to early B-stars (30000 K)
c                                   to        mid B-        (23000 K)
c                                   to        A-supergiants (10000 K)

c Rolf Kudritzki, Feb. 5, 1999
c..modified for use with TYCHO by Dave Arnett, April 21, 1999

c PROGRAM TO CALCULATE MASS-LOSS RATE AND TERMINAL VELOCITY
c OF HOT STARS WINDS USING MODIFIED ANALYTICAL RADIATION DRIVEN
c WIND APPROACH BY KUDRITZKI ET AL. 1989
c AND EMPIRICAL (DEPTH INDEPENDENT) FORCE MULTIPLIERS ADJUSTED TO
c REPRODUCE OBSERVED WIND MOMENTUM - LUMINOSITY RELATIONSHIP AND
c TERMINAL VELOCITIES


c input XLOGL =  log L/L_sun
c       STMASS =  M/M_sun
c       TEFF =  T_eff
c       YPS =  N(He)/N(H)
c output XMDCM =  M_dot in g/s
c        VINF =  V_inf in cm/s
c        RST =  R/R_sun
c        GA =  Eddington Gamma

c theory does not work for T_eff < 7500 K

      IF ( TEFF .LT. 7500.0d0 ) THEN
        WRITE(*,*)' T_eff < 7500 K : theory not applicable'
        WRITE(*,*)' extension to lower T_eff requires'
        WRITE(*,*)' different approach'
        STOP
      ENDIF
c luminosity
      STLUM = 10.0d0**XLOGL
      TTT   = TEFF/5770.0d0
      TTT   = TTT*TTT*TTT*TTT
c radius in R_sun
      RST   = SQRT(STLUM/TTT)
c radius in centimeter
      RSTCM = RST*6.9599d10

c XIHE =  number of electrons provided by helium nucleus
c       we choose a very simple approach here
c       needs update using real photospheric ionization

      IF (TEFF .GE. 30000.0d0) THEN
c O_stars
         XIHE = 2.0d0
         GO TO 101
      ENDIF
      IF (TEFF .GE. 10000.0d0) THEN
c B-supergianst
         XIHE = 1.0d0
         GO TO 101
      ENDIF
c A-supergiants
      XIHE = 0.0d0

c force multiplier parameters
c     !!!!! this will need further refinement !!!

 101  BETA = 1.0d0
      IF (TEFF .GE. 30000.0d0) THEN
c O-stars
      ALPHA = 0.665d0
      DELTA = 0.05d0
      XK0   = 0.10d0
      GO TO 102
      ENDIF
      IF (TEFF .GE. 23000.0d0) THEN
c early B_supergiants
      ALPHA = 0.625d0
      DELTA = 0.075d0
      XK0   = 0.095d0
      GO TO 102
      ENDIF
      IF (TEFF .GE. 10000.0d0) THEN
c mid B_supergiants
      ALPHA = 0.575d0
      DELTA = 0.075d0
      XK0   = 0.05d0
      GO TO 102
      ENDIF

c A-supergiants
      ALPHA = 0.500d0
      DELTA = 0.002d0
      XK0   = 0.135d0

c redifining force multiplier k according to metallicity

 102  XMETFAC = ZMET**(1.0d0-ALPHA)
      XK      = XK0*XMETFAC
      DALPHA  = 1.0d0
      ALEFF   = ALPHA*DALPHA
      ALEFAC  = ALEFF/(1.0d0-ALEFF)

c electron scattering coefficient

      SIG = SIGMAE(YPS,XIHE)
c Eddington Gamma
      GA = EDGAMMA(STLUM,STMASS,SIG)
      IF (GA .GT. 1.0d0) THEN
        WRITE(*,*)' from RADWIND: we are above the Eddington limit'
        WRITE(*,*)' Gamma  = ',GA
c        STOP
      ENDIF

c escape velocity in km/s. NOTE: it includes correction for
c electron scatt.

      VESCKM = STMASS*(1.0d0-GA)/RST
      VESCKM = 617.6d0*SQRT(VESCKM)
c in cm/s
      VESC = VESCKM*1.0d5
      VESCQ = VESC*VESC
      CA = VESCQ*RSTCM/2.0d0
      VIERPI = 4.0d0*3.1415d0
      CNW = 1.5033d24*SIG/VIERPI/RSTCM/RSTCM
      CNW11 = CNW*1.0d-11
c isothermal sound speed in cm/s
      VSOUND = TEFF*(2.0d0+(1.0d0+XIHE)*YPS)/(1.0d0+4.0d0*YPS)
      VSOUND = 9.085d3*SQRT(VSOUND)
      VSSQ = VSOUND*VSOUND
      VSQUAR = VSSQ/VESCQ
c in km/s
      VSOUKM = VSOUND/1.0d5
c proton thermal velocity in cm/s
      VPROT = 1.2848d4*SQRT(TEFF)

c u_c = R/r_crit. critical point coordinate
c approximation according to Kudritzki et al. 89

      UCRIT = UC(ALPHA,DELTA,VSQUAR)

c velocity at critical point
C Kudritzki et al., 89 approximations for v_c
C            used as starting approximations
      VC = VAUCR(VSOUND,VESC,ALPHA)
      VCSQ = VC*VC

c now calculation of mass-loss rate
C modified Kudritzki et al., 89 approximations for M_dot

      BFAC  = 1.0d0-VSSQ/VCSQ
      CAFAC = 1.0d0-4.0d0*VSQUAR/UCRIT
      ASC   = CA*CAFAC
      YCRIT = ALPHA*ASC/(1.0d0-ALPHA)
      YCRIT = YCRIT/BFAC
      IAPPR = 1
c modified CAK approximation
      XMDCA =  XCAKNEW(XK,ALPHA,ALEFF,STLUM,SIG,CA,CAFAC,BFAC,VPROT,
     &  VIERPI)
c modified Kudritzki et al., 1989
      XMDCM = XFCNEW(XMDCA,UCRIT,YCRIT,VC,CNW11,ALPHA,DELTA,RSTCM,
     1 IAPPR)

C modified Kudritzki et al., 89 approximations for v_inf

      VFAC = CAFAC/BFAC
      XX = XINTAP(ALPHA,BETA,DELTA,UCRIT)
      XX = XX+VCSQ/ALEFAC/VESCQ/VFAC
C v_inf in cm/s
      VINF = VESC*SQRT(XX)
      VINF = VINF*SQRT(ALEFAC)*SQRT(VFAC)
      RETURN
      END


