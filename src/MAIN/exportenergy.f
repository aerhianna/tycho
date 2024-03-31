c
c
c The following subroutine calculates the nuclear reaction rates with
c   special attention to detail required for calculating solar neutrino
c   fluxes. The neutrino fluxes are evaluated also in this subroutine.
c
c
c All physical quantities are cgs. The reaction rates are in 
c  gm^{-1}s^{-1}.
c
c  Last change: 
c   I made some purely cosmetic changes on 1/20/96.
c
c
c All numbers in this subroutine have been calculated by John Bahcall 
c  so that they agree with the modern numbers in Neutrino 
c  Astrophysics (1989) or much more recent results, as indicated in 
c  comments lines.
c
c
c The original structure of this  subroutine ia based upon the  Yale 
c  energy generation  subroutine.  However, I have made  changes and 
c  improvements in the physics, in the coding, in the input 
c  parameters, and in the reactions that are taken into account.
c  In practice, the only feature that was preserved from the original
c  subroutine is the interface to the external program. Also, I have 
c  added calculations of the neutrino fluxes, which is the reason
c  that I began this work. 
c
c  I include many comment statements in order to make the program  
c  transparent and simple to revise.  
c
c Some statements have been added to simplify the process of revising 
c  the subroutine as improved data become available.  For example, I have 
c  inserted explicitly the standard cross section factors in a way that 
c  it is easy to revise:
c
c  SStandard(I) = 1.0 when the Ith cross section has the standard value.
c  
c  The standard values of SStandard(I) used here are those that given in
c   Table 1 of Bahcall and Pinsonneault (1992).  Published in the 
c   Reviews of Modern Physics, 64,885, 1992.
c
c The nuclear energy release to the star from each reaction 
c  is taken from Bahcall and Ulrich, RMP 60, 297 (1988) and takes account 
c  accurately of neutrino energy loss. 
c
c I have added a new section at the end of the subroutine that 
c  calculates the solar neutrino fluxes at the Earth.  These fluxes
c  are in the units of cm^-2 sec^-2 per gm. To get the flux
c  from a shell, multiply by the mass of the shell in units of
c  grams.  The fluxes are in a common block, Fluxes.  I also calculate
c  the fictional neutrino fluxes associated with the He3 + He3 and with
c  the He3 + He4 rections; these fictional fluxes are useful diagnostics
c  of the solar model.
c
c Weakscreening is a parameter passed in the Flux common block. To obtain 
c  the Graboske et al. and Salpeter standard results, use: 
c  weakscreening = 0.03.  To investigate the effect of always using 
c  weak screening, use a value for weakscreening greater than unity, 
c  e. g., 30.  
c
c  The value of SStandard(17) for hep is taken from Carlson et al (1991)
c  Phys. Rev. C 44, 619. It corresponds to an S sub 0 =
c  1.3 E-20 keV-barns, a factor of 0.1625  smaller than indicated by the 
c  older measurements and calculations used in Neutrino Astrophysics.
c
c *****************************************
c Changing nuclear reaction cross sections.
c *****************************************
c
c Nuclear cross sections can be changed simply by inserting new numbers
c  for the Data values of SStandard(I), which are the ratios of the 
c  desired cross section factors to the values given in Tables 3.2 and 
c  3.4 of Neutrino Astrophysics (1989). If the value given in Neutrino
c  Astrophysics is used, then SStandard(I) = 1.0 .  To increase the cross
c  section for reaction K by a factor of two compared to the standard 
c  value, set SStandard(K) = 2.0 .  To determine which value of I goes 
c  with which reaction, see the section just below.
c   
c The energy derivatives enter in a form in which they are divided by 
c  the absolute values of the cross sections at zero energy.  Thus if 
c  the shape of the cross section extrapolation is unchanged and only
c  the intercept of S(E) at zero energy is changed, then no correction
c  need be made for the derivatives.  They are automatically scaled 
c  correctly.  The exact way that the derivatives eneter the 
c  calculations is described in the section labeled ``Defining the 
c  Q(I)' that is presented below.
c
c ************************************
c Identifying the reactions.
c ************************************
c
c  The value of J denotes which of the reactions the coefficients
c  refer to:  
c  J = 1, pp; J = 2, He3+He3; J = 3, He3+ He4; J =4, P + C12;  J = 5, 
c   p+C13; J = 6. p + N14; J = 7, p + O16.
c
c  Reactions J = 8, 13 are not relevant for the solar interior; they are
c   holdovers from the earlier Yale code.  
c  Reaction 14 is pep; reaction 15 is Be7 electron capture; reaction 16 
c   is  Be7 proton capture; reaction 17 is the hep reaction.
c
c   Reactions 14-17 were not included in the Yale energy generation
c   code, but they are most of the story for the solar neutrino problem. 
c
c   The branching of the N15 + p reactions is treated in a series of 
c    separate statements following the calculation of the Be7 + p 
c    reaction. See the definitions of F3 and F4.  If the cross-section
c    factors of the N15 + p reactions are revised, then the numerical 
c    coefficients must be changed in the definition of C12alpha and 
c    O16gamma.
c
c For Q1(I), ...,Q(5(I), I = 8 corresponds to the Be7 +  p reaction.
c  This assignment for I = 8 is only valid for the listed Q's and not
c  for other arrays in the program.
c
c IU is the shell number.
c
c The user should change the subroutine variable list to correspond to 
c  those quantities that the user wishes to extract.
c Note that the neutrino fluxes and the weakscreening parameter are 
c  currently passed via a Common statement.
c Two examples of the Subroutine labeling are given.
c  First example:
       Subroutine Energy(EPP1,EPP2,EPP3,ECN,E3AL,SUM1,DL,TL,
     $   X,Y,XHe3,XC12,XC13,XN14,
     $   XO16,jnbfkt,IU)
c  Second example:
c      SUBROUTINE Energy(EPP1,EPP2,EPP3,ECN,E3AL,PEP,PET,SUM1,DL,
c     *TL,PDT,PDP,X,Y,Z,XHE3,XC12,XC13,XN14,XN15,XO16,XO17,
c     *XO18,XH2,XLI6,XLI7,XBE9,IU,HR1,HR2,HR3,HR4,HR5,HR6,HR7,
c     *HR8,HR9,HR10,HR11,HR12,HR13,HF1,HF2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CTLIM/ATIME(11),TCUT(5),TSCUT,TENV0,TENV1,TENV,TGCUT,LNEUTR
      COMMON/CONST1/ CLN,CLNI,C4PI,C4PIL,C4PI3L,CC13,CC23
      Common/Fluxes/Flux(10),weakscreening
      DIMENSION HR1(1000),HR2(1000),HR3(1000),HR4(1000),HR5(1000),
     1 HR6(1000),HR7(1000),HR8(1000),HR9(1000),HR10(1000),HR11(1000),
     2 HR12(1000),HR13(1000),HF1(1000),HF2(1000)
      DIMENSION XFRAC(13),RATE(13),DRATRO(13),DRATT(13),UTOT(13),
     1 DSCR (13),DSCT(13),DG(13),ZPRD(13),Z53(13),Z43(13),Z23(13),
     2 Z86(13),Q1(8),Q2(8),Q3(8),Q4(8),Q5(8),Q6(7),Q7(7),Q8(7),EG(50)
      DIMENSION ANUC(13), ZNUC(13), SStandard(17)

      real*8 jnbfkt(13),jnbrate
c
      character*10 creac(17)
      character*5  cnuc(13)

      data cnuc / 'n','p','d','t','he3','he4','c12','c13','n14',
     1     'o16','o18','ne20','mg24'/
c
c ***************************
c ANUC are atomic mass units.  
c ***************************
c
c They (ANUC(I) ) include the electron masses as well as the nuclear 
c  masses.  
c  The scale is the mass of C12 divided by 12 or 931.49432 MeV,
c  which is 1.6605402 times 10^{-24} gm. The values are obtained
c  by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
c  adding to this the atomic mass number, A.  The value for Be7, which
c  is used implicitly in this subroutine is, 7.016930.
c
      DATA ANUC/1.008665,1.007825,2.014102,3.01605,3.01603,4.002603,12.,
     1   13.00335,14.00307,15.99491,17.99916,19.99244,23.98504/,ZNUC/0.,
     2   1.,1.,1.,2.,2.,6.,6.,7.,8.,8.,10.,12./,NELEM/13/
c 
c The isotopes are neutron, H1, D, H3, He3, He4, C12, C13, N14, O16, O18, 
c  Ne20, Mg24, respectively. All of these numbers were checked.
c Nelem is the number of isotopes included.
c
c *****************
c Defining the Q(I)
c *****************
c **************************************************************************
c The quantities Q1(J), Q2(J), ...,Q5(J) are the terms in Equation 3.14 of
c  Neutrino Astrophysics and in EQuation 53 of Fowler, Caughlan, and 
c  Zimmerman (1967) eq. 53, in both cases multiplied by T sub 9 ^(-2/3).
c  The reactions corresponding to each J are listed above, under:
c   Identifying the rections.
c  For this set of parameters, and only for this set of parameters,
c  J = 8 corresponds to the Be7 + p reaction.
c **************************************************************************
c
c General expression for Q's:
c
c T_9^(-2/3)[S_eff/S(0)] = 
c   [T_9^(-2/3) + Q1(I)T_9^(-1/3) + Q2(I) + Q3(I)T_9^(1/3) + 
c    Q4(I)T_9^(2/3) + Q(5)T_9 ]
c
c By comparison with Equation 3.14 of Neutrino Astrophysics, we see that 
c
c  Q1 = (5/(12*tau))*T_9^(-1/3) .
c  Q2 = (S'/S)(E_0))*T_9^(-2/3) .
c  Q3 = (S'/S)(35/36)(k*10^9 K)
c  Q4 = (S''/2S)(E_0^2)(T_9^(-4/3)
c  Q5 = (89/72)(S''/S)(E_0)(kT)(T_9^(-5/3)
c
c Each of the Q's is independent of temperature (T), as can be seen from
c  Equations 3.10 and 3.11 .
c
c All of the values of the Q1, ...., Q5 have been recalculated, using 
c  where needed nuclear cross sections given in Tables 3.2 and 3.4 of Neutrino 
c  Astrophysics.
c
c ******************************************************************
c Q6 is the coefficient of the temperature term in the definition of
c  tau, equation 3.10 of Neutrino Astrophysics. 
c  tau = -q6*(T sub 9 to the (-1/3) power ).
c ******************************************************************
c
c Slight changes have been made in the previous values of Q6 to make 
c  the data more accurate.
c Note that q6 is negative.
c
c *************************************************
c Q7 is the constant in front of the reaction rate. 
c *************************************************
c
c The general relation is:  
c  Q7 = 70.62860 + ((ln(Z0*Z1/A))/3) -ln(A0*A1) + ln(S sub 0)
c       -ln(1 + delta_01) 
c  Here S sub 0 is the cross section factor in units of keV-barns.
c  The numerical values used here are taken from Tables 3.4 and 3.2
c   of Nuclear Astrophysics.
c  The quantity delta_01 is non-zero (equal to unity) only when the 
c   two reacting nuclei, 0 and 1, are identical.
c
c Q8 reflects a term in the exponential that occurs in the rates, the
c  term being proportional to e^( constant*T_9^2). See Harris, et al. 
c  1983, Annual Rev. Astron. Astrophys.21, 165 (1983) for the meaning
c  of this obscure term.
c
c The values of Q1(I),...Q7(I) given below were obtained using the data
c  in Neutrino Astrophysics with the aid of auxilary computer codes
c  that generated accurate evaluations. The numbers given here are 
c  in many cases completely different from the values in the original
c  Yale code.
c
c NRXNS is the number of reactions being tracked.
c
      DATA Q1/0.1231,.0339,.0326,.0304,.0303,.0274,.02494,.040565/,
     1  Q2/1.079,-.0616,-.2115,.6645,.9601,-.7788,-1.173, -0.3640/,
     2  Q3/.9304,-.01464,-.04809,.14155,.2041,-.1491,-.205,-0.1034/,
     3  Q4/0.,0.,0.,3.627,1.394,.2612,.734,0.0/,
     4  Q5/0.,0.,0.,1.9647,.753,.12717,.3261,0.0/,
     5  Q6/-3.3796,-12.272,-12.822,-13.6869,-13.7142,-15.2247,-16.6888/,
     6  Q7/20.8951,76.5995,67.8028,69.1290,70.3799,69.8507,70.8002/,
     7  Q8/0.,0.,0.,1.0504,0.7874,0.090155,0./,
     8  NRXNS/13/,C6/2.302585/,
c
c December 27, 1993: cross section factors (includes improved 
c  calculation of p-p reaction, with better matrix nuclear element and 
c  vacuum polarization, by Kamionkowski and Bahcall (ApJ, Jan.10, 1994)
c  and vacuum polarization effects on the 3He-3He, 3He-4He, 7Be-p, and 
c  14N-p reactions by Kamionkowski and Bahcall (1994,Phys. Rev. C,49).
c  Use these to  get Bahcall-Pinsonneault 1995, Rev. Mod. Phys. Oct.
c  This is the latest version as of 7/1995.
c
     9  SStandard/0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0,1.0,1.0,
     $  1.0,1.0,1.0,1.0,1.0,0.92088,0.1625/
c
c
c May 6, 1992: cross section factors (used for Bahcall-Pinsonneault, Rev. 
c  Mod. Phys. 64, 885, 1992).
c
c     9  SStandard/0.9828,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
c     $  1.0,1.0,1.0,1.0,0.9218,0.1625/
c
c
c Cross-Section factors from Neutrino Astrophysics: For Checking only.
c
c     9  SStandard/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
c     $  1.0,1.0,1.0,1.0,1.0,1.0/
c
c *********************************************************************
c The values of SStandard(I) are to be changed from unity if the cross 
c  section factors are not the ones given in Neutrino Astrophysics,
c  Table 3.2 and Table 3.4 .
c
c *********************************************************************
c 
c ZPRD is used in screening calculations. It is the product of the 
c  charges of the interacting ions. ZPRD was checked. Z86 is used
c  in calculating intermediate screening and is defined by Graboske et 
c  al, Ap. J. 181, page 465 (1973), in Table 4.  Z86 was checked and 
c  some numerical values were made slightly more accurate.  Z53, Z43,
c  and Z23 are also defined in Table 4 (see above).  Since they are
c  only used in strong screening, the values of Z53, Z43, and Z23
c  were not checked.
c
      DATA ZPRD/1.,4.,4.,6.,6.,7.,8.,12.,16.,12.,14.,12.,36./,Z53/1.175,
     1   3.73,3.73,4.804,4.804,5.385,5.941,9.014,11.24,9.014,10.15,
     2   9.104,23.28/,Z43/0.52,1.31,1.31,1.488,1.488,1.61,1.721,2.577,
     3   3.025,2.577,2.81,2.577,5.668/,Z23/-0.413,-0.655,-0.655,-0.643,
     4   -0.643,-0.659,-0.673,-0.889,-0.946,-0.889,-0.92,-0.889,-1.36/,
     5   Z86/1.630,5.917,5.917,8.302,8.302,9.520,10.716,16.192,20.978,
     6   16.192,18.606,16.192,45.6635/,
     7   C21/5.240358E-8/
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i = 1, 17
         creac(i) = '          '
      enddo
      creac(1) = 'p+p'
      creac(2) = 'he3+he3'
      creac(3) = 'he3+he4'
      creac(4) = 'p+c12'
      creac(5) = 'p+c13'
      creac(6) = 'p+n14'
      creac(7) = 'p+o16'
ccccccccccccccccccccccccccccccccccccccccccccccc

c
c Define next the fractional abundances by mass of the important 
c  isotopes.
c
c X, Y, Z, XHe3,..., XBE9 are the mass fractions of the isotopes.
c  For the Sun, can ignore the abundance of neutron, H2, H3, Ne20,
c   and Mg24, which are, respectively, XFRAC(I) for I = 1,3,4,12,13.
c
      XFRAC(1) = 0.0
      XFRAC(2) = X
      XFRAC(3) = 0.0
      XFRAC(4) = 0.0
      XFRAC(5) = XHE3
      XFRAC(6) = Y
      XFRAC(7) = XC12
      XFRAC(8) = XC13
      XFRAC(9) = XN14
      XFRAC(10) = XO16
      XFRAC(11) = XO18
      XFRAC(12) = 0.0
      XFRAC(13) = 0.0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,'(/a30)')'EXPORTENERGY: '
      write(*,'(a5,a5,3a12)')'i','nuc','Z','A','X'
      do i = 1, 13
         write(*,'(i5,2x,a5,1p3e12.3)')i,cnuc(i),
     1        znuc(i),anuc(i),xfrac(i)
      enddo
      write(*,'(2(a10,1pe12.4))')'rho',10.0d0**dl,'T',10.0d0**tl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c *******************************************************************
c Begin calculation of screening correction.
c *******************************************************************
c
c  The basic references are Salpeter, Australian Journal of Physics, 
c  Vol. 7, 373 (1954). The formula for weak screening that is being 
c  programmed is Equation (25) of this paper. The other important 
c  references are: DeWitt, Graboske, and Cooper, Ap. J. 181, 439 (1973)
c  and Graboske et al., Ap. J. 181, 457 (1973). The values of EMU and
c  ZET are essential for computing weak screening; the value of AMU is
c  used in an non-essential way in this computation. XTR is used in 
c  computing intermediate screening.
c
c
c ANU is one over the mean molecular weight of the ions, mu sub I .
c EMU is one over the electron mean molecular weight, mu sub e.
c  EMU is used here as the name for the second part of the zeta function
c  in the Salpeter expression for weak screening.
c XTR is used later in the intermediate screening calculation.  The average 
c  of the quantity Z**(3b -1) is equal to XTR/AMU.  
c ZET is the first part of the Salpeter screening zeta variable.
c
c mu = sum over I of [X(I)/A(I)].
c mu sub e = sum over I [ Z(I)*X(I)/A(I)].
c 
c
      AMU = 0.
      EMU = 0.
      XTR = 0.
      ZET = 0.
      DO 10 I = 1,NELEM
         TRM = XFRAC(I)/ANUC(I)
         AMU = AMU+TRM
         EMU = EMU+TRM*ZNUC(I)
         XTR = XTR+TRM*ZNUC(I)**1.58
         ZET = ZET+TRM*ZNUC(I)**2
   10 CONTINUE
c
c
c DL and DT are the the log10 of the density and temperature. 
c  The unit of temperature is 10^9 K and the unit of density is
c  gm per cm^3 .
c PDT and PDP are the derivatives of the density with respect to 
c  temperature and density.
c
c DD = log rho to the base 10.
c c6 = ln10.  C6 is conversion between log10 and ln.
c Convert density to unlogged form.
c RWE = rho/(mu sub e), i. e., the number of electrons divided by 
c  Avogadro's number.
c
      RWE = ( DEXP(C6*DL) )*EMU
c
c RWE is used later in computing the screening correction.

cccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,'(5(a7,1pe14.5))')'Y',amu,'Y*Z',emu,'xtr',xtr,
     1     'Y*ZZ',zet,'rho*Ye',rwe
ccccccccccccccccccccccccccccccccccccccccccccccccc

c
      DD = DL
c
c Set rates equal to zero for the Log_10(T) < 6.0.
c 
      IF(TL.LE.6.0) THEN
         EN = -20.
         PEP = 0.
         PET = 0.
         DO 20 I = 1,NRXNS
            EG(I) = 0.
            RATE(I) = 0.
            DG(I) = 0.
   20    CONTINUE
         GO TO 200
      ENDIF
c
c T9P13 is the temperature in units of 10^9 degrees K to the Plus 1/3 
c  power.  Minus is denoted by M.  Here T9 is the temperature in units
c  of 10^9 K, converted from the log_10 (T) and Rho is the density in 
c  cgs units.
c
      RHO=DEXP(C6*DD)
      T9 = DEXP(C6*(TL - 9.0D0))
      T9P13 = T9**0.33333333333333D0
      T9P23 = T9P13**2
      T9M13=1./T9P13
      T9M23=T9M13**2
      T9M1=1./T9
      T9M2=T9M1**2
      T9M12=1./DSQRT(T9)
      T9M32=T9M1*T9M12
c
c ***********************************
c f prime/f
c ***********************************
c
c The next piece of code computes first the Fermi energy divided by kT, where
c  PFMC2 is the Fermi momentum divided by mass of electron times c, all
c  squared and EFMKT is the Fermi energy divided by kT.  The quantity 
c  FPRF is the ratio of f prime to f in Salpeter's screening correction.
c  The value of FPRF is determined in the intermediate case by an 
c  interpolation formula depending upon the degree of degeneracy, as 
c  measured by DEGD = Log10(E_F/kT). The numerical values for the fit
c  were taken from Salpeter's original paper, Figure 1. The only changes 
c  in this part of the subroutine were the correction of the error in 
c  the definition of RWE (see above) and refinements of the coefficients
c  in the expressions for PFMC2 and EFMKT.
c
      PFMC2=1.017677E-4*RWE**0.6666667
      EFMKT=5.92986*T9M1*(DSQRT(1.+PFMC2)-1.)
      IF(EFMKT.LE.1.E-2) THEN
         FPRF=1.0
      ELSE
         DEGD=DLOG10(EFMKT)
         IF(DEGD.GE.1.5) THEN
            FPRF=0.0
         ELSE
            FPRF=0.75793-0.54621*DEGD-0.30964*DEGD**2+0.12535*DEGD**3+
     *      0.1203*DEGD**4-0.012857*DEGD**5-0.014768*DEGD**6
         ENDIF
      ENDIF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,'(5(a7,1pe14.5))')'pf/mc2',pfmc2,'ef/kT',efmkt,
     1     'fpr/f',fprf
cccccccccccccccccccccccccccccccccccccccccccc

c 
c End of calculation of Fermi energy divided by kT and the interpolation
c  formula for (f prime/f), which appears in Salpeter's screening 
c  formula.
c
c Now we get to the computation of weak screening (see also page 61 of 
c  Neutrino Astrophysics, which gives only a simplified formula) and 
c  also the more complicated intermediate and strong screening cases 
c  (see references to Ap. J. 181 above for intermediate
c  and strong screening), especially page 465.
c
c XXL is used in all the screening formulae.  This quantity is the 
c  is the function called Lambda sub zero by Graboske et al. The
c  quantity XXL**0.86 = XXL8 is used in calculating intermediate 
c  screening.
c Zcurl is the quantity defined by equation (4) of DeWitt et al; it
c  is their Z with a curly symbol on its top. Zcurl is used in weak
c  and in intermediate screening and was first defined by Salpeter.
c  Zcurl is the same as Salpeter's zeta except for the factor of 1/AMU.
c
c (Z sub 1 times Z sub 2)*XXL*Zcurl gives the weak screening factor,
c  the same as Salpeter or as Equation (19) of Graboske et al.
c
c Z bar is the same as Z bar of De Witt et al and is the average charge
c  of the ions.  It is equal to EMU/AMU.  The quantity (Z bar)**0.28 =
c  Z28 occurs in the computation of intermediate screening.
c
c XXL6 is used for computing strong strong screening.  
c
c The notation used here is explained in large part by
c  equation (4) of Dewitt et al., Ap. J. 181, page 439.  
c The final expression for weak screening is exactly equal to Salpeter's
c  formula, which includes a degeneracy correction. The more general
c  expressions are given in Table 4 and equation (19) of Graboske et al.
c
      XXL=5.9426E-6*T9M32*DSQRT(RHO*AMU)
      XXL6=XXL**0.666667
      XXL8=XXL**0.86
      ZCURL=DSQRT((ZET+FPRF*EMU)/AMU)
      ZBAR=EMU/AMU
      Z58=ZCURL**0.58
      Z28=ZBAR**0.28
      Z33=ZBAR**0.333333
      TM1=XXL*ZCURL

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,'(5(a7,1pe14.5))')'xxl',xxl,'YZZ+YZ',zcurl,
     1        'tm1',tm1,'wflg',weakscreening
ccccccccccccccccccccccccccccccccccc

c 
c Compute screening for each of the reactions.
c
      DO 30 I=1,NRXNS
         UWK=TM1*ZPRD(I)
         IF(UWK.LE.weakscreening) THEN
c
c Weakscreening is a numerical parameter passed in the Flux common 
c  block. To obtain the Graboske et al. and Salpeter standard results,
c  use: weakscreening = 0.03.  For the standard solar model, this is the
c  value that should be adopted. To investigate the effect of always using 
c  weak screening, use a large value for weakscreening, e. g., 30.  As
c  long as weakscreening is assumed to be bigger than one, the program
c  will always calculate for the Sun with the weak screening 
c  approximation.
c 
c Utot is the final screening correction which appears in the
c  rate expression as: exp(Utot) .
c DSCR is the logarithmic derivative with respect to density of the
c  screening correction, d log (e^(U_tot)) /d log rho .
c DSCT is the logarithmic derivative of the screening with respect to T,
c d log (e^(U_tot)) /d log T .  The limiting formulae given in the first
c option are obvious since in the weak limit Salpeter's formula shows 
c that u is proportional to the square root of rho times T^(-3/2).
c
            UTOT(I)=UWK
            DSCR(I)=0.5*UWK
            DSCT(I)=-1.5*UWK
         ELSE
            UINT=0.38*XXL8*XTR*Z86(I)/(AMU*Z58*Z28)
            IF(UWK.LE.2.) THEN
               UTOT(I)=UINT
               DSCR(I)=0.43*UINT
               DSCT(I)=-1.29*UINT
            ELSE
               USTR=0.624*Z33*XXL6*(Z53(I)+0.316*Z33*Z43(I)+0.737*
     *         Z23(I)/(ZBAR*XXL6))
               IF(USTR.LT.UINT.OR.UWK.GE.5.) THEN
                  UTOT(I)=USTR
                  DSCR(I)=0.208*Z33*(Z53(I)+0.316*Z33*Z43(I))*XXL6
                  DSCT(I)=-3.*DSCR(I)
               ELSE
                  UTOT(I)=UINT
                  DSCR(I)=0.43*UINT
                  DSCT(I)=-1.29*UINT
               ENDIF
            ENDIF
         ENDIF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if( i .le. 7 )then
            write(*,'(a12,5(a7,1pe14.5))')creac(i),
     1           'Z1*Z2',zprd(i),'Uweak',uwk,'utot',utot(i),
     2           'f0',exp(utot(i))
         endif
cccccccccccccccccccccccccccccc

   30 CONTINUE

c
c ****************************************************************
c End of screening calculation. Weak and Intermediate screening forms
c  are given correctly.  Strong screening was not checked because it is
c  not relevant for the Sun.
c ****************************************************************
c
      NZ=1
      IF(X.EQ.0.0) THEN
         F1=0.
         F2=0.
         F3=0.
         F4=0.
         GOTO 50
      ENDIF
      NZ=8
c
c **************************************************************
c  Calculate reaction rates for the three principal rections of 
c   the pp chain: pp, He3+He3, He3 +He4, and the four proton
c   burning reactions on C12: C13, N14, and O16.
c **************************************************************
c
c R1 is (T sub 9)^(-3/2) times (S sub eff)/(S sub 0). The correct
c  expression for (S sub eff)/(S sub 0) is given in equation 3.14 of
c  Neutrino Astrophysics.  The numerical form that is used is equation
c  52 of Fowler, Cauhlan, and Zimmerman, Vol. 5, 1967.
c
c Rate(I) is the rate of the different reactions per second per gram,
c  except that the mass fractions are omitted at this point and put in
c  later.
c DRATT is logarithmic derivative of rate with respect to temperature,
c  d log Rate divided by d log T, log to base 10.
c DRATRO is the logarithmic derivative of the rate with respect to 
c  density, d log Rate/d log rho, log to base 10.
c
      DO 40 I=1,7
         R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
         RATE(I)=RHO*R1*DEXP(Q6(I)*T9M13+Q7(I)+(Q8(I)*T9)**2+UTOT(I))
         Rate(I) = Rate(I)*SStandard(I)
         IF(RATE(I).LT.1.E-30) THEN
            RATE(I)=0.
            DRATT(I)=0.
         ELSE
            DRATRO(I)=1.+DSCR(I)
            DRATT(I)=DSCT(I)-(Q6(I)*T9M13+(2.*T9M23+Q1(I)*T9M13-Q3(I)*
     *      T9P13-2.*Q4(I)*T9P23-3.*Q5(I)*T9)/R1)/3.+2.*(Q8(I)*T9)**2
         ENDIF
   40 CONTINUE
c
c ***************************************************************
c End of calculation of reaction rates for first 7 reactions.
c ***************************************************************
c
c ***********************************************
c Add a description of Be7 burning
c ***********************************************
c
c Calculate the burning of Be7 by protons (which produces the most 
c  experimentally accessible solar neutrinos) and the burning of Be7
c  by electron capture (the dominant process). 
c
c The electron capture rate in sec^(-1) is given by equation (3.18) of
c  Neutrino Astrophysics.  Can omit the factor of rho (in cgs units)
c  which also appears in the Be7 + proton rate. The Be7 + proton capture
c  rate is given by Table 3.2 and Equation (3.12).  Note that 1 over mu
c  sub e is equal to EMU.
c
c Will use the notation Be7 + proton = Be7proton and 
c  Be7 + e = Be7electron.
c
c F1 is the fraction of the Be7 that is burned by electron capture.
c F2 is the fraction of the Be7 that is burned by proton capture.
c F3 is the fraction of the N14 that is burned by p, alpha reaction.
c F4 is the fraction of the N15 that is burned by p, gamma reaction.
c See Table 21 of Bahcall and Ulrich (1988), Rev. Mod. Phys. 60.
c
         Be7electron = (1.752E-10)*T9M12*(1.0 + 0.004*(1000.*T9 - 16.))
         Be7electron = Be7electron*EMU*SStandard(15)
         Temp3 = (-10.26202*T9M13)
         Be7proton = (3.126571E+5)*X*SStandard(16)*Dexp(temp3)
c
c Include for Be7proton the T9M23 factor and all corrections proportional to 
c  Q1,...,Q5 from Equation 3.14 of Neutrino Astrophysics. These 
c  corrections are defined earlier in this subroutine.
c
         QRBe7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
     $          + Q4(8)*T9P23 + Q5(8)*T9
         Be7proton = Be7proton*QRBe7
c
c Calculate the screening correction for Be7 + p reaction.  Use weak and
c  intermediate screening formulae.
c
         ZPRDBe7p = 4.0
         Z86Be7p = 5.7790
         UWK = TM1*ZPRDBe7p
         IF(UWK.LE.weakscreening) THEN
            UTOTBe7p = UWK
         ELSE
            UINT = 0.38*XXL8*XTR*Z86Be7p/(AMU*Z58*Z28)
            UTOTBe7p = UINT
          End if
c
c End of calculation of screening correction for Be7 + p reaction.
c
         Be7proton = Be7proton*DEXP(UTOTBe7p)
c
c Multiply rates by factor of rho/[(atomic mass unit)*A(Be7)] to get 
c  in units of gm^{-1}.  Call factor cAMUBe7.
c
         cAMUBe7 = rho*8.582295E+22
         Be7proton = cAMUBe7*Be7proton
         Be7electron = cAMUBe7*Be7electron
c
         write(*,'(2(a20,1pe12.3))')'Be7+p',Be7proton,
     1        'Be7+e',Be7electron
ccccccccccccccccccccccccccccccccccccccccc

c
c End of multiplication inserted November 6, 1990.
c
         F1 = Be7electron/(Be7electron + Be7proton)
         F2 = Be7proton/(Be7electron + Be7proton)
c
c *********************************************************************
c End of calculation of crucial Be7 electron capture and proton capture
c  rates and their ratio.
c *********************************************************************
c
c ***************************
c N15 + p branching.
c ***************************
c
c The following statements compute the branching of N15(p,alpha)C12
c  and N15(p,gamma)O16. These statements replace outdated statements
c  in the Yale code. The CNO cross-section factors are from Table 3.4
c  of Neutrino Astrophysics.  The ratio of the reactions depends only 
c  upon the effective zero energy S-factor, which is S0(zero energy) 
c  times the combination of temperature and S-factor derivatives/S0
c  that was used previously as R1 in the rate calculations. The 
c  numerical coefficients that appear in the rate were represented 
c  by the Q1(J),...Q5(J) for the other reactions.
c  The Qvalues for the N15 reactions have been computed separately.
c
      O16gamma = T9M23 + 0.0272969*T9M13 + 0.205379 + 0.0392714*T9P13
     $          + 5.999033*T9P23 + 2.91700256*T9
c 
c Multiply by the value of S0 in kev-b.
c
      O16gamma = O16gamma*64.
c 
c Multiply by the value of S0 in kev-b.
c
      C12alpha = T9M23 + 0.0272969*T9M13 + 1.97164 + 0.3770055*T9P13
     $          + 13.65933*T9P23 + 6.641790*T9
      C12alpha = C12alpha*78000.         
      F3 = C12alpha/(C12alpha + O16gamma)
      F4 = 1. - F3
c
c End of new routine for the branching of N15 + p .
c

   50 DO 60 I=NZ,NRXNS
         RATE(I)=0.
         DRATRO(I)=0.
         DRATT(I)=0.
   60 CONTINUE
c
c  I have eliminated at this point the calculation of the alpha-capture rates
c   and the C12 + C12 rate.  These rates are certainly negligible in
c   the Sun and in general probably below 25 million degrees. Thus 
c   reactions 8-13 are omitted from this version of the subroutine.
c
c *******************
c EG(I)
c *******************
c
c Multiply the rates per gram, Rate(I), by the abundances of the 
c  reacting species by mass, to get the total rates per gram, EG.
c Only reactions 1-7 are important in the Sun.  The first three, I = 1,3
c  are the most important since they determine the rate of the pp cycle.
c  The next three, I = 4-6, are significant for the relatively rare CNO
c  cycle. The rate for reaction 7 is a correction to the CNO cycle.
c  The reactions, I = 8-13, are not important in the Sun and are 
c  omitted here.
c
  100 EG(1)=RATE(1)*X*X
      EG(2)=RATE(2)*XHE3*XHE3
      EG(3)=RATE(3)*XHE3*Y
      EG(4)=RATE(4)*X*XC12
      EG(5)=RATE(5)*X*XC13
      EG(6)=RATE(6)*X*XN14
      EG(7)=RATE(7)*X*XO16
c
c Set rates of reactions 8-13 equal to zero.
c
       Do Izero = 8,13
         EG(Izero) = 0.0
       End do    
c
c The following statements are left in the form of comments so that the 
c  rates of I = 8-13 are easily identified.  The statements are left 
c  over from the original Yale code.
c
c      EG(8)=RATE(8)*Y*XC13
c      EG(9)=RATE(9)*Y*X016
c      EG(10)=RATE(10)*Y*XC12
c      EG(11)=RATE(11)*Y*XN14
c      EG(12)=RATE(12)*Y**3
c      EG(13)=RATE(13)*XC12*XC12
c 
c ******************************************************************
c
c ****************************************
c Energy generation.
c ****************************************
c
c Calculate energy generation by multiplying rates per gram per sec by 
c  the energy release.  The energies are taken from Table 21 of Bahcall and 
c  Ulrich (1988), Rev. Mod. Phys. 60, 297. This table is based upon a careful
c  calculation of the average amount of energy loss by neutrinos for 
c  each reaction. The numbers for the C12 + p reaction sequence and the
c  C13 + p reaction are broken down separately for this 
c  subroutine.
c
c The final numbers are in erg per gm per second.
c
c Define the constant to convert MeV's to ergs. The numbers that appear 
c  are in MeV so they can be easily identified.
c
      convert = 1.602177E-6
c
c The multiplying constants below are in MeV.
c
      DG(1)=EG(1)*6.664*convert
      DG(2)=EG(2)*12.860*convert
      DG(3)=EG(3)*(1.586+F1*17.394+F2*11.499)*convert
      DG(4)=EG(4)*3.457372*convert
      DG(5)=EG(5)*7.550628*convert
      DG(6)=EG(6)*(9.054+F3*4.966+F4*12.128)*convert
      DG(7)=EG(7)*3.553*convert
c
c *******************************************************************
c End of calculation of energy release.
c *******************************************************************
c
c Set to zero all of the alpha-particle and C12-burning reactions.
c
      DRATRO(8) = 0.0
      DRATT(8) = 0.0
      DG(8) = 0.0
      DRATRO(9) = 0.0
      DRATT(9) = 0.0
      DG(9) = 0.0
      DRATRO(10) = 0.0
      DRATT(10) = 0.0
      DG(10) = 0.0
      DRATRO(11) = 0.0
      DRATT(11) = 0.0
      DG(11) = 0.0
      DRATRO(12) = 0.0
      DRATT(12) = 0.0
      DG(12) = 0.0
      DRATRO(13) = 0.0
      DRATT(13) = 0.0
      DG(13) = 0.0
c
c End of xeroing out of reactions 8-13.
c
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      DO 110 I=1,NRXNS
C 
c *******************************************************************
c Sum of the total energy generation in ergs per grm per second with 
c derivatives with respect to density and to temperature.
c *******************************************************************
c
c Sum1 = sum of all energy generation rates. Note that the burning of
c  Be7 is included in DG(3) above.
c Sum2 = Sum over I of DG(I)* [d log Rate(I) / d log rho ].
c Sum3 = Sum over I of DG(I)* [d log Rate(I) / d log T ].
c
         SUM1=SUM1+DG(I)
         SUM2=SUM2+DG(I)*DRATRO(I)
         SUM3=SUM3+DG(I)*DRATT(I)
  110 CONTINUE
      IF(SUM1.LE.1.E-12) THEN
         EN=-20.
         PEP=0.
         PET=0.
         DO 120 I=1,NRXNS
            EG(I)=0.
  120    CONTINUE
      ELSE
c
c ******************************************************
c Global quantities that are returned by the subroutine.
c ******************************************************
c
         EN=DLOG10(SUM1)
         PEPD=SUM2/SUM1
         PEP=PDP*PEPD
         PET=SUM3/SUM1+PDT*PEPD
      ENDIF
c
c PDP = d Log rho/ d log P; PDT = d log rho/ d log T.
c
c *****************************************************
c End of computation of the global quantities.
c *****************************************************
c
      DO 130 I=1,NRXNS
         IF(RATE(I).LE.1.E-5) RATE(I) = 0.0
  130 CONTINUE
c
c ******************************************************
c Rates per 10^9 years per atomic mass unit: HRk(IU)
c ******************************************************
c
c HR1, ..., HR13 are the rates of the individual reactions.
c  The interpretation of which reaction goes with which symbol can be 
c  made easily by looking at the definitions of the EG(I)'s.
c  The abundances are updated in subroutine Kemcom using these matrices.
c
c C21 is the product of (10^9 years/1 second)*(1 atomic mass unit/1 
c  gram). I have used here sidereal year in converting to seconds.
c
  200 HR1(IU)=RATE(1)*C21
      HR2(IU)=RATE(2)*C21
      HR3(IU)=RATE(3)*C21
      HR4(IU)=RATE(4)*C21
      HR5(IU)=RATE(5)*C21
      HR6(IU)=RATE(6)*C21
      HR7(IU)=RATE(7)*C21
      HR8(IU)=RATE(8)*C21
      HR9(IU)=RATE(9)*C21
      HR10(IU)=RATE(10)*C21
      HR11(IU)=RATE(11)*C21
      HR12(IU)=RATE(12)*C21
      HR13(IU)=RATE(13)*C21
      HF1(IU)=F3
      HF2(IU)=F1
c
c ****************************************
c End of computation of HRk(IU).
c ****************************************
c
c ****************************************
c Calculating the total energy generation.
c ****************************************
c 
c The summation of the energies is given in Table 21 of Neutrino
c  Astrophysics.
c
c
c EPP1 includes the energy generated by the pp reaction, by the H2 + p
c  reaction, and by the He3 + He3 reaction.  See Table 21 of Neutrino
c  Astrophysics.
c
      EPP1 = DG(1)+DG(2)
c
c EPP2 includes the energy generated by the He3 + He4 reaction and by 
c  the burning of Be7 through electron capture.
c
      EPP2 = EG(3)*(1.586 + F1*17.394)*convert
c
c EPP3 includes the energy generated by the He3 + He4 reaction and by 
c  the burning of Be7 through proton capture.
c
      EPP3 = EG(3)*(1.586 + F2*11.499)*convert
c
c ECN is the energy generated through the CNO cycle.
c
      ECN=DG(4)+DG(5)+DG(6)+DG(7)
c
c E3AL is the energy generated through the alpha-burning reactions and
c  is negligible for the Sun.
c
      E3AL=DG(8)+DG(9)+DG(10)+DG(11)+DG(12)
c
c     WRITE (9,999) X,ECN,Y,E3AL,T,DD,IU
c999  FORMAT(1X,' ENGE ',6(1PE11.3),I3)
c
c PEP and PET are the derivatives of the total energy generation rate 
c  with respect to pressure and temperature.
c

      write(*,'(3(a10,0pf12.5))')'epp1',epp1,'epp2',epp2,'epp3',epp3,
     1     'ecn',ecn,'etot',sum1

      do i = 1,7
         write(*,'(i5,2x,a10,1pe12.3)')i,creac(i),dg(i)
      enddo


c ****************************************************************
c Calculation of Neutrino Fluxes
c ****************************************************************
c
c This part of the subroutine calculates the neutrinos fluxes in 
c  number per gram per square centimeter per second at the Earth's surface
c  (assuming nothing happens to the neutrinos after they are created).
c See Tables 3.1 and 3.3 of Neutrinos Astrophysics or equations 6.1-6.8
c  for the reactions. The order of the reactions is the same as in 
c  equations 6.1-6.8 .
c
c Define 4*pi*(AU)**2 .
c
      fourpiAU2 = 2.812295E+27            
c
c Flux of pp neutrinos.
c
      Flux(1) = Eg(1)/fourpiAU2
c
c Flux of pep neutrinos. Use Equation 3.17 of Neutrino Astrophysics.
c
      Flux(2) = (3.4848E-6)*RWE*T9M12*(1.0 + 20.*T9)*Eg(1)
      Flux(2) = Flux(2)*SStandard(14)/fourpiAU2
c
c Flux of hep neutrinos.  Use Equation 3.12 directly.
c
      Q6hep = -6.1399
c
c Q6 is the negative of the coefficient of T9M13 in tau, equation 3.10.
c
      Flux(3) = (1.71724E+11)*Rho*T9M23*DEXP(Q6hep*T9M13)
c
c The derivatives of the cross section factor are not known and are 
c  taken to be zero.  The only term from equation 3.14 that survives
c  is 5/(12*tau).
c
      Flux(3) = (1.0 + 0.067862*T9P13)*SStandard(17)*Flux(3)
c
c
c Calculate weak or intermediate screening for hep neutrinos. 
c
         ZPRDHe3p = 2.0
         Z86He3p = 3.08687
         UWK = TM1*ZPRDHe3p
         IF(UWK.LE.weakscreening) THEN
            UTOTHe3p = UWK
         ELSE
            UINT = 0.38*XXL8*XTR*Z86He3p/(AMU*Z58*Z28)
            UTOTHe3p = UINT
          End if
c
c End of calculation of screening correction for He3 + p reaction.
c
      Flux(3) = Flux(3)*Dexp(UTOTHe3p)
      Flux(3) = Flux(3)*X*XHe3/fourpiAU2
c
c Compute Be7massfraction. This is not required for the neutrino
c  fluxes since Be7 is always in equilibrium with the slower production
c  rate of He3 + He4.  However, it is of interest in some applications
c  to know the Be7 mass fraction, so I compute it here and it can be
c  extracted with a common statement if desired.
c
      Be7massfraction = EG(3)/(Be7proton + Be7electron)
c
c End of November 6, 1990  addition.
c
c
c Flux of Be7 neutrinos.
c
      Flux(4) = Eg(3)*F1/fourpiAU2
c
c Flux of B8 neutrinos.
c
      Flux(5) = Eg(3)*F2/fourpiAU2
c
c Flux of N13 neutrinos.
c
      Flux(6) = Eg(4)/fourpiAU2
c
c Flux of O15 neutrinos.
c
      Flux(7) = Eg(6)/fourpiAU2
c
c Flux of F17 neutrinos.
c
      Flux(8) = Eg(7)/fourpiAU2
c
c Flux of fictional He3 + He3 neutrinos.
c
      Flux(9) = Eg(2)/fourpiAU2
c
c Flux of fictional He3 + He4 neutrinos.
c
      Flux(10) = Eg(3)/fourpiAU2
c
c
c End of Neutrino Flux routine.
c

      avogadro = 6.02214d23
      write(*,'(a5,9a13)')'i','reaction','screening',
     1     'rho NA sigv','jnb*A','noscr jnb*A','fkt','fkt/jnb*A','flow'
      do i = 1, 7
c..Yale definition uses mass fractions, so A0*A1 factor is in rate
         if( i .eq. 1 )then
            jnbrate = rate(i)/avogadro*anuc(2)*anuc(2)
         elseif( i .eq. 2 )then
            jnbrate = rate(i)/avogadro*anuc(5)*anuc(5)
         elseif( i .eq. 3 )then
            jnbrate = rate(i)/avogadro*anuc(6)*anuc(5)
         elseif( i .eq. 4 )then
            jnbrate = rate(i)/avogadro*anuc(7)*anuc(2)
         elseif( i .eq. 5 )then
            jnbrate = rate(i)/avogadro*anuc(8)*anuc(2)
         elseif( i .eq. 6 )then
            jnbrate = rate(i)/avogadro*anuc(9)*anuc(2)
         elseif( i .eq. 7 )then
            jnbrate = rate(i)/avogadro*anuc(10)*anuc(2)
         else
            jnbrate = rate(i)/avogadro
         endif
c..factor out screening to get raw NA sigma v
c         jnbrate = jnbrate * dexp( -utot(I) )
         write(*,'(i5,2x,a11,1p8e13.3)')i,creac(i),
     1        dexp( utot(i) ),rate(i)/avogadro, 
     2        jnbrate, jnbrate*dexp( -utot(I) ), jnbfkt(i),
     3        jnbfkt(i)/jnbrate,
     4        eg(i)/avogadro
      enddo

      write(*,'(a30/)')'leaving exportenergy.f'

      RETURN
      END
c **********************************************************************
c Comments included down to here.
c ***********************************************************************
