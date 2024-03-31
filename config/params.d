'.....................................................................'
' TYCHO:                                                              '
'.....................................................................'
' Sequence identifier(keep this spacing).....' 'prefix'  'A1'
' newnet (revise network,0=no)...............' 'newnet'    0
' nbin (flag for binary=0,else formatted)....' 'nbin'      1
'.....................................................................'
' Stoptime:0=no,1.436d17=sun,1.9724e16=hyad..' 'stime'     4.0d17
' Stop: species gone at center(i=1-522,0=no).' 'ixstop'    0
' Stop: species at center below (Xi).........' 'xstop'     1.0d-6
' Number of timesteps for this run...........' 'll'        50000
' Pause (0=no,>0=each time step,<0=never)....' 'ipause'    0
' Cpu time in seconds for this run (0=ignore)' 'runtmx'    0.0
'.....................................................................'
' Flag to make solar model for test..........' 'nsoltest'  0
' Flag to freeze abundances(0=change)........' 'mset'      0	
' nopac (-1=thomson,0=opal2,1=opal1).........' 'nopac'     1
' scale opacity (0=no, else yes).............' 'nkscale'   0
' scale thoul diffusion (0=no,else = factor).' 'fthoul'    1.0d0
' (-1=complete ionization,else=compute)......' 'noion'     0
' (0=radiation pressure,else=no).............' 'norad'     0
' (0=H2 molecule,else=no)....................' 'nomole'    0
' (0=H2 coulomb plasma,else=no)..............' 'nocoul'    0
' nburn (-1=no,0=opath,1=explicit)...........' 'nburn'     1
' non-nuclear neutrino flag(-1=no,0=on)......' 'neutro'    0
' flag to force Bahcall rates(0=no)..........' 'jnb'       0
' flag to force OPAL EOS (0=no)..............' 'nopaleos'  1
' monitor iteration at k=istep(0=no effect)..' 'istep'     0
' izams (0=no,1=stop at zero age main seq.)..' 'izams'     0
' maximum envelope mass/total..0.010e0.......' 'fmxenv'    0.01d0
' minimum envelope mass/total...1e-3.........' 'fmnenv'    1.0d-3
' minimum temperature in envelope.....6e4....' 'tmenv'     6.0d4
' fraction for envelope derivatives..(.02)...' 'epsilon'   0.2d-1
' network: fractional change.................' 'delchi'    0.1d0
' network: minimum abundance used in dtime...' 'chimin'    1.0d-4
' network: fractional increase in dtime......' 'fdtn'      1.414d0
' network: fractional of nucleons switching..' 'fdysum'    4.0d-1
' network: subcycles allowed.................' 'ncymax'    200
' network: write if subcycles .ge.this value.' 'ncytest'   295
' rotation: no rotation=0, 1 and 2=yes.......' 'mrot'      2
' conv. velocity integration 0=yes (cmix.f)..' 'nsweep'    0
'.....................................................................'
' Debugging (0=normal,1=step,2=iteration)....' 'nbug'      0
' Stop at nup failed updates.................' 'nup'       150
' GR or Newtonian gravity.(0=GR,1=newton)....' 'newt'      0
' Constant outer radius (=1; ignore=0).......' 'nouter'    0
' mode(0=ex.hd;1=hs).........................' 'mode'      1
' Edit write interval........................' 'l2'        90000
' Model dump interval........................' 'l3'        100
' Number of iterations allowed...............' 'iter'      50
' Convection mode(-1=no,0=dS,1=Schw,2=L+mix).' 'modec'     3
' Mixing length/pressure scale height.1.643..' 'alphaml'   1.6d0
' MLT geometric parameter(relative to B-V)...' 'uuml'      1.0d0
' Overshoot parameter (1/ml).................' 'hmlfak'    0.0d0
' Mixing (0=PAY mix+diff;1=mix;2=diff;3=none)' 'mixmode'   2
' Surface mode(0=surf.phot.;1=mov/phot.;2=en)' 'modes'     2
' Zoning mode(0=binary;1=relax;2=smooth).....' 'modez'     0
' time step (seconds) (0=old)................' 'dth1'      0.0d6
' iteration residual (fractional)............' 'resid'     1.0d-5
' fractional temperature change.....0.01.....' 'cdelt'     0.01
' fractional volume change..........0.03.....' 'cdelv'     0.03
' fractional abundance change.......0.05.....' 'cdeln'     0.05
' zones allowed (0=no rezone,-k=k equal dms).' 'ktot'      3900
' force single zone to be added..............' 'kforce'    0
' time steps per rezone......................' 'l4'        4
' d log rho/d zone.(decade per zone rho 0.2).' 'dlnv'      0.2
' fractional mass in outer zone......1.0d-4..' 'xmmin'     1.0d-4
' biggest fractional zone mass.....0.01d0....' 'dmmax'     0.01d0 
' rezone for shell (edge of flame s5/s5max)..' 'facte'     1.0e-2 
' vline (*0.06)..............................' 'vline'     3.0d-1
' pline (*0.06)..............................' 'pline'     3.0d-1
' 0.1*(radius increase)......................' 'drmax'     1.0d-1
' smooth zone masses in rezone(yes=1)........' 'ismoo'     0
' map envelope onto Henyey grid (0=no).......' 'mapenv'    0
'.....................................................................'
' accretion rate (Msol/year).....2.5d-8......' 'peryear'   2.5d-9
' flag (0=no,1=acc,<0=loss),needs modes=2....' 'mloss'     -5
' loss: -1=only BSG,-2=reimers,-3=bloecker (see advec.f)..............'
'.......-4=schroeder, -5=vink.........................................'
' Pulsational mass loss 2, WR 1, both 3......' 'altloss'   0
'.....................................................................'
' include blackowen mass loss? (0=n,1=y).....' 'bomloss'   1
'.....................................................................'
' graphics(default=0=yes,else=none)..........' 'igraf'     0
' device1 for pgplot(/xwin,/cps,/ps)..panels.' 'device1'   '1/xwin'
' pgplot x-axis (0=r,1=logr,2=m,3=k).........' 'ixflag'    2
' pgplot panel 4(0=dT,1=nablas,2=gh).........' 'nabflg'    1
' pgplot panel 2 (color:0=frac(m),1=max X)...' 'nxflg'     1
' pgplot panels 234 ( low value x-axis)......' 'fbot'      0.
' pgplot panels 234 (high value x-axis)......' 'ftop'      0.
' powers of 10 in abundance..................' 'xlogmin'  -6.0
' powers of 10 in abundance..................' 'xlogmax'   0.0
' velocity scale for convection (cm/s).......' 'cvsc'      5.0d5
' device2 for pgplot(/xwin,/cps,/ps)..HR.....' 'device2'   '2/xwin'
' mimimum log Te for HR plot.................' 'gtmin'     3.4
' maximum log Te for HR plot.................' 'gtmax'     4.6
' minimum log L  for HR plot.................' 'glmin'     -1.0
' maximum log L  for HR plot.................' 'glmax'     3.0
' device3 for pgplot(/xwin,/cps,/ps)..cv.....' 'device3'   '3/xwin'
' mimimum log x for cv plot..................' 'pg3xmin'   0.
' maximum log x for cv plot..................' 'pg3xmax'   0.
' minimum log y  for cv plot.................' 'pg3ymin'   0.
' maximum log y  for cv plot.................' 'pg3ymax'   0.
' x variable for cv plot(0=time,1=model).....' 'ixpg3'     1
' y variable for cv plot(0=mass,1=logr)......' 'iypg3'     0
' device4 for pgplot(/xwin,/cps,/ps)..db.....' 'device4'   '7/xwin'
' mimimum x for db plot......................' 'pg4xmin'   100.
' maximum x for db plot......................' 'pg4xmax'   200.
' minimum y  for db plot.....................' 'pg4ymin'   -2.0
' maximum y  for db plot.....................' 'pg4ymax'   2.0
' x variable for db plot(0=zone k)...........' 'ixpg4'     0
' y variable(0-8=dr,nab,vc,s,H,He,ic,adv,dif)' 'iypg4'     0
