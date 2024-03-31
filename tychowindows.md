# Tycho Windows

## Main Window

### Section 1:

First line basic settings information

Green for normal operation. Will turn red when iteration does not converge but usually recovers
on its own.

time(s) - stellar age in seconds

dt(s) - length of timestep in seconds

log Te - log(Teff)

log L - log(L/L⊙)

vrot - equatorial rotational velocity (cm/s)

alphaml - not used in this version

Model - total number of converged models for this track

steps - total number of timesteps since starting this instance of the code

kk - number of zones in the calculation

it - number of iterations for the Henyey integration to converge

J53 - angular momentum in units of 1053 g cm2 s-1

nadd - number of zones to be added in next rezoning

(T,V,isotope,dR,dt) - variable and zone limiting the timestep

(T,V) - variable and zone determining the number of iterations

wtest - convergence of Henyey integration

echk - diagnostic of total energy erg/g including binding energy

eta(2) - neutron excess in innermost zone

ndel - number of zones to be deleted in next rezoning

M/sol - total mass in M⊙

Men/sol - mass in envelope part of integration in M⊙

dM/dt - mass loss rate in M⊙/y

He - mass fraction of He at surface

netadd - total number of zones to be added in next rezone

R(cm) - stellar radius in cm

T(c) - central temperature in K

rho(c) - central density in g cm-3

T(kk) - temperature at outermost zone of the interior Henyey integration in K

metals - mass fraction of metals at surface

mixmode - flag for treatment of convective boundaries

knu - zone of peak neutrino cooling

epnu - rate of neutrino cooling at knu in erg g-1 s-1

knuc - zone of peak nuclear energy generation

epnuc - rate of nuclear energy generation at knuc in erg g-1 s-1

Lnuc - total luminosity from nuclear energy generation in L/L⊙

ncytot - number of subcycles required by nuclear network

### Section 2:

log(T)/log(rho) plot

Green lines mark boundaries between different size nuclear networks

Dark gray lines mark boundaries between different equations of state

Curve is interior structure profile of star. Red line is envelope part of integration. Cross hatches
are interior part of integration with color corresponding to most abundant species in zone.

red vertical line - location where number of iterations (it) set.

white vertical line - location where timestep (dt) is set

green vertical line - location where nuclear abundance changing most rapidly labeled with
fastest changing isotope

### Section 3:

blue hatched line - nuclear energy generation normalized to peak zone

white line - log temperature normalized to peak value

gray line - integrated luminosity from nuclear burning and neutrinos

magenta line - luminosity carried by enthalpy of convection

green line - integrated luminosity from all sources including mechanical work L/L⊙

yellow hatches - convective velocity (km/s). Changes to red for envelope part of Henyey
integration.

orange line u(r) - radial expansion or contraction velocity

blue line - sound speed, often invisible at default scale

purple line v(rot) - rotational velocity (km/s). Does not display if zero.

Left vertical axis ticks for L in L/L⊙.

right vertical axis ticks for velocities in km/s.

### Section 4:

Abundances - log mass fraction. Colors correspond to labels. Species that do not change with
mass coordinate not shown.

### Section 5:

green line (nab) - actual nabla of star (dlnT/dlnP)

yellow line (nad) - adiabatic nabla

red line (nrad) - nabla established by radiative diffusion

orange line (doux) - Ledoux mean molecular weight gradient

dark blue line (rich) - Richardson number, usually invisible at default scale

light blue line (log dm) - log of zone mass/total mass

arrows - indicates where zones will be added/removed in next rezone

Horizontal axis for sections 3-5 can be 
set to r, log r, m, or zone number k in params.d
Plot limits can be set in params.d.

## HR Window:

HR diagram of evolutionary track. Limits set in params.d. If a calculation is stopped and
restarted old timesteps will appear in blue, timesteps since restart in yellow.

## Convection Window:

Extent of convection in star shown in red. If a calculation is stopped and restarted old
timesteps will appear in grayscale.

Zone of peak nuclear burning shown in blue. Secondary peak shown in green. If a calculation is
stopped and restarted old timesteps will appear in yellow and orange, respectively.

Axes can be set in params.d to M or logR on Y axis, t or step number on x axis.

Can be reproduced by analysis program cvplot, which uses data in file cv.(prefix).