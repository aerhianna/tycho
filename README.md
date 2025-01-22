# Getting Started with Tycho

## System Requirements:

The code runs on any *NIX-based system. For most of you this means Linux or macOS. It
requires fortran and c compilers.
 - Most Linux distributions come with gcc and gfortran.
 - If you are using macOS you will need to install
   - Xcode from the app store
   - the Xcode command line tools by entering sudo xcode-select —install on the command
  line
   - XQuartz https://www.xquartz.org
   - Mac users will also need to download and install a fortran compiler. See the instructions below for installing with Homebrew.

**Linux (WSL)**
```
$ sudo apt install gcc gfortran make libx11-dev --assume-yes  
```
> `--assume-yes` is optional. It skips the prompt asking you to confirm the install

Verify that `gcc` installed correctly
```
$ gcc --version

gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
Copyright (C) 2021 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

Verify that `gfortran` installed correctly
```
$ gfortran --version

GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
Copyright (C) 2021 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

Verify that `make` installed correctly
```
$ make --version

GNU Make 4.3
Built for x86_64-pc-linux-gnu
Copyright (C) 1988-2020 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
```

Ensure that key scripts have executable permissions
```bash
cd src/PGPLOT
sudo chmod +x makemake makedoc makehtml makehelp maketex
```

**MacOS**

If you are using macOS you will need to install
   - Xcode from the app store
   - the Xcode command line tools by entering sudo xcode-select —install on the command
  line
   - XQuartz https://www.xquartz.org
   - Mac users will also need to download and install a fortran compiler. 

Two security settings will need to be checked.
1) Under System Settings -> Privacy & Security set Allow Applications downloaded from App Store and identified developers
2) Under System Settings -> Privacy & Security -> Developer Tools allow the Terminal to run software locally that does not meet the system's security policy.
By default macOS does not allow unsigned executables to run as a
security measure. If you do not see the options above, most likely because you are running an old version of macOS,
enter sudo spctl —master-disable at the command line.

You can Install the [Homebrew package manager](https://brew.sh/) and use it to install gcc and gfortran. 

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Use Homebrew to install the necessary compilers to build Tycho on your Mac

```
brew install gcc gfortran make
```

Ensure `gcc` installed correctly

```
% gcc --version

Apple clang version 15.0.0 (clang-1500.1.0.2.5)
Target: x86_64-apple-darwin22.6.0
Thread model: posix
InstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin
```

Ensure `gfortran` installed correctly

```
% gfortran --version

GNU Fortran (Homebrew GCC 13.2.0) 13.2.0
Copyright (C) 2023 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

Ensure `make` installed correctly

```
% make --version

GNU Make 3.81
Copyright (C) 2006  Free Software Foundation, Inc.
This is free software; see the source for copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

This program built for i386-apple-darwin11.3.0
```

## Directory Structure:

After unpacking the tarball you will see a root directory called tycho-spring21/ *(NOTE: This repo is the unpacked version)*. Beneath this are
subdirectories. The ones you are most likely to use are:

```bash
├── burn/ # Utilities for creating reaction networks and doing single zone calculations
├── cmd_files/ # Templates for all the *.in files used by analysis and setup programs
├── run/ # Working space for the code
├── src/ # Source code
└── test/ # Compiled binaries and testing routines
```
## Compilation and Setup:

TYCHO uses GNUmake to compile the code, build libraries, and create executables.

1. Set up the makefile for your computer. The following command should grab the correct template for your operating system.
   
   ```bash
   cp src/makefile.env.$OSTYPE src/makefile.env
   ```

    <details>
        <summary>Example `makefile.env`</summary>

    - An example of where you might need to make changes is below.

    ```bash
    # for OSX
    # Change LIBPG and LIBX11 to match the location of your Tycho installation.
    LIBPG = /Users/payoung/tycho-1-14/src/PGPLOT/lib/libpgplot.a
    LIBX11 = /opt/X11/lib/libX11.dylib

    # NOTE: LIBAPP and LIBFOUND may need to be uncommented for older versions of macOS.
    #LIBFOUND = /System/Library/Frameworks/Foundation.framework/Versions/
    Current/Foundation -lcc_dynamic
    #LIBAPP = /System/Library/Frameworks/AppKit.framework/Versions/
    Current/AppKit -lcc_dynamic
    LIBALL = /Users/payoung/tycho-clean/src/PGPLOT/lib/libpgplot.a \
    /opt/X11/lib/libX11.dylib \
    # /System/Library/Frameworks/Foundation.framework/Versions/
    Current/Foundation -lgcc \
    # /System/Library/Frameworks/AppKit.framework/Versions/Current/
    AppKit -lgcc
    INCX11 = -I/opt/X11/include
    # for Linux Intel with gfortran
    #OPT = -g3 -O2 -C -Wall
    # for powerpc64
    #OPT = -g3 -O2 -C -Wall -mcpu=powerpc64
    #OPTHI = -O3 -finline-functions -funroll-loops
    #FORT = gfortran
    #OPT = -O2 -g
    PROF = -pg
    # for macOS with gfortran
    FORT = gfortran
    CC = gcc
    #OPT = -O3 -g -msse -msse2
    #OPT = -O2 -g -mcmodel=medium -m64
    OPT = -O3 -march=native -g -m64 -fallow-argument-mismatch
    MIPS =
    FFLAGS = $(OPT) $(MIPS)
    CFLAGS =
    ```

    </details>

1. You will need to add the PGPLOT_DIR and PG_FONT variables to your .tcshrc or .bashrc
or .profile file. Here is an example. An example .tcshrc file is available on Canvas. Here is a
resource for converting between bash and tcsh: https://hyperpolyglot.org/unix-shells.

```shell
# tcsh:
setenv TYCHO_HOME ${HOME}/git/tycho
setenv PGPLOT_DIR ${TYCHO_HOME}/src/PGPLOT/lib/
setenv PG_FONT ${TYCHO_HOME}/src/PGPLOT/lib/grfont.dat
```
```bash
# bash:
export TYCHO_HOME="$HOME/git/tycho"
export PGPLOT_DIR="$TYCHO_HOME/src/PGPLOT/lib/"
export PGPLOT_FONT="$TYCHO_HOME/src/PGPLOT/lib/grfont.dat"
```

Make sure to restart your shell for the settings to take effect (e.g. `exec $SHELL`)

3. Compile the code. 
   - You need to be in `src/`.
   - Type `make clean` to remove any old compiled subroutines. The command make all will make tycho8 and most of the setup/analysis programs. To make an individual program use `make tycho8` or the appropriate name (`make hrt`, `make gennuc`, etc.). These will end up in `test/`
4. The simplest way to organize different projects with different sets of models is to create a
new directory for each project.
   - The `run/` directory is an example. You can copy `run/` or  make a new directory for your problem.
   - If you do the latter, copy the executables from `test/` to this directory as well as a `params.*` file from `run/` to `params.d`.
   - If you created a new project directory from scratch, copy `init-data` (found in `bin/`) to your problem directory and run it. 
   - Copy an initial model from `run/` to the working directory. 
5. You are now ready to move on to the steps we will discuss in class.

## Tycho Data

Tycho produces data is several ways:
 - An online log in a terminal window
 - Online graphics in 3 PGPLOT windows
 - Model files: Regular ”dumps" of snapshots of interior structure properties (which are re-
startable) with the format ?? followed by a 5 digit number (i.e. ab00000).
 - Time sequences for position in the HR diagram (hr.??) and for convective status (cv.??)
   - Here "??" represents a two-character symbol identifying the computational sequence, and is user defined in the control file "params.d".

Control file: params.d

## Programs for model modification

genex: This auxillary program takes a model (in file "old.model") and maps it into "new.model",
with optional changes:

 - mass and radius scaling
 - abundance scaling
 - solid body rotation

Control file: genex.in

## Analysis Suite

In additon to the hydrodynamics/evolution code Tycho, a suite of analysis programs are
included. These use PGPLOT graphics and allow rapid investigation of both snapshots and
time sequences. The most commonly used are:

*Snapshots (use model files):*

genplot (workhorse graphics of model files, plot just about any variable you can think of for the
stellar interior)

gennuc (abundance profiles in stellar interior)

genrate (details of reactions in individual zones)

*Time sequences:*

hrt (surface quantities, including HR diagram, radius, abundances, mass loss rates, etc. Uses
hr.?? files)

cvplot (evolution of convective structure, uses cv.?? files)

The control files for all these programs are named according to the program, plus a postfix:
“.in”, i.e. genplot.in.

## Creating initial models

In the /run directory you will find a variety of stellar models. Models at the beginning of the pre-
MS will be labeled ??00000. Any of these can be scaled to create an initial model of arbitrary
properties. To find the mass of a model, consult the header. Entry 4 will be mass in solar
masses. Entry 5 will be Z (mass fraction of heavy elements).


TYCHO 8.00 a0 <mark>1.00</mark> 522 <mark>1.5E-02</mark> 1.60 .709 2 Thu Dec 15 11:54:15 2016

 - copy starting model to old.model (i.e. `cp ab00000 old.model`)
 - edit genex.in (described below) (i.e. `nano genex.in`)
 - run `genex` to create new.model (i.e. `./genex`)
 - copy new.model to imodel (i.e. `cp new.model imodel`)

```
'.....................................................................
'
' GENEX:
'
'.....................................................................
'
' Scaling(0=M,R;1=sol scale;2=OPAL)..........' 'ifm' 0
' lumin.ne.0..reevaluate luminosity=radiative' 'lumin' 0
' OPAL opacity type (0=type2,1=type1)........' 'nopac' 1
' OPAL EOS (0=no,1=yes)......................' 'nopaleos' 1
'..............................................MASS and RADIUS
scaling'
' ifm=0:Fractional mass scaling (0=ignore).' 'fm' 1.056d0
' ifm=0:Fractional radius scaling (0=ignore).' 'fr' 1.0d0
'.........................................SYSTEMATIC ABUNDANCE
SCALING'
' ifm=1:metallicity scale relative to solar..' 'zscale' 0.0d-0
' ifm=1:production ratio He4/z...............' 'hetoz' 2.1d0
'......................................PARAMETERIZED ABUNDANCE
SCALING'
' ifm=2:mass fraction for metallicity........' 'zpop' 0.01881d0
' ifm=2:mass fraction for He4................' 'xhe' 0.2676
' ifm=2:mass fraction for extra C12 (OPAL)...' 'xc12' 0.00
' ifm=2:mass fraction for extra O16 (OPAL)...' 'xo16' 0.00
'.....................................................................
'
' Solid body rotation if omegbar .gt. 0/sec..' 'omegbar' 0.0d-7
'.....................................................................
'
' Mixing mode flag...........................' 'mixmode' 2
‘.....................................................................
'
```

genex will scale mass and radius (set ifm=0) or abundances relative to solar (ifm=1). Mass and
radius are scaled relative to old.model (i.e., if old.model is a 2 Msol star, you would set fm=2.0
to make a 4 Msol star). Scaling radius by the same amount will usually be close enough for the
new model to converge in Tycho. Scaling from the nearest existing model is usually desirable.

Abundances are scaled relative to a reference composition, by default the solar abundances in
Lodders (2010). So for 1.5 Zsol, you would set ifm=1 and zscale = 1.5. The production of He
with metals (hetoz) is observationally determined and shouldn’t be changed. Compositions with
non-solar abundance ratios can also be used. We will discuss this if you choose to work on
this for a project.

You will not use ifm=2. You can also create a rotating model by setting omegbar to a non-zero
value. omegbar is rotational frequency in radian/sec.

## Running Tycho

After compiling Tycho, copy /test/tycho8 to your working directory (i.e. in /tycho-spring21, cp
test/tycho8 run/). To run, type ./tycho8 or tycho8 if ./ is in your path. The same applies to other
programs such as genex.

Parameters controlling various code options are set in the file params.d. An annotated example
is in [paramsannotated.md](./paramsannotated.md). The code will always consult this file when started.

Tycho will always read initial conditions from imodel. To restart a run that you have stopped,
copy the last output model file to imodel (i.e. cp ab00100 imodel). You can restart from any
model, but be aware that later models with the same prefix will not be overwritten, so you may
want to change the prefix or delete later models. The hr and cv files will be overwritten from the
restarting point, so you should not need to worry about them.

To stop Tycho gracefully, in a terminal while in your working directory, type `touch starlock` and
enter. This will end the run and create a model file from which you can restart.

A basic explanation of the runtime graphics is contained in [tychowindows.md](./tychowindows.md). Keys to the
model, hr, and cv file formats will be provided in case you want to write your own analysis
code.

## FAQ

### Error: `Couldn't find program "pgxwin_server"`

Here's what the full error message typically looks like:
```
PGPLOT /xw: Couldn't find program "pgxwin_server" in the directory named
PGPLOT /xw: in your PGPLOT_DIR environment variable, or in any directory
PGPLOT /xw: listed in your PATH environment variable.
Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO IEEE_OVERFLOW_FLAG IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
STOP  pgopen error
```

**Potential Solutions**

Ensure your PGPLOT_DIR variable is set correctly
```bash
$ echo $PGPLOT_DIR
<PATH TO YOUR PGPLOT DIR>  # If blank, your PGPLOT_DIR isn't being set at all
$ ls $PGPLOT_DIR | grep pgxwin_server
pgxwin_server  # If blank, your PGPLOT_DIR is pointing to the wrong folder.
```

It's also possible that `pgxwin_server` doesn't have executable permissions. Try adding those.
```bash
sudo chmod +x $PGPLOT_DIR/pgxwin_server
```

### Error: `Timed out waiting for program PGPLOT/lib//pgxwin_server to start`

Here's what the full error message typically looks like:
```
[xcb] Unknown sequence number while appending request
[xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
[xcb] Aborting, sorry about that.
pgxwin_server: ../../src/xcb_io.c:157: append_pending_request: Assertion `!xcb_xlib_unknown_seq_number' failed.
PGPLOT /xw: Waiting for $TYCHO_HOME/src/PGPLOT/lib//pgxwin_server to start (timeout in 7 seconds).
PGPLOT /xw: Timed out waiting for program $TYCHO_HOME/src/PGPLOT/lib//pgxwin_server to start
Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO IEEE_OVERFLOW_FLAG IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
STOP  pgopen error
```

**Potential Solutions**

Honestly not sure exactly how I fixed this one, but here's what I gathered from my shell `history`

I followed the instructions for [verifying WSLg](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps)
```
sudo apt install gnome-text-editor -y
gnome-text-editor ~/.bashrc
```

I then updated my system
```
sudo apt update
sudo apt full-upgrade
```

Finally, I did another `make clean tycho8`. The new build launched shortly thereafter.
