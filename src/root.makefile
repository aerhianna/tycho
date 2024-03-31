#
#
#	|------------------------------------------------|
#	|  Makefile for restructured TYCHO:              |
#	|  Casey A. Meakin, U.Chicago/U.Arizona (2007)   |
#	|------------------------------------------------|
#
#
#------------------------------------------------ SETUP DIRECTORIES
r_main     = ${rootdir}/MAIN
r_aux      = ${rootdir}/AUX
r_lin      = ${rootdir}/LINSOLVE
r_inc      = ${rootdir}/INC
r_program  = ${rootdir}/PROGRAM
r_dimen    = ${rootdir}/DIMEN

l_main     = ${localdir}/MAIN
l_aux      = ${localdir}/AUX
l_lin      = ${localdir}/LINSOLVE
l_inc      = ${localdir}/INC
l_program  = ${localdir}/PROGRAM
l_dimen    = ${localdir}/DIMEN

l_pgplot   = ${rootdir}/PGPLOT

#  LOCAL OBJECTS

f77_l_main = $(wildcard ${l_main}/*.f)
obj_l_main = $(patsubst %.f,%.o,$(f77_l_main))

f77_l_aux  = $(wildcard ${l_aux}/*.f)
gcc_l_auxc = $(wildcard ${l_aux}/*.c)
obj_l_aux  = $(patsubst %.f,%.o,$(f77_l_aux))
obj_l_auxc = $(patsubst %.c,%.o,$(gcc_l_auxc))

f77_l_lin = $(wildcard ${l_lin}/*.f)
obj_l_lin = $(patsubst %.f,%.o,$(f77_l_lin))

lsource = ${f77_l_main} ${f77_l_program} ${f77_l_aux} ${gcc_l_auxc} ${f77_l_lin} 
lobjects = ${obj_l_main} ${obj_l_aux} ${obj_l_auxc} ${obj_l_lin}

# added wda 4-29-09
#LIBPG = -L$(localdir)/PGPLOT/lib -lpgplot
########################################
#  ROOT OBJECTS

f77_r_main = $(wildcard ${r_main}/*.f)
obj_r_main = $(patsubst %.f,%.o,$(f77_r_main))

f77_r_aux  = $(wildcard ${r_aux}/*.f)
gcc_r_auxc = $(wildcard ${r_aux}/*.c)
obj_r_aux  = $(patsubst %.f,%.o,$(f77_r_aux))
obj_r_auxc = $(patsubst %.c,%.o,$(gcc_r_auxc))

f77_r_lin = $(wildcard ${r_lin}/*.f)
obj_r_lin = $(patsubst %.f,%.o,$(f77_r_lin))

rsource  = ${f77_r_main} ${f77_r_program} ${f77_r_aux} ${gcc_r_auxc} ${f77_r_lin}
robjects = ${obj_r_main} ${obj_r_aux} ${obj_r_auxc} ${obj_r_lin}

#  INCLUDES
r_inc_files = $(wildcard ${r_inc}/*)
l_inc_files = $(wildcard ${l_inc}/*)
dimen_files = $(wildcard ${l_dimen}/*)

#  DEPEND
${robjects}:${r_inc_files}
${lobjects}:${l_inc_files}

