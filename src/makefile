###############################################################
# Makefile for TYCHO: 
#
# revised 8-27-02 (darnett)
#
# revised Sept. 23, 2007 (c.a.meakin) 
# (1) to query for site specific makefile
# (2) to handle run directory specific compilation
#     using local source files
# revised Jan. 10, 2008 (c.a.meakin)
# (1) 
#
#
#
################################################################


#  SETUP DIRECTORIES

localdir = ./
rootdir  = ./

include ${rootdir}/makefile.env


#  SETUP OBJECTS/DEPENDENCIES

include ${rootdir}/root.makefile
LIBALL += $(LIBPG)

#  EXPLICIT/IMPLICIT RULES
# (generic f77)
.f.o:
	$(FORT) $(FFLAGS) -I${l_inc} -I${r_inc} -c $< -o $@
.c.o:
	$(CC) -c $< -o $@


# PROGRAM COMPILATION/LOAD

${localdir}/libtycho.a: ${robjects}
	ar cr ${localdir}/libtycho.a ${robjects}

# build PGPLOT library
tychopgplot:
	cd ${l_pgplot} && make 

${localdir}/libprogram.a: ${localdir}/libtycho.a ${lobjects}        
	cp ${localdir}/libtycho.a ${localdir}/libprogram.a
	ar r libprogram.a ${lobjects} 


tycho8:	tychopgplot ${localdir}/libprogram.a ${l_program}/tycho.f 
	$(FORT) $(FFLAGS) ${l_program}/tycho.f -o tycho8 -L./ -lprogram -I${l_inc} $(LIBALL)
	mv tycho8 ../test

hrplot: ${localdir}/libprogram.a ${l_program}/hrplot.f
	$(FORT) $(FFLAGS) ${l_program}/hrplot.f -o hrplot -L./ -lprogram -I${l_inc} $(LIBALL)
	mv hrplot ../test

cvplot: ${localdir}/libprogram.a ${l_program}/cvplot.f
	$(FORT) $(FFLAGS) ${l_program}/cvplot.f -o cvplot -L./ -lprogram -I${l_inc} $(LIBALL)
	mv cvplot ../test

genplot: ${localdir}/libprogram.a ${l_program}/genplot.f
	$(FORT) $(FFLAGS) ${l_program}/genplot.f -o genplot -L./ -lprogram -I${l_inc} $(LIBALL)
	mv genplot ../test

gennuc: ${localdir}/libprogram.a ${l_program}/gennuc.f
	$(FORT) $(FFLAGS) ${l_program}/gennuc.f -o gennuc -L./ -lprogram -I${l_inc} $(LIBALL)
	mv gennuc ../test

genrate: ${localdir}/libprogram.a ${l_program}/genrate.f
	$(FORT) $(FFLAGS) ${l_program}/genrate.f -o genrate -L./ -lprogram -I${l_inc} $(LIBALL)
	mv genrate ../test

genex: ${localdir}/libprogram.a ${l_program}/genex.f
	$(FORT) $(FFLAGS) ${l_program}/genex.f -o genex -L./ -lprogram -I${l_inc} $(LIBALL)
	mv genex ../test

omit: ${localdir}/libprogram.a ${l_program}/omit.f
	$(FORT) $(FFLAGS) ${l_program}/omit.f -o omit -L./ -lprogram -I${l_inc} $(LIBALL)
	mv omit ../test

ratios: ${localdir}/libprogram.a ${l_program}/ratios.f
	$(FORT) $(FFLAGS) ${l_program}/ratios.f -o ratios -L./ -lprogram -I${l_inc} $(LIBALL)
	mv ratios ../test

hrt:    ${localdir}/libprogram.a ${l_program}/hrt.f
	$(FORT) $(FFLAGS) ${l_program}/hrt.f -o hrt -L./ -lprogram -I${l_inc} $(LIBALL)
	mv hrt ../test

testenvel:    ${localdir}/libprogram.a ${l_program}/testenvel.f
	$(FORT) $(FFLAGS) ${l_program}/testenvel.f -o testenvel -L./ -lprogram -I${l_inc} $(LIBALL)
	mv testenvel ../test



all: 
	make tycho8 hrplot hrt cvplot genplot gennuc genrate genex omit ratios

clean:
	-rm -f ${robjects} ${lobjects} ${other_obs} ${localdir}/lib*.a tycho8
	-rm -f ${l_pgplot}/lib/*
	-rm -f pgplot_env*






