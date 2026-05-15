localdir = ./
rootdir  = ../

SITE=$(shell /bin/hostname)
MAKEFILESITE=${rootdir}/SITES/$(SITE)/makefile.$(SITE)
include ${MAKEFILESITE}

all: lib/libpgplot.a lib/pgxwin_server

lib/libpgplot.a:
	echo $(SITE)
	echo $(FORT)
	echo $(CC)
	echo $(FFLAGS)
	echo $(LIBX11)

	/usr/bin/env FC="$(FORT)" CC="$(CC)" FFLAGS="$(FFLAGS)" XLIBS="$(LIBX11)" \
	XINCS="$(INCX11)" /bin/sh ./make_pgplot_tycho.sh ./
	make -C lib

lib/pgxwin_server: lib/libpgplot.a
	make -C lib
