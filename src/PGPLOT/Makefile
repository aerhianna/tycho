localdir = ./
rootdir  = ../

include ${rootdir}/makefile.env

all: lib/libpgplot.a lib/pgxwin_server

lib/libpgplot.a:
	@echo "###################################################################"
	@echo [DEBUG][EXECUTING] src/PGPLOT/Makefile :: lib/libpgplot.a
	@echo [DEBUG] SITE=$(SITE)
	@echo [DEBUG] FORT=$(FORT)
	@echo [DEBUG] CC=$(CC)
	@echo [DEBUG] FFLAGS=$(FFLAGS)
	@echo [DEBUG] LIBX11=$(LIBX11)

	/usr/bin/env FC="$(FORT)" CC="$(CC)" FFLAGS="$(FFLAGS)" XLIBS="$(LIBX11)" \
	XINCS="$(INCX11)" /bin/sh ./make_pgplot_tycho.sh ./
	make -C lib

lib/pgxwin_server: lib/libpgplot.a
	make -C lib
