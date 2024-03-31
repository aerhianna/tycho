#!/bin/sh

echo "###################################################################"
echo [DEBUG][EXECUTING] src/PGPLOT/make_pgplot_tycho.sh
echo "Attempting to build PGPLOT library. Copyright Tim Pearson / Caltech"
PWD=`pwd`

echo [DEBUG] FC=$FC
echo [DEBUG] CC=$CC
echo [DEBUG] FFLAGS=$FFLAGS
echo [DEBUG] XLIBS=$XLIBS
echo [DEBUG] XINCS=$XINCS

#XLIBS=\"$XLIBS
#echo $XLIBS

cd $1
cd sys_linux

cat tycho.conf.in|sed s,'@FC@',$FC,g|sed s,'@CC@',$CC,g|sed s,'@XLIBS@',"$XLIBS",g| \
     sed s,'@XINCS@',"$XINCS",g|sed s,'@FFLAGS@',"$FFLAGS",g \
    > tycho.conf

mkdir --parents ../lib
cd ../lib
rm *.o *.a
cp ../drivers.list .
../makemake $PWD/../ linux tycho

PGPLOT_DIR=$PWD
PGPLOT_FONT=$PWD/grfont.dat
echo "#!/bin/sh" > pgplot_env.sh
echo "export PGPLOT_DIR="$PGPLOT_DIR >> pgplot_env.sh
echo "export PGPLOT_FONT="$PGPLOT_FONT >> pgplot_env.sh

echo "#!/bin/csh" > pgplot_env.csh
echo "setenv PGPLOT_DIR "$PGPLOT_DIR >> pgplot_env.csh
echo "setenv PGPLOT_FONT "$PGPLOT_FONT >> pgplot_env.csh

cp pgplot_env.sh ../../
cp pgplot_env.csh ../../
