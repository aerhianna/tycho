MYENV=myenv
[ -n "$1" ] && MYENV="$1"


THISDIR="$(dirname $0)"  # The path to the directory containing this script

cd $THISDIR

echo "[INFO] Downloading PGPLOT source code"
./download_pgplot.sh

echo "[INFO] Creating a makefile.env file for your OS ($OSTYPE)"
cp src/makefile.env.$OSTYPE src/makefile.env

echo "[INFO] Compiling Tycho and related programs"
mkdir -p test
cd src && make all && cd ../

echo "[INFO] Creating environment directory $MYENV"
mkdir -p $MYENV

echo "[INFO] Copying built programs to $MYENV"
cp test/cvplot $MYENV
cp test/genex $MYENV
cp test/gennuc $MYENV
cp test/genplot $MYENV
cp test/genrate $MYENV
cp test/hrplot $MYENV
cp test/hrt $MYENV
cp test/omit $MYENV
cp test/ratios $MYENV
cp test/tycho8 $MYENV

echo "[INFO] Copying data files to $MYENV"
cp data/* $MYENV
mv $MYENV/EOSdata $MYENV/EOSdata_H-He

echo "[INFO] Copying config files to $MYENV"
cp config/* $MYENV

echo "[INFO] Copying a sample base model to $MYENV"
cp models/a000000 $MYENV/imodel

echo "[INFO] Adding a .gitignore file to $MYENV to exclude auto-added environment files from source control"
echo "# Ignore all auto-added environment files" > $MYENV/.gitignore
ls $MYENV >> $MYENV/.gitignore