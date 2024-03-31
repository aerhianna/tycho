MYENV=myenv
[ -n "$1" ] && MYENV="$1"


THISDIR="$(dirname $0)"  # The path to the directory containing this script

cd $THISDIR

# Download PGPLOT
./download_pgplot.sh

# Configure `makefile.env`
cp src/makefile.env.$OSTYPE src/makefile.env

# Compile Tycho and related programs
mkdir --parents test
cd src && make all && cd ../

# Create environment folder
mkdir --parents $MYENV

# Copy in built programs (requires `cd src && make all` first)
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

# Copy in data files
cp data/* $MYENV
mv $MYENV/EOSdata $MYENV/EOSdata_H-He

# Copy in config files
cp config/* $MYENV

# Copy in a base model file
cp models/a000000 $MYENV/imodel

# Create a .gitignore files
echo "# Ignore all auto-added environment files" > $MYENV/.gitignore
ls $MYENV >> $MYENV/.gitignore