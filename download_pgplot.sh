# Download + Extract to src/PGPLOT
pgplotdir=src/PGPLOT
[ ! -f pgplot.tar.gz ] && curl ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz --output pgplot.tar.gz
tar -xzf pgplot.tar.gz
mkdir --parents $pgplotdir
cp --recursive --force pgplot/* $pgplotdir

# Change compiler to `gfortran`
linux_g77_conf="$pgplotdir/sys_linux/g77_gcc.conf"
sed --in-place "s/FCOMPL=\"g77\"/FCOMPL=\"gfortran\"/g" $linux_g77_conf
# Enable PostScript and Xwin displays
drivers_list="$pgplotdir/drivers.list"
sed --in-place "s/! PSDRIV/  PSDRIV/g" $drivers_list
sed --in-place "s/! XWDRIV/  XWDRIV/g" $drivers_list

# Cleanup
rm -r pgplot