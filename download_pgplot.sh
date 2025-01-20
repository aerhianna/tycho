#Tycho comes with its own distribution of PGPLOT. This files is retained for reference only.
# Download + Extract to src/PGPLOT
pgplotdir=src/PGPLOT
[ ! -f pgplot.tar.gz ] && curl ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz --output pgplot.tar.gz
echo "[INFO] Extracting PGPLOT source code to $pgplotdir"
tar -xzf pgplot.tar.gz
mkdir -p $pgplotdir  # -p for creating parent directories if they don't exist
cp -R -f pgplot/* $pgplotdir  # -R for recursive, -f to overwrite

echo "[INFO] Customizing PGPLOT for usage with TYCHO"
# Change compiler to `gfortran`
linux_g77_conf="$pgplotdir/sys_linux/g77_gcc.conf"
sed -i .bak "s/FCOMPL=\"g77\"/FCOMPL=\"gfortran\"/g" $linux_g77_conf  # -i .bak for in-place editing
# Enable PostScript and Xwin displays
drivers_list="$pgplotdir/drivers.list"
sed -i .bak "s/! PSDRIV/  PSDRIV/g" $drivers_list
sed -i .bak "s/! XWDRIV/  XWDRIV/g" $drivers_list

# Cleanup
echo "[INFO] Removing temporary files"
rm -r pgplot
