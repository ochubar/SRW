#!/bin/bash

# Compile fftw-2.15 and srwlpy.so.
# Maksim Rakitin
# 2016-03-24

#---> Dirs and files:
root_dir=$PWD
gcc_dir=$root_dir/cpp/gcc
py_dir=$root_dir/py
ext_dir=$root_dir/ext_lib
fftw_dir=$ext_dir/fftw-2.1.5
log_fftw=$root_dir/log_fftw.log
log_srw=$root_dir/log_srw.log

echo -e "\nWelcome to SRW compiling script!\n" 

#---> Compile FFTW2:
echo "    Compiling FFTW2..."
cd $ext_dir
if [ -d "$fftw_dir" ]; then
    rm -rf $fftw_dir
fi

tar -zxf fftw-2.1.5.tar.gz
cd $fftw_dir

echo "        Configuring FFTW2..."
./configure --enable-float --with-pic --prefix=$fftw_dir/tmp >> $log_fftw 2>&1

echo "        Updating CFLAGS..."
sed 's/^CFLAGS = /CFLAGS = -fPIC /' -i $fftw_dir/Makefile

echo "        Making FFTW2..."
make >> $log_fftw 2>&1 && \
make install >> $log_fftw 2>&1 && \
cp fftw/.libs/libfftw.a $ext_dir

echo "    Compilation of FFTW2 finished."
echo ""

#---> Compile srwlpy.so:
echo "    Compiling SRW..."
cd $root_dir
make all > $log_srw 2>&1
echo "    Compilation of SRW finished."
echo ""

#---> Clean:
echo "    Cleaning..."
rm -rf $fftw_dir
rm -rf $ext_dir/libfftw.a $gcc_dir/libsrw.a $gcc_dir/srwlpy.so $py_dir/build/
echo ""

echo "Finished successfully!"
echo ""

