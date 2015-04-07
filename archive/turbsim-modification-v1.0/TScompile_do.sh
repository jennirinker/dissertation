#!/bin/bash
#
# Usage:
# 		$ ./TScompile_do.sh <directory-name>
#
# Linux bash script to compile TurbSim:
#		1) Compiles TurbSim, cleans up
#		2) Makes TurbSim executable
#		3) Moves file into /turbsim/
#
# Jenni Rinker, Duke University/NWTC
# 24-Mar-2015

# get dirname to compile, check it exists, change dirs
dname=$1
echo "  Compiling directory $dname"
if [ -d $dname ]; then 
   cd $dname
else
   echo "  ERROR: Directory does not exist"
   exit
fi


echo "  Compiling TurbSim..."

# compile turbsim, make it exectuable
make turbsim
echo "  Cleaning..."
make clean
echo "  Making executable..."
chmod +x TurbSim

# rearrange /turbsim/, relocate to upper directory
mv -f TurbSim turbsim/

echo "  Script complete."