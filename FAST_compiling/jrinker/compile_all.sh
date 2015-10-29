#!/bin/bash

if [ -d build_structure ] ; then rm -rf build_structure ; fi

tar -zxvf src/build_structure.tar.gz
tar -zxvf src/AD_v13.00.02a-bjj.tar.gz -C build_structure/aero/
tar -zxvf src/FAST_v7.02.00d-bjj.tar.gz -C build_structure/fast/
tar -zxvf src/NWTC_Lib_v1.07.00b-mlb.tar.gz -C build_structure/nwtc/
tar -zxvf src/Crunch_v3.02.00c-mlb.tar.gz -C build_structure/crunch/
unzip src/InflowWind_v1.02.00b-adp.exe -d build_structure/inflow
unzip src/MBC_v1.00.00a-gsb.exe -d build_structure/mbc
unzip src/MCrunch_v1.00.00ab-gjh.exe -d build_structure/mcrunch
unzip src/TurbSim_v1.06.00.exe -d build_structure/turbsim

echo
echo cd build_structure
sleep 2
cd build_structure

echo
sleep 2
echo Patching 
patch --dry-run -p1 -i all.patch

echo
sleep 2
echo Patching 
patch -p1 -i all.patch

CURMATLABPATH="\/opt\/matlab\/r2012b\/"
NEWMATLABPATH="\/afs\/ies.auc.dk\/sw\/pack\/matlab-r2011a\/"
NEWMATLABPATH="\/afs\/ies.auc.dk\/sw\/pack\/matlab-r2012a\/opt\/MATLAB\/R2012a\/"

sed -i.bak "s/$CURMATLABPATH/$NEWMATLABPATH/" Makefile
echo
echo
grep ^MATLAB_DIR Makefile
sleep 2

echo
echo make turbsim
sleep 2
make turbsim

echo
echo make clean
sleep 2
make clean

echo
echo make crunch
sleep 2
make crunch

echo
echo make clean
sleep 2
make clean

echo
echo make ...
sleep 2
make

echo
echo make test_cases
sleep 2
make test_cases

echo done ...