#
# Makefile to compile TurbSim v1.06.00
# on 64-bit Linux with gfortran 4.8.2
#
# Jenni Rinker, Duke University/NWTC
# 19-feb-2015
#
# Compile TurbSim with "$ make turbsim"
#

FC = gfortran
export FFLAGS = -fdefault-real-8 -ffree-line-length-none -fPIC -O2 -mtune=native

NWTC_LIB = nwtc
TURBSIM_LIB = turbsim

export BASEDIR = $(PWD)
export OUTDIR = $(PWD)/build

ALL_OBJ = $(OUTDIR)/*.o

turbsim: turb_lib fix_build_dir
	$(FC) $(ALL_OBJ) $(OUTDIR)/native/*.o -ldl -o TurbSim

# Needed to compile SysSubs.mod
fix_build_dir:
	@touch $(OUTDIR)/SysGnuLinux.o
	@mkdir -p $(OUTDIR)/native
	@mv $(OUTDIR)/SysGnuLinux.o $(OUTDIR)/native/SysGnu.o

turb_lib: nwtc_lib
	$(MAKE) -C $(TURBSIM_LIB)/Source/RNG
	$(MAKE) -C $(TURBSIM_LIB)/Source/NetLib/FFTPACK
	$(MAKE) -C $(TURBSIM_LIB)/Source/NetLib/LAPACK/BLAS
	$(MAKE) -C $(TURBSIM_LIB)/Source/NetLib/LAPACK
	$(MAKE) -C $(TURBSIM_LIB)/Source

nwtc_lib:
	$(MAKE) -C $(NWTC_LIB)/source

.PHONY: clean veryclean

clean:
	$(MAKE) clean -C $(NWTC_LIB)/source; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/RNG; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/NetLib/FFTPACK; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/NetLib/LAPACK/BLAS; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/NetLib/LAPACK
veryclean: clean
	rm -rf build include nwtc
