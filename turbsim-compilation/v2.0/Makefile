
FC = gfortran
export FFLAGS = -fdefault-real-8 -ffree-line-length-none -fPIC -O2 -mtune=native

TURBSIM_LIB = trunk

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

turb_lib:
	$(MAKE) -C $(TURBSIM_LIB)/Source/dependencies/NWTC_Library
	$(MAKE) -C $(TURBSIM_LIB)/Source/dependencies/RNG
	$(MAKE) -C $(TURBSIM_LIB)/Source/dependencies/NetLib
	$(MAKE) -C $(TURBSIM_LIB)/Source

.PHONY: clean veryclean

clean:
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/dependencies/NWTC_Library; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/dependencies/RNG; \
	$(MAKE) clean -C $(TURBSIM_LIB)/Source/dependencies/NetLib; \
#veryclean: clean
#	rm -rf build include nwtc
