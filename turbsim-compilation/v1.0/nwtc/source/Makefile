FC = gfortran

SRC = DoubPrec.f90
SRC += ModMesh.f90
SRC += SysGnuLinux.f90
SRC += NWTC_IO.f90
SRC += NWTC_Num.f90
SRC += NWTC_Aero.f90
SRC += NWTC_Library.f90

OBJ = $(SRC:.f90=.o)

ifdef BASEDIR
	INCLUDE = $(BASEDIR)/include
	MODDIR = $(BASEDIR)/include
else
	MODDIR = ./
endif

ifndef OUTDIR
	OUTDIR = ./
endif

all: checkdir $(SRC:.f90=.d) $(addprefix $(OUTDIR)/,$(OBJ))

$(OUTDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ -J$(MODDIR)	-I$(INCLUDE)

%.d: %.f90
	@set -e; rm -f $@; \
	$(FC) -E -cpp -MM $< -MF $@.$$$$ > /dev/null; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

checkdir:
	@mkdir -p $(OUTDIR); \
	mkdir -p $(MODDIR);

-include $(SRC:.f90=.d)
      
.PHONY: clean checkdir
clean:
	rm -rf *.d *.s $(OUTDIR)/*.o
