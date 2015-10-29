FC = gfortran

SRC = SharedInflowDefs.f90
SRC += HHWind.f90
SRC += FFWind.f90
SRC += FDWind.f90
SRC += CTWind.f90
SRC += UserWind.f90
SRC += HAWCWind.f90
SRC += InflowWindMod.f90

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
	rm -rf *.d *.s *.o
