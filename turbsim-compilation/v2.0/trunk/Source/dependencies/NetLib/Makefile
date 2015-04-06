FC = gfortran

SRC_F = fftpack4.1.f
SRC   = NWTC_FFTPACK.f90
SRC  += NWTC_LAPACK.f90

OBJ = $(SRC:.f90=.o)
OBJ_F = $(SRC_F:.f=.o)

ifdef BASEDIR
	INCLUDE = $(BASEDIR)/include
	MODDIR = $(BASEDIR)/include
else
	MODDIR = ./
endif

ifndef OUTDIR
	OUTDIR = ./
endif

all: checkdir $(SRC:.f90=.d) $(addprefix $(OUTDIR)/,$(OBJ)) $(SRC_F:.f=.d) $(addprefix $(OUTDIR)/,$(OBJ_F))

$(OUTDIR)/%.o: %.f*
	$(FC) $(FFLAGS) -c $< -o $@ -J$(MODDIR)	-I$(INCLUDE)

%.d: %.f*
	@set -e; rm -f $@; \
	$(FC) -E -cpp -MM $< -MF $@.$$$$ > /dev/null; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

checkdir:
	@mkdir -p $(OUTDIR); \
	mkdir -p $(MODDIR);

-include $(SRC:.f90=.d)
-include $(SRC_F:.f=.d)

.PHONY: clean checkdir
clean:
	rm -rf *.d *.s *.o
