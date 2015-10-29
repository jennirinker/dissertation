FC = gfortran

SRC = Crunch_Mods.f90
SRC += Crunch_Subs.f90
SRC += Crunch_IO.f90
SRC += Parser.f90
SRC += AzimAver.f90
SRC += CalcChan.f90
SRC += ExtEvent.f90
SRC += ExtrValu.f90
SRC += Filter.f90
SRC += LoadRose.f90
SRC += MoveAver.f90
SRC += PeakList.f90
SRC += PMF.f90
SRC += Rainflow.f90
SRC += Stats.f90
SRC += SetProg.f90
SRC += Crunch.f90

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
