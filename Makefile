MF90      = mpif90

#FFLAGS = -g -Wall -fbounds-check -O -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace \
            -fno-range-check 

FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check  


MPI       = -DMPI=1 

BINDIR    = ./bin
TARGET    = $(BINDIR)/mc3deconv
OBJS      = src/mc3decon.o src/params.o src/pt.o src/mcmc.o \
            src/init.o src/mt19937.o src/read_data.o \
	    src/calclogPPD.o src/conv.o src/temp.o src/output.o \
            src/pt_akuhara.o

all: $(TARGET) 

$(TARGET): $(OBJS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(MF90) $(FFLAGS) $(MPI) $^ -o $@ 

src/pt_akuhara.o: params.mod
src/mc3decon.o: params.mod
src/mcmc.o:     params.mod mt19937.mod
src/init.o:     params.mod mt19937.mod
src/read_data.o: params.mod
src/calclogPPD.o: params.mod
src/temp.o: mt19937.mod
src/output.o: params.mod

clean: 
	-rm -f *.mod src/*.o

#------------------------------------------------------------
# Pattern rule
#------------------------------------------------------------
$(OBJS): %.o: %.f90
	$(MF90) $(FFLAGS) -c $< -o $*.o 
%.mod: src/%.f90 src/%.o
	@:

