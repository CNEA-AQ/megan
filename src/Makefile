.SUFFIXES: .o .f90

FC  =gfortran
LIBS=-L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm
INC =-I/usr/include

FCFLAGS=-O2 -ffree-line-length-none #-Wunused -Wall#
LDFLAGS=-O2

OBJS=prep_megan.o megan_voc.o megan_nox.o megan_driver.o #bdsnp.o
#OBJS=megan_voc.o megan_nox.o megan_driver.o #bdsnp.o
EXE =../exe/megan_v3.3.exe

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
clean:
	rm -f *.o *.mod
