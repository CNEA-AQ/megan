.SUFFIXES: .o .f90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
INC  = -I/usr/include

FCFLAGS = -O2 -ffree-line-length-none -Wunused 
LDFLAGS = -O2
EXE  = prep_megan.exe
OBJS = prep_megan.o 

prepmegan4cmaq: $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)

%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@

clean:
	rm -f *.o *.mod
