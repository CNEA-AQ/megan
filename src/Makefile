.SUFFIXES: .o .f90

FC  =gfortran
LIBS=-L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm
INC =-I/usr/include

FCFLAGS=-O2 -ffree-line-length-none #-Wall#-Wunused 
LDFLAGS=-O2

OBJS=module_proj.o meganv32.o main_megan.o #module_utils.o 
EXE =meganv3.2.exe

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
clean:
	rm -f *.o *.mod                             
