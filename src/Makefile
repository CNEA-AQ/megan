.SUFFIXES: .o .f90

FC  =ifx
#NETCDF= /glade/u/apps/casper/24.12/spack/opt/spack/netcdf/4.9.2/openmpi/5.0.6/oneapi/2024.2.1/tlt6
#NETCDF_LIB =/glade/u/apps/casper/24.12/spack/opt/spack/netcdf/4.9.2/openmpi/5.0.6/oneapi/2024.2.1/tlt6/lib
#NETCDF_INC =/glade/u/apps/casper/24.12/spack/opt/spack/netcdf/4.9.2/openmpi/5.0.6/oneapi/2024.2.1/tlt6/include
LIBS=-L${NETCDF_LIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm
INC =-I${NETCDF_INC}

#FCFLAGS=-O2 -ffree-line-length-none -w #-Wunused #-Wall 
#LDFLAGS=-O2
#DEBUG = -g -traceback -check all -fpe0 -O0
DEBUG = -traceback
FCFLAGS= $(DEBUG) -O2 -extend-source #-Wunused #-Wall 
LDFLAGS= $(DEBUG) -O2

OBJS= constants_module.o misc_definitions_module.o bio_types.o area_mapper_grw.o \
      datetime_module.o prep_megan.o megan_voc.o megan_nox.o megan_driver.o #bdsnp.o
EXE =../exe/megan_v3.3.exe

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
clean:
	rm -f *.o *.mod
