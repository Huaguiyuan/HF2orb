# ulimit -s unlimited
#source /opt/intel/fc/10.0.023/bin/ifortvars.sh
#source /opt/intel/idb/10.0.023/bin/idbvars.sh

#FC	= /usr/local/server/intel/intel_fc_90/bin/ifort
FC	= ifort 
FC = /opt/intel_fce_80/bin/ifort
FC = gfortran
# gfortran:
#DEBUG = -std=f2003 -g -fimplicit-none -Wall -pedantic 
#Intel:
OPTIM	= -O3
DEBUG   = -check all -traceback -g  -O0  
#gfortran:
OPTIM	= -O3
DEBUG  =  -g -O0 

OPTIONS = $(OPTIM)
OPTIONS = $(DEBUG)
#OPTIONS = 
#Stuttgart and UT, for Graz comment -i-static:
LLIBS   =  -llapack -lblas -i-static
#LLIBS   = -lblas -llapack -lblas -lg2c -i-static

OBJS    = diag_mod.o nrtype.o module.o coodinate.o akw.o ham.o prob.o spin_structure.o  random.o  ## add to the list with the .o
MODS    = nrtype.mod diag_mod.mod table_site.mod coodinate.mod akw.mod ham.mod prob.mod spin_structure.mod  random.mod 


main: main.o  $(MODS) Makefile
	$(FC) -o main main.o  $(OBJS) $(LLIBS)

main.o: main.f90 Makefile diag_mod.o module.o
	$(FC) $(OPTIONS) -c main.f90

diag_mod.o diag_mod.mod : diag_mod.f90 nrtype.mod Makefile  
	$(FC) $(OPTIONS) -c diag_mod.f90

module.o table_site.mod : module.f90  Makefile  
	$(FC) $(OPTIONS) -c module.f90

coodinate.o coodinate.mod : coodinate.f90  Makefile   
	$(FC) $(OPTIONS) -c coodinate.f90

akw.o akw.mod : akw.f90  Makefile   
	$(FC) $(OPTIONS) -c akw.f90

ham.o ham.mod : ham.f90  Makefile   
	$(FC) $(OPTIONS) -c ham.f90

prob.o prob.mod : prob.f90  Makefile  
	$(FC) $(OPTIONS) -c prob.f90

spin_structure.o spin_structure.mod : spin_structure.f90  Makefile   
	$(FC) $(OPTIONS) -c spin_structure.f90

random.o random.mod : random.f  Makefile   
	$(FC) $(OPTIONS) -c random.f



nrtype.o nrtype.mod : nrtype.f90 Makefile
	$(FC) $(OPTIONS) -c nrtype.f90





clean : 
	rm *.o *.mod main
