# Makefile for compiling output modules of Atomsk

SOU=$(wildcard out_*.f90)
OOU=$(SOU:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

output: $(OOU)
	$(FORTRAN) $(FFLAGS) -cpp -o ..$(SEP)$(OBJ)$(SEP)writeout.o -c writeout.f90
