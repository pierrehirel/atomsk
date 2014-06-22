# Makefile for compiling input modules of atomsk

SIN=$(wildcard in_*.f90)
OIN=$(SIN:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

input: $(OIN)
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)readin.o  -c readin.f90