# Makefile for compiling options modules of atomsk

SOP=$(wildcard opt_*.f90)
OOP=$(SOP:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

options: $(OOP)
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)options.o  -c options.f90