# Makefile for compiling cbindings modules of Atomsk

SCMP=$(wildcard cbindings_*.f90)
OCMP=$(SCMP:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -cpp -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

cbindings: $(OCMP)
