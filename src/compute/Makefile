# Makefile for compiling compute modules of Atomsk

SCMP=$(wildcard compute_*.f90)
OCMP=$(SCMP:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -cpp -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

compute: $(OCMP)
