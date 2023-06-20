# Makefile for compiling options modules of atomsk

SOP=disloc_iso.f90   \
    disloc_aniso.f90 \
    disloc_loop.f90  \
    $(wildcard opt_*.f90)
OOP=$(SOP:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -cpp -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

options: $(OOP)
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)options.o  -c options.f90
