# Makefile for compiling common modules of atomsk

# Modes modules
# Note: "mode_interactive.f90" and "mode_list.f90" depend on "mode_normal.f90"
#      and must NOT be placed in this list
SMO=mode_difference.f90     \
    mode_normal.f90         \
    mode_create.f90         \
    mode_merge.f90          \
    1ia_dlpoly_history.f90  \
    1ia_qe_out.f90          \
    1ia_lmc.f90             \
    1ia_xsf.f90             \
    1ia_xyz.f90             \
    mode_unwrap.f90         \
    mode_allinone.f90       \
    mode_epola.f90          \
    mode_electricdipoles.f90\
    mode_rdf.f90            \
    mode_polycrystal.f90    \
    mode_interpolate.f90    \
    mode_nye.f90            \
    mode_average.f90        \
    mode_density.f90        \
    mode_centrosym.f90

# Object files
OMO=$(SMO:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)$@  -c $<

modes: $(OMO)
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)mode_list.o -c mode_list.f90
	$(FORTRAN) $(FFLAGS) -o ..$(SEP)$(OBJ)$(SEP)modes.o -c modes.f90
	$(FORTRAN) $(FFLAGS) -cpp -o ..$(SEP)$(OBJ)$(SEP)mode_interactive.o -c mode_interactive.f90
