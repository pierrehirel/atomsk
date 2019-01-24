# Makefile for compiling common modules of atomsk

# List of source files (order is very important here)
# Localization (l10n) modules: additional languages can be added here
L10N=messages_FR.f90 messages_DE.f90

SGE=comv.f90              \
    constants.f90         \
    math.f90              \
    functions.f90         \
    atoms.f90             \
    subroutines.f90       \
    sort.f90              \
    display_messages.f90  \
    messages_EN.f90       \
    $(L10N)               \
    messages.f90          \
    expreval.f90          \
    resize.f90            \
    spacegroups.f90       \
    neighbors.f90         \
    files.f90             \
    symops.f90            \
    guess_format.f90      \
    read_cla.f90          \
    determine_H.f90       \
    readconf.f90          \
    qepw_ibrav.f90        \
    in_stl.f90

# Object files
OGE=$(SGE:.f90=.o)

VPATH=..$(SEP)$(OBJ)
#Default rule to make targets
# NOTE: compilation of "spacegroups.f90" with GNU Fortran (gfortran)
#      and with the flag "-ftree-pre" (activated with -O2) causes a compiler error.
#      Therefore "-fno-tree-pre" is used to disable it.
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -fno-tree-pre -cpp -o ..$(SEP)$(OBJ)$(SEP)$@  -c $< 

general: $(OGE)
