
# General Makefile for compiling Atomsk.

#####################################################################
#
#              THIS MAKEFILE SHOULD NOT BE MODIFIED !!
#
#####################################################################
# (except for debugging/development purposes)
# If you want to use your own compiler/flags/options
# it is recommended to create a new Makefile defining the
# necessary variables and calling this one. See e.g. the
# files "Makefile.ifort" or "Makefile.g95" or "Makefile.windows".

# REQUIREMENTS:
# - Atomsk uses some intrinsic functions of the Fortran 2008 standard,
#   therefore the compiler must be compatible with this standard.
#   If you use GNU Fortran (gfortran), version 9 or greater is required.
#

# USAGE:
#    make atomsk
# Or, to use parallel compilation:
#    make -j3 atomsk
# If you wish to use another makefile, e.g. Intel compiler:
#    make -f Makefile.ifort atomsk

export
# Version of the program
VER="master"

# Name of the final binary
BIN=atomsk
# Name of final library (should be libatomsk.ext)
ifndef ATOMSKLIB
  ATOMSKLIB=libatomsk.so
endif

# Subdirectory where objects will be stored during compilation
OBJ=OBJ

# Name of the first Makefile called during current run
MF=$(word 1, $(MAKEFILE_LIST) )

# Define folder separator <<< ARCHITECTURE-DEPENDENT
ifdef ComSpec
  # Microsoft(R) Windows(R) uses the backslash
  SEP=\\
  MKDIR=if not exist $(OBJ) mkdir
else
  # UNIX or Linux use the forward slash
  SEP=/
  MKDIR=mkdir -p
endif

mode=static

# COMPILER AND OPTIONS
# if compiler was not defined before, the following defaults are used
ifndef FORTRAN
  # GNU Fortran as the default Fortran compiler
  # NOTE: use gfortran version 9 or greater, or any compiler compatible with Fortran 2008.
  FORTRAN=gfortran
endif
  
# Comment the following line to disable parallelization with OpenMP
OPENMP=-DOPENMP -fopenmp

# -O2 gives good optimization/performance
# -fno-backslash prevents backslash characters to be misunderstood
# -I is used to specify folders where the compiler should look for *.mod files
# -J is used to specify the folder where the compiler will store the *.mod files
  
LDFLAGS+=-L.
OPTIONS=

# Add -fpic to compile Atomsk as a shared library
ifeq ($(shared),yes)
  OPTIONS+=-fPIC
endif
  
# Un-comment the following line to use 8-bits integers by default
# (number of atoms up to 9x10^18 instead of 2x10^9)
ifeq ($(bigint),yes)
  OPTIONS+=-fdefault-integer-8
endif
  
# DEBUGGING FLAGS
ifeq ($(debug),yes)
  OPTIONS+=-fbacktrace -g -Wall -Og #-fbounds-check
endif

modes compute aux caux input options output cbindings 110n general objdir : FFLAGS=-O2 $(OPENMP) -fno-backslash -I..$(SEP)$(OBJ) -J..$(SEP)$(OBJ) $(OPTIONS)
all atomsk lib clib : FFLAGS=-O2 $(OPENMP) -fno-backslash -I$(OBJ) $(OPTIONS)

ifndef LAPACK
  # Path to the LAPACK library
  # If it is installed system-wide just use "-llapack"
  LAPACK=-llapack
endif

# DEFINITION OF PATHS:
# (for UNIX/Linux systems only)
# (do NOT insert a / at the end of the paths)
ifndef INSTPATH
  INSTPATH=/usr/local
endif
# BINPATH: the binary will be copied there
ifndef BINPATH
  BINPATH=$(INSTPATH)/bin
endif
# DOCUMENTATION PATH: the content of the "doc" directory
# will be copied in the DOCPATH/atomsk folder
ifndef DOCPATH
  DOCPATH=$(INSTPATH)/share/doc
endif
# MAN PATH: the file man/atomsk.1.gz
# will be copied in the MANPATH folder
ifndef MPATH
  MPATH=$(INSTPATH)/share/man/man1
endif
ifndef CONFPATH
  CONFPATH=$(INSTPATH)/etc
endif

# Name for the tar.gz archive (do NOT add extension like .tar.gz)
ARCNAME=atomsk_$(VER)

# Define folders where objects are found
VPATH=$(OBJ)
# Define folders where source files are found
vpath %.f90 include:input:options:output:modes

# Define default rule to make targets
%.o : %.f90
	$(FORTRAN) $(FFLAGS) -o $(OBJ)$(SEP)$@  -c $<


# TARGETS

# Default target is to display help
.DEFAULT: help
help: message
	@echo "  Available targets for make:"
        ifeq ($(MF),Makefile)
	  @echo "- make atomsk       compile the program using ${FORTRAN}"
	  @echo "- make help         display this message"
        else
	  @echo "- make -f ${MF} atomsk"
	  @echo "                    compile the program using ${FORTRAN}"
	  @echo "- make -f ${MF} help"
	  @echo "                    display this message"
        endif
	@echo "- make install      install the program in ${BINPATH}"
	@echo "                    (requires SuperUser rights)"
	@echo "- make clean        remove compiled files from current directory"
	@echo "- make uninstall    uninstall the program from ${BINPATH}"
	@echo "                    (requires SuperUser rights)"
	@echo "- make tar          create a tar.gz archive of the source code"
	@echo "                    (requires tar and gzip programs)"
	@echo ""

message:
	@echo "       o---o     ___________"
	@echo "      o---o|     A T O M S K"
	@echo "      |   |o"
	@echo "      o---o      Version ${VER}"
	@echo ""

# All = compile the code
all: atomsk


# Compile the program
# Note that the preprocessor must be called (-cpp)
atomsk: message modes compute general
	$(FORTRAN) $(FFLAGS) -cpp -o $(BIN) $(OBJ)$(SEP)*.o atomsk.f90 $(LDFLAGS) $(LAPACK)
	@echo ""
	@echo "    \o/ Compilation was successful!"
	@echo ""
	@echo "    <i> To install Atomsk system-wide, you may now run:"
	@echo "          sudo make install"
	@echo ""

# Compile the program as a library
# Note that the preprocessor must be called (-cpp)
lib: message compute modes
	$(FORTRAN) $(FFLAGS) -cpp -shared -o $(ATOMSKLIB) $(OBJ)$(SEP)*.o $(LAPACK)
	$(FORTRAN) $(FFLAGS) -cpp -o $(BIN) atomsk.f90 $(LDFLAGS) -latomsk $(LAPACK)
	@echo ""
	@echo "    \o/ Atomsk was compiled as a shared library file, 'atomsk.so'."
	@echo ""

clib: message compute modes cbindings
	$(FORTRAN) $(FFLAGS) -cpp -shared -o $(ATOMSKLIB) $(OBJ)$(SEP)*.o $(LAPACK)
	$(FORTRAN) $(FFLAGS) -cpp -o $(BIN) atomsk.f90 $(LDFLAGS) -latomsk $(LAPACK)
	@echo ""
	@echo "    \o/ Atomsk was compiled as a shared library file with cbindings, 'atomsk.so'."
	@echo ""

# Compile modes modules
modes: aux compute FORCE
	$(MAKE) $(MFLAGS) -C modes

# Compile compute modules
compute: general FORCE
	$(MAKE) $(MFLAGS) -C compute

# Compile input, options and output - this can be run in parallel
aux: general
	@echo "..> Compiling input, options and output modules..."
	$(MAKE) $(MFLAGS) caux

# Compile input, options and output - this can be run in parallel
caux: input options output

# Compile input modules
input: FORCE
	$(MAKE) $(MFLAGS) -C input

# Compile options modules
options: FORCE
	$(MAKE) $(MFLAGS) -C options

# Compile output modules
output: FORCE
	$(MAKE) $(MFLAGS) -C output

# Compile output modules
cbindings: FORCE
	$(MAKE) $(MFLAGS) -C cbindings

# Compile common modules ("include") - NOTE: this MUST NOT be run in parallel (-j1)
l10n: general
	$(MAKE) $(MFLAGS) -j1 -C include

# Compile common modules ("include") - NOTE: this MUST NOT be run in parallel (-j1)
general: objdir
	$(MAKE) $(MFLAGS) -j1 -C include

# Create OBJ directory (the - sign indicates to ignore any error)
objdir:
	-$(MKDIR) $(OBJ)

# Dummy target to force the compilation of input, options, output
FORCE:


# Clean current directory
.PHONY: clean
clean:
	@rm -rf $(OBJ) input$(SEP)*~ *.mod *.o .$(SEP)*$(SEP)*.o .$(SEP)*$(SEP)*.mod *.x *.exe *~ *$(SEP)*~ $(BIN)
	@echo ">>> All binary and temporary files were deleted."


# Create an archive of the package (only on UNIX and Linux systems)
# NOTE: Source files are modified with sed to match target version defined above
tar:
	@sed -i "/CHARACTER(LEN=24),PARAMETER:: version =/ c\CHARACTER(LEN=24),PARAMETER:: version = '${VER}'" include/globalvar.f90
	@sed -i "/<p>&copy;/ c\<p>&copy; Pierre Hirel 2010 - Version ${VER}</p>" ../doc/en/index.html
	@sed -i "/<p>&copy;/ c\<p>&copy; Pierre Hirel 2010 - Version ${VER}</p>" ../doc/fr/index.html
	@tar --exclude-vcs -czf ../$(ARCNAME).tar.gz                            \
	          ../doc/*/*.html ../doc/*/*.css   ../doc/*/*.svg               \
	          ../doc/*/*.png  ../doc/index.html                             \
	          ../etc/atomsk.conf                                            \
	          ../examples/                                                  \
	          ../man/atomsk  ../src/Makefile                                \
	          ../src/atomsk.f90 ../src/include/*.f90 ../src/input/*.f90     \
	          ../src/options/*.f90 ../src/output/*.f90 ../src/modes/*.f90   \
	          ../src/compute/*.f90 ../src/Makefile.* ../src/*/Makefile      \
	          ../tools/*  ../README ../LICENSE ../CHANGELOG
	@echo ">>> The program was successfully archived in ../$(ARCNAME).tar.gz"


# Install the program (only on UNIX and Linux systems)
install:
	@echo ""
	@echo "  >>> Installation of Atomsk version $(VER):"
	@echo "  ..> The program was successfuly installed in $(BINPATH)"
	@echo "      To run it, enter 'atomsk' in a terminal."
	@mkdir -p $(DOCPATH)/atomsk
	@rm -rf $(DOCPATH)/atomsk/*
	@cp -rf ../doc/* $(DOCPATH)/atomsk/
	@chmod -R a+r $(DOCPATH)/atomsk/
	@echo "  ..> The html documentation was installed. You may read it by entering the"
	@echo "      following address in your Web browser: ${DOCPATH}/atomsk/index.html"
	@mkdir -p $(MPATH)
	@gzip -c ../man/atomsk >$(MPATH)/atomsk.1.gz
	@echo "  ..> To read the manual, type 'man atomsk'."
	@chmod +x ../tools/*.sh
	@mkdir -p $(BINPATH)
	@cp -f ../tools/*.sh $(BINPATH)
	@mkdir -p $(CONFPATH)
	@cp -rf ../etc/atomsk.conf $(CONFPATH)
	@mkdir -p $(BINPATH)
	@cp -rf $(BIN) $(BINPATH)
	@echo ""


# Uninstall the program from the system (only on UNIX and Linux systems)
uninstall:
	@echo ""
	@echo "  >>> Removing Atomsk and its documentation from $(INSTPATH)"
	@rm -rf $(BINPATH)/$(BIN)
	@rm -rf $(DOCPATH)/atomsk
	@rm -rf $(MPATH)/atomsk.1.gz
	@echo "      Uninstallation successfull."
	@echo "      NOTE: configuration file(s) 'atomsk.conf' may still exist in folders /etc/ and /home/.config/"
	@echo ""
