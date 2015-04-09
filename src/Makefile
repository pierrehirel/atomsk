# Makefile for compiling atomsk on a Unix/Linux machine.
# Tested and works with GNU make and gfortran >4.4

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

# USAGE:
#    make atomsk
# Or, to use parallel compilation:
#    make -j 3 atomsk

export

# Version of the program
VER=0.8.1

# Name of the final binary
BIN=atomsk

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

# COMPILER AND OPTIONS
# if compiler was not defined before, the following defaults are used
ifndef FORTRAN
  # GNU Fortran as the default Fortran compiler
  # NOTE: use gfortran version 4.4 or greater. Versions prior to gfortran 4.4 may cause problems
  FORTRAN=gfortran

  # -O2 gives good optimization/performance
  # -fno-backslash prevents backslash characters to be misunderstood
  # -I is used to specify folders where the compiler should look for *.mod files
  # -J is used to specify the folder where the compiler will store the *.mod files
  FFLAGS=-O2 -fno-backslash -I..$(SEP)$(OBJ) -J..$(SEP)$(OBJ)
  
  # DEBUGGING FLAGS
  #FFLAGS+=-fbacktrace -fbounds-check -g -Wall
endif

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

# Name for the tar.gz archive (do NOT add extension)
ARCNAME=atomsk_b$(VER)

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
help:
	@echo "       o---o                     O"
	@echo "      o---o|   _____________     |"
	@echo "      |   |o    A T O M S K     ,0;-o"
	@echo "      o---o                   o'   'O"
	@echo "Available targets for make:"
        ifeq ($(MF),Makefile)
	  @echo "- make help         display this message"
	  @echo "- make atomsk       compile the program using ${FORTRAN}"
        else
	  @echo "- make -f ${MF} help"
	  @echo "                    display this message"
	  @echo "- make -f ${MF} atomsk"
	  @echo "                    compile the program using ${FORTRAN}"
        endif
	@echo "- make install      install the program in ${BINPATH}"
	@echo "                    (requires SuperUser rights)"
	@echo "- make clean        remove compiled files from current directory"
	@echo "- make uninstall    uninstall the program from ${BINPATH}"
	@echo "                    (requires SuperUser rights)"
	@echo "- make tar          create a tar.gz archive of the source code"
	@echo "                    (requires tar and gzip programs)"
	@echo ""


# All = compile the code
all: atomsk


# Compile the program
# Note that the preprocessor must be called (-cpp)
atomsk: modes
	$(FORTRAN) $(FFLAGS) -cpp -o $(BIN) $(OBJ)$(SEP)*.o atomsk.f90 -I $(OBJ) $(LAPACK)
	@echo "\o/ Compilation was successful!"

# Compile modes modules
modes: aux
	$(MAKE) $(MFLAGS) -C modes

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

# Compile common modules - this MUST NOT be run in parallel
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
tar:
	@tar -czf ../$(ARCNAME).tar.gz                                          \
	          ../doc/*/*.html ../doc/*/*.css   ../doc/*/*.svg               \
	          ../doc/*/*.png  ../doc/index.html                             \
	          ../etc/atomsk.conf                                            \
	          ../examples/                                                  \
	          ../man/atomsk                                                 \
	          ../src/atomsk.f90 ../src/include/*.f90 ../src/input/*.f90     \
	          ../src/options/*.f90 ../src/output/*.f90 ../src/modes/*.f90   \
	          ../src/Makefile ../src/Makefile.* ../src/*/Makefile           \
	          ../tools/*                                                    \
	          ../README ../LICENSE ../CHANGELOG --exclude-vcs
	@echo ">>> The program was successfully archived in ../$(ARCNAME).tar.gz"


# Install the program (only on UNIX and Linux systems)
install:
	@cp -rf ../etc/atomsk.conf /etc/
	@mkdir -p $(BINPATH)
	@cp -rf $(BIN) $(BINPATH)
	@echo ">>> Program was successfuly installed in $(BINPATH)"
	@mkdir -p $(DOCPATH)/atomsk
	@rm -rf $(DOCPATH)/atomsk/*
	@cp -rf ../doc/* $(DOCPATH)/atomsk/
	@chmod -R a+r $(DOCPATH)/atomsk/
	@echo ">>> Html documentation is available from $(DOCPATH)/atomsk/index.html"
	@mkdir -p $(MPATH)
	@gzip -c ../man/atomsk >$(MPATH)/atomsk.1.gz
	@echo ">>> The man page is available from $(DOCPATH)/atomsk/"
	@chmod +x ../tools/*.sh
	@cp -f ../tools/*.sh $(BINPATH)


# Uninstall the program from the system (only on UNIX and Linux systems)
uninstall:
	@rm -rf $(BINPATH)/$(BIN)
	@rm -rf $(DOCPATH)/atomsk
	@rm -rf $(MPATH)/atomsk.1.gz
	@echo ">>> Program atomsk was uninstalled from $(INSTPATH)."
