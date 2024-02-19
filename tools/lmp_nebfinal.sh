#!/bin/bash

#USAGE:
# lmp_nebfinal.sh <file1.lmp> <file2.lmp> ...

#DESCRIPTION:
# This script converts a LAMMPS data file (.lmp) into
# a format suitable to serve as final image in a NEB calculation.

if [[ -z "$@" ]]; then
  printf "Convert LAMMPS data files into neb final image. \n"
  printf "Usage: lmp_nebfinal.sh <file1.lmp> [<file2.lmp> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    printf ">>> Converting $f into final NEB image..."
    if [ -e $f ]; then
      #Check that it is a LAMMPS file
      islmp=$(grep "atom types" $f | wc -l)
      if [ $islmp = 1 ]; then
        # Generate random file name
        of="$(mktemp ./XXXXXXXX.tmp)"
        # Write number of atoms in first line
        grep atoms$ $f | awk '{print $1}' > ${of}
        #In lines after "Atoms", copy field 1 (atom-id), and last 3 fields (x y z)
        awk '/^Atoms /,EOF { if(NF>3) print $1 "\t" $(NF-2) "\t" $(NF-1) "\t" $NF }' $f >> ${of}
        # Replace original file
        mv ${of} ${f}
        printf " Done.\n"
      else
        printf " Not a LAMMPS data file, skipping.\n"
      fi
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
