#!/bin/bash

#USAGE:
# lmp_ortho2tri.sh <file1.lmp> <file2.lmp> ...

#DESCRIPTION:
# This script converts a LAMMPS data file (.lmp)
# from orthogonal to triclinic, by inserting a line
# "0.0000000       0.00000000       0.00000000  xy xz yz"
# if necessary.

if [[ -z "$@" ]]; then
  printf "Add a line 'xy xz yz' into LAMMPS data files (*.lmp). \n"
  printf "Usage: lmp_ortho2tri.sh <file1.lmp> [<file2.lmp> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    printf ">>> Converting $f to triclinic..."
    if [ -e $f ]; then
      #Check that it is a LAMMPS file
      islmp=$(grep "atom types" $f | wc -l)
      if [ $islmp = 1 ]; then
        #Check if it is already triclinic
        if [ "$(grep 'xy xz yz' $f | wc -l)" -ne "0" ] ; then
          printf " Already triclinic, skipping.\n"
        else
          sed -i '/zhi/a\  0.0000000    0.00000000    0.00000000  xy xz yz' $f
          printf " Done.\n"
        fi
      else
        printf " Not a LAMMPS file, skipping.\n"
      fi
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
