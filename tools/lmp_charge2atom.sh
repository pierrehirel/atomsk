#!/bin/bash

#USAGE:
# lmp_charge2atom.sh <file1.lmp> <file2.lmp> ...

#DESCRIPTION:
# This script keeps only the first two and
# the last three columns of data from
# one or many files. It is intended
# for changing a LAMMPS data file from
# the atom_style "charge" into "atomic".

if [[ -z "$@" ]]; then
  printf "Keep only first two and last three columns of LAMMPS data files (*.lmp) to comply to 'atom_style atom'. \n"
  printf "Usage: lmp_charge2atom.sh <file1.lmp> [<file2.lmp> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    if [ -e $f ]; then
      printf ">>> Converting $f for 'atom_style atom'..."
      #Check that it is a LAMMPS file
      islmp=$(grep "atom types" $f | wc -l)
      if [ $islmp = 1 ]; then
        # If the keyword "Atoms" is not followed by "# atomic" add it
        charge=$(grep "Atoms" $f | grep "atomic" | wc -l)
        if [ $charge = 0 ]; then
          sed -i '/Atoms/ c\Atoms  # atomic' $f
        fi
        # Keep only columns 1, 2, and three last columns
        awk '{if(NR>9 && NF>4 && $1>0 && $2>0)
                {print $1 "\t" $2 "\t" $(NF-2) "\t" $(NF-1) "\t" $NF}
              else
                {print}
             }' $f >/tmp/temp.lmp
        mv -f /tmp/temp.lmp $f
        printf " Done.\n"
      else
        printf " Not a LAMMPS data file, skipping.\n"
      fi
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
