#!/bin/bash

#USAGE:
# lmp_atom2charge.sh <file1.lmp> <file2.lmp> ...

#DESCRIPTION:
# This script keeps the first two columns
# of data, adds "0.0" as a third column, and
# copies the last three columns of data,
# to one or many files. It is intended
# for changing a LAMMPS data file from
# the atom_style "atomic" into "charge".

if [[ -z "$@" ]]; then
  printf "Add a column containing zeros to LAMMPS data files (*.lmp) to comply to 'atom_style charge'. \n"
  printf "Usage: lmp_atom2charge.sh <file1.lmp> [<file2.lmp> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    printf ">>> Converting $f for 'atom_style charge'..."
    if [ -e $f ]; then
      #Check that it is a LAMMPS file
      islmp=$(grep "atom types" $f | wc -l)
      if [ $islmp = 1 ]; then
        # Generate random file name
        of="$(mktemp ./XXXXXXXX.tmp)"
        # If the keyword "Atoms" is not followed by "# charge", add it
        atomic=$(grep "Atoms" $f | grep "charge" | wc -l)
        if [ $atomic = 0 ]; then
          sed -i '/Atoms/ c\Atoms  # charge' $f
        fi
        nrAtom=$(awk '/^Atom/{print NR}' $f)
        # Add a column of zero after keyword "Atom"
        awk -v n="${nrAtom}" '{if(NR>n && NF>4 && $1>0 && $2>0)
                {print $1 "\t" $2 "\t0.0\t" $(NF-2) "\t" $(NF-1) "\t" $NF}
              else
                {print}
             }' $f > ${of}
        mv -f ${of} $f
        printf " Done.\n"
      else
        printf " Not a LAMMPS data file, skipping.\n"
      fi
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
