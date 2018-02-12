#!/bin/bash
set -e

#USAGE:
# qepw_bohr.sh <file1.pw> [<file2.pw> ...]

#DESCRIPTION:
# This script sets the lines starting with "ATOMIC_POSITIONS"
# and "CELL_PARAMETERS" to "bohr" in one or several input file(s)
# for Quantum Espresso/pwscf.

if [[ -z "$@" ]]; then
  printf "Set lines 'ATOMIC_POSITIONS' and 'CELL_PARAMETER' to 'bohr' in Quantum Espresso PWscf files (*.pw). \n"
  printf "Usage: qepw_bohr.sh <file1.pw> [<file2.pw> ...] \n"
else
  FILES="$@"
  for f in $FILES; do
    printf ">>> Replacing 'angstroms' with 'bohr' in $f ..."
    if [ -e $f ]; then
      sed -i "/CELL_PARAMETERS/ c\CELL PARAMETERS bohr" $f
      sed -i "/ATOMIC_POSITIONS/ c\ATOMIC_POSITIONS bohr" $f
      printf " Done.\n"
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
