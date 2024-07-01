#!/bin/bash

#USAGE:
# dat_rm0.sh <file1.dat> <file2.dat> ...

#DESCRIPTION:
# This script reads "dat" file(s) containing
# six columns of real numbers, and removes lines
# where columns 4 to 6 only contain zeros.

if [[ -z "$@" ]]; then
  printf "Remove lines where columns 4-6 contain only zeros in data files (*.dat). \n"
  printf "Usage: dat_rm0.sh <file1.dat> [<file2.dat> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    printf ">>> Removing lines with zero vectors in $f ..."
    if [ -e $f ]; then
      # Generate random file name
      of="$(mktemp ./XXXXXXXX.tmp)"
      awk '{if($4!=0.0 && $5!=0.0 && $6!=0.0) print $0 }' $f >> ${of}
      mv -f ${of} $f
      printf " Done.\n"
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
