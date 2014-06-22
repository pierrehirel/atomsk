#!/bin/bash

#USAGE:
# dat_rm0.sh <file1.dat> <file2.dat> ...

#DESCRIPTION:
# This script reads "dat" file(s) containing
# six columns of real numbers, and removes lines
# where columns 4 to 6 only contain zeros.

if [[ -z "$@" ]]; then
  printf "Usage: dat_rm0.sh <file1.dat> [<file2.dat> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    printf ">>> Removing lines with zero vectors in $f ..."
    if [ -e $f ]; then
      awk '{if($4!=0.0 && $5!=0.0 && $6!=0.0) print $0 }' $f >>/tmp/temp.dat
      mv -f /tmp/temp.dat $f
      printf " Done.\n"
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
