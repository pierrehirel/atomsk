#!/bin/bash
set -e

#USAGE:
# dat_mulvec.sh <factor> <file1.dat> <file2.dat> ...

#DESCRIPTION:
# This script reads "dat" file(s) containing
# six columns of real numbers, and multiplies
# columns 4 to 6 by the given <factor>.
# The purpose is to scale the vectors in dat files
# containing atom positions followed by a vector
# (x y z vx vy vz).

factor=0.0

if [[ -z "${@:2}" ]]; then
  printf "Usage: dat_mulvec.sh <factor> <file1.dat> [<file2.dat> ...] \n"
else
  factor=$1
  FILES="${@:2}"
  for f in $FILES
  do
    printf ">>> Multiplying vectors by $factor in $f ..."
    if [ -e $f ]; then
      awk -v fac=$factor '{print $1 " " $2 " " $3 " " fac*$4 " " fac*$5 " " fac*$6}' $f >>/tmp/temp.dat
      mv -f /tmp/temp.dat $f
      printf " Done.\n"
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
