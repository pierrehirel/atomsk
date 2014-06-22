#!/bin/bash
set -e

#USAGE:
# cfg_setA.sh <factor> <file1.cfg> [<file2.cfg> ...]

#DESCRIPTION:
# This script sets the scaling factor A to the
# given <factor> in one or many CFG files.
# This can be useful when Atomeye complains
# about "ATOM_COORDINATION_MAX = 24 exceeded".

if [[ -z "${@:2}" ]]; then
  printf "Usage: cfg_setA.sh <factor> <file1.cfg> [<file2.cfg> ...] \n"
else
  factor=$1
  isreal="$(echo "$factor" | egrep --null "^[[:digit:]]+(\.[[:digit:]]+)*$" )"
  if [[ -z "$isreal" ]]; then
    isreal="false"
  fi
  if [[ "$isreal" ]]; then
    FILES="${@:2}"
    for f in $FILES; do
      printf ">>> Setting A=$factor in $f ..."
      if [ -e $f ]; then
        if [ "$(head -c 6 $f)" != "Number" ]; then
          printf " Not a CFG file, skipping.\n"
        else
          sed -i "/A =/ c\A = $factor Angstrom (basic length-scale)" $f
          printf " Done.\n"
        fi
      else
        printf " File doesn't exist, skipping.\n"
      fi
    done
  else
    echo "X!X ERROR: incorrect scaling factor $factor."
  fi
fi
