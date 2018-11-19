#!/bin/bash

#USAGE:
# dat2csv.sh <file1.dat> <file2.dat> ...

#DESCRIPTION:
# This script reads "dat" file(s) containing columns
# of data, and introduces commas to separate columns
# so as to comply to Comma-Separated Values (CSV) format.

if [[ -z "$@" ]]; then
  printf "Convert space-separated data files (*.dat) into comma-separated values format (*.csv). \n"
  printf "Usage: dat2csv.sh <file1.dat> [<file2.dat> ...] \n"
else
  FILES="$@"
  for f in $FILES
  do
    printf ">>> Converting $f to CSV format ..."
    if [ -e $f ]; then
      # Set name of output file
      of="$(basename "$f" .dat).csv"
      # Print all fields, separated with a comma, into output file
      awk 'BEGIN {OFS=","}; {$1=$1 ; print $0}' $f > $of
      printf " Done.\n"
    else
      printf " File doesn't exist, skipping.\n"
    fi
  done
fi
