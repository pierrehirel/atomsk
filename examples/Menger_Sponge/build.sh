#!/bin/bash

rm -f unitcell.xsf sponge_*.cfg

# Create unit cell of aluminum
# NOTE: you may replace that by any lattice of your liking
atomsk --create fcc 4.02 Al unitcell.xsf

# Create grid file
# A Menger sponge is created by cutting each face of a cube
# in nine equal parts, and removing the central part.
# We create a grid containing nine elements, where only
# the central element is 1 and will be selected
echo "0 0 0" >  grid_sponge.txt
echo "0 1 0" >> grid_sponge.txt
echo "0 0 0" >> grid_sponge.txt

# Initial input file name
f1="unitcell.xsf"
# Initial output file name
f2="sponge_0.cfg"
# Initial duplication factor (for unit cell)
n=2

# Loop to generate Menger sponge
# The same grid is applied to all faces of the cube.
# You may play with the max. value of i, but know
# that the number of atoms increases EXPONENTIALLY with imax!
# Sponges of successive sizes are saved in files named "sponge_$i.cfg"
for ((i=1;i<=4;i++)) ; do

  atomsk $f1                         \
         -dup $n $n $n               \
         -select grid grid_sponge.txt \
         -rmatom select -swap X Z    \
         -select grid grid_sponge.txt \
         -rmatom select -swap Z Y    \
         -select grid grid_sponge.txt \
         -rmatom select              \
         $f2

  f1=$f2
  f2="sponge_${i}.cfg"
  n=3

done

