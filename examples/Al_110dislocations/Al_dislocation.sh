#!/bin/bash

# This script uses the atomsk program to build
# a supercell of aluminium and then insert two
# 1/2[110> edge dislocations of opposite Burgers vectors.
# The dislocations are constructed by insertion of a half-plane.

rm -f al_disloc*

# Lattice constant of fcc Al (A)
a0=4.02

# Burgers vector of 1/2[110] dislocation
b=$(echo "0.5*$a0*sqrt(2.0)" | bc -l)

atomsk --create fcc $a0 Al orient [0-11] [100] [011] \
       -duplicate 40 30 2 \
       -dislocation 0.251*BOX 0.251*BOX edge_add z y -$b 0.33 \
       -dislocation 0.751*BOX 0.751*BOX edge_add z y $b 0.33 \
       al_disloc.xsf cfg

#Note: dislocation stresses can be visualized with Atomeye
