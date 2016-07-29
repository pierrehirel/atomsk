#!/bin/bash

# 1/2[110](1-11) edge dislocation in fcc aluminium
# built by superimposing two supercells.
# The top supercell has one more atomic plane than the
# bottom supercell, and both have the same length.
# X = direction of glide (and parallel to Burgers vector)
# Y = normal to glide plane
# Z = direction of dislocation line
# NOTE: this system is designed to be periodic along X and Z, but not along Y.
# NOTE2: the system constructed here is *not* relaxed and does not contain a dislocation!
#      The dislocation will only appear after relaxation with a proper atomistic potential.

rm -f Al_unitcell.xsf Al_top.* Al_bottom.* Al.*

# Lattice constant of fcc Al
a=4.02

# Define box size
boxxT=60               # Number of unit lattices along X in top half crystal
let boxxB=$boxxT-1     # Number of unit lattices along X in bottom half crystal
boxy=8                 # Height of top (or bottom) crystal along Y
boxz=1                 # Number of unit lattices along dislocation line

# Define applied deformation
s1=$(echo "-0.5/$boxxT" | bc -l)
s2=$(echo "0.5/$boxxB" | bc -l)

# Build unit cell of fcc Al oriented
atomsk --create fcc $a Al orient [110] [1-11] [1-1-2] Al_unitcell.xsf

# Build top half crystal: compressed along X
atomsk Al_unitcell.xsf -dup $boxxT $boxy $boxz -deform x $s1 0.0 Al_top.xsf
# Build bottom half crystal: elongated by $a/2 along X
atomsk Al_unitcell.xsf -dup $boxxB $boxy $boxz -deform x $s2 0.0 Al_bottom.xsf

# System with dislocation: merge top and bottom along Y
atomsk --merge y 2 Al_bottom.xsf Al_top.xsf Al.cfg lmp exyz

