#!/bin/bash

# This script illustrates how to generate a cylindrical
# polycrystal of copper (Cu) using Atomsk.

# The radius of the cylinder (A)
R=100

# Size of the box (slightly larger than cylinder)
# Here the box is cubic, but you can use different sizes along X,Y,Z
d=$(echo "2*$R + 20" | bc -l)

# First, generate the unit cell of copper
atomsk --create fcc 3.61 Cu Cu_unitcell.xsf

# Second, write parameters for the polycrystal into a file
# Here we ask for 15 grains with random positions and orientations
echo "box $d $d $d" > poly.txt
echo "random 15" >> poly.txt

# Finally, run Atomsk with the mode "--polycrystal".
# Use option "-wrap" to make sure that all atoms are inside the box.
# Use option "-select" to select atoms that are out of the cylinder,
# then option "-rmatom" to remove selected atoms.
atomsk --polycrystal Cu_unitcell.xsf poly.txt \
       -wrap \
       -select out cylinder Z 0.5*box 0.5*box $R \
       -rmatom select \
       Cu_polycrystal.cfg

# The final file "Cu_polycrystal.cfg" can be visualized with Atomeye.
