#!/bin/bash

# This script uses the atomsk program to build
# a polycrystal of silicon.

rm -f Si_*

# 1: use mode "--create" to build a unit cell of silicon with diamond structure
atomsk --create diamond 5.431 Si Si_unitcell.xsf

# 2: write the parameters for the Voronoi construction
#    of the polycrystal into a file. Here we want a cubic supercell
#    of size 160x160x160 A^3, and 14 grains with random
#    positions and random crystallographic orientations
echo "box 160 160 160" > Si_voronoi.txt
echo "random 14" >> Si_voronoi.txt

# 3: generate the polycrystal
atomsk --polycrystal Si_unitcell.xsf Si_voronoi.txt Si_polycrystal.cfg -wrap

# The final system "Si_polycrystal.cfg" can be visualized with Atomeye.
