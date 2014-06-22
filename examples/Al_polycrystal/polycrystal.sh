#!/bin/bash

# This script uses the atomsk program to build
# a polycrystal of fcc aluminium.

rm -f Al_*

# 1: use mode "--create" to build a unit cell of fcc aluminium
atomsk --create fcc 4.02 Al Al_unitcell.xsf

# 2: write the parameters for the Voronoi construction
#    of the polycrystal into a file. Here we want a cubic supercell
#    of size 160x160x160 A^3, and 14 grains with random
#    positions and random crystallographic orientations
echo "box 160 160 160" > Al_voronoi.txt
echo "random 14" >> Al_voronoi.txt

# 3: generate the polycrystal
atomsk --polycrystal Al_unitcell.xsf Al_voronoi.txt Al_polycrystal.cfg -wrap

# The final system "Al_polycrystal.cfg" can be visualized with Atomeye.
