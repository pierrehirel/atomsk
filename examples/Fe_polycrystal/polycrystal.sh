#!/bin/bash

# This script uses the atomsk program to build
# a polycrystal of bcc iron.

rm -f Al_*

# 1: use mode "--create" to build a unit cell of bcc iron
atomsk --create bcc 2.9 Fe Fe_unitcell.xsf

# 2: write the parameters for the Voronoi construction
#    of the polycrystal into a file. Here we want a cubic supercell
#    of size 160x160x160 A^3, and 12 grains with random
#    positions and random crystallographic orientations
echo "box 160 160 160" > Fe_voronoi.txt
echo "random 14" >> Fe_voronoi.txt

# 3: generate the polycrystal
atomsk --polycrystal Fe_unitcell.xsf Fe_voronoi.txt Fe_polycrystal.cfg -wrap

# The final system "Fe_polycrystal.cfg" can be visualized with Atomeye.
