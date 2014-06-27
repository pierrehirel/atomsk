#!/bin/bash

# This script uses the atomsk program to build
# a 2-D polycrystal of rocksalt magnesium oxide (MgO).

rm -f MgO_*

# 1: use mode "--create" to build a unit cell of rocksalt MgO
atomsk --create rs 4.5 Mg O MgO_unitcell.xsf

# 2: write the parameters for the 2-D Voronoi construction
#    of the polycrystal into a file. We ask for a cubic supercell
#    of size 200x200x1 A^3, and 8 grains with random
#    positions and random crystallographic orientations.
#    Since the box is smaller than the unit cell along Z,
#    atomsk will automatically set it to the unit cell dimension.
echo "box 200 200 1" > MgO_voronoi.txt
echo "random 8" >> MgO_voronoi.txt

# 3: generate the polycrystal
atomsk --polycrystal MgO_unitcell.xsf MgO_voronoi.txt MgO_polycrystal.cfg

# The final system "MgO_polycrystal.cfg" can be visualized with Atomeye.
