#!/bin/bash

rm -f al_supercell.*

# This script shows how to use atomsk to build a
# supercell of fcc aluminium with two surface steps.
# First the unit cell with the crystallographic
# orientation X=[0-11], Y=[100], Z=[011] is created,
# then it is duplicated to form a 40x20x4 supercell,
# then atoms in boxes near the top surface are removed.
# The final result is written to XYZ, XSF and CFG formats
# for visualization, and to the LAMMPS data format.

atomsk --create fcc 4.02 Al orient [0-11] [100] [011] \
       -duplicate 40 20 4 \
       -select in box 0 BOX-4 0  20 BOX BOX \
       -cut above -INF X \
       -select in box BOX-20 BOX-4 0  BOX BOX BOX \
       -cut above -INF X \
       al_supercell.xyz cfg xsf lmp
