#!/bin/bash

# This script shows how to use atomsk to build
# a supercell of silicon (diamond structure) with
# a surface step

rm -f Si.*

# First the unit cell is created with the crystallographic
# orientation X=[-110], Y=[001], Z=[110], then it is expanded
# to form a supercell. Atoms close to the top surface are selected,
# and selected atoms that are above the (-111) plane located 70 A
# away from the origin are deleted, so a surface step remains.
# The final result is output to XYZ, CFG and XSF formats for
# visualization, and to LAMMPS data format
atomsk --create diamond 5.431 Si orient [-110] [001] [110] -e 20 14 3 -cut above BOX-1.5 y -select above BOX-12.5 y -cut above 70 [-111] xyz cfg xsf lmp
