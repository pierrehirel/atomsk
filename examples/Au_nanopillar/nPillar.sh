#!/bin/bash

# This script shows how to use atomsk to build
# a nanopillar of gold in one command-line.

rm -f nanopillar*

# The mode "--create" is used to generate the unit cell of fcc gold.
# Then it is expanded (-e) to form a 40x40x40 supercell.
# Then all the atoms outside of the cylinder placed at the middle of
# the XY plane and with a diameter of 50.0 A are selected (-select)
# and deleted (-cut). Note that a slab of gold of thickness dz=10 A
# is kept at the bottom of the pillar.
# The resulting system is output to XSF and CFG formats.

atomsk --create fcc 4.08 Au -e 40 40 40 -select out cylinder z box/2 box/2 50.0 -cut above 10.0 z nanopillar.xsf cfg

