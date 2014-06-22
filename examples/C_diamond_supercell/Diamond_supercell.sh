#!/bin/bash

# This script shows how to use atomsk to build
# a supercell of carbon diamond

# First we create the unit cell, then we expand it
# and we output to XYZ and CFG for visualization
atomsk --create diamond 3.567 C -e 4 4 1 C_diamond_supercell xsf cfg
