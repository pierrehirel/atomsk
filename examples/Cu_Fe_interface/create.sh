#!/bin/bash

# This example shows how to build an interface
# between fcc copper and bcc iron using atomsk.

rm -f interface.* *.atsk

# First, create the bottom crystal of fcc Cu
atomsk --create fcc 3.61 Cu orient [11-2] [111] [1-10] -expand 16 13 10 bottom_cu.atsk

# Second, create the top crystal of bcc Fe
atomsk --create bcc 2.87 Fe orient [1-12] [110] [-111] -expand 10 20 10 top_fe.atsk

# Finally, merge the two systems and output to XYZ, XSF, CFG and LAMMPS data format
atomsk --merge z 2 bottom_cu.atsk top_fe.atsk interface.xyz xsf cfg lmp

# Remove temporary files
rm -f *.atsk

