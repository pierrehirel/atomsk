#!/bin/bash

rm -f Loop.cfg Loop.lmp

# Create a unit cell of bcc iron
# Duplicate it to form a 40x40x40 supercell
# Rotate the crystal so that X=[121], Y=[-101], Z=[1-11]
# Introduce a dislocation loop in the middle of the box, with b=1/2[1-11]
# Rotate the crystal back to its original orientation
# Write the final result into Loop.cfg and Loop.lmp
atomsk --create bcc 2.87 Fe  \
       -dup 40 40 40 \
       -orient [100] [010] [001] [121] [-101] [1-11] \
       -disloc loop 0.51*box 0.51*box 0.51*box Y 40 0 0 2.4855 0.33 \
       -alignx \
       Loop.cfg lmp

# The final system "Loop.cfg" may be visualized with Atomeye or OVITO.
