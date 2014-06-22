#!/bin/bash

### This script uses atomsk (and some bash scripting)
### to create a screw 1/2 <111> dislocation in iron
### using anisotropic elasticity. In addition a cylinder
### is defined around the dislocation core; atoms inside
### this cylinder are mobile, those outside are fixed.

rm -f Fe_dislo*

### Define a few variables:
### Lattice constant (A)
alat=2.87
### Dimensions of unit cell along X and Y
### (this should match the crystal orientation below)
uX=$(echo "$alat*sqrt(6.0)" | bc -l)
uY=$(echo "$alat*sqrt(2.0)" | bc -l)
### Length of Burgers vector 1/2 [1-11]
b=$(echo "$alat*0.5*sqrt(3.0)" | bc -l)
### Number of unit cells along X, Y, Z
### This can be changed to make the system bigger or smaller
eX=20
eY=30
eZ=1
### Position of the dislocation (also the center of cylinder)
posX=$(echo "$uX*$eX*0.5"|bc)
posY=$(echo "$uY*(0.5*$eY+0.666666666)"|bc)
### Radius (in A) of cylinder around dislocation core
radius=50.0

### Run atomsk!
### Summary of the command:
### Create the oriented unit cell (mode --create)
### Expand it to make a supercell (option -e)
### Read material properties from "Fe_prop.txt" (option -prop)
###     (elastic tensor + system orientation)
### Build 1/2 <111> screw dislocation, line along Z (option -disloc)
###     (elastic tensor was defined before, so
###      anisotropic elasticity will be used)
### Fix atoms at the boundaries (options select and fix)
### Output to XSF, CFG and BOP formats

atomsk --create bcc $alat Fe orient [121] [-101] [1-11] -e $eX $eY $eZ -prop Fe_prop.txt -disloc $posX $posY screw z y $b 0.0 -select out cylinder z $posX $posY $radius -fix all above -100.0 x Fe_dislo.xsf cfg bop

### Note: dislocation stresses and fixed atoms can be visualized
###       with Atomeye: see auxiliary properties in "Fe_dislo.cfg"

