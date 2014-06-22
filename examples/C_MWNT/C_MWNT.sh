#!/bin/bash

## This script builds a carbon multi-wall nanotube (MWNT)
## using atomsk. When creating a nanotube, its axis
## is along Z and the nanotube is centered at (0,0).

rm -f *.atsk mwnt.*

## Create (8,0) CNT
atomsk --create nanotube 2.6 8 0 C cnt1.atsk

## Create (16,0) CNT
atomsk --create nanotube 2.6 16 0 C cnt2.atsk

## Merge the two nanotubes, and reapeat it 4 times along its axis (i.e. Z axis).
## Note that the biggest NT is read first so that its cell encloses the two nanotubes
atomsk --merge 2 cnt2.atsk cnt1.atsk mwnt.xsf xyz cfg -e 1 1 4

## Remove temporary files
rm -f *.atsk
