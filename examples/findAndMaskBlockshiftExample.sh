#!/bin/bash

#take a look at exampleAlignment.nex at taxon 'O. sativai AA', coordinate 1394
#for a typical block-shift

#read the nexus alignment, locate block-shifts, output a summary of them
#and write a "masked" alignment to inputAndOutput/exampleAlignment.nex.masked
../scripts/findAndMaskBlockshifts.py inputAndOutput/exampleAlignment.nex
