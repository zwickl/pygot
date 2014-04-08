#!/bin/bash

#take a look at exampleAlignment.nex at taxon 'O. sativai AA', coordinate 1394
#for a typical block-shift

#read the nexus alignment, locate block-shifts, output a summary of them
#and write a "masked" alignment to inputAndOutput/exampleAlignment.nex.masked.
#but do NOT mask the outgroup O._brachyantha_FF, which is divergent enough to cause
#likely false positive block-shifts detection
../scripts/findAndMaskBlockshifts.py --ignore-patterns "brachyantha" inputAndOutput/exampleAlignment.nex

