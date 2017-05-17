#!/bin/bash

#indicate the location of the script files, or just set blank if they are in your path already
#SCRIPTDIR=""
SCRIPTDIR=../scripts/

#pipe a list of bootstrap treefiles to dendropyScoreTriples.py script, read a quartet file and summarize 
#the trees, outputting tabular summary files into the treeset1 and treeset2 directories
ls inputAndOutput/treeset1/*.blink.00*.boot.tre | $SCRIPTDIR/dendropyScoreTriples.py -o inputAndOutput/treeset1 --quartet-file quartetList
if [ "$?" -ne "0" ];then echo "Problem running examples!"; exit; fi

ls inputAndOutput/treeset2/*.blink.00*.boot.tre | $SCRIPTDIR/dendropyScoreTriples.py -o inputAndOutput/treeset2 --quartet-file quartetList
if [ "$?" -ne "0" ];then echo "Problem running examples!"; exit; fi

mkdir -p figures

#pass the corresponding quartet output files from each of the two output directories to the plotting
#script, which will create pdf files for each quartet
for QUART in nivara.punctata.barthii.glaberrimaM.dat  sativaj.punctata.barthii.glaberrimaM.dat
do
    #plot cumulative figures for each treatment
    #$SCRIPTDIR/cumulativeAndPieFigure.py -i inputAndOutput/treeset1/$QUART inputAndOutput/treeset2/$QUART -o figures/$QUART.cumulative.pdf --columns 2 \
    $SCRIPTDIR/cumulativeAndPieFigure.py -i inputAndOutput/treeset1/$QUART -o figures/$QUART.cumulative.pdf --columns 2 \
        --titles "Treatment 1" "Treatment 2"
    if [ "$?" -ne "0" ];then echo "Problem running examples!"; exit; fi

    #plot flux diagrams of comparison of pairs of treatments
    ../scripts/calculateFlux.py inputAndOutput/treeset1/$QUART inputAndOutput/treeset2/$QUART | ../scripts/plotFlux.py -o figures/$QUART.flux.pdf
    if [ "$?" -ne "0" ];then echo "Problem running examples!"; exit; fi
done

