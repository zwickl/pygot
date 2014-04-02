#!/usr/bin/env python
import sys
import re
from argparse import ArgumentParser
from itertools import izip
from pygot.quartets import triplet_support_details, locus_support_details

#use argparse module to parse commandline input
parser = ArgumentParser(description='Compare the triplet resolutions obtained for a set of loci using different treatments')

parser.add_argument('--min-support', type=float, default=0.5,
                    help='minimum support to be considered resolved, 0.5 <= min_val <= 1.0, default 0.5')

parser.add_argument('--full-summary', action='store_true', default=False,
                    help='output more verbose summary that just a 4x4 matrix of triplet flux')

parser.add_argument('infile', nargs='+', type=str)

parser.add_argument('-o', '--output', type=str, default=None, 
                    help='file to output summaries to (default is stdout)')

options = parser.parse_args()

if len(options.infile) == 1:
    #deprecated
    #input in the old format was a pasting together of the following columns of the dendropyScoreTriples.py output for each treatment
    #../../../../../../prank.genes/analyses/boot.garli/runs//aligned.blink.00000.00000.11T.noDupes.6356C.boot.tre    0.9350  0.0150  0.0000  0.9517  0.0317  0.0167  2834944 4098535 38790   3564993
    lines = [ line.strip().split() for line in open(options.infile[0], 'rbU') if not 'file' in line ]
    filenames = [ (line[0], line[11]) for line in lines ]

    #get a tuple of lists 
    tuples = [ ([float(el) for el in line[1:4]], [float(el) for el in line[12:15]] ) for line in lines ]
else:
    first = triplet_support_details(filename=options.infile[0])
    sec = triplet_support_details(filename=options.infile[1])

    #this assumes that the ordering and content of the loci is the same in both files, which it might not be in some cases.  
    #Would need work in that case, as in commonLinePaste.py
    filenames = [ (loc1.treefile, loc2.treefile) for loc1, loc2 in izip(first, sec) ]
    supportTuples = [ (loc1.support_value_list, loc2.support_value_list) for loc1, loc2 in izip(first, sec) ]
   
#this will allow boot support values to be used, with support <= min-support being unresolved
#this will convert from the tuple of potentially floating point support values read from the file to a tuple of ints, with 1 
#indicating the best tree
#add value to all of the elements, and convert all to ints, this just makes anything with sufficient support to truncate to 1
#HOWEVER, having exactly 0.5 for two triples would give a tie, and should really be uninformative, so need to check for that
binaryTuples = []
adder = 1.0 - options.min_support
for tset1, tset2 in supportTuples:
    first = tuple(int(el + adder) for el in tset1)
    sec = tuple(int(el + adder) for el in tset2)
    if sum(first) > 1:
        first = (0, 0, 0)
    if sum(sec) > 1:
        sec = (0, 0, 0)

    binaryTuples.append((first, sec))

#throughout here an index of 0 means unresolved
tupleToIndex = { ('0', '0', '0'):0, ('1', '0', '0'):1, ('0', '1', '0'):2, ('0', '0', '1'):3, (0, 0, 0):0, (1, 0, 0):1, (0, 1, 0):2, (0, 0, 1):3}

counts = [[0, 0, 0, 0] for _ in xrange(2)]
changeMatrix = [[0, 0, 0, 0] for _ in xrange(4)]

summaryFile = open(options.output, 'w') if options.output else sys.stdout

for num, (tup1, tup2) in enumerate(binaryTuples):
    try:
        first = tupleToIndex[tup1]
    except KeyError:
        sys.exit('problem with tup1: %s' % str(tup1))
    try:
        sec = tupleToIndex[tup2]
    except KeyError:
        sys.exit('problem with tup2: %s' % str(tup2))
    #these are the number of times each triple resolution or unresoved is supported by each
    counts[0][first] += 1
    counts[1][sec] += 1

    #this is the number of times that a first triplet resolution changes to a sec triplet resolution in the second treatment
    changeMatrix[first][sec] += 1
    if options.full_summary:
        summaryFile.write('%s\t%s\t%d\t%d\n' % (filenames[num][0], filenames[num][1], first, sec))

tot1, tot2 = sum(counts[0]), sum(counts[1])

if tot1 != tot2:
    sys.exit("something doesn't add up!: %d, %d" % (tot1, tot2))

if options.full_summary:
    summaryFile.write('counts:\n')
    for t in xrange(4):
        summaryFile.write("%2d %4d %4d %4.2f %4.2f\n" % ( t, counts[0][t], counts[1][t], counts[0][t] / float(tot1), counts[1][t] / float(tot2) ))
    summaryFile.write("tot%4d %4d\n\n" % (tot1, tot2))

    summaryFile.write('flux matrix:\n')

summaryFile.write('%s\n' % '\n'.join(['\t'.join([ str(el) for el in row]) for row in changeMatrix]))

