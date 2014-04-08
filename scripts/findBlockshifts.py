#!/usr/bin/env python
import sys
import re
from collections import Counter
from argparse import ArgumentParser

from pygot.utils import extract_sequences_and_stuff_from_nexus_file
from pygot.utils import proportion_type


#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

parser = ArgumentParser(description='look for regions taxon by taxon that are poorly aligned')

parser.add_argument('--output-full-window-range', action='store_true', default=False,
                    help='output coordinates of the actual windows that trip the mismatch criterion, rather than starting at first mismatch and ending at last within the window (default false)')

parser.add_argument('-d', '--debug', action='store_true', default=False,
                    help='output a bunch of extra crap')

parser.add_argument('-mtc', '--min-tax-in-column', type=int, default=4,
                    help='min number of taxa present in column to calculate measure (default 4)')

parser.add_argument('-l', '--region-length', type=int, default=10,
                    help='length of poorly aligned run to screen for (default 10)')

parser.add_argument('-c', '--consensus-proportion', type=proportion_type, default=0.5,
                    help='proportion of non-ambiguous bases in a column that must have the same state for the column to be considered (must be GREATER THAN this value, default 0.5)')

parser.add_argument('-m', '--mismatch-proportion', type=proportion_type, default=0.5,
                    help='min proportion of region that does not match consensus to be considered poorly aligned (default 0.5)')

parser.add_argument('-o', '--outputfile', type=str, default=None, 
                    help='file to write output to (default stdout)')

#variable number of arguments
parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of NEXUS filenames to search (none for stdin)')

#now process the command line
options = parser.parse_args()


ofile = open(options.outputfile, 'w') if options.outputfile else sys.stdout


sys.stderr.write('Algorithm details:\n')
sys.stderr.write('\tRegion length of window tested: %d\n' % options.region_length)
sys.stderr.write('\tMinimum taxa to consider column: %d\n' % options.min_tax_in_column)
sys.stderr.write('\tRequired proportion of column conserved\n\t\t(proportion must be > than this value): %d\n' % options.consensus_proportion)
sys.stderr.write('\tRequired proportion of sequence within\n\t\twindow not matching consensus\n\t\t(proportion must be >= than this value): %.2f\n' % options.mismatch_proportion)
if not options.output_full_window_range:
    sys.stderr.write('\tIdentified regions reported from first mismatch base to last\n')
else:
    sys.stderr.write('\tIdentified regions reported by start to end of window coordinates\n')


for nfile in options.filenames:
    ofile.write('%s\n' % nfile)
    
    taxToSequenceDict = {}
    beginningLinesInNexus, endLinesInNexus = [], []
    extract_sequences_and_stuff_from_nexus_file(nfile, taxToSequenceDict, beginningLinesInNexus, endLinesInNexus)
    if not taxToSequenceDict:
        sys.exit('no sequences in alignment file %s?' % nfile)
    taxAndSequenceList = sorted([ [k, v.lower()] for k, v in taxToSequenceDict.items() ])
    
    seqLen = len(taxAndSequenceList[0][1])
    cols = []
    consensuses = ['-'] * seqLen
    skipStates = ['-', 'n', 'N', '?']
    #find the consensus for each column
    for index in xrange(seqLen):
        #grab non-ambiguous bases in the column
        noGaps = [t[1][index] for t in taxAndSequenceList if t[1][index] not in skipStates]
        numNonAmbig = len(noGaps)
        if options.debug:
            sys.stderr.write('%d\t%s\t%d\n' % (index, noGaps, numNonAmbig))
        
        #if there are enough non-ambig to pay attention to the column
        #otherwise the default consensus of all columns was already set to - above
        if numNonAmbig >= options.min_tax_in_column:
            #this gives a list of tuples of counts, with each being like ('a', 2), from most to least common
            c = Counter(noGaps)
            #this gives a list of length 1 with the most common tuple.  If there is a tie,
            #it returns an arbitrary one, but required proportion of > 0.5 should be required anyway
            state, count = c.most_common(1)[0]
            if options.debug:
                sys.stderr.write('%d\t%s\t%d\n' % (index, state, count))
            if count > (options.consensus_proportion * numNonAmbig):
                consensuses[index] = state
    
    reqMismatchNum = int(options.region_length * options.mismatch_proportion)
    for tax, seq in taxAndSequenceList:
        mismatchList = [ 0 if seq[index] == consensuses[index] or consensuses[index] == '-' or seq[index] in skipStates else 1 for index in xrange(seqLen) ]
        #standardize the taxon names
        #otax = re.sub('[ _]', '.', tax)
        otax = re.sub('[ ]', '_', tax)
        ofile.write('%s' % otax)
        foundAt = -1
        for startIndex in xrange(seqLen - options.region_length):
            window = mismatchList[startIndex:startIndex+options.region_length]
            mismatchCount = sum(window)
            if options.debug:
                sys.stderr.write('%s\t%d\t%d\t%s\n' % (tax, startIndex, mismatchCount, window))
            if foundAt < 0 and mismatchCount >= reqMismatchNum:
                if not options.output_full_window_range:
                    startOffset = window.index(1)
                    foundAt = startIndex + startOffset
                else:
                    foundAt = startIndex
                ofile.write(' %d' % (foundAt + 1))
            if foundAt >= 0 and (mismatchCount < reqMismatchNum or startIndex == (seqLen - options.region_length - 1)):
                if not options.output_full_window_range:
                    end = startIndex + options.region_length - list(reversed(window)).index(1) - 1
                else:
                    end = startIndex + options.region_length

                ofile.write('-%d' % (end + 1))
                foundAt = -1
        ofile.write('\n')

