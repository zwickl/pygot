#!/usr/bin/env python
import re
import shlex
import os
from argparse import ArgumentParser

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment

from pygot.biopyutils import *
#from oryzautils import *
from pygot.utils import extract_sequences_and_stuff_from_nexus_file


#use argparse module to parse commandline input
parser = ArgumentParser()

parser.add_argument('-m', '--masterfile', type=str, default=None, required=True, 
                    help='file containing list of gffs, toplevels, names, etc. (required)')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='series of alignment files')

parser.add_argument('-na', '--no-alignments', default=False, action='store_true', 
                    help='don\'t actually write alignment files, just the code strings')

#now process the command line
options = parser.parse_args()

#file with lines containing short taxon identifiers, sequence files and gff files for each taxon
allTaxonInfo = get_taxon_genomic_information_dict(options.masterfile, readToplevels=False, usePickle=True)

for nfile in options.filenames:
    print nfile
    taxToSequenceDict = {}
    beginningLinesInNexus = []
    endLinesInNexus = []
    extract_sequences_and_stuff_from_nexus_file(nfile, taxToSequenceDict, beginningLinesInNexus, endLinesInNexus)
    print len(taxToSequenceDict), len(beginningLinesInNexus), len(endLinesInNexus)
                   
    parsedAlignment = extract_all_information_for_seqs_in_alignments(nfile)

    convertedSequences = {}
    intronsStrippedSequences = {}
    exonsStrippedSequences = {}
    parsed = parsedAlignment[0][1]
    for taxon in oryza.taxon_names:
        shortName = oryza.long_name_to_short(taxon)
        if shortName in allTaxonInfo:
            this = [x[1] for x in parsed if x[0] == taxon]
            if len(this) > 0:
                thisParsedSeqInfo = this[0]
                feat = allTaxonInfo[shortName].gff_feature_dict[re.sub('FGT', 'FG', thisParsedSeqInfo.name)]
                #D'OH! problem here with satj, in that start of gene is not start of coding due to UTRs being in exons
                #THIS WOULD NEED TO BE CHANGED IF UTRs WERE INCLUDED IN SATIVA SEQS AGAIN
                #adjust_feature_coords([feat], -int_feature_location(feat)[0])
                if feat.strand == -1:
                    adjust_feature_coords([feat], -find_cds_end_coordinate(feat))
                else:
                    adjust_feature_coords([feat], -find_cds_start_coordinate(feat))

                exons = get_features_by_name(feat, "CDS")
                exonPos = []
                for exon in exons:
                    r = int_feature_location(exon)
                    exonPos.extend(xrange(r[0], r[1]))
                if feat.strand == -1:
                    #this will "flip" the positions, since the first base of the gene on the neg strand
                    #is the one with the largest coord
                    n = max(exonPos)
                    exonPos = [ n - pos for pos in exonPos ]
                exonPos.sort()

                rawPos = 0
                seq = taxToSequenceDict[taxon]

                newSeq = [ char.upper() for char in seq ]
                seqLen = len(newSeq)
               
                #this is to mask off a terminal stop codon, which CDS sequences aligned as AA's won't have
                #there is probably a better algorithm for this
                stops = [ ['T', 'A', 'G'], ['T', 'A', 'A'], ['T', 'G', 'A'] ]
                bases = [ b for b in newSeq if b != '-' ]
                endCodon = bases[-3:]
                codonPos = 2
                if endCodon in stops:
                    for pos in xrange(len(newSeq) - 1, -1, -1):
                        if newSeq[pos] == endCodon[codonPos]:
                            newSeq[pos] = 'N'
                            codonPos -= 1
                            if codonPos == -1:
                                break

                #noIntronSeq = copy.deepcopy(newSeq)
                #noExonSeq = copy.deepcopy(newSeq)
                noIntronSeq = list(newSeq)
                noExonSeq = list(newSeq)
                for alignedPos, char in enumerate(seq):
                    if char.isalpha():
                        if rawPos not in exonPos:
                            newSeq[alignedPos] = char.lower()
                            noIntronSeq[alignedPos] = 'N'
                        else:
                            noExonSeq[alignedPos] = 'N'

                    if char not in ['-', '?']:
                        rawPos += 1

                convertedSequences[taxon] = newSeq
                intronsStrippedSequences[taxon] = noIntronSeq
                exonsStrippedSequences[taxon] = noExonSeq

    #write an alignment that is identical except in the case of intron and exon characters
    #write an alignment that only contains exon characters, with characters annotated as introns changed to N's
    #write an alignment that only contains intron characters, with characters annotated as exons changed to N's
    if not options.no_alignments:
        for (sequences, directory) in [ (convertedSequences, "alteredCases"), (intronsStrippedSequences, "intronCharsStripped"), (exonsStrippedSequences, "exonCharsStripped") ]:
            if not os.path.exists(directory):
                os.makedirs(directory)

            f = directory + '/' + 'aligned.blink.' + parsedAlignment[0][0] + '.nex'
            outnex = open(f, 'w')
            outnex.write('%s' % ''.join(beginningLinesInNexus[:-1]))
            outnex.write('[this can be used to view upper and lower case as different characters]\n[format datatype=standard respectcase missing=? gap=- symbols="a c g t n A C G T N";]\n')
            outnex.write('%s' % ''.join(beginningLinesInNexus[-1]))

            for taxon in oryza.taxon_names:
                if taxon in sequences:
                    outnex.write("\'%s\'\t%s\n" % (taxon, ''.join(sequences[taxon])))

            outnex.write('%s' % ''.join(endLinesInNexus))
        
        #now write versions will all missing columns removed
        for (sequences, directory) in [ (intronsStrippedSequences, "intronCharsStripped.collapsed"), (exonsStrippedSequences, "exonCharsStripped.collapsed") ]:
            if not os.path.exists(directory):
                os.makedirs(directory)

            collapsedSequences = dict((k, []) for k in sequences.keys())
            for site in xrange(seqLen):
                thisCol = []
                for tax, seq in sequences.items():
                    thisCol.append(seq[site])
                if len(''.join(thisCol).translate(None, 'nN?-')):
                    for tax, seq in sequences.items():
                        collapsedSequences[tax].append(seq[site])

            newLen = len(collapsedSequences.values()[0])
            lenStr = 'nchar=%d;' % newLen

            f = directory + '/' + 'aligned.blink.' + parsedAlignment[0][0] + '.nex'
            outnex = open(f, 'w')
            #outnex.write('%s' % ''.join(beginningLinesInNexus[:-1]))
            outnex.write('%s' % re.sub('nchar.*;', lenStr, ''.join(beginningLinesInNexus[:-1])))
            outnex.write('[this can be used to view upper and lower case as different characters]\n[format datatype=standard respectcase missing=? gap=- symbols="a c g t n A C G T N";]\n')
            outnex.write('%s' % ''.join(beginningLinesInNexus[-1]))

            for taxon in oryza.taxon_names:
                if taxon in sequences:
                    #outnex.write("\'%s\'\t%s\n" % (taxon, ''.join(sequences[taxon])))
                    outnex.write("\'%s\'\t%s\n" % (taxon, ''.join(collapsedSequences[taxon])))

            outnex.write('%s' % ''.join(endLinesInNexus))

    countColumns = True
    if countColumns:
        align = MultipleSeqAlignment( [SeqRecord(Seq(''.join(convertedSequences[taxon]), generic_dna), id=taxon) for taxon in oryza.taxon_names if taxon in convertedSequences] )
        empty_align = MultipleSeqAlignment( [SeqRecord(Seq('', generic_dna), id=taxon) for taxon in oryza.taxon_names if taxon in convertedSequences] )

        mixedG = []
        mixedNG = []
        intronG = []
        intronNG = []
        exonG = []
        exonNG = []
        pureG = []

        codeString = []
        informativeCounts = []
        #codeTranslation will be indexed by codeTranslation[bool(gaps)][bool(upper)][bool(lower)]
        #X shouldn't be possible
        codeTranslation = [['X', 'I'], ['E', 'M']], [['G', 'J'], ['F', 'N']]
        '''
        G = full gap

        E, F = full exon (no gap, gap)
        I, J = full intron (no gap, gap)
        M, N = mixed intron/exon (no gap, gap)
        '''

        for colnum in xrange(len(align[0])):
            #this should pull out individual columns
            colStr = align[:, colnum]
            lower = re.search("[a-z]", colStr)
            upper = re.search("[A-Z]", colStr)
            gaps = re.search("-", colStr)
            
            lowerCountNoNs = len(re.findall("[a-mo-z]", colStr))
            upperCountNoNs = len(re.findall("[A-MO-Z]", colStr))
            informativeCount = lowerCountNoNs + upperCountNoNs

            if not (lower or upper or gaps):
                raise RuntimeError('WTF is up with seq column %d?: %s' % (colnum, colStr))
            codeString.append(codeTranslation[bool(gaps)][bool(upper)][bool(lower)])
            informativeCounts.append(informativeCount)
            if gaps:
                if lower and upper:
                    mixedG.append(colnum)
                elif lower:
                    intronG.append(colnum)
                elif upper:
                    exonG.append(colnum)
                else:
                    pureG.append(colnum)
            else:
                if lower and upper:
                    mixedNG.append(colnum)
                elif lower:
                    intronNG.append(colnum)
                elif upper:
                    exonNG.append(colnum)
        
        mixedIntronsAndExons = sorted(mixedG + mixedNG)
        intronsOnly = sorted(intronG + intronNG)
        exonsOnly = sorted(exonG + exonNG)
        noFullIntrons = sorted(mixedIntronsAndExons + exonsOnly)
        noFullExons = sorted(mixedIntronsAndExons + intronsOnly)
        noGaps = sorted(mixedNG + exonNG + intronNG)
        someGaps = sorted(mixedG + exonG + intronG)

        codeStringFilename = 'aligned.blink.%s.dat' % extract_core_filename(nfile)
        directory = 'codeStrings'
        if not os.path.exists(directory):
            os.makedirs(directory)
        with open(directory + '/' + codeStringFilename, 'w') as strFile:
            #strFile.write('%s\n' % '\n'.join(codeString))
            strFile.write('%s\n' % '\n'.join(['%s\t%s\t%d' % (code, count, site + 1) for (code, count, site) in zip(codeString, informativeCounts, xrange(len(codeString)))]))

    if not options.noAlignments:
        outnex.write('[\nmixed %d\nmixedG %d\nmixedNG %d\nintron %d\nintronG %d\nintronNG %d\nexon %d\nexonG %d\nexonNG %d\nnoGaps %d\nsomeGaps %d\n]\n' % 
            len(mixedIntronsAndExons), 
            len(mixedG), 
            len(mixedNG), 
            len(intronsOnly), 
            len(intronG), 
            len(intronNG), 
            len(exonsOnly), 
            len(exonG), 
            len(exonNG), 
            len(noGaps), 
            len(someGaps), 
            )

        newAligns = [ 
            (mixedIntronsAndExons, 'mixedIntronsAndExons'),
            (intronsOnly, 'intronsOnly'),
            (exonsOnly, 'exonsOnly'),
            (noFullIntrons, 'noFullIntrons'), 
            (noFullExons, 'noFullExons'),
            (noGaps, 'noGaps'),
            (someGaps, 'someGaps')
            ]

        outputFilename = 'aligned.blink.%s.nex' % extract_core_filename(nfile)

        for (sites, directory) in newAligns:
            if not os.path.exists(directory):
                os.makedirs(directory)
            thisAln = copy.deepcopy(empty_align)
            for site in sites:
                thisAln += align[:, site:site + 1]
            if len(thisAln[0].seq):
                fname = directory + '/' + outputFilename
                #outfile = open(fname, 'w')
                #AlignIO.write(thisAln, outfile, 'nexus')
                #outfile.close()
                #I think that this will work and automatically close the file, since we need to open and read it below
                AlignIO.write(thisAln, fname, 'nexus')
               
                #this is a little silly, but to add in this one line we need to read and rewrite the whole fie
                outfile = open(fname, 'w')
                for l in open(fname, 'r'):
                    if 'format' in l:
                        outfile.write('[this can be used to view upper and lower case as different characters]\n[format datatype=standard respectcase missing=? gap=- symbols="a c g t n A C G T N";]\n')
                    outfile.write(l)
                outfile.write("[columns:\n%s\n]\n" % ' '.join([str(c + 1) for c in sites]))
                outfile.close()
            
        outnex.close()
