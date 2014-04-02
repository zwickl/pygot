#!/usr/bin/env python
import sys
import os
import re
import collections
import glob
import itertools
from argparse import ArgumentParser
from pygot.dendroutils import CustomTreeList
from pygot import generate_quartet_lists_from_file
from pygot.utils import extract_all_information_for_seqs_in_alignments

try:
    from dendropy import TreeList, Tree, treesplit, treecalc, TaxonSet
except ImportError:
    sys.exit('Problem importing dendropy modules.  This package is a prerequisite to run this script')


def filter_out_strings_by_pattern(toFilter, patterns):
    '''used to ignore some treefiles in initial list'''
    if not patterns:
        return toFilter
    filtering = []
    for name in toFilter:
        for patt in patterns:
            if search(patt, name):
                break
        else:
            filtering.append(name)
    return filtering


def dendropy_score_triples(quartets, treefiles, oryza_names=False, write_trees=False, output_dir='./', filter_file=None, gene_order_file=None):
    #allow override of alignment/treefile lists if passed as second argument
    filterPatterns = []
    if filter_file:
        filterPatterns = [ line.strip() for line in open(filter_file, 'rU') ]

    if os.path.isfile(treefiles):
        tfile = open(treefiles, 'rU')
    else:
        tfile = sys.stdin
    treefileList = []
    for line in tfile:
        unfiltered = glob.glob(line.strip())
        if filterPatterns:
            filtered = filter_out_strings_by_pattern(unfiltered, filterPatterns)
            if len(filtered) > 0:
                treefileList.extend(filtered)
        else:
            treefileList.extend(unfiltered)

    if not treefileList:
        sys.exit("No trees found! Check paths in treefileList")

    #bizarrely, glob does not necessarily list the files in alphnumeric order.  It seems to on 
    #the mac, but it looks totally arbitrary on linux. It doesn't really matter here, but sort
    #for consistency
    treefileList.sort()

    allQuartets = generate_quartet_lists_from_file()

    #dump quartet names if desired
    #for requiredLabels in allQuartets:
    #    print '.'.join([ '%s' % (label.split()[1] if len(label.split()) > 1 and label.split()[0][:2] in ['O.', 'L.'] else re.sub(' ', '_', label)) for label in requiredLabels]) + ' \\'

    #this will get the coordinates of all seqs in all alignments into a dictionary, but ignore if it isn't there
    #end-users won't use this
    alignmentInfo = {}
    if os.path.isfile('alignmentList'):
        for line in open('alignmentList', 'rb'):
            unfiltered = glob.glob(line.strip())
            filtered = filter_out_strings_by_pattern(unfiltered, filterPatterns)
            if filtered:
                alignmentInfo.update( [ (al[2].short_filename, al) for al in extract_all_information_for_seqs_in_alignments(filtered) ] )
    #this gives a CoordinateSet object for each alignment, indexed by the short alignment filename
    alignmentDict = dict( [ [ al[2].short_filename, al[2] ] for al in alignmentInfo.values() ] )

    #again, not for end users
    if gene_order_file:
        #this will get the gene order of all seqs into a dictionary
        splitMap = [ line.split() for line in open(gene_order_file, 'rb') ]
        #dict keys are seq names, values are numerical order 
        geneOrderDict = dict( (line[1], line[0]) for line in splitMap )
        if len(splitMap[0]) > 2:
            geneCoordDict = dict( (line[1], line[2]) for line in splitMap )
    else:
        geneOrderDict = None
        geneCoordDict = None

    if write_trees:
        #matchingTreeDict = dict( (tuple(quart), (CustomTreeList(), CustomTreeList(), CustomTreeList(), CustomTreeList())) for quart in allQuartets )
        matchingTreeDict = dict( (tuple(quart), ([], [], [], [])) for quart in allQuartets )

    #read trees in sets - dendropy can be super memory intensive if loading tons of trees at once
    chunkStart = 0
    chunkSize = 20
    formats = ['nexus', 'newick']
    while chunkStart < len(treefileList):
        treefileSubList = treefileList[ chunkStart:min(chunkStart + chunkSize, len(treefileList)) ]
        sys.stdout.write('reading %d treefiles ...\n' % chunkSize)
        
        #using my derived Treelist is necessary below when testing bipartitions
        #make a treelist for each treefile, which might contain only one tree
        for form in formats:
            try:
                if chunkStart == 0:
                    sys.stdout.write('attempting to read tree files in %s format ...\n' % form)
                allTreesPerTreefile = [ CustomTreeList.get_from_path(f, form) for f in treefileSubList ]
                #now we know what the format is, and don't need to try on the next chunk
                formats = [form]
                break
            except:
                pass

        '''
        #could do pickle writing of trees like this:
        out = open('pickledTrees.pkl', 'wb')
        cPickle.dump(allTreesPerTreefile, out, protocol=1)
        out.close()
        
        #and reading like this, but it doesn't gain a whole lot, maybe 10 or 20%
        inf = open('pickledTrees.pkl', 'rb')
        allTreesPerTreefile = cPickle.load(inf)
        '''

        if not allTreesPerTreefile:
            sys.exit("read no trees!")

        #loop over the desired quartets to test
        sys.stdout.write("processing ...\n")
        for requiredLabels in allQuartets:
            #outfilename = '.'.join([ '%s' % label.split()[1] for label in requiredLabels]) + '.dat'
            #outfilename = '.'.join([ '%s' % (label.split()[1] if len(label.split()) > 1 and label.split()[0][:3] in ['O. ', 'L. '] else label) for label in requiredLabels]) + '.dat'
            #outfilename = ''
            if output_dir:
                #outfilename = os.path.normpath(output_dir)
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                elif not os.path.isdir(output_dir):
                    raise IOError

            outfilename = os.path.join(output_dir, '.'.join([ '%s' % (label.split()[1] if len(label.split()) > 1 and label.split()[0][:2] in ['O.', 'L.'] else re.sub(' ', '_', label)) for label in requiredLabels]) + '.dat')
            #outfilename += '.'.join([ '%s' % (label.split()[1] if len(label.split()) > 1 and label.split()[0][:3] in ['O. ', 'L. '] else re.sub(' ', '_', label)) for label in requiredLabels]) + '.dat'
            fout = open(outfilename, 'w') if chunkStart == 0 else open(outfilename, 'a')
          
            '''
            possibleGroups = combine_components_and_uniqueify(requiredLabels)
            possibleSplits = []
            for s in possibleGroups:
                thisSet = []
                for item in s:
                    if len(item) == 2 and (isinstance(item, list) or isinstance(item, tuple)):
                        thisSet.append(list(item))
                if thisSet:
                    possibleSplits.append(thisSet)
            print possibleSplits
            print len(possibleSplits)
            '''

            numFound = 0
            firstOfChunk = True
            for (trees, treefile) in itertools.izip(allTreesPerTreefile, treefileSubList):
                #verify that the correct taxa are present
                if trees.taxon_set.has_taxa(labels=requiredLabels):
                    numFound += 1
                    #masked = True
                    #if masked:
                    if len(requiredLabels) == 4:
                        #this is the newer scheme, passing the mask instead of actually pruning the tree
                        #faster for quartets
                        theseTrees = trees
                        
                        if len(requiredLabels) == 4:
                            possibleSplits = [ [requiredLabels[0], requiredLabels[1]], [requiredLabels[0], requiredLabels[2]], [requiredLabels[0], requiredLabels[3]] ]
                        else:
                            #this would generate all possible splits, if one wanted to test for certain splits within a larger tree
                            possibleSplits = [ list(sp) for sp in itertools.combinations(requiredLabels, 2) ]
                        
                        backboneSplitMask = theseTrees.taxon_set.get_taxa_bitmask(labels=requiredLabels)

                        possibleSplitMasks = []
                        #these sets will have only one bipartition each currently, but the frequency function can handle more
                        for splitSet in possibleSplits:
                            possibleSplitMasks.append([theseTrees.taxon_set.get_taxa_bitmask(labels=splitSet)])
                        
                        freqs = [ (theseTrees.masked_frequency_of_splitlist(split_bitmask=splitMask, mask=backboneSplitMask), num) for (num, splitMask) in enumerate(possibleSplitMasks) ]
                        maxFreq = max(freqs)
                        if write_trees:
                            #results = [ (theseTrees.masked_frequency_of_splitlist(split_bitmask=splitMask, returnMatches=True, mask=backboneSplitMask), num) for (num, splitMask) in enumerate(possibleSplitMasks) ]
                            #the indexing is screwed up here, as results is a list of ((freq, trees), num) tuples
                            #freqs = [ (res[0][0], res[1]) for res in results ]
                            cons = theseTrees.consensus(min_freq=0.0)
                            #cons.label = treefile
                            if maxFreq[0] > 0.5:
                                #matchingTreeDict[tuple(requiredLabels)][maxFreq[1]].append(cons)
                                matchingTreeDict[tuple(requiredLabels)][maxFreq[1]].append((treefile, cons.as_newick_string()))
                            else:
                                #matchingTreeDict[tuple(requiredLabels)][3].append(cons)
                                matchingTreeDict[tuple(requiredLabels)][3].append((treefile, cons.as_newick_string()))
                        
                        if chunkStart == 0 and firstOfChunk:
                            fout.write('file\t%s\n' % ('\t'.join(repr(sp) for sp in possibleSplits)))
                        firstOfChunk = False
                        
                    else:
                        #this is the only way that I've worked out to do this with > 4 taxa
                        #it also generates and tests partially resolved trees

                        #this is complicated, and somewhat annoying.  MUST create new taxon_set object for new TreeSet, because otherwise 
                        #Taxon objects will be copied, and may have already had things like their split bitmasks set based on the full tree
                        #Couldn't figure out how to have them reindexed when taxa are removed.
                        #This creates a TreeList with a new TaxonSet and copies the trees over.  Then it infers the proper taxon_set from 
                        #the first tree and assigns it to the TreeList taxon_set.  Then that taxon_set is asigned back to all of the trees
                        taxonSet = TaxonSet()
                        theseTrees2 = CustomTreeList(trees, taxon_set=taxonSet)
                        #you can pass update_splits=True to t.retain_taxa_with_labels, but it borks things, I think because the taxon_set
                        #has to be updated first
                        for t in theseTrees2:
                            t.retain_taxa_with_labels(labels=requiredLabels)
                        theseTrees2.taxon_set = theseTrees2[0].infer_taxa()
                        theseTrees2.reindex_subcomponent_taxa()
                        for t in theseTrees2:
                            treesplit.encode_splits(t)
                        
                        possibleTrees = CustomTreeList(taxon_set=theseTrees2.taxon_set)
                        possibleTrees.generate_all_trees_for_taxon_list(requiredLabels, min_bipartitions=0, max_bipartitions=5, criterion=None)

                        freqs = [ (theseTrees2.frequency_of_identical_trees(targetTree), num) for num, targetTree in enumerate(possibleTrees) ]
                        maxFreq = max(freqs)
                    
                        if chunkStart == 0 and firstOfChunk:
                            fout.write('file\t%s\n' % ('\t'.join(t.as_newick_string() for t in possibleTrees)))
                        firstOfChunk = False
                    

                    fout.write("%s\t" % treefile)

                    if len(trees) > 1:
                        fout.write('%s' % '\t'.join(['%.4f' % f[0] for f in freqs]))
                        #the proportions may not sum to 1 (if bootstraps), so divide any unassigned proportion among the triples and
                        #output these new rescaled columns.  This will mainly be good for barycentric plots, where the props must sum to 1.0
                        s = sum([f[0] for f in freqs])
                        ambigProportion = (1.0 - s) / len(freqs)
                        fout.write('\t%s\t' % '\t'.join(['%.4f' % (max(f[0] + ambigProportion, 0.0)) for f in freqs]))
                    else:
                        #if there is only one tree (i.e., not bootstrap or bayesian samples) output as ints, add 
                        #small value to ensure that a floating point of 1.0ish truncates to 1 as int     
                        fout.write('%s\t' % '\t'.join(['%d' % int(f[0] + 0.1) for f in freqs]))
                        #output this twice, to match the behavior when doing bootstrap/bayes, since there there is another set of columns 
                        #with the support rescaled to sum to 1.0
                        fout.write('%s\t' % '\t'.join(['%d' % int(f[0] + 0.1) for f in freqs]))
                    
                    if oryza_names:
                        from dzutils import extract_core_filename
                        shortTreeFile = extract_core_filename(treefile)
                    else:
                        shortTreeFile = treefile

                    if alignmentDict:
                        fout.write('%s' % alignmentDict[shortTreeFile].row_output(taxon_labels=requiredLabels))
                    else:
                        fout.write('0\t0\t0\t0\t')

                    #output gene order info too
                    if not geneOrderDict:
                        fout.write('0\t0\t0\t0\t')
                    elif geneOrderDict and alignmentInfo:
                        #DEBUG
                        try:
                            #this is ridiculous
                            #
                            info = alignmentInfo[shortTreeFile]
                            idict = dict([i for i in info[1]])
                            fout.write('%s\t' % '\t'.join([geneOrderDict[re.sub('[.][0-9]+$', '', idict[lab].name)] for lab in requiredLabels]))
                        except StandardError:
                            fout.write('0\t0\t0\t0\t')
                    
                    #don't think that I had any reason for having this as 1% cutoff.  The resolution output isn't really used for anything anyway
                    #if maxFreq[0] > 0.01:
                    if maxFreq[0] > 0.0:
                        #fout.write('(%s,%s)\n' % (possibleSplits[ maxFreq[1] ][0], possibleSplits[ maxFreq[1] ][1]))
                        #fout.write('%s\n' % (repr(possibleSplits[ maxFreq[1] ])))
                        if len(requiredLabels) == 4:
                            fout.write('(%s,%s)\n' % (possibleSplits[maxFreq[1]][0], possibleSplits[maxFreq[1]][1]))
                        else:
                            fout.write('%.4f:%s\n' % (maxFreq[0], possibleTrees[maxFreq[1]].as_newick_string()))
                    else:
                        fout.write('(noRes)\n')
        
        fout.close()
        chunkStart += chunkSize

    for requiredLabels in allQuartets:
        outfilename = os.path.join(output_dir, '.'.join([ '%s' % (label.split()[1] if len(label.split()) > 1 and label.split()[0][:2] in ['O.', 'L.'] else re.sub(' ', '_', label)) for label in requiredLabels]) + '.dat')
     
        #outfilename = ''
        '''
        for label in requiredLabels:
            splabel = label.split()
            outfilename += '%s.' % (splabel[1] if splabel[0][:3] in ['O. ', 'L. '] else splabel[0]) 
        '''
        #if output_dir:
        #    outfilename = os.path.normpath(output_dir)
        #outfilename += '.'.join([ '%s' % (label.split()[1] if len(label.split()) > 1 and label.split()[0][:3] in ['O. ', 'L. '] else re.sub(' ', '_', label)) for label in requiredLabels]) + '.dat'
        if os.stat(outfilename).st_size == 0:
            #this is kind of cheesy, but it's easier to just create the potential output file above and
            #then remove it if nothing was written to it
            sys.stdout('no trees found containing quartet %s found\n' % requiredLabels)
            os.remove(outfilename)
        else:
            if write_trees:
                treeLists = matchingTreeDict[tuple(requiredLabels)]
                for num, quartTrees in enumerate(treeLists[:3]):
                    tempTreeList = CustomTreeList.get_from_string(';'.join(tree[1] for tree in treeLists[num]), 'newick')
                    for tnum, tree in enumerate(tempTreeList):
                        tree.label = treeLists[num][tnum][0]
                    #quartTrees.write(open('.'.join([outfilename, str(num), 'tre']), 'w'), schema='nexus')
                    #this is hardcoded for oryza names at the moment
                    stream = open(re.sub('dat', requiredLabels[0].split()[1] + '-' + requiredLabels[num+1].split()[1] + '.tre', outfilename), 'w')
                    tempTreeList.write(stream, schema='nexus')
                    #stream.write('\n'.join(treeLists[num]) + '\n')
                    #quartTrees.write(stream, schema='nexus')
                   
                tempTreeList = CustomTreeList.get_from_string(';'.join(tree[1] for tree in treeLists[3]), 'newick')
                for tnum, tree in enumerate(tempTreeList):
                    tree.label = treeLists[3][tnum][0]
                stream = open(re.sub('dat', 'U.tre', outfilename), 'w')
                tempTreeList.write(stream, schema='nexus')
                #stream.write('\n'.join(treeLists[3]) + '\n')
                #treeLists[3].write(stream, schema='nexus')


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    #use argparse module to parse commandline input
    parser = ArgumentParser(description='Read a quartet file and use dendropy to read and summarize a set of trees')

    parser.add_argument('-o', '--output-dir', type=str, default='./',
                        help='Directory in which to write files summarizing each quartet. Default is current dir.')

    parser.add_argument('--filter-file', type=str, default=None,
                        help='Optional file containing regex or literal patterns that cause matching treefiles to be ignored')

    parser.add_argument('--gene-order-file', type=str, default=None,
                        help='Optional file containing gene order information')

    default_quartet_file = 'quartetList'
    parser.add_argument('--quartet-file', type=str, default=default_quartet_file,
                        help='File containing list of quartets to summarize (see docs for format). Default is %s.' % default_quartet_file)

    default_treelist_file = 'treefileList'
    parser.add_argument('--treelist-file', type=str, default=default_treelist_file,
                        help='File containing a list of relative paths to tree files (or glob patterns of tree files) to read and summarize. \
                        A list of treefiles may also be piped in on the command line.  If nothing is piped, default filename is %s.' % default_treelist_file)

    parser.add_argument('--oryza-names', action='store_true', default=False,
                        help='Do some name manipulations specific to Oryza data')

    parser.add_argument('--write-trees', action='store_true', default=False,
                        help='output treefiles for each treefile-quartet resolution combination')

    options = parser.parse_args()

    dendropy_score_triples(quartets=options.quartet_file, 
                            treefiles=options.treelist_file,
                            oryza_names=options.oryza_names,
                            write_trees=options.write_trees,
                            output_dir=options.output_dir,
                            filter_file=options.filter_file,
                            gene_order_file=options.gene_order_file)


