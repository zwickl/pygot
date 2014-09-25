#!/usr/bin/env python
import sys
import re
from os.path import expandvars
import copy
import itertools

import Bio
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
import BCBio
from pygot.utils import read_from_file_or_pickle

#my extensions and functions for working with biopython objects


class TaxonGenomicInformation:
    '''Stores correspondence between a set of sequences, a gff file referring to those 
    sequences and a corresponding toplevel assembly.  Some of these may be omitted.
    Important members:
    name            - taxon name passed to __init__
    self.short_name - first 7 characters of taxon name passed to __init__
    toplevel_record_dict   - dictionary of toplevel SeqRecord.id to toplevel SeqRecords
                        empty if toplevel_filename not passed to __init__
    toplevel_record_list   - list of toplevel SeqRecords in toplevel_record_dict
                        empty if toplevel_filename not passed to __init__
    
    seq_dict        - dictionary of individual SeqRecords (presumably genes or mRNA)
                        existing dict can be passed to __init__, or generated/added to
                        with a passed seq_filename
                        empty if seq_dict or seq_filename not passed to __init__

    if a gff_filename is passed to __init__, several things happen:
        -gff_seqrecord_list is a list filled with SeqRecords for each feature referenced in the gff
        -gff_feature_dict is a dict of either feature Aliases (preferably) or id to SeqFeatures
            that are held by the SeqRecords in gff_seqrecord_list
        -toplevel_record_list is updated to have gff features added to the SeqRecords
        -seq_dict is updated with features from gff
    '''
    
    def __init__(self, name, seq_dict=None, seq_filename=None, toplevel_filename=None, gff_filename=None, usePickle=False):
        '''
        print '#############'
        print name
        print seq_dict
        print seq_filename
        print toplevel_filename
        print gff_filename
        print '#############'
        '''
        self.name = name
        self.short_name = name[0:7]
        
        if toplevel_filename is not None:
            toplevel_filename = expandvars(toplevel_filename)
            if not usePickle:
                #it can be handy to have a dict of the toplevel seq(s) recs which may just be a single chrom
                self.toplevel_record_dict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open(toplevel_filename), "fasta"))
                #pull the toplevel reqs back out as a list of seq recs
                self.toplevel_record_list = [it[1] for it in self.toplevel_record_dict.iteritems()]
            else:
                #this is slightly faster reading the pickle than re-parsing
                self.toplevel_record_list = read_from_file_or_pickle(toplevel_filename, toplevel_filename + '.list.pickle', Bio.SeqIO.parse, "fasta")
                if not isinstance(self.toplevel_record_list, list):
                    #it may be an iterator or generator rather than list
                    self.toplevel_record_list = [ it for it in self.toplevel_record_list ]
                self.toplevel_record_dict = Bio.SeqIO.to_dict(self.toplevel_record_list)

            '''
            print '1#####################'
            for top in self.toplevel_record_list:
                if top.features:
                    print len(self.toplevel_record_list)
                    print type(self.toplevel_record_list)
                    print top
                    print len(top.features)
                    print top.features
                    break
            print '/1#####################'
            '''
        else:
            self.toplevel_record_dict = dict()
            self.toplevel_record_list = []
        
        #if filenames are passed in, do the necessary reading and add them to any dicts/lists that might have been passed
        if seq_dict is None:
            seq_dict = dict()
        if seq_filename is not None:
            seq_filename = expandvars(seq_filename)
            if not usePickle:
                self.seq_dict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open(seq_filename), "fasta"))
            else:
                #haven't actually tested pickle here
                self.seq_dict = Bio.SeqIO.to_dict(read_from_file_or_pickle(seq_filename, Bio.SeqIO.parse, "fasta"))
        else:
            self.seq_dict = dict()
        self.seq_dict.update(seq_dict)

        if gff_filename is not None:
            #pickling the gffs doesn't seem to help - pickle file is much larger than actual gff and slower to read
            gff_filename = expandvars(gff_filename)
            #this will asign features in the gff to the toplevel seqs - the toplevel_record_dict will be unaltered
            #print len(self.toplevel_record_dict.items())
            #I think that we can avoid reading this if we don't have toplevels, but it might bork some script
            if self.toplevel_record_dict:
                self.toplevel_record_list = [feat for feat in BCBio.GFF.parse(gff_filename, base_dict=self.toplevel_record_dict)]
                #now assign back to the dict
                self.toplevel_record_dict = dict([ (top.id, top) for top in self.toplevel_record_list ])
            '''
            print '2#####################'
            for top in self.toplevel_record_list:
                if top.features:
                    print len(self.toplevel_record_list)
                    print type(self.toplevel_record_list)
                    print top
                    print len(top.features)
                    print top.features
                    break
            print '/2#####################'
            '''
            #print type(self.toplevel)
            #print len(self.toplevel)
            #print "SELF.TOPLEVEL"
            #print self.toplevel[0]
            #print "/SELF.TOPLEVEL"
            self.gff_seqrecord_list = [rec for rec in BCBio.GFF.parse(gff_filename)]
            '''
            print '3#####################'
            print len(self.gff_seqrecord_list)
            for rec in self.gff_seqrecord_list:
                if rec.features:
                    print rec
                    print len(rec.features)
                    print rec.features
                    break
            if len(self.gff_seqrecord_list) > 1:
                print self.gff_seqrecord_list[1]
            print '/3#####################'
            '''
            self.seq_dict = Bio.SeqIO.to_dict([feat for feat in BCBio.GFF.parse(gff_filename, base_dict=self.seq_dict)])

        else:
            self.gff_seqrecord_list = []
            self.seq_dict = dict()

        if toplevel_filename is not None and gff_filename is not None:
            self.feature_to_toplevel_record_dict = dict()
            for top in self.toplevel_record_list:
                for feat in top.features:
                    if 'Alias' in feat.qualifiers:
                        self.feature_to_toplevel_record_dict[feat.qualifiers['Alias'][0]] = top
                    else:
                        self.feature_to_toplevel_record_dict[feat.id] = top

        self.gff_feature_dict = {}
        if self.gff_seqrecord_list:
            for rec in self.gff_seqrecord_list:
                if rec.features[0].type.lower() in ['chromosome', 'contig', 'scaffold']:
                    #trying to verify when this is happening
                    #exit('DEBUG: FIRST RECORD IS %s' % rec.features[0].type)
                    toIter = rec.features[1:]
                else:
                    toIter = rec.features
                for g in toIter:
                    if 'Alias' in g.qualifiers:
                        self.gff_feature_dict[g.qualifiers['Alias'][0]] = g
                    else:
                        self.gff_feature_dict[g.id] = g

    def output(self):
        print '%s\t%d features\t%d dictFeatures\t%d sequences\t%d toplevel records' % (self.name, 
            len(self.gff_seqrecord_list), len(self.gff_feature_dict.keys()), len(self.seq_dict.keys()), len(self.toplevel_record_list))

    def toplevel_record_for_feature(self, feat):
        if isinstance(feat, str):
            if feat in self.feature_to_toplevel_record_dict:
                return self.feature_to_toplevel_record_dict[feat]
            else:
                raise ValueError('toplevel for feature named %s not found!' % feat)
        else:
            name = feat.qualifiers['Alias'][0] if 'Alias' in feat.qualifiers else feat.id
            if name in self.feature_to_toplevel_record_dict:
                return self.feature_to_toplevel_record_dict[name]
            else:
                raise ValueError('toplevel for feature named %s not found!' % name)


def instantiate_taxon_genomic_information(taxon, gff_filename, usePickle, toplevel_filename=None ):
    return TaxonGenomicInformation(taxon, gff_filename=gff_filename, toplevel_filename=toplevel_filename, usePickle=False)


def get_taxon_genomic_information_dict(source, report=True, readToplevels=True, usePickle=False, useSMP=False):
    #file with lines containing short taxon identifiers, sequence files and gff files for 
    #each taxon
    #like this (on one line)
    #ObartAA	blah	/Users/zwickl/Desktop/GarliDEV/experiments/productionOryza2/gramene34_split/gffs/bartAA.fullWithFixes.gff \
        #/Users/zwickl/Desktop/GarliDEV/experiments/productionOryza2/gramene34_split/toplevels/Oryza_barthii-toplevel-20110818.fa
    if isinstance(source, list):
        masterFilenames = [ s.strip.split() for s in source ] if len(source[0]) == 1 else source
    else:
        masterFilenames = [ line.strip().split() for line in open(source, 'rb') if len(line.strip()) > 0 ]

    allTaxonInfo = {}
    if useSMP:
        import pp
        job_server = pp.Server(2)
        sys.stderr.write("Starting parallel python with %d workers\n" % job_server.get_ncpus())
        if not readToplevels:
            taxon = masterFilenames[0]
            funcs = [ job_server.submit(instantiate_taxon_genomic_information, (taxon[0], taxon[2], usePickle), (TaxonGenomicInformation, expandvars), ("dzbiopython", "BCBio.GFF", "Bio.SeqIO")) for taxon in masterFilenames[:] ]
        else:
            funcs = [ job_server.submit(instantiate_taxon_genomic_information(taxon[0], gff_filename=taxon[2], toplevel_filename=taxon[3], usePickle=usePickle)) for taxon in masterFilenames ]
        exit('CHECK parallel usage in get_taxon_genomic_information_dict')
        #can this guarantee that the funcs from the list comp is in the same order as masterFilenames?
        job_server.wait()
        for t, i in itertools.izip([ t[0] for t in masterFilenames ], funcs):
            allTaxonInfo[t] = i()
        sys.stderr.write("%d done\n" % len(allTaxonInfo))
    else:
        for taxon in masterFilenames:
            #ended up not using sequence files, just getting everything from toplevels
            if not readToplevels:
                allTaxonInfo[ taxon[0] ] = TaxonGenomicInformation(taxon[0], gff_filename=taxon[2], usePickle=usePickle)
            else:
                allTaxonInfo[ taxon[0] ] = TaxonGenomicInformation(taxon[0], gff_filename=taxon[2], toplevel_filename=taxon[3], usePickle=usePickle)
            if report:
                allTaxonInfo[ taxon[0] ].output()
    return allTaxonInfo


def parse_feature_name(feature, errorIsFatal=True):
    '''Return what I'm treating as the name of the SeqFeature, which is 
    stored as one of the arbitrarily named qualifiers, named differently
    for IRGSP and OGE gff's
    
    5/21/13 - Switched this to prefer Name over Alias, which fixes an issue
    with MAKER annotations and hopefully won't bork anything else.
    6/4/13 - Also added ID as a potential feature name
    '''

    #print feature.qualifiers
    for qual in ['Name', 'Alias', 'ID']:
        if qual in feature.qualifiers:
            return feature.qualifiers[qual][0]
       
    #print feature
    if errorIsFatal:
        raise Exception('unable to parse a feature name!:')
    else:
        sys.stderr.write('unable to parse a name!')


def int_feature_location(feat):
    '''Deals with the annoying fact that hoops must be jumped through to make FeatureLocations into numbers
    NOTE: returns location coords in standard biopython format, with start counting from zero, and
    end being one PAST the last base'''
    return (feat.location.start.position, feat.location.end.position)


def flattened_subfeature_iterator(feature, reverse=False):
    '''generator to essentially get flattened series of features-subfeatures-subsubfeatures, etc.
    in pre-order. Given
    F   SF  SSF
    1   A   a
            b
        B   c
    2   A   a
    etc. would give
    1AabBc2Aa
    '''
    #print feature.id, feature.type
    yield feature
    if reverse:
        for subfeat in feature.sub_features[::-1]:
            for subsub in flattened_subfeature_iterator(subfeat, reverse=True):
                yield subsub
    else:
        for subfeat in feature.sub_features:
            for subsub in flattened_subfeature_iterator(subfeat):
                yield subsub


def find_cds_start_coordinate(feature):
    '''this will find the "start" of the cds of a gene, which will be the rightmost
    base of a - strand gene, and the first of a + strand gene.  Either way the 
    first CDS listed is the first in the gene, I think.  It needs to deal with the
    fact that the CDS features may be a variable number of layers down from the feature
    passed in.
    NOTE: what the coords mean is VERY misleading.  For a plus strand start is the index
    of the start, and end is one past the index of the last base
    For minus, start is one to the right of the beginning base (reading left to right),
    and end is the index of the leftmost base.
    So, extacting [start:end] will work properly for plus strand,
    and [end:start] will work properly for minus'''

    strandIndex = 0 if feature.strand == 1 else 1
    for feat in flattened_subfeature_iterator(feature):
        if feat.type.lower() == 'cds':
            return int_feature_location(feat)[strandIndex]


def find_cds_end_coordinate(feature):
    '''see find_cds_start_coordinate notes'''
    strandIndex = 1 if feature.strand == 1 else 0
    for feat in flattened_subfeature_iterator(feature, reverse=True):
        if feat.type.lower() == 'cds':
            return int_feature_location(feature)[strandIndex]


def extract_seqrecord_between_outer_cds(rec, ifeat):
    '''This will grab everything between the first and
    last cds of a gene, mainly as a way to chop off any
    UTRs.  This does NOT preperly set the features of the
    returned SeqRecord.
    '''
    if ifeat.sub_features[0].type.lower() == 'mrna':
        if len(ifeat.sub_features) > 1 and ifeat.sub_features[1].type.lower() == 'mrna':
            raise ValueError('Multiple mRNA features found! Pass only one.')
        feat = copy.deepcopy(ifeat.sub_features[0])
    else:
        feat = copy.deepcopy(ifeat)
    
    if feat.strand == -1:
        feat.sub_features.sort(key=lambda x: x.location.start, reverse=True)
    else:
        feat.sub_features.sort(key=lambda x: x.location.start)
    
    '''
    print 'ORIGINAL FEAT'
    print feat
    print feat.location
    print '/ORIGINAL FEAT'
    '''
    s, e = find_cds_start_coordinate(feat), find_cds_end_coordinate(feat)
    print 'start, end', s, e
    tempFeat = copy.deepcopy(feat)
    tempFeat.sub_features = []
    #print
    #print tempFeat.extract(rec)
    #sorting of the cds is necessary for the extraction to work right
    
    tempFeat.location = FeatureLocation(s, e) if feat.strand == 1 else FeatureLocation(e, s)
    
    '''
    print 'TEMP FEAT'
    print tempFeat
    print tempFeat.location
    print '/TEMP FEAT'
    '''
    extracted = tempFeat.extract(rec)
    extracted.name, extracted.id = parse_feature_name(feat), parse_feature_name(feat)

    #DEBUG
    tempFeat.sub_features = collect_features_within_boundaries(feat, min([s, e]), max([s, e]))
    extracted.features = [tempFeat]

    return extracted


def collect_features_within_boundaries(feature, start, end, relativeIndeces=False):
    '''This will pull out the outermost cds fully or partially included in the range,
    and all other intervening features
    this doesn't care about strand, and start should always be < end
    NOTE:the start and end indeces are zero offset
    '''

    #I think that len(feature) actually does a sum of the lengths of all of the subfeatures
    #so, if both exons and cds are included it is longer than the actual length
    realFeatureLength = int_feature_location(feature)[1] - int_feature_location(feature)[0]
    
    '''
    print 'num subfeat', len(feature.sub_features)
    print 'feature loc:', feature.location
    #print 'my feature loc', int_feature_location(feature)
    print 'strand', feature.strand
    print 'feature length:', realFeatureLength
    print 'boundaries:', start, end
    '''

    if start > end:
        #be forgiving here, allowing for reversed coordinates
        if feature.strand == -1:
            print 'flipping start and end coords in collect_features_within_boundaries'
            start, end = end, start
        else:
            raise ValueError('start > end in collect_features_within_boundaries?')

    if relativeIndeces:
        if feature.strand == -1:
            #alignment was of the sequences in sense direction, which didn't know about
            #strand, so need to flip
            start, end = realFeatureLength - end, realFeatureLength - start
        start += int_feature_location(feature)[0]
        end += int_feature_location(feature)[0]
        print 'adjusted boundaries:', start, end

    if feature.type.lower() == 'cds':
        raise ValueError('pass a feature above CDS to find_nearest_cds_boundaries')

    if feature.sub_features[0].type.lower() == 'mrna':
        feature = feature.sub_features[0]

    collectedSubfeatures = []
    #print start, end
    cdsStart = sys.maxint
    cdsEnd = -1
    for sub in feature.sub_features:
        if sub.type.lower() == 'cds':
            fstart, fend = int_feature_location(sub)
            if fend > start and fstart < cdsStart:
                cdsStart = fstart
            if fstart < end and fend > cdsEnd:
                cdsEnd = fend
    #this should grab exons that correspond the cds as well as everything else in between
    for sub in feature.sub_features:
        fstart, fend = int_feature_location(sub)
        if fstart >= cdsStart and fend <= cdsEnd:
            collectedSubfeatures.append(sub)
    if len(collectedSubfeatures) == 0:
        #exit('no features collected?')
        print 'NO CDS FEATURES WITHIN BOUNDARIES: %d, %d!' % (start, end)
    return collectedSubfeatures


def get_first_cds(feature):
    '''this will find the "start" of the cds of a gene, which will be the rightmost
    base of a - strand gene, and the first of a + strand gene.  Either way the 
    first CDS listed is the first in the gene, I think.  It needs to deal with the
    fact that the CDS features may be a variable number of layers down from the feature
    passed in.'''

    strandIndex = 0 if feature.strand == 1 else 1
    for feat in flattened_subfeature_iterator(feature):
        if feat.type.lower() == 'cds':
            return feat
    '''
    if feature.type.lower() == 'cds':
        return feature
    for sub in feature.sub_features:
        if sub.type.lower() == 'cds':
            return sub
        for subsub in sub.sub_features:
            if subsub.type.lower() == 'cds':
                return subsub
            for subsubsub in subsub.sub_features:
                if subsubsub.type.lower() == 'cds':
                    return subsubsub
    '''
    print feature
    raise ValueError('no cds found!')


def sort_feature_list_by_id(recList):
    recList.sort(key=lambda rec: rec.features[0].qualifiers['ID'])


def sort_feature_list(recList):
    '''allow the passed object to be either a list of SeqRecords, which will be sorted based on their
    first feature name, or a single SeqRecord, with a list of features that is to be sorted.
    Prefer the qualifier 'Alias', which is in IRGSP and properly maintains ordering there
    '''
    if isinstance(recList, list):
        qual = 'ID'
        if 'Alias' in recList[0].features[0].qualifiers:
            qual = 'Alias'
        try:
            recList.sort(key=lambda rec: rec.features[0].qualifiers[qual])
        except KeyError:
            sys.stderr.write('ERROR qualifier %s not found\n' % qual)
            for feat in recList.features:
                if qual not in feat.qualifiers:
                    sys.stderr.write('feature:\n%s' % feat)
            exit(1)

    elif isinstance(recList, SeqRecord):
        qual = 'ID'
        if recList.features and 'Alias' in recList.features[0].qualifiers:
            qual = 'Alias'
        try:
            recList.features.sort(key=lambda feat: feat.qualifiers[qual])
        except KeyError:
            sys.stderr.write('ERROR qualifier %s not found\n' % qual)
            for feat in recList.features:
                if qual not in feat.qualifiers:
                    sys.stderr.write('feature:\n%s' % feat)
            exit(1)
    else:
        exit("what is recList?")


def sort_feature_list_by_coordinate(recList):
    '''allow the passed object to be either a list of SeqRecords, which will be sorted based 
    on the start coord of their first feature, or a single SeqRecord, with a list of features that is to be sorted.
    '''
    #recList is a list of SeqRecords, try to sort the SeqRecords by the start coord of their first feature
    #SeqRecords themselves don't have coords
    if isinstance(recList, list):
        try:
            recList.sort(key=lambda rec: rec.features[0].location.start.position)
        except KeyError:
            sys.stderr.write('ERROR could not sort by coordinate')
            for feat in recList.features:
                sys.stderr.write('feature:\n%s' % feat)
            exit(1)
        #for each SeqRecord in the list, want to sort its list of features too.
        #this will hit the second isinstance here
        for rec in recList:
            sort_feature_list_by_coordinate(rec)

    elif isinstance(recList, SeqRecord):
        try:
            recList.features.sort(key=lambda feat: feat.location.start.position)
        except KeyError:
            sys.stderr.write('ERROR could not sort by coordinate')
            for feat in recList.features:
                sys.stderr.write('feature:\n%s' % feat)
            exit(1)
    #reclist is a feature, want to sort it's subfeatures
    elif isinstance(recList, Bio.SeqFeature):
        try:
            print recList
            exit()
            recList.sub_features.sort(key=lambda feat: feat.location.start.position)
        except KeyError:
            sys.stderr.write('ERROR could not sort by coordinate')
            for feat in recList.sub_features:
                sys.stderr.write('feature:\n%s' % feat)
            exit(1)
    else:
        exit("what is recList?")


def make_gene_order_map(geneOrderFilename):
    '''read the indicated file, which should be a simple file with columns indicating gene number and gene name'''
    mapFile = open(geneOrderFilename, 'rb')
    splitMap = [ line.strip().split() for line in mapFile ]
    mapDict = dict([ (line[1], int(line[0])) for line in splitMap ])
    '''
    mapDict = {}
    for line in splitMap:
        mapDict[line[1]] = int(line[0])
    '''
    return mapDict


def adjust_feature_coords(features, delta):
    '''Shift all feature and subfeature coords by delta'''
    for feature in features:
        start, end = feature.location.start.position + delta, feature.location.end.position + delta
        feature.location = FeatureLocation(start, end)
        adjust_feature_coords(feature.sub_features, delta)


def remove_features(features, namesToRemove):
    '''Remove any features with the specified types'''
    featsToDelete = []
    for feature in features:
        for rem in namesToRemove:
            if feature.type.lower() == rem.lower():
                featsToDelete.append(feature)
    for f in featsToDelete:
        features.remove(f)
    #now remove sub_features from anything not alreay removed
    for feature in features:
        remove_features(feature.sub_features, namesToRemove)


def get_features_by_name(feature, name):
    name = name.lower()
    if feature.type.lower() == name:
        return [feature]
    featsToReturn = []
    for feat in flattened_subfeature_iterator(feature):
        if feat.type.lower() == name:
            featsToReturn.append(feat)

    '''
    for sub in feature.sub_features:
        if sub.type.lower() == name:
            featsToReturn.append(sub)
        for subsub in sub.sub_features:
            if subsub.type.lower() == name:
                featsToReturn.append(subsub)
            for subsubsub in subsub.sub_features:
                if subsubsub.type.lower() == name:
                    featsToReturn.append(subsubsub)
    '''
    return featsToReturn


def remove_stop_codon_from_feature(feature):
    '''this just chops off the last three bases of the gene, adjusting gene, mRNA, CDS and exon locations
    it is a simpler version of adjust_for_out_of_phase_cds'''
    assert feature.type == 'gene'

    start, end = int_feature_location(feature)
    adjust = (end, end - 3) if feature.strand == 1 else (start, start + 3)
    #easiest thing to do here will just be to find all mention of a given coordinate and adjust it,
    #which will work for cds, exon, mRNA, gene, etc.
    toSearch = [ feature ]
    for sub in feature.sub_features:
        toSearch.append(sub)
        for subsub in sub.sub_features:
            toSearch.append(subsub)
            for subsubsub in subsub.sub_features:
                toSearch.append(subsubsub)

    for s in toSearch:
        start, end = int_feature_location(s)
        if start == adjust[0]:
            s.location = FeatureLocation(adjust[1], end)
        if end == adjust[0]:
            s.location = FeatureLocation(start, adjust[1])


def append_string_to_feature_names(feature, appStr):
    '''add _appStr to the contents of the following qualifiers, if they exist'''
    for field in ['ID', 'Alias', 'Name', 'Parent']:
        if field in feature.qualifiers:
            feature.qualifiers[field][0] += '_%s' % str(appStr)

    for sub in feature.sub_features:
        append_string_to_feature_names(sub, appStr)


def substitute_feature_names(feature, oldName, newName):
    #print "#######"
    #print feature
    #print oldName
    #print newName
    #print "#######"
    '''replace any appearances of oldName with newName, in id's, qualifiers, etc'''
    if feature.id is not None:
        feature.id = re.sub(oldName, newName, feature.id)
    for field in ['ID', 'Alias', 'Name', 'Parent']:
        if field in feature.qualifiers:
            feature.qualifiers[field][0] = re.sub(oldName, newName, feature.qualifiers[field][0])

    for sub in feature.sub_features:
        substitute_feature_names(sub, oldName, newName)

    if 'Parent' in feature.qualifiers and len(feature.qualifiers['Parent']) > 1:
        exit('MULTIPLE PARENTS?')


