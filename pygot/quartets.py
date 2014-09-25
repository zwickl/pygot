#!/usr/bin/env python
import sys
import random


class locus_support_details(object):
    '''The basic details of the results of a particular locus (alignment) tree search for a particular triplet + outgroup of 
    taxa.  i.e., one line from the old table format
    can pass the individual elements that would be extracted from a line like
    ../../..//mafft.cds/analyses/boot.garli/runs//aligned.blink.00025.00013.11T.noDupes.2751C.boot.tre   0.5900  0.0000  0.0450  0.7117  0.1217  0.1667  0   0   0   0   0   0   0   0   (O. sativaj AA,O. punctata BB)
    or just pass the line and have the __init__ parse it
    '''
    def __init__(self, treefile=None, support_value_list=None, taxon_locus_coordinate_list=None, taxon_locus_order_list=None, locus_info=None, support_threshold=0.5):

        if locus_info:
            if isinstance(locus_info, list):
                self.locus_string = '\t'.join(locus_info)
                split_info = locus_info
            else:
                self.locus_string = locus_info.strip()
                #need to split on tabs here because some columns can have spaces within
                split_info = locus_info.strip().split('\t')
            treefile = treefile or split_info[0]
            support_value_list = support_value_list or split_info[1:4]
            taxon_locus_coordinate_list = taxon_locus_coordinate_list or split_info[7:11]
            taxon_locus_order_list = taxon_locus_order_list or split_info[11:15]

        self.treefile = treefile
        self.support_threshold = support_threshold
        
        self.support_value_list = [ float(val) for val in support_value_list ]
        support_sum = sum( self.support_value_list )
        if support_sum:
            self.rescaled_support_value_list = [ val / support_sum for val in self.support_value_list if support_sum ]
        else:
            self.rescaled_support_value_list = [ 1.0 / 3.0 for _ in self.support_value_list ]

        self.taxon_locus_coordinate_list = [ int(val) for val in taxon_locus_coordinate_list ]
        self.taxon_locus_order_list = [ int(val) for val in taxon_locus_order_list ]

        for num, val in enumerate(self.support_value_list):
            if val > self.support_threshold:
                self.supported_tree_num = num
                break
        else:
            self.supported_tree_num = -1

    def __repr__(self):
        if self.locus_string:
            return self.locus_string
        return '%s\t%s\t%s\t%s\t%s' % (self.treefile, '\t'.join(self.support_value_list), '\t'.join(self.rescaled_support_value_list), '\t'.join(self.taxon_locus_coordinate_list), '\t'.join(self.taxon_locus_order_list))

    def __hash__(self):
        return hash(self.treefile)

    def __eq__(self, other):
        return other and self.treefile == other.treefile

    def __ne__(self, other):
        return not self.__eq__(other)


class triplet_support_details(list):
    '''A list of the results of many loci (locus_support_details) allowing various summaries and manipulations.
    Derived from type list, so list operations can be done on self.'''

    def __init__(self, locus_list=None, filename=None):
        if filename:
            if locus_list:
                sys.exit('pass either a locus_list or filename to triplet_support_details.__init__')
            inLines = [ line.strip() for line in open(filename, 'rb') if not 'file' in line]
            inLines = [ line.split('\t') for line in inLines if len(line) ]

            locus_list = [ locus_support_details(line[0], line[1:4], line[7:11], line[11:15]) for line in inLines ]
        self.extend(locus_list or []) 

    def __getslice__(self, i, j):
        #this is some esoteric stuff needed for slicing in addition to __getitem__
        return triplet_support_details(self[max(0, i):max(0, j):])
        #return slice(max(0, i), max(0, j))

    '''
    def __getitem__(self, key):
        if isinstance(key, slice):
            return triplet_support_details( [self[i] for i in xrange(*key.indices(len(self)))] )
        else:
            return self[key]
    '''

    def coordinate_sort(self, taxon=0):
        self.sort(key=lambda loc: loc.taxon_locus_coordinate_list[taxon])

    def order_sort(self, taxon=0):
        self.sort(key=lambda loc: loc.taxon_locus_order_list[taxon])

    def tree_num_list(self):
        return [ loc.supported_tree_num for loc in self ]

    def coordinate_list(self, taxon=0):
        return [ loc.taxon_locus_coordinate_list[taxon] for loc in self ]

    def treefile_list(self):
        return [ loc.treefile for loc in self ]

    def locus_distance_list(self, taxon=0):
        dists = []
        coord = self[0].taxon_locus_coordinate_list[taxon] 
        for loc in self[1:]:
            dists.append(loc.taxon_locus_coordinate_list[taxon] - coord)
            coord = loc.taxon_locus_coordinate_list[taxon]
        return dists

    def informative_tree_num_list(self):
        return [ loc.supported_tree_num for loc in self if loc.supported_tree_num != -1 ]

    def extract_informative_loci(self):
        return triplet_support_details( [ loc for loc in self if loc.supported_tree_num != -1 ] )

    def extract_loci_with_function(self, func):
        return triplet_support_details( [ loc for loc in self if func(loc) ] )

    def extract_loci_in_coordinate_window(self, start, end, taxon):
        ret = triplet_support_details( [ loc for loc in self if start <= loc.taxon_locus_coordinate_list[taxon] < end ] )
        ret.coordinate_sort()
        return ret

    def split_into_runs(self, informative_only=True):
        runs = []
        
        if informative_only:
            loci = self.extract_informative_loci()
        else:
            loci = self
        
        match_start = 0
        topology = self[0].supported_tree_num
        for num, locus in enumerate(self[1:], 1):
            if locus.supported_tree_num != topology or num == len(self):
                runs.append(triplet_support_details(self[match_start:num]))
                match_start = num
            topology = locus.supported_tree_num
        
        return runs

    def shuffle(self):
        random.shuffle(self)

    def supported_tree_counts(self):
        return [ len([loc for loc in self if loc.supported_tree_num == tnum]) for tnum in [ 0, 1, 2, -1 ] ]

    def supported_tree_proportions(self):
        if self:
            return [ supp / float(len(self)) for supp in self.supported_tree_counts() ]
        else:
            return [ 0.0, 0.0, 0.0 ]

    def windowed_tree_counts(self, window_size, window_stride, taxon=0):
        largest_coordinate = max([locus.taxon_locus_coordinate_list[taxon] for locus in self])
        return [ self.extract_loci_in_coordinate_window(start, start + window_size, taxon).supported_tree_counts() for start in xrange(1, largest_coordinate + 1, window_stride) ]

    def windowed_tree_proportions(self, window_size, window_stride, taxon=0):
        largest_coordinate = max([locus.taxon_locus_coordinate_list[taxon] for locus in self])
        return [ self.extract_loci_in_coordinate_window(start, start + window_size, taxon).supported_tree_proportions() for start in xrange(1, largest_coordinate + 1, window_stride) ]

    def count_adjacent_tree_matches(self):
        return sum( [self[num].supported_tree_num == self[num + 1].supported_tree_num for num in xrange(len(self) - 1)] )

    def match_distance_tuples(self, taxon=0):
        return [ (self[num].supported_tree_num == self[num + 1].supported_tree_num, self[num + 1].taxon_locus_coordinate_list[taxon] - self[num].taxon_locus_coordinate_list[taxon]) for num in xrange(len(self) - 1)] 

    def max_coordinate(self, taxon=0):
        return max(locus.taxon_locus_coordinate_list[taxon] for locus in self)

    def support_mean(self, topology=0, informative_only=True, only_if_topology_preferred=False):
        if informative_only:
            loci = self.extract_informative_loci()
        else:
            loci = self
        if only_if_topology_preferred:
            sublist = [locus.support_value_list[topology] for locus in loci if locus.supported_tree_num == topology]
            return sum(sublist) / float(len(sublist))
        else:
            return sum([locus.support_value_list[topology] for locus in loci]) / float(len(loci))

    def support_sum(self, topology=0, informative_only=True, only_if_topology_preferred=False):
        if informative_only:
            loci = self.extract_informative_loci()
        else:
            loci = self
        if only_if_topology_preferred:
            sublist = [locus.support_value_list[topology] for locus in loci if locus.supported_tree_num == topology]
            return sum(sublist)
        else:
            return sum([locus.support_value_list[topology] for locus in loci])

    def __repr__(self):
        outstr = ''
        for locus in self:
            outstr += '%s\n' % locus
        return outstr


def generate_quartet_lists_from_file(filename='quartetList'):
    '''Generates the possible quartets of taxa from specifically formatted taxon strings.

    Can be one per file line, i.e., 
    one, three, four, five
    two, three, four, five
    Or, you can use colons to indicate groups of taxa to make quartets combinatorically,
    i.e. this would be equivalent to the two lines above:
    one:two, three, four, five
    '''

    allQuartets = set()
    comboQuartets = []

    for line in open(filename, 'rb'):
        #ignore blank lines and "comments" starting with #
        if len(line.strip()) > 0 and line[0] != '#':
            if ':' in line:
                comboQuartets.append([tax.strip() for tax in line.split(',')])
            else:
                allQuartets |= set([tuple([tax.strip() for tax in line.split(',')])])

    quartsFromCombos = []
    for combo in comboQuartets:
        thisCombo = []
        for subset in combo:
            thisSubset = []
            for el in subset.split(':'):
                if thisCombo:
                    thisSubset.extend([ p + [el] for p in thisCombo ])
                else:
                    thisSubset.append([el])
            thisCombo = thisSubset
        quartsFromCombos.extend(thisCombo)
    allQuartets |= set([tuple(q) for q in quartsFromCombos])

    #allQuartets was initially a set of tuples, to ensure that duplicate quarts weren't kept, now
    #convert back to list of lists
    return [list(q) for q in allQuartets]


