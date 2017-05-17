#!/usr/bin/env python
import copy
import dendropy

#this deals with changes in DendroPy 4
from dendropy.utility import bitprocessing
'''
try:
    from dendropy.calculate import treesplit
except ImportError:
    from dendropy import treesplit
'''

def compat_encode_bipartitions(tree, **kwargs):
    '''Convenience function dealing with different ways of encoding splits in DP4'''
    if hasattr(tree, "encode_bipartitions"):
        if 'delete_outdegree_one' in kwargs:
            val = kwargs.pop('delete_outdegree_one')
            kwargs['collapse_unrooted_basal_bifurcation'] = val
        tree.encode_bipartitions(**kwargs)
    elif not hasattr(tree, "bipartition_encoding") or not tree.bipartition_encoding:
        tree.encode_bipartitions(**kwargs)

class CustomTree(dendropy.Tree):
    '''Override of denodropy.Tree, which redefines equality as an RF distance of 0.'''

    def __eq__(self, other):
        return dendropy.treecalc.symmetric_difference(self, other) == 0
    
    def __hash__(self):
        h = hash(''.join(sorted(self.as_newick_string())))
        return h


class CustomTreeList(dendropy.TreeList):
    '''dendropy.TreeList with a number of functions overriden to tweak functionality
    '''

    def __contains__(self, item):
        '''Overridden function to allow basic use of 'in' keyword.
        
        Identity defined by RF distance of 0.
        if treeX in treelistY: 
            treeX.something()
        NOT very efficient
        NOTE: Other cheaper pre-filtering steps could be added (e.g.
        compare number of nodes) but these are dangerous because of 
        how Dendropy handles rooting of trees
        '''
        for t in self:
            if dendropy.treecalc.symmetric_difference(t, item) == 0:
                return True
        return False

    def frequency_of_identical_trees(self, targetTree):
        '''Return the proportion of trees in self that match targetTree. 
        
        Identity defined by RF distance of 0.
        See NOTE in __contains__ regarding efficiency
        '''
        count = 0
        for tree in self:
            if dendropy.treecalc.symmetric_difference(tree, targetTree) == 0:
                count += 1

        return float(count) / len(self)

    #def masked_frequency_of_split(self, **kwargs):
    def masked_frequency_of_bipartition(self, **kwargs):
        """Adaptation of dendropy.TreeList.frequency_of_bipartition that takes a taxon mask. 
        
        This allows identifying splits on a subset of taxa within a larger tree without
        pruning any tree structures, which is much slower.

        Given a split or bipartition specified as:

            - a split bitmask given the keyword 'split_bitmask'
            - a list of `Taxon` objects given with the keyword `taxa`
            - a list of taxon labels given with the keyword `labels`
            - a list of oids given with the keyword `oids`

        this function returns the proportion of trees in self in which the 
        split is found.
        """
        partialMask = kwargs["mask"] if "mask" in kwargs else self.taxon_namespace.all_taxa_bitmask()

        if "split_bitmask" in kwargs:
            targetSplit = kwargs["split_bitmask"]
        else:
            targetSplit = self.taxon_namespace.get_taxa_bitmask(**kwargs)
            k = kwargs.values()[0]
            if bitprocessing.num_set_bits(targetSplit) != len(k):
                raise IndexError('Not all taxa could be mapped to split (%s): %s' 
                    % (self.taxon_namespace.split_bitmask_string(targetSplit), k))
        found = 0
        total = 0
        for tree in self:
            tree.compat_encode_bipartitions()
            total += 1
            compSplit = (~targetSplit & partialMask)
            #for test_split in tree.split_edges:
            for test_split in tree.reference_tree.bipartition_encoding:
                if not treesplit.is_compatible(test_split, targetSplit, partialMask):
                    break
                masked_test = (test_split & partialMask)
                if targetSplit == masked_test or compSplit == masked_test:
                    found += 1
                    break

        return float(found) / total

    def masked_frequency_of_splitlist(self, returnMatches=False, **kwargs):
        """As masked_frequency_of_split, but counts trees that contain a list of splits.

        Given a LIST of splits or bipartitions specified as:

            - a split bitmask given the keyword 'split_bitmask'
            - a list of `Taxon` objects given with the keyword `taxa`
            - a list of taxon labels given with the keyword `labels`
            - a list of oids given with the keyword `oids`

        this function returns the proportion of trees in self
        in which all of the splits are found.
        NOTE: This is not that useful in some cases here you call it sucessively with
        different numbers of splits and expect the freqs to add up to 1.0
        """
        #if returnMatches is requested, return matching trees
        matches = []

        partialMask = kwargs["mask"] if "mask" in kwargs else self.taxon_namespace.all_taxa_bitmask()

        if "split_bitmask" in kwargs:
            targetSplits = kwargs["split_bitmask"]
        else:
            split = self.taxon_namespace.get_taxa_bitmask(**kwargs)
            k = kwargs.values()[0]
            if treesplit.count_bits(split) != len(k):
                raise IndexError('Not all taxa could be mapped to split (%s): %s' 
                    % (self.taxon_namespace.split_bitmask_string(split), k))

        found = 0
        total = 0
        for tnum, tree in enumerate(self):
            compat_encode_bipartitions(tree)
            total += 1
            matchedSplits = 0
            incompatible = False
            #work through the required splits
            for num, targetSplit in enumerate(targetSplits):
                compSplit = (~targetSplit & partialMask)
                #work through the splits in this tree
                for test_split in tree.bipartition_encoding:
                    #mask out unimportant taxa
                    masked_test = (test_split.split_bitmask & partialMask)
                    #don't need to test anything if masked_test is empty 
                    #(i.e., no taxa in partialMask appear on opposite sides of test_split
                    if masked_test:
                        #if not treesplit.is_compatible(test_split.split_bitmask, targetSplit, partialMask):
                        if not dendropy.Bipartition.is_compatible_bitmasks(test_split.split_bitmask, targetSplit, partialMask):
                            incompatible = True
                            break
                        elif targetSplit == masked_test or compSplit == masked_test:
                            matchedSplits += 1
                            break
                if incompatible:
                    break
            if not incompatible and matchedSplits == len(targetSplits):
                found += 1
                if returnMatches:
                    matches.append(copy.deepcopy(tree))
        if returnMatches:
            return float(found) / total, matches
        else:
            return float(found) / total

    def generate_all_trees_for_taxon_list(self, taxon_list, min_bipartitions=None, max_bipartitions=None, criterion=None):
        '''Will call functions to generate newick strings representing all possible trees for taxon set.  
        
        Work must be done here to make that list fully unique allowing for differences in tree rotation.  
        Can pass min and max bipartitions to control resolvedness of trees, or omit to only generate fully resolved.
        This is impractically slow for > 6 taxa.
        >>> TL = MyTreeList()
        >>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd'])
        >>> len(TL)
        3
        >>> TL = MyTreeList()
        >>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e'])
        >>> len(TL)
        15

        #don't think that this works for trees with polytomies
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e'], min_bipartitions=0, max_bipartitions=2)
        #>>> len(TL)
        #26
        
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e', 'f'])
        #>>> len(TL)
        105
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e', 'f', 'g'])
        #>>> len(TL)
        945
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
        #>>> len(TL)
        10395
        '''
        ntax = len(taxon_list)

        #default to fully resolved trees
        min_bipartitions = max(min_bipartitions, 0) if min_bipartitions else ntax - 3
        max_bipartitions = max(max_bipartitions, ntax - 3) if max_bipartitions else ntax - 3

        componentLists = combine_components_and_uniqueify(taxon_list, min_components=ntax - max_bipartitions, max_components=ntax - min_bipartitions, criterion=criterion)

        fullString = ''
        for componentList in componentLists:
            #the component list is a list of tuples, so actually the repr is exactly newick representation
            fullString += repr(componentList) + ';'
        self.read_from_string(fullString, 'newick')

        #print len(self)
        #for t in self:
        #    print t.as_newick_string()
        newList = TreeList(self, taxon_namespace=self.taxon_namespace)
        #print len(newList)
        self[:] = []
        for num, tr in enumerate(newList):
            #if tr not in TreeList(self[num+1:]):
            if tr not in self:
                self.append(tr)
                #print tr.as_newick_string()

    def as_python_source(self, tree_list_name=None, tree_list_args=None, oids=False):
        """As dendropy.TreeList.as_python_source, but source instantiates a CustomTreeList
        
        Returns string that will rebuild this tree list in Python.
        """
        p = []

        if tree_list_name is None:
            tree_list_name = "tree_list_%s" % id(self)

        if self.label is not None:
            label = "'" + self.label + "'"
        else:
            label = "None"
        if oids:
            oid_str = ', oid="%s"' % self.oid
        else:
            oid_str = ""
        if tree_list_args is None:
            tree_list_args = ""
        else:
            tree_list_args = ", " + tree_list_args
        p.append("%s = CustomTreeList(label=%s%s%s)" 
            % (tree_list_name,
               label,
               oid_str,
               tree_list_args))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        taxon_map = {}
        for taxon in self.taxon_namespace:
            tobj_name = taxon_obj_namer(taxon)
            if taxon.label is not None:
                label = "'" + taxon.label + "'"
            else:
                label = "None"
            if oids:
                oid_str = ', oid="%s"' % taxon.oid
            else:
                oid_str = ""
            p.append("%s = %s.taxon_namespace.require_taxon(label=%s%s)" 
                % (tobj_name,
                   tree_list_name,
                   label,
                   oid_str))
            taxon_map[taxon] = tobj_name

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for tree in self:
            tree_obj_name = "tree_%s" % id(tree)
            if tree.label is not None:
                label = "'" + tree.label + "'"
            else:
                label = "None"
            if oids:
                oid_str = ', oid="%s"' % tree.oid
            else:
                oid_str = ""
            p.append("%s = dendropy.Tree(label=%s, taxon_namespace=%s.taxon_namespace%s)" 
                % (tree_obj_name,
                   label,
                   tree_list_name,
                   oid_str))
            p.append("%s.append(%s, reindex_taxa=False)" % (tree_list_name, tree_obj_name))
            if oids:
                p.append("%s.seed_node.oid = '%s'" % (tree_obj_name, tree.seed_node.oid))
            for node in tree.preorder_node_iter():
                for child in node.child_nodes():
                    if node is tree.seed_node:
                        nn = "%s.seed_node" % tree_obj_name
                    else:
                        nn = node_obj_namer(node)
                    if child.label is not None:
                        label = "'" + child.label + "'"
                    else:
                        label = "None"
                    if child.taxon is not None:
                        ct = taxon_obj_namer(child.taxon)
                    else:
                        ct = "None"
                    if oids:
                        oid_str = ', oid="%s"' % child.oid
                    else:
                        oid_str = ""
                    p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                            (node_obj_namer(child),
                             nn,
                             label,
                             ct,
                             child.edge.length,
                             oid_str))
                    if oids:
                        p.append('%s.edge.oid = "%s"' % (node_obj_namer(child), child.edge.oid))

        return "\n".join(p)

    '''
    #haven't finished implmenting this yet
    def resolve_polytomies(self, source_tree, update_splits=False, rng=None):
        """
        Arbitrarily resolve polytomies using 0-length splits.

        If `rng` is an object with a sample() method then the polytomy will be
            resolved by sequentially adding (generating all tree topologies
            equiprobably
            rng.sample() should behave like random.sample()
        If `rng` is not passed in, then polytomy is broken deterministically by
            repeatedly joining pairs of children.
        """
        polytomies = []
        indeces = []
        for node in self.postorder_node_iter():
            nchild = len(node.child_nodes())
            if nchild > 2:
                polytomies.append(node)
                indeces.append(range(nchild))

        for node in polytomies:
            children = node.child_nodes()
            nc = len(children)
            if nc > 2:
                if rng:
                    to_attach = children[2:]
                    for child in to_attach:
                        node.remove_child(child)
                    attachment_points = children[:2] + [node]
                    while len(to_attach) > 0:
                        next_child = to_attach.pop()
                        next_sib = rng.sample(attachment_points, 1)[0]
                        next_attachment = Node()
                        p = next_sib.parent_node
                        p.add_child(next_attachment)
                        p.remove_child(next_sib)
                        next_attachment.add_child(next_sib)
                        next_attachment.add_child(next_child)
                        attachment_points.append(next_attachment)
                else:
                    while len(children) > 2:
                        nn1 = Node()
                        nn1.edge.length = 0
                        c1 = children[0]
                        c2 = children[1]
                        node.remove_child(c1)
                        node.remove_child(c2)
                        nn1.add_child(c1)
                        nn1.add_child(c2)
                        node.add_child(nn1)
                        children = node.child_nodes()
        if update_splits:
            self.update_splits()
    '''


