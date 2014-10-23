pygot
=====

PhYloGenOmic Tools : Python scripts for processing and summarizing phylogenomic analyses

Initial functionality will include scripts used in 
Zwickl DJ, Stein JC, Wing RA, Ware D, Sanderson MJ. 2014. Disentangling methodological 
    and biological sources of gene tree discordance on oryza (poaceae) chromosome 3. 
    Systematic Biology. 63(5):645-659.

Note that it will be difficult to understand what this all means without reading the 
above paper.

Currently contains functionality to:
quartet summaries:
    summarize sets of bootstrap trees (required for the following)
    plot cumulative support diagrams
    plot flux diagrams
alignment block-shifts:
    detect block-shifts
    mask block-shifts

Documentation and scripts are being improved, so email me if you want to use this and
have any trouble or questions
-Derrick Zwickl
zwickl@email.arizona.edu

##############################

INSTALLATION:

This could use some work and testing yet, but the typical python install procedure of
python setup.py install 
should work (you might need to add sudo to the start of that line).  

If you get errors with that command, you might need to update your python "setuptools" 
package, i.e.:
easy_install setuptools
(again, you might need to add sudo to the start of that line)

Or you can just run things in the example directories, as demonstrated in the sample 
scripts.  (although you'll need to install the below dependencies yourself)

If you install it globally you may need to put the script install location into your 
PATH environment variable, or copy the script files (in the scripts directory) 
somewhere else in your path.

Dependencies:
-The excellent Dendropy library:

Sukumaran J, Holder MT. 2010. DendroPy: a Python library for phylogenetic computing. 
    Bioinformatics. 26(12):1569-1571.
http://pythonhosted.org/DendroPy/

is required for all quartet based summaries. The pygot setup procedure above will 
attempt to install Dendropy for you if necessary.
I highly recommend Dendropy for your phylogenetic computing needs!

-The matplotlib library (http://matplotlib.org/) is required for making any figures.  
This can be non-trivial to get installed, and downloading an installer package for 
your system generally seems to work better than using easy_install or pip.  The pygot 
setup script will not attempt to install this for you.

-The block-shift related functions have no external dependencies.

##############################

QUARTET BASED SUMMARIES, CUMULATIVE SUPPORT DIAGRAMS AND FLUX DIAGRAMS

See masterScript.sh in examples directory if you'd like to jump in without reading the 
following ...

Terminology (as used in the paper):
locus       - A set of sequences considered an orthologous cluster.  Note that not all 
              loci need to have the same number or set of taxa.
treatment   - Arbitrary protocol for progressing from the sequences at a locus to a 
              set of bootstrap trees or samples from a Bayesian posterior distribution.  
              Includes at least alignment and tree inference.  Treatments could differ 
              from one another in the alignment algorithm used, whether alignments are 
              cleaned, tree inference method, evolutionary model, etc.
alignment   - The result of applying a particular treatment to a particular locus.
quartet     - Set of four taxa used as basis of tree summaries.  Can also be considered 
              a "rooted triplet" if you are willing to choose one as an outgroup.  Three 
              possible resolutions of the quartet (or triplet) exist, which a fourth 
              possibility being an unresolved polytomy.

Input: 
    - samples of bootstrap or posterior trees infered for each of a set of alignments. 
      Note that bootstrap or posterior majority rule consensus trees are NOT sufficient.  
      The individual tree samples are needed.
    - specified quartets

Minimal workflow is this:
1. Make alignments for a particular set of taxa
2. Obtain bootstrap or posterior tree samples for those alignments
3. Specify quartets of taxa in a file, in a particular format
4. Run the dendropyScoreTriples.py script to analyze the topological relationships 
    of the quartets of taxa embedded within your sampled trees
5. Analyze the output of dendropyScoreTriples.py, either directly looking at the 
    frequency of each of the triplet resolutions, or making a "cumulative support 
    diagram" using cumulativeAndPieFigure.py, a much richer summary

If you used multiple treatments for a single set of loci, those treatments can 
additionally be compared by means of a "flux diagram", made with plotFlux.py.

Note that there are a large number of command line arguments to the plotting 
scripts that can customize exactly how the figures look.

#############################

DETECTION OF BLOCK-SHIFT ARTIFACTS

Block shifts are defined as stretches of an alignment in which one (or a few) 
taxon differs from generally conserved columns.  None of the columns individually 
are out of the ordinary, but taken as a whole the region is generally clearly 
misaligned.  Because of this typical alignment "cleaning" methods (GBLOCKS, 
Aliscore, etc.) don't pick up on the error.  This is not a typical or realistic 
example, but a block shift might look like this:

T1 ACGTACGTACGT
T2 ACGTACGTACGT
T3 ACGTACGTACGT
T4 ACGAAAAAAA--

The script simply scans along each sequence in a NON-INTERLEAVED nexus alignment 
file, and reports on regions that appear to be block-shifts in any of the taxa.  
False positive block-shifts can be sometimes identified, especially in distantly 
related taxa.  Take a look at the identified regions of your alignments to get a 
feeling for whether identified regions are real or not.  The details of the 
detection algorithm can be changed with command line flags.

This would be the output for the above example, indicating that bases 4-10 in 
taxon t4 are identified as a block-shift:

alignment.nex
T1
T2
T3
T4 4-10

The masked alignment output would look like this:
T1 acgtacgtacgt
T2 acgtacgtacgt
T3 acgtacgtacgt
T4 acgNNNNNNN--

#############################
