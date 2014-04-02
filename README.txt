pygot
=====

PhYloGenOmic Tools : Python scripts for processing and summarizing phylogenomic analyses

Initial functionality will include scripts used in 
Zwickl D.J., Stein J.C., Wing R.A., Ware D. and Sanderson M.J. 2014. Disentangling Methodological 
    and Biological Sources of Gene Tree Discordance on Oryza (Poaceae) Chromosome 3. Syst. Biol. In press. 

Note that it will be difficult to understand what this all means without reading the above paper.

Currently contains functionality to:
summarize sets of bootstrap trees (required for the following)
plot cumulative support diagrams
plot flux diagrams

Documentation and scripts are being improved, so email me if you want to use this and have any trouble
or questions
-Derrick Zwickl
zwickl@email.arizona.edu

##############################

See masterScript.sh in examples directory if you'd like to jump in without reading the following ...


Terminology (as used in the paper):
locus       - A set of sequences considered an orthologous cluster.  Note that not all loci need to 
              have the same number or set of taxa.
treatment   - Arbitrary protocol for progressing from the sequences at a locus to a set of bootstrap
              trees or samples from a Bayesian posterior distribution.  Includes at least alignment 
              and tree inference.  Treatments could differ from one another in the alignment algorithm 
              used, whether alignments are cleaned, tree inference method, evolutionary model, etc.
alignment   - The result of applying a particular treatment to a particular locus.
quartet     - Set of four taxa used as basis of tree summaries.  Can also be considered a "rooted triplet"
              if you are willing to choose one as an outgroup.  Three possible resolutions of the 
              quartet (or triplet) exist, which a fourth possibility being an unresolved polytomy.

Input: 
    - samples of bootstrap or posterior trees infered for each of a set of alignments. Note that
      bootstrap or posterior majority rule consensus trees are NOT sufficient.  The individual tree
      samples are needed.
    - specified quartets

Minimal workflow is this:
1. Make alignments for a particular set of taxa
2. Obtain bootstrap or posterior tree samples for those alignments
3. Specify quartets of taxa in a file, in a particular format
4. Run the dendropyScoreTriples.py script to analyze the topological relationships of the quartets
    of taxa embedded within your sampled trees
5. Analyze the output of dendropyScoreTriples.py, either directly looking at the frequency of 
    each of the triplet resolutions, or making a "cumulative support diagram" using 
    cumulativeAndPieFigure.py, a much richer summary


If you used multiple treatments for a single set of loci, those treatments can additionally be 
compared by means of a "flux diagram", made with plotFlux.py.

