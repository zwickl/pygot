#!/usr/bin/env python
import sys
from os import path
import re
from itertools import izip, cycle, chain, repeat, islice
import argparse 
from math import ceil

try:
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
except ImportError:
    sys.exit('matplotlib needs to be installed to use this module')

from pygot.utils import flatten_array, linspace, find_shortest_unique_leading_substrings


def padded(toPad, n, fillvalue=None):
    '''return an iterator of n elements, either slicing toPad or padding it
    >>> [ e for e in padded((1, 2, 3), 2) ]
    [1, 2]
    >>> [ e for e in padded((1, 2, 3), 4, 0) ]
    [1, 2, 3, 0]
    '''
    padNum = n - len(toPad)
    if padNum <= 0:
        return islice(toPad, n)
    return chain(toPad, repeat(fillvalue, padNum))


def cycle_n_elements(toCycle, n):
    '''return an iterator of length n, drawing items from iterable toCycle, cycling if n > len(toCycle)
    If n < len(toCycle), identical to padded(toCycle, n))
    >>> [ e for e in cycle_n_elements([1, 2, 3], 4) ]
    [1, 2, 3, 1]
    >>> [ e for e in cycle_n_elements([1, 2, 3], 2) ]
    [1, 2]
    '''
    return iter(el for el, _ in izip(cycle(toCycle), xrange(n)))


def cycle_n_elements_m_times(toCycle, n, m):
    '''this is confusing - repeat the following m times, chaining into a single iterable:
    draw n items from toCycle, cycling if n > len(toCycle)
    (see cycle_n_elements)
    >>> [ el for el in cycle_n_elements_m_times([1, 2, 3], 4, 2) ]
    [1, 2, 3, 1, 1, 2, 3, 1]
    '''
    return chain.from_iterable(cycle_n_elements(toCycle, n) for __ in xrange(m))


def read_files_and_split_columns_as_strings(filenames, skipRows=None, ignorePatts=None, allowMissing=False):
    '''This is very much like functions from the csv module, but reads multiple files and allows 
    mechanisms for skipping rows based on a regex or skipping a specified number of leading 
    rows'''
    data = []
    compiledPatts = []

    if isinstance(filenames, str):
        filenames = [filenames]

    if ignorePatts:
        if isinstance(ignorePatts, str):
            ignorePatts = [ignorePatts]
        compiledPatts = [ re.compile(patt) for patt in ignorePatts ]
    
    for f in filenames:
        if not path.exists(f) and allowMissing:
            sys.stderr.write('WARNING: skipping missing file %s\n' % f)
        else:
            theseData = []
            fileData = [ l.strip() for l in open(f, 'rb') ]
            for line in fileData[skipRows:]:
                for cpat in compiledPatts:
                    if cpat.search(line):
                        break
                #the freaky for-else construct ...
                else:
                    theseData.append(line.split())
            data.append(theseData)
    return data


def translate_annotation_code(c):
    #code to indicate intron/exon regions
    '''
    G = full gap

    E, F = full exon (no gap, gap)
    I, J = full intron (no gap, gap)
    M, N = mixed intron/exon (no gap, gap)
    '''
    translate = {'E':5, 'F':5, 'I':-5, 'J':-5, 'M':0, 'N':0}
    return translate[c]


def translate_annotation_code_string(lines, column):
    translate = {'E':5, 'F':5, 'I':-5, 'J':-5, 'M':0, 'N':0}
    return [ translate[line[column]] for line in lines ]


def string_to_float_or_zero(val):
    try:
        return float(val)
    except ValueError:
        return 0.0


def path_to_plot_title(filename, sep='\n'):
    '''
    ../mafft.cds/
    ../mafft.cds.allBadAlignAA/
    ../mafft.cds.gblocks/
    ../mafft.genes/
    ../mafft.genes.gblocks/
    ../mafft.genes.intronsStripped/
    ../mafft.genes.noFullIntrons/
    ../mafft.genes.stripNs/
    
    all
    11T

    boot
    ml
    '''

    filename = filename.lower()

    plotTitle = ''
    if 'mafft' in filename:
        plotTitle += 'MAFFT'
    elif 'prank' in filename:
        plotTitle += 'PRANK'
    elif 'muscle' in filename:
        plotTitle += 'MUSCLE'

    plotTitle += sep

    if 'cds.allbadalignaa' in filename:
        plotTitle += 'CDS-BADCOLAA '
    elif 'cds.allbadalign' in filename:
        plotTitle += 'CDS-BADCOL '
    elif 'cds.maskedbadalign.onlyaa' in filename:
        plotTitle += 'CDS-MASKEDAA '
    elif 'cds.maskedbadalign.nobrach' in filename:
        plotTitle += 'CDS-MASKEDNOBRACH '
    elif 'cds.maskedbadalign' in filename:
        plotTitle += 'CDS-MASKED '
    elif 'cds.gblocks.maskedbadalign' in filename:
        plotTitle += 'CDS+GBLOCKS%s-MASKED' % sep
    elif 'cds.gblocks' in filename:
        plotTitle += 'CDS+GBLOCKS'
    elif 'cds' in filename:
        plotTitle += 'CDS'
    elif 'genes.gblocks.maskedbadalign' in filename:
        plotTitle += 'GENE+GBLOCKS%s-MASKED' % sep
    elif 'genes.gblocks' in filename:
        plotTitle += 'GENE+GBLOCKS'
    elif 'genes.intronsstripped' in filename:
        #plotTitle += 'GENE-INTRON'
        plotTitle += 'MASKED-INTRON'
    elif 'genes.nofullintrons' in filename:
        plotTitle += 'GENE-FULLI'
    elif 'genes.allbadalignAA' in filename:
        plotTitle += 'GENE-BADCOLAA '
    elif 'genes.allbadalign' in filename:
        plotTitle += 'GENE-BADCOL '
    elif 'genes.maskedbadalign.onlyaa' in filename:
        plotTitle += 'GENE-MASKEDAA '
    elif 'genes.maskedbadalign.nobrach' in filename:
        plotTitle += 'GENE-MASKEDNOBRACH '
    elif 'genes.maskedbadalign' in filename:
        plotTitle += 'GENE-MASKED '
    elif 'genes.stripns' in filename:
        plotTitle += 'GENE-STRIPNs '
    elif 'genes' in filename:
        plotTitle += 'GENE'
    else:
        plotTitle += 'UNKNOWN'

    if 'ssr' in filename:
        plotTitle += 'SSR'

    return plotTitle


def simpler_path_to_plot_title(filename, sep='\n'):
    '''
    This just simplies the names a bit.  maskedBadAlign.noBrach is just ALMASK
    and other almasks are unknown.  Also, is +GBLOCKS and +ALMASK
    ../mafft.cds/
    ../mafft.cds.allBadAlignAA/
    ../mafft.cds.gblocks/
    ../mafft.genes/
    ../mafft.genes.gblocks/
    ../mafft.genes.intronsStripped/
    ../mafft.genes.noFullIntrons/
    ../mafft.genes.stripNs/
    
    all
    11T

    boot
    ml
    '''

    filename = filename.lower()

    plotTitle = ''
    if 'mafft' in filename:
        plotTitle += 'MAFFT'
    elif 'prank' in filename:
        plotTitle += 'PRANK'
    elif 'muscle' in filename:
        plotTitle += 'MUSCLE'

    plotTitle += sep

    if 'cds.allbadalignaa' in filename:
        plotTitle += 'CDS-BADCOLAA '
    elif 'cds.allbadalign' in filename:
        plotTitle += 'CDS-BADCOL '
    elif 'cds.gblocks.maskedbadalign' in filename:
        plotTitle += 'CDS+GBLOCKS%s+ALMASK' % sep
    elif 'cds.maskedbadalign.nobrach' in filename:
        plotTitle += 'CDS+ALMASK '
    elif 'cds.maskedbadalign ' in filename:
        plotTitle += 'CDS+UNKNOWN '
    elif 'cds.gblocks' in filename:
        plotTitle += 'CDS+GBLOCKS'
    elif 'cds.trimal' in filename:
        plotTitle += 'CDS+TRIMAL'
    elif 'cds' in filename:
        plotTitle += 'CDS'
    elif 'genes.gblocks.maskedbadalign' in filename:
        plotTitle += 'GENE+GBLOCKS%s+ALMASK' % sep
    elif 'genes.gblocks' in filename:
        plotTitle += 'GENE+GBLOCKS'
    elif 'genes.trimal' in filename:
        plotTitle += 'GENE+TRIMAL'
    elif 'genes.intronsstripped' in filename:
        #plotTitle += 'GENE-INTRON'
        plotTitle += 'MASKED-INTRON'
    elif 'genes.nofullintrons' in filename:
        plotTitle += 'GENE-FULLI'
    elif 'genes.allbadalignAA' in filename:
        plotTitle += 'GENE-BADCOLAA '
    elif 'genes.allbadalign' in filename:
        plotTitle += 'GENE-BADCOL '
    elif 'genes.maskedbadalign.nobrach' in filename:
        plotTitle += 'GENE+ALMASK '
    elif 'genes.maskedbadalign' in filename:
        plotTitle += 'GENE+UNKNOWN '
    elif 'genes.stripns' in filename:
        plotTitle += 'GENE-STRIPNs '
    elif 'genes' in filename:
        plotTitle += 'GENE'
    else:
        plotTitle += 'UNKNOWN'

    return plotTitle


def simpler_path_to_plot_title_nuc_default_cds(filename, sep='\n'):
    '''
    As simpler_path_to_plot_title, but makes the nucleotide level CDS
    alignment the default (i.e., MAFFT CDS), while the protein level is
    distinguished (MAFFT CDS-P)

    This just simplies the names a bit.  maskedBadAlign.noBrach is just ALMASK
    and other almasks are unknown.  Also, is +GBLOCKS and +ALMASK

    EDIT: changing to +GBMASK, +ALIMASK and +BLSMASK for masking methods
          changing CDS-P to CDSPROT

    ../mafft.cds                            ../mafft.cds.maskedBadAlign             ../mafft.genes                          ../mafft.genes.maskedBadAlign
    ../mafft.cds.allBadAlignAA              ../mafft.cds.maskedBadAlign.noBrach     ../mafft.genes.allBadAlign              ../mafft.genes.maskedBadAlign.noBrach
    ../mafft.cds.gblocks                    ../mafft.cds.maskedBadAlign.onlyAA      ../mafft.genes.gblocks                  ../mafft.genes.maskedBadAlign.onlyAA
    ../mafft.cds.gblocks.maskedBadAlign     ../mafft.cds.nuc                        ../mafft.genes.gblocks.maskedBadAlign   ../mafft.genes.noFullIntrons
    ../mafft.cds.gblocks.nuc                ../mafft.cds.nuc.maskedBadAlign.noBrach ../mafft.genes.intronsStripped          ../mafft.genes.stripNs
    
    '''

    filename = filename.lower()

    plotTitle = ''
    if 'mafft' in filename:
        plotTitle += 'MAFFT'
    elif 'prank' in filename:
        plotTitle += 'PRANK'
    elif 'muscle' in filename:
        plotTitle += 'MUSCLE'

    plotTitle += sep

    if 'cds.allbadalignaa' in filename:
        plotTitle += 'CDSPROT-BADCOLAA '
    elif 'cds.allbadalign' in filename:
        plotTitle += 'CDSPROT-BADCOL '
    elif 'cds.gblocks.maskedbadalign' in filename:
        plotTitle += 'CDSPROT+GBMASK%s+BLSMASK' % sep
    elif 'cds.nuc.maskedbadalign.nobrach' in filename:
        plotTitle += 'CDS+BLSMASK '
    elif 'cds.maskedbadalign.nobrach' in filename:
        plotTitle += 'CDSPROT+BLSMASK '
    elif 'cds.maskedbadalign ' in filename:
        plotTitle += 'CDSPROT+UNKNOWN '
    elif 'cds.gblocks.nuc' in filename:
        plotTitle += 'CDS+GBMASK'
    elif 'cds.nuc.trimal' in filename:
        plotTitle += 'CDS+TRIMMASK'
    elif 'cds.nuc.aliscore' in filename:
        plotTitle += 'CDS+ALIMASK'
    elif 'cds.gblocks' in filename:
        plotTitle += 'CDSPROT+GBMASK'
    elif 'cds.trimal' in filename:
        plotTitle += 'CDSPROT+TRIMMASK'
    elif 'cds.aliscore' in filename:
        plotTitle += 'CDSPROT+ALIMASK'
    elif 'cds.nuc' in filename:
        plotTitle += 'CDS'
    elif 'cds' in filename:
        plotTitle += 'CDSPROT'
    elif 'genes.gblocks.maskedbadalign' in filename:
        plotTitle += 'GENE+GBMASK%s+BLSMASK' % sep
    elif 'genes.gblocks' in filename:
        plotTitle += 'GENE+GBMASK'
    elif 'genes.trimal' in filename:
        plotTitle += 'GENE+TRIMMASK'
    elif 'genes.aliscore' in filename:
        plotTitle += 'GENE+ALIMASK'
    elif 'genes.intronsstripped' in filename:
        plotTitle += 'MASKED-INTRON'
    elif 'genes.nofullintrons' in filename:
        plotTitle += 'GENE-FULLI'
    elif 'genes.allbadalignAA' in filename:
        plotTitle += 'GENE-BADCOLAA '
    elif 'genes.allbadalign' in filename:
        plotTitle += 'GENE-BADCOL '
    elif 'genes.maskedbadalign.nobrach' in filename:
        plotTitle += 'GENE+BLSMASK '
    elif 'genes.maskedbadalign' in filename:
        plotTitle += 'GENE+UNKNOWN '
    elif 'genes.stripns' in filename:
        plotTitle += 'GENE-STRIPNs '
    elif 'genes' in filename:
        plotTitle += 'GENE'
    else:
        plotTitle += 'UNKNOWN'

    return plotTitle


def path_to_plot_title_random_supermatrix(filename, sep='\n'):
    plotTitle = path_to_plot_title(filename, sep)

    if 'correctaanoindica' in filename:
        plotTitle += '%scorrect AA-indica' % sep
    elif 'correctallnoindica' in filename:
        plotTitle += '%scorrect All-indica' % sep
    elif 'correctaa' in filename:
        plotTitle += '%scorrect AA' % sep
    elif 'correctall' in filename:
        plotTitle += '%scorrect All' % sep
    elif 'correctBBCC' in filename:
        plotTitle += '%scorrect AA-BB-CC' % sep

    return plotTitle


def filename_to_paren_trees(s, oryza=False):
    if oryza and isinstance(s, str):
        #ORYZA SPECIFIC STUFF - probably better for most users to pass in
        #the paren tree strings with the --paren-labels option
        sp = s.split('/')
        sp = sp[-1]
        sp = sp.split('.')
        
        nameMap = {'sativaj': 'J', 'sativai': 'I'}
        first = nameMap.get(sp[0].lower(), sp[0][0].upper())

        if 'punc' in sp[1].lower() and not 'brach' in s.lower():
            sec = nameMap.get(sp[2].lower(), sp[2][0].upper())
            third = nameMap.get(sp[3].lower(), sp[3][0].upper())
            paren = [ '(%s,(%s,%s))' % (first, sec, third), 
                    '((%s,%s),%s)' % (first, sec, third), 
                    '((%s,%s),%s)' % (first, third, sec),
                    '(%s,%s,%s)' % (first, sec, third) ]
        else:
            sec = nameMap.get(sp[1].lower(), sp[1][0].upper())
            third = nameMap.get(sp[2].lower(), sp[2][0].upper())
            paren = [ '((%s,%s),%s)' % (first, sec, third), 
                    '((%s,%s),%s)' % (first, third, sec), 
                    '(%s,(%s,%s))' % (first, sec, third),
                    '(%s,%s,%s)' % (first, sec, third) ]
    else:
        if isinstance(s, file):
            first, sec, third, fourth = '1', '2', '3', '4'
        else:
            sp = s.split('/')
            sp = sp[-1]
            sp = sp.split('.')
            
            first, sec, third, fourth = find_shortest_unique_leading_substrings(sp[:4])

        paren = [ '((%s,%s),(%s,%s))' % (first, sec, third, fourth), 
                '((%s,%s),(%s,%s))' % (first, third, sec, fourth), 
                '((%s,%s),(%s,%s))' % (first, fourth, sec, third),
                '(%s,%s,%s,%s)' % (first, sec, third, fourth) ]

    return paren


def prepare_plot_kwargs(kwargDict, optionList, fontscaler=1.0):
    '''
    parse and prepare kwargs in various formats, add new options,
    This is mainly for use with the PlottingArgumentParser defined below
    and would be called like this:
    
        superTitleKwargs = prepare_kwargs(opt.additionalSuperTitleKwargs,
                [('fontsize', opt.superTitleFontsize)],
                fontscaler=scaler)
        plt.suptitle(opt.superTitle, **superTitleKwargs)

    where 
    -the first argument is any extra arbitrary kwargs that might have been
        passed on the command line
    -then a list of things to set specifically.  Defaults for those should 
        have been set when the parser was instantiated, but can be overridden
    
    kwargDict - can just be a list of strings, i.e. ['blah=2', 'poo=4']
        or a list of tuples [('blah', 2), ('poo', 4)]
        or an actual dictionary {'blah':2, 'poo':4}
    optionList - list of (key, value) tuples, added to returned dict ONLY if
        not already in it
    fontscaler - value to rescale fonts by, if desired - deprecated and may not work
    returns outKwargs - dict of kwargs to pass to other functions
    
    >>> prepare_plot_kwargs(['blah=2', 'poo=4'], [('fontsize', 10)])
    {'blah': 2.0, 'fontsize': 10.0, 'poo': 4.0}
    >>> prepare_plot_kwargs([('blah', 2), ('poo', 4)], [('fontsize', 10)])
    {'blah': 2.0, 'fontsize': 10.0, 'poo': 4.0}
    >>> prepare_plot_kwargs({'blah':2, 'poo':4}, [('fontsize', 10)])
    {'blah': 2.0, 'fontsize': 10.0, 'poo': 4.0}
    '''

    outKwargs = {}
    if kwargDict:
        error = False
        if isinstance(kwargDict, dict):
            outKwargs = dict(kwargDict)
        elif isinstance(kwargDict, list) or isinstance(kwargDict, set) or isinstance(kwargDict, tuple):
            #There are shorter and easier ways to do the below, but I'm not sure that they guarantee
            #that things are added to the dictionary in list order, which is needed so that later 
            #ones override earlier
            outKwargs = {}
            if isinstance(kwargDict[0], str):
                for key, val in (kw.split('=') for kw in kwargDict):
                    outKwargs[key] = val
            elif len(kwargDict[0]) == 2:
                for key, val in kwargDict:
                    outKwargs[key] = val
            else:
                error = True
        else:
            error = True
        if error:
            sys.exit('kwargDict must be a dictionary, or a list, set or tuple containing strings, \
                    or a list, set or tuple containing strings containing 2-tuples')
   
    #this just ensures that if a kwarg is specified at the command line for something, 
    #it trumps any actual argparse option for it
    for arg, setting in optionList:
        outKwargs.setdefault(arg, setting)

    for key, val in outKwargs.iteritems():
        #try to convert values to numbers where possible
        try:
            outKwargs[key] = float(val)
        except ValueError:
            pass
        except TypeError:
            pass
        if val == 'False':
            outKwargs[key] = False
        elif val == 'True':
            outKwargs[key] = True
        
        #this doesn't really work as intended
        if 'fontsize' in key or 'labelsize' in key:
            outKwargs[key] = outKwargs[key] * fontscaler

    return outKwargs


class ArgparseActionAppendToDefault(argparse.Action):
    '''This is in a way related to the above prepare_plot_kwargs function.
    Normally defaults can be set on argparse options, but will be overridden if the 
    argument appears on the command line.  This will allow arguments passed on the
    command line to simply be appended to the default list.  This would mainly be 
    used for kwargs specified on the command line and a PlottingArgumentParser
    instantiated with some kwargs set as defaults. Because the values that would
    come from the commandline appear later, they should trump earlier ones in the
    prepare_plot_kwargs function.
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        #print '%r %r %r' % (self.dest, self.default, values)
        if not hasattr(self, 'default'):
            raise ValueError('only makes sense to call AppendToDefaultArgparseAction \
                    when default value to argument is defined')
        if not isinstance(self.default, list):
            raise ValueError('only makes sense to call AppendToDefaultArgparseAction \
                    when defaults are in a list')
        if isinstance(values, str):
            values = values.split()

        setattr(namespace, self.dest, self.default + values)


class PlottingHistogramArgumentParser(argparse.ArgumentParser):
    '''A generic argparse parser that comes prepackaged with a bunch of common arguments
    for matplotlib/pyplot preloaded.
    Default values for the created parser can be passed in as keyword arguments.  Pass
    False to remove a built in option.'''

    def __init__(self, **kwargs):
        '''These are the defaults which can be overwridden with keyword arguments passed
        to the constructor.
        '''
        defaultSubplotKwargs = kwargs.pop('defaultSubplotKwargs', [])
        defaultHistKwargs = kwargs.pop('defaultHistKwargs', [])
        defaultNumBins = kwargs.pop('defaultNumBins', 20)
        defaultBarFaceColor = kwargs.pop('defaultBarFaceColor', 'gray')
        
        if defaultNumBins is not False:
            self.add_argument('-n', '--num-bins', dest='numBins', type=int, default=defaultNumBins,
                                help='number of evenly spaced histogram bins (default %s)' % str(defaultNumBins))

        if defaultBarFaceColor is not False:
            self.add_argument('-c', '--bar-color', dest='barFaceColor', type=str, default=defaultBarFaceColor,
                                help='color of histogram bars (default %s)' % str(defaultBarFaceColor))

 
class PlottingArgumentParser(argparse.ArgumentParser):
    '''
    A generic argparse parser that comes prepackaged with a bunch of common arguments
    for matplotlib/pyplot preloaded.
    
    Default values for the created parser can be passed in as keyword arguments.  Pass
    a values of 'SUPPRESS' to remove a built in option.

    e.g.
    
    my_parser = PlottingArgumentParser(defaultColors='SUPPRESS', defaultMarkers='o*x')
    
    creates a parser with the defaults as defined below, except the default for
    the markers is set to 'o*x' and colors is completely removed as a command line
    option.
    '''

    def __init__(self, 
             prog=None,
             usage=None,
             description=None,
             epilog=None,
             version=None,
             parents=[],
             formatter_class=argparse.ArgumentDefaultsHelpFormatter,
             prefix_chars='-',
             fromfile_prefix_chars=None,
             argument_default=None,
             conflict_handler='error',
             add_help=True,
             **kwargs):
        #the default named kwargs above are just those from the normal ArgumentParser, 
        #except for the formatter_class.  These are the defaults which can be overwridden 
        #with keyword arguments passed to the constructor.
        
        #markers "x+|_" were not working in the defaultMarkers, for some reason
        option_defaults = {
                'defaultSubplotKwargs': [],
                #plot defaults
                'defaultGreys': ['0.0'],
                'defaultColors': None,
                'defaultMarkerSizes': [12.0],
                'defaultMarkers': 'Dos*^p',
                'defaultLineWidth': [3.0],
                'defaultLineStyle': ['-', '--', ':'],
                'defaultCycleStylesWithinFiles': True,
                'defaultPlotKwargs': [],
                        
                #axis defaults
                'defaultYrange': None,
                'defaultXrange': None,
                'defaultAxisTickFontsize': 16.0,
                'defaultAxisLabelFontsize': 32.0,
                        
                'defaultXTickFunction': None,
                'defaultTickKwargs': [],
                'defaultXTickKwargs': [],
                'defaultYTickKwargs': [],
                        
                'defaultXTickLabelKwargs': ['weight=bold'],
                'defaultYTickLabelKwargs': ['weight=bold'],
                'defaultXtick_label_location': 'm',
                'defaultYtick_label_location': 'm',
                'defaultDrawAxesAtZero': None,
                'defaultDrawVerticalGrid': False,
                'defaultDrawHorizontalGrid': False,
                        
                #axis labels
                'defaultXlabel': None,
                'defaultYlabel': None,
                'defaultXlabelLocation': 's',
                'defaultYlabelLocation': 's',
                'defaultXlabelKwargs': [],
                'defaultYlabelKwargs': [],
                        
                #title defaults
                'defaultTitleFunc': None,
                'defaultTitles': None,
                'defaultTitleFontsize': 18.0,
                'defaultSuperTitle': None,
                'defaultSuperTitleFontsize': 32.0,
                'defaultTitleHalign': 'center',
                'defaultTitleValign': 'top',
                'defaultTitleKwargs': [],
                'defaultSuperTitleKwargs': [],
                        
                'defaultDZaxKwargs': ['leftw=4', 'lefts=solid', 'rightw=4', 'rights=solid', 
                        'topw=4', 'tops=solid', 'bottomw=4', 'bottoms=solid'],

                'defaultNcol': 2,
                'defaultOnePlot': None,
                'defaultSubplotPerFunction': None,
                'defaultInput': True,
                'defaultOutput': True,

                'defaultMissingOk': None,
                'defaultDataColumnFunc': None,

                'defaultHatches': None,

                'defaultSkipRows': 0,
                'defaultRowIgnorePatterns': [],

                'defaultSeriesNames': None,

                'defaultHistogram': None,
                'defaultHistogramKwargs': [],
                'defaultHistogramBins': 20,
                       
                'defaultLegendTextKwargs': [],
                        
                'defaultBarGraph': None,
                'defaultTkGui': False
                }

        #this just pulls out and removes any kwargs matching the keys for the defaults defined above, 
        #and then passes any remaining on to the base __init__
        reduced_kwargs = dict(kwargs)
        for key, val in kwargs.iteritems():
            if key[0:7] == "default":
                if key not in option_defaults:
                    sys.exit('unknown option %s passed to PlottingArgumentParser.__init__' % key)
                else:
                    if val == 'SUPPRESS':
                        option_defaults.pop(key)
                    else:
                        option_defaults[key] = val
                    reduced_kwargs.pop(key)

        #since the default kwargs are now explictly included in the __init__ and explictly passed 
        #to the ArgumentParser __init__, reduced_kwargs should actually be empty here

        super(PlottingArgumentParser, self).__init__(prog=prog, usage=usage, description=description, epilog=epilog, version=version, parents=parents, formatter_class=formatter_class, prefix_chars=prefix_chars, fromfile_prefix_chars=fromfile_prefix_chars, argument_default=argument_default, conflict_handler=conflict_handler, add_help=add_help, **reduced_kwargs)
        
        if 'defaultSubplotKwargs' in option_defaults:
            self.add_argument('--subplot-kwargs', nargs='*', type=str, default=option_defaults['defaultSubplotKwargs'], action=ArgparseActionAppendToDefault, 
                    help='kwargs to be passed on to matplotlib Figure.subplots_adjust function')
        
        if 'defaultDZaxKwargs' in option_defaults:
            self.add_argument('--dz-axis-kwargs', nargs='*', type=str, default=option_defaults['defaultDZaxKwargs'], action=ArgparseActionAppendToDefault,
                    help='settings for axis appearance - interpreted by ME not MPL. w=weight, s=style, examples: rights=bold topw=5.')
        ################
        #argument for plotting of lines/points
        styleArgs = self.add_argument_group('ARGUMENTS FOR POINT AND LINE STYLES')

        if 'defaultGreys' in option_defaults:
            styleArgs.add_argument('-gv', '--grey-values', nargs='*', type=str, default=option_defaults['defaultGreys'],
                                help='floating point values in range 0.0 - 1.0 (black to white) indicating grey value cycle of lines')

        if 'defaultColors' in option_defaults:
            styleArgs.add_argument('-cv', '--color-values', nargs='*', type=str, default=option_defaults['defaultColors'],
                                help='valid matplotlib color specs indicating color cycle of lines')
        
        if 'defaultMarkers' in option_defaults:
            styleArgs.add_argument('-m', '--markers', type=str, default=option_defaults['defaultMarkers'],
                                help='single string with marker designations that will be cycled through for series')

        if 'defaultMarkerSizes' in option_defaults:
            styleArgs.add_argument('-ms', '--marker-sizes', nargs='*', type=float, default=option_defaults['defaultMarkerSizes'],
                                help='floating point values indicating size of markers in points for series 1 2 3, or a single value to apply to all')

        if 'defaultLineWidth' in option_defaults:
            styleArgs.add_argument('-lw', '--line-width', nargs='*', type=float, default=option_defaults['defaultLineWidth'],
                                help='point size of lines to cycle through')
 
        if 'defaultLineStyle' in option_defaults:
            styleArgs.add_argument('-ls', '--line-style', nargs='*', type=str, default=option_defaults['defaultLineStyle'],
                                help='styles of lines to cycle through')
 
        if 'defaultCycleStylesWithinFiles' in option_defaults:
            self.add_argument('--style-per-file', default=option_defaults['defaultCycleStylesWithinFiles'], action='store_false',
                                help='use the same line/point styles for all series within a given file, rather \
                                than cycling through styles for each series _within_ a file')

        if 'defaultPlotKwargs' in option_defaults:
            self.add_argument('--plot-kwargs', nargs='*', type=str, default=option_defaults['defaultPlotKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to matplotlib axes.plot function')
        
        ############
        #axis labels
        axisLabelArgs = self.add_argument_group('ARGUMENTS FOR AXIS LABELING')
        if 'defaultXlabel' in option_defaults:
            axisLabelArgs.add_argument('-xl', '--x-label', type=str, default=option_defaults['defaultXlabel'],
                                help='label for x axis')

            if 'defaultXlabelLocation' in option_defaults:
                axisLabelArgs.add_argument('-xll', '--x-label-location', type=str, default=option_defaults['defaultXlabelLocation'], choices=['s', 'm', 'n'],
                                    help='where xlabels should appear in multiplot: s (single centered), m (one per plot on margins), n (none)')
        
            if 'defaultXlabelKwargs' in option_defaults:
                axisLabelArgs.add_argument('--x-label-kwargs', nargs='*', type=str, default=option_defaults['defaultXlabelKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to pyplot.xlabel function')
        
        if 'defaultYlabel' in option_defaults:
            axisLabelArgs.add_argument('-yl', '--y-label', type=str, default=option_defaults['defaultYlabel'],
                                help='label for y axis')

            if 'defaultYlabelLocation' in option_defaults:
                axisLabelArgs.add_argument('-yll', '--y-label-location', type=str, default=option_defaults['defaultYlabelLocation'], choices=['s', 'm', 'n'],
                                    help='where ylabels should appear in multiplot: s (single centered), m (one per plot on margins), n (none)')
        
            if 'defaultYlabelKwargs' in option_defaults:
                axisLabelArgs.add_argument('--y-label-kwargs', nargs='*', type=str, default=option_defaults['defaultYlabelKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to pyplot.ylabel function')
        
        if 'defaultXlabel' in option_defaults or 'defaultYlabel' in option_defaults:
            if 'defaultAxisLabelFontsize' in option_defaults:
                axisLabelArgs.add_argument('-als', '--axis-label-fontsize', type=float, default=option_defaults['defaultAxisLabelFontsize'],
                                    help='font size of axis labels')
            
        ##########
        #axis range and tics
        axisArgs = self.add_argument_group('ARGUMENTS FOR AXIS RANGES AND TICKS')
        
        if 'defaultYrange' in option_defaults:
            axisArgs.add_argument('-yr', '--y-range', nargs='*', type=float, default=option_defaults['defaultYrange'],
                                help='min and max values on y axis pass one min max pair to have it shared, or a series of min max min max for each subplot (default is auto determined by pyplot from the data)')

        if 'defaultXrange' in option_defaults:
            axisArgs.add_argument('-xr', '--x-range', nargs='*', type=float, default=option_defaults['defaultXrange'],
                                help='min and max values on x axis pass one min max pair to have it shared, or a series of min max min max for each subplot (default is auto determined by pyplot from the data)')

        if 'defaultAxisTickFontsize' in option_defaults:
            axisArgs.add_argument('-ats', '--axis-tick-fontsize', type=float, default=option_defaults['defaultAxisTickFontsize'],
                                help='font size of graph tick labels')

        if 'defaultXTickFunction' in option_defaults:
            axisArgs.add_argument('--x-tick-function', type=str, default=option_defaults['defaultXTickFunction'],
                                help='optional lambda or named function to get x-axis tick labels from input file')

        if 'defaultTickKwargs' in option_defaults:
            axisArgs.add_argument('--tick-kwargs', nargs='*', type=str, default=option_defaults['defaultTickKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to matplotlib axes.tick_params function, applied to both axes')
        
        if 'defaultXTickKwargs' in option_defaults:
            axisArgs.add_argument('--x-tick-kwargs', nargs='*', type=str, default=option_defaults['defaultXTickKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to matplotlib axes.tick_params function, applied to only x axis')
        
        if 'defaultYTickKwargs' in option_defaults:
            axisArgs.add_argument('--y-tick-kwargs', nargs='*', type=str, default=option_defaults['defaultYTickKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to matplotlib axes.tick_params function, applied to only y axis')
        
        if 'defaultXTickLabelKwargs' in option_defaults:
            axisArgs.add_argument('--x-tick-label-kwargs', nargs='*', type=str, default=option_defaults['defaultXTickLabelKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to axes.set_xticklabels')

        if 'defaultYTickLabelKwargs' in option_defaults:
            axisArgs.add_argument('--y-tick-label-kwargs', nargs='*', type=str, default=option_defaults['defaultYTickLabelKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to axes.set_yticklabels')

        if 'defaultDrawAxesAtZero' in option_defaults:
            axisArgs.add_argument('--draw-axes-at-zero', default=option_defaults['defaultDrawAxesAtZero'], action='store_true', 
                                help='whether to draw horizontal and vertical lines through 0, 0')
       
        if 'defaultDrawVerticalGrid' in option_defaults:
            axisArgs.add_argument('--draw-vertical-grid', default=option_defaults['defaultDrawVerticalGrid'], action='store_true', 
                                help='whether to draw vertical lines from top to bottom axis at major tick points')
       
        if 'defaultDrawHorizontalGrid' in option_defaults:
            axisArgs.add_argument('--draw-horizontal-grid', default=option_defaults['defaultDrawHorizontalGrid'], action='store_true', 
                                help='whether to draw horizontal lines from left to right axis at major tick points')
       
        if 'defaultXtick_label_location' in option_defaults:
            axisArgs.add_argument('--x-tick-label-location', default=option_defaults['defaultXtick_label_location'], choices=['m', 'e', 'n'],
                    help='whether to draw X tick labels only on (m)arginal plots, (n)one, or on (e)very plot.')
       
        if 'defaultYtick_label_location' in option_defaults:
            axisArgs.add_argument('--y-tick-label-location', default=option_defaults['defaultYtick_label_location'], choices=['m', 'e', 'n'],
                    help='whether to draw Y tick labels only on (m)arginal plots, (n)one, or on (e)very plot.')

        #########
        #titles
        titleArgs = self.add_argument_group('ARGUMENTS FOR PLOT AND SUBPLOT TITLES')
        if 'defaultTitleFunc' in option_defaults or 'defaultTitle' in option_defaults:
            titleType = titleArgs.add_mutually_exclusive_group()
            titleType.add_argument('--titles', nargs='*', type=str, default=option_defaults['defaultTitles'],
                                help='plot titles, HACKY!, must supply one per SERIES, although if multiple series per plot then later will supercede earlier')
            
            titleType.add_argument('--title-func', type=str, default=option_defaults['defaultTitleFunc'],
                                help='function used to map arbitrary strings (datafile path names) to plot names')
            
            if 'defaultTitleHalign' in option_defaults:
                titleArgs.add_argument('-tha', '--title-horiz-align', type=str, default=option_defaults['defaultTitleHalign'], choices=['left', 'right', 'center'],
                                    help='horizontal alignment of title')

            if 'defaultTitleValign' in option_defaults:
                titleArgs.add_argument('-tva', '--title-vert-align', type=str, default=option_defaults['defaultTitleValign'], choices=['top', 'bottom', 'center', 'baseline'],
                                    help='vertical alignment of title')
        
            if 'defaultTitleKwargs' in option_defaults:
                titleArgs.add_argument('--title-kwargs', nargs='*', type=str, default=option_defaults['defaultTitleKwargs'], action=ArgparseActionAppendToDefault,
                                    help='kwargs to be passed on to matplotlib axes.setTitle function')

            if 'defaultTitleFontsize' in option_defaults:
                titleArgs.add_argument('-ts', '--title-fontsize', type=float, default=option_defaults['defaultTitleFontsize'],
                                    help='font size of title')

        if 'defaultSuperTitle' in option_defaults:
            titleArgs.add_argument('-st', '--super-title', type=str, default=option_defaults['defaultSuperTitle'],
                                help='plot super title on multiplot')

            if 'defaultSuperTitleFontsize' in option_defaults:
                titleArgs.add_argument('-stf', '--super-title-fontsize', type=float, default=option_defaults['defaultSuperTitleFontsize'],
                                    help='font size of super title of multiplot')

            if 'defaultSuperTitleKwargs' in option_defaults:
                titleArgs.add_argument('--super-title-kwargs', nargs='*', type=str, default=option_defaults['defaultSuperTitleKwargs'], action=ArgparseActionAppendToDefault,
                                    help='kwargs to be passed on to pyplot.suptitle function')
       
        ##############
        dataArgs = self.add_argument_group('ARGUMENTS FOR PARSING AND MANIPULATING DATA')
        
        if 'defaultInput' in option_defaults:
            #dataArgs.add_argument('-i', '--input', dest='inFiles', nargs='*', type=FileType('rb'), default=None, 
            #                    help='Series of intput files with data to be plotted')
            dataArgs.add_argument('-i', '--input', dest='inFiles', nargs='*', type=str, default=None, 
                                help='Series of intput files with data to be plotted')

        if 'defaultDataColumnFunc' in option_defaults:
            dataArgs.add_argument('--data-column-func', nargs='*', type=str, default=option_defaults['defaultDataColumnFunc'],
                                help='functions to select or convert column(s) in datafiles to a plotable series for pyplot.plot. \
                                To plot a single series of only x values, something like \'lambda rows: ([float(r[1]) - float(r[3]) for r in rows])\'\
                                will work.  To plot x-y points, \'lambda rows: ([float(r[1]) for r in rows], [float(r[3])/float(r[2]) for r in rows])\'\
                                Multiple functions can be used, in which case each is applied to each datafile.')
        
        if 'defaultSkipRows' in option_defaults:
            dataArgs.add_argument('--skip-rows', type=int, default=option_defaults['defaultSkipRows'],
                                help='number of rows to skip at start of data file')
        
        if 'defaultRowIgnorePatterns' in option_defaults:
            dataArgs.add_argument('--row-ignore-patterns', nargs='*', type=str, default=option_defaults['defaultRowIgnorePatterns'],
                                help='lines in the data file matching these regex patterns are ignored')
        
        if 'defaultMissingOk' in option_defaults:
            dataArgs.add_argument('--missing-ok', action='store_true', default=False, 
                                help='allow some input files to be missing')

        ###############
        plotType = self.add_mutually_exclusive_group()
        if 'defaultSubplotPerFunction' in option_defaults:
            plotType.add_argument('--subplot-per-function', default=False, action='store_true',
                                help='plot each data column function on a new subplot, rather than using one subplot per _file_')
        
        if 'defaultOnePlot' in option_defaults:
            plotType.add_argument('--one-plot', default=False, action='store_true',
                                help='plot all data from all files and data functions on one axis, rather than using one file per subplot')
        
        if 'defaultNcol' in option_defaults:        
            self.add_argument('-c', '--columns', type=int, default=option_defaults['defaultNcol'],
                                help='number of figures to place per row in a multiplot')

        if 'defaultOutput' in option_defaults:
            #self.add_argument('-o', '--output', dest='outFile', type=FileType('w'), default=sys.stdout, 
            #                    help='file to write output to (default stdout or display plot to screen)')
            self.add_argument('-o', '--output', dest='outFile', type=str, default=None, 
                                help='file to write output to (default stdout or display plot to screen)')

        #this is specialized and not generally useful (but used for pies)
        if 'defaultHatches' in option_defaults:
            self.add_argument('-p', '--patterns', nargs='*', type=str, default=option_defaults['defaultHatches'], 
                                help='string with four characters representing patterns (hatching) for the four wedges \
                                out of \"/ \\ | - + x X o O . *\" and blank')

        if 'defaultSeriesNames' in option_defaults:
            self.add_argument('--series-names', nargs='*', type=str, default=None, 
                                help='names to give each plotted series, which will trigger the creation of a legend')

        if 'defaultLegendTextKwargs' in option_defaults:
            self.add_argument('--legend-text-kwargs', nargs='*', type=str, default=option_defaults['defaultLegendTextKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to the legend function, if --series-names are specified')

        if 'defaultHistogram' in option_defaults:
            histogramArgs = self.add_argument_group('ARGUMENTS FOR PLOTTING OF HISTOGRAMS')

            histogramArgs.add_argument('--histogram', default=False, action='store_true',
                                help='plot a histogram rather than line plot, disregarding or redefining some options (difference from --bar-graph is that histogram automatically bins data)')

            histogramArgs.add_argument('--histogram-bins', default=option_defaults['defaultHistogramBins'], type=int,
                                help='the number of bins to calculate, if this is a histogram plot')

            if 'defaultHistogramKwargs' in option_defaults:
                histogramArgs.add_argument('--histogram-kwargs', nargs='*', type=str, default=option_defaults['defaultHistogramKwargs'], action=ArgparseActionAppendToDefault,
                                help='kwargs to be passed on to hist function, if --histogram is used')

            histogramArgs.add_argument('--normalize-histogram', default=False, action='store_true',
                                help='normalilze histogram series, so that multiple series are comparable')

        if 'defaultBarGraph' in option_defaults:
            barGraphArgs = self.add_argument_group('ARGUMENTS FOR PLOTTING BAR GRAPHS')
            
            barGraphArgs.add_argument('--bar-graph', default=False, action='store_true',
                                help='plot a bar graph rather than line plot, disregarding or redefining some options (difference from --histogram is that histogram automatically bins data)')

        if 'defaultTkGui' in option_defaults:
            self.add_argument('--gui', dest='useGui', default=False, action='store_true',
                                help='Start a tkinter-based GUI for plotting - experimental, and may not work')


def make_figure_and_subplots(nrows, ncols, sharex=True, sharey=True):
    '''Not overly helpful wrapper to the pyplot.subplots function'''
    
    #squeeze=False ensures that the axes are returned in a 2-d array, even if only 1x1, which  
    #standardizes things
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey, squeeze=False)
    axes = flatten_array(axes)

    #if don't do this, aspect ratio of subplots won't be same as a single plot would be, 
    #but this also means that fonts and points will be effectively scaled down 
    oldSize = fig.get_size_inches()
    fig.set_size_inches(oldSize[0] * ncols, oldSize[1] * nrows)
    return fig, axes


def prepare_all_kwargs(opt):
    #allKwargDict = collections.defaultdict(lambda: {})
    allKwargDict = {}
    
    allKwargDict['superTitleKwargs'] = prepare_plot_kwargs(opt.super_title_kwargs,
        [('fontsize', opt.super_title_fontsize)])
    
    allKwargDict['subplotKwargs'] = prepare_plot_kwargs(opt.subplot_kwargs, [])
    
    allKwargDict['xLabelKwargs'], allKwargDict['yLabelKwargs'] = [ prepare_plot_kwargs(kw, [('size', opt.axis_label_fontsize)]) for kw in [opt.x_label_kwargs, opt.y_label_kwargs ] ]
    
    #various kwargs
    allKwargDict['titleKwargs'] = prepare_plot_kwargs(opt.title_kwargs, 
            [('fontsize', opt.title_fontsize), ('horizontalalignment', opt.title_horiz_align), ('verticalalignment', opt.title_vert_align)])

    allKwargDict['tickKwargs'] = prepare_plot_kwargs(opt.tick_kwargs,
            [('labelsize', opt.axis_tick_fontsize)])
    
    #these can only be combined like this because no specific values are being passed in
    others = [ ('xTickKwargs', 'x_tick_kwargs'), ('yTickKwargs', 'y_tick_kwargs'), ('xTickLabelKwargs', 'x_tick_label_kwargs' ), ('yTickLabelKwargs', 'y_tick_label_kwargs'), ('plotKwargs', 'plot_kwargs'), ('histogramKwargs', 'histogram_kwargs'), ('legendTextKwargs', 'legend_text_kwargs') ]
    for var, option in others:
        if hasattr(opt, option):
            allKwargDict[var] = prepare_plot_kwargs(opt.__dict__[option], []) 
        else:
            allKwargDict[var] = []
    
    return allKwargDict
        

def full_plot_routine(opt, fileData):
    '''This does a whole bunch of stuff related to making subplots, plotting data, cycling through styles,
    etc.  It is directed by the options returned from the parse_args call to a PlottingArgumentParser.
    
    There are a number of ways to line up plots/files/series, etc.
    A=(axis,alternatively called subplot), F=file, S=series(i.e., data function)
    1. Everything goes onto one axis. (opt.onePlot = True)
        A   A   A   A
        F   F   F2  F2
        S   S2  S   S2
    2. Everything goes on its own axis. (opt.subplotPerFunction = True)
        A1  A2  A3  A4
        F   F   F2  F2
        S   S2  S   S2
    3. Each file has its own axis, possibly with multiple series (default)
        A1  A1  A2  A2
        F   F   F2  F2
        S   S2  S   S2
 
    A number of important quantities:
        These don't depend on which of the options above are used:
            numFiles  =  what it sounds like 
            numFuncs  =  number of data funcs _per file_
            numSeries =  numFiles * numFuncs
        
        These do depend on options 1-3 
            numUsedSubplots = the actual number of axes that series are put on
                Option 1: numUsedSubplots = 1
                Option 2: numUsedSubplots = numSeries
                Option 3: numUsedSubplots = numFiles
            
            totalSubplots = may be greater than numUsed subplots, for example
                if a 2 x 2 matrix of subplots are created but only 3 are used.
                len(axes) will be equal to this.

    '''

    numFiles = len(fileData)
    numFuncs = len(opt.data_column_func)
    numSeries = numFiles * numFuncs

    if opt.subplot_per_function:
        numUsedSubplots = numSeries
    elif opt.one_plot:
        numUsedSubplots = 1
    else:
        numUsedSubplots = numFiles

    #################
    #now the plotting

    #make a multiplot (i.e., call pyplot.subplots) even if there is only one plot, to standardize below code
    multiplot = True
    nrows = int(ceil(numUsedSubplots / float(opt.columns)))
    ncols = min(opt.columns, numUsedSubplots)
    totalSubplots = nrows * ncols
    sharex = sharey = True
    if opt.x_tick_label_location == 'e':
        sharex = False
    if opt.y_tick_label_location == 'e':
        sharey = False
    fig, axes = make_figure_and_subplots(nrows, ncols, sharex, sharey)

    assert(len(axes) == totalSubplots)

    allKwargDict = prepare_all_kwargs(opt)

    fig.subplots_adjust(**allKwargDict['subplotKwargs'])
    
    if opt.super_title:
        plt.suptitle(opt.super_title, **allKwargDict['superTitleKwargs'])

    #if this isn't explicitly set it is figured out by pyplot and applied equally to all subplots
    if opt.y_range:
        plt.ylim(opt.y_range)

    if opt.x_range:
        plt.xlim(opt.x_range)

    #make a single centered label spanning individual plots
    if opt.x_label and opt.x_label_location == 's':
        fig.text(allKwargDict['xLabelKwargs'].pop('x'), allKwargDict['xLabelKwargs'].pop('y'), opt.x_label, **allKwargDict['xLabelKwargs'])
    if opt.y_label and opt.y_label_location == 's':
        fig.text(allKwargDict['yLabelKwargs'].pop('x'), allKwargDict['yLabelKwargs'].pop('y'), opt.y_label, **allKwargDict['yLabelKwargs'])

    #########################
    def make_plot_iterator():
        '''It will be easiest to set up a bunch of lists of length numSeries, and then iterate through
        these to get the proper association between data, functions, axes, etc.
        This is a closure defined within the outer function
        '''
        if opt.subplot_per_function:
            #the # of axes should already be the # of plots, and will match up to d1f1 d1f2 d2f1 d2f2 etc
            #just slice off any axes that won't be used
            axesToIterate = islice(axes, numSeries)
        elif opt.one_plot:
            #axes is a list of length 1, so repeat it numSeries times
            assert(len(axes) == 1)
            axesToIterate = axes * numSeries 
        else:
            #a single axis for each datafile, which will be equivalent to subplotPerFunction  unless there are multiple functions per file 
            axesToIterate = chain.from_iterable([ [ax] * numFuncs for ax in islice(axes, numFiles)])

        #this will be equivalent - the infiles are only used for titles anyway
        dataToIterate = chain.from_iterable( [dat] * numFuncs for dat in fileData )
        inFilesToIterate = chain.from_iterable( [ [inf] * numFuncs for inf in opt.inFiles ] )

        #functions to iterate should just cycle through the functions for each file
        functionsToIterate = cycle_n_elements_m_times( opt.data_column_func, numFuncs, numFiles )

        #style stuff
        if opt.one_plot or opt.subplot_per_function:
            #each series (whether from different files or not) gets a sucessive style within the one plot - thus, generally just keep cycling
            markersToIterate = cycle_n_elements(opt.markers, numSeries)
            markerSizesToIterate = cycle_n_elements(opt.marker_sizes, numSeries)
            lineWidthsToIterate = cycle_n_elements(opt.line_width, numSeries)
            lineStylesToIterate = cycle_n_elements(opt.line_style, numSeries)
            colorValuesToIterate = cycle_n_elements(opt.color_values if opt.color_values else opt.grey_values, numSeries)
            hatchesToIterate = cycle_n_elements(opt.patterns or [''], numSeries)
            seriesNamesToIterate = padded(opt.series_names or [], numSeries, '_nolegend_')

        elif opt.style_per_file:
            #this is the default really, sucessive styles within subplots
            ''' This won't work, I believe
            markersToIterate = markerSizesToIterate = lineWidthsToIterate = lineStylesToIterate = []
            for it, source in [(markersToIterate, opt.markers), (markerSizesToIterate, opt.marker_sizes), 
                    (lineWidthsToIterate, opt.line_width), (lineStylesToIterate, opt.line_style)]:
                it = cycle_n_elements_m_times(source, numFuncs, numFiles)
            '''
            markersToIterate = cycle_n_elements_m_times(opt.markers, numFuncs, numFiles)
            markerSizesToIterate = cycle_n_elements_m_times(opt.marker_sizes, numFuncs, numFiles)
            lineWidthsToIterate = cycle_n_elements_m_times(opt.line_width, numFuncs, numFiles)
            lineStylesToIterate = cycle_n_elements_m_times(opt.line_style, numFuncs, numFiles)
            
            colorValuesToIterate = cycle_n_elements_m_times(opt.color_values or opt.grey_values, numFuncs, numFiles)
            hatchesToIterate = cycle_n_elements_m_times(opt.patterns or [''], numSeries, numFiles)
            seriesNamesToIterate = cycle_n_elements_m_times(padded(opt.series_names or [], numSeries, '_nolegend_'), numSeries, numFiles)

        else:
            #use a single style for all series in a single file
            markersToIterate = chain.from_iterable( [ [mark] * numFuncs for mark in cycle_n_elements(opt.markers, numFiles) ] )
            markerSizesToIterate = chain.from_iterable( [ [size] * numFuncs for size in cycle_n_elements(opt.marker_sizes, numFiles) ] )
            lineWidthsToIterate = chain.from_iterable( [ [width] * numFuncs for width in cycle_n_elements(opt.line_width, numFiles) ] )
            lineStylesToIterate = chain.from_iterable( [ [style] * numFuncs for style in cycle_n_elements(opt.line_style, numFiles) ] )
            colorValuesToIterate = chain.from_iterable( [ [color] * numFuncs for color in cycle_n_elements(opt.color_values or opt.grey_values, numFiles) ] )
            hatchesToIterate = chain.from_iterable( [ [hatch] * numFuncs for hatch in cycle_n_elements(opt.patterns or [''], numFiles) ] )
            seriesNamesToIterate = padded(opt.series_names or [], numSeries, '_nolegend_')  
        
        '''
        print len([ d for d in axesToIterate]), 
        print len([ d for d in dataToIterate]),
        print len([ d for d in functionsToIterate])
        print opt.markers, 
        print [ m for m in markersToIterate ]
        print opt.marker_sizes,
        print [ m for m in markerSizesToIterate ]
        print opt.line_width,
        print [ m for m in lineWidthsToIterate ]
        print opt.line_style,
        print [ m for m in lineStylesToIterate ]
        print opt.grey_values,
        print [ m for m in colorValuesToIterate ]
        print [ m for m in hatchesToIterate ]
        print opt.series_names,
        print [ m for m in seriesNamesToIterate ]
        exit()
        '''

        return izip(
                inFilesToIterate,
                dataToIterate,
                functionsToIterate,
                markersToIterate,
                markerSizesToIterate,
                lineWidthsToIterate,
                lineStylesToIterate,
                colorValuesToIterate,
                hatchesToIterate,
                seriesNamesToIterate,
                axesToIterate)
    #####################

    #cycle through markers, marker sizes, colors, axes, etc., plotting one series per loop
    #the izip object returned by make_plot_iterator will return tuples until the end of the shortest is reached
    #subplot is a matplotlib.axes instance
    for num, (inFile, series, function, marker, markerSize, lineWidth, lineStyle, color, hatch, seriesName, subplot) in enumerate(make_plot_iterator()):
        #turn off extra ticks
        subplot.get_xaxis().tick_bottom()
        subplot.get_yaxis().tick_left()
        
        #if location option is m(arginal), add a normal axis label for plots on left or bottom margin
        if opt.y_label:
            if opt.y_label_location == 'm' or not multiplot:
                if num % ncols == 0:
                    subplot.set_ylabel(opt.y_label, **allKwargDict['yLabelKwargs'])
        if opt.x_label:
            if opt.x_label_location == 'm' or not multiplot:
                if num >= totalSubplots - ncols:
                    subplot.set_xlabel(opt.x_label, **allKwargDict['xLabelKwargs'])

        #set the title - this could get set multiple times for a single file/plot with multiple funcs
        if opt.title_func and opt.title_func != 'None':
            if isinstance(opt.title_func, str):
                subplot.set_title(eval(opt.title_func)(inFile, sep=' '), **allKwargDict['titleKwargs'])
            else:
                subplot.set_title(opt.title_func(inFile, sep=' '), **allKwargDict['titleKwargs'])
        elif opt.titles:
            titleNum = num if len(opt.titles) > 1 else 0
            subplot.set_title(opt.titles[titleNum], **allKwargDict['titleKwargs'])
     
        #this is a little goofy, as it will be re-set for each subplot, despite being shared
        if opt.x_tick_function:
            tickFunc = eval(opt.x_tick_function)
            #note that the tick_params calls below mainly affect the ticks themselves, rather
            #than the tick _labels_.  Although, the tick_params do set the label fontsizes.
            #This pyplot.xticks call sets the location and labels of the ticks
            #but any text property kwargs can also be passed
            xTickData = tickFunc(series)
            if len(xTickData) == 2:
                #this is assuming that the tickData returned from the function is two lists, 
                #one specifying tick locations and the other the labels
                #plt.xticks(*tickData, rotation=90, weight='bold')
                plt.xticks(*xTickData, **allKwargDict['xTickLabelKwargs'])
            else:
                #this is treating the tickData just as labels, to be placed at 0, 1, 2, etc.
                #not that useful
                plt.xticks(range(len(series)), xTickData, **allKwargDict['xTickLabelKwargs'])
                #plt.xticks(range(len(series)), tickData, rotation=90, weight='bold')
                #plt.xticks(range(len(series)), tickFunc(series))
        else:
            subplot.set_xticklabels(subplot.get_xticks(), **allKwargDict['xTickLabelKwargs'])

        ##########################################################################################
        #UPDATE the getting of the tick locs and setting them back as the labels has problems.
        #If the yrange is being auto determined by the data series and changes, or is the default 
        #0-1 and then changes to something else the labels will get totally borked.

        if not opt.y_range:
            sys.exit('Enter an explicit y range for now.  See code')

        #this is all a ridiculous hack.  
        #turns out that the cleanest way to set the xticklabel kwargs is by calling subplot.get_yticks()
        #and passing that to subplot.set_yticklabels.  The problem is that what is returned are float numbers
        #and passing them directly to subplot.set_yticklabels makes the tick labels have a decimal point, even
        #if they are integers.  So, this checks whether the numbers are actually all ints, and if so ensures that
        #that are printed as such.
        def equal_float(a, b):
            return abs(a - b) <= sys.float_info.epsilon

        ytick_locs = subplot.get_yticks()
        for num, loc in enumerate(ytick_locs):
            if not equal_float(loc, float(int(loc))):
                break
        else:
            labels_as_ints = [ str(int(loc)) for loc in ytick_locs ]
            ytick_locs = labels_as_ints

        subplot.set_yticklabels(ytick_locs, **allKwargDict['yTickLabelKwargs'])
        #subplot.set_yticklabels(subplot.get_yticks(), **allKwargDict['yTickLabelKwargs'])
        ###########################################################################################

        #change tick properties
        if allKwargDict['tickKwargs']:
            subplot.tick_params(**allKwargDict['tickKwargs'])
        if allKwargDict['xTickKwargs']:
            subplot.tick_params(axis='x', **allKwargDict['xTickKwargs'])
        if allKwargDict['yTickKwargs']:
            subplot.tick_params(axis='y', **allKwargDict['yTickKwargs'])
        
        #do the actual data evaluation and plotting
        #print series
        dataFunc = eval(function)
        toPlot = dataFunc(series)
        print toPlot
        '''to plot just x values, lambda looks like this:
        lambda rows: [float(r[2]) for r in rows]
        i.e., output is [x1, x2, ...]
        for x,y, like this:
        lambda rows: ([float(r[2]) for r in rows], [float(r[3]) for r in rows])
        but, the plot function wants [x1, x2, ...], [y1, y2, ...] as the first
        two arguments, hence the * to remove the tuple containing the two lists,
        which I think is required for the lambda'''
        #if isinstance(toPlot[0], list):
        if False:
            if (not hasattr(opt, 'histogram') or not opt.histogram) and (not hasattr(opt, 'bar_graph') or not opt.bar_graph):
                subplot.plot(*dataFunc(series), 
                        marker=marker, 
                        markersize=markerSize, 
                        linewidth=lineWidth,
                        linestyle=lineStyle,
                        label=seriesName,
                        color=color,
                        mfc=None,
                        mec=color,
                        **allKwargDict['plotKwargs'])
            else:
                if hasattr(opt, 'histogram') and opt.histogram:
                    bins = linspace(opt.x_range[0] if opt.x_range else 0.0, opt.x_range[1] if opt.x_range else 1.0, opt.histogram_bins)
                    outN, outBins, patches = subplot.hist(*dataFunc(series), 
                            bins=bins,
                            facecolor=color,
                            normed=opt.normalize_histogram,
                            **allKwargDict['histogramKwargs'])
                else:
                    patches = subplot.bar(*dataFunc(series), 
                            color=color, label=seriesName,
                            **allKwargDict['histogramKwargs'])
                    if opt.series_names:
                        subplot.legend(frameon=False, prop=font_manager.FontProperties(**allKwargDict['legendTextKwargs']))
                for pat in patches:
                    pat.set_hatch(hatch)

        else:
            if (not hasattr(opt, 'histogram') or not opt.histogram) and (not hasattr(opt, 'bar_graph') or not opt.bar_graph):
                subplot.plot(dataFunc(series), 
                        marker=marker, 
                        markersize=markerSize, 
                        linewidth=lineWidth, 
                        label=seriesName,
                        color=color,
                        **allKwargDict['plotKwargs'])
            else:
                if hasattr(opt, 'histogram') and opt.histogram:
                    bins = linspace(opt.x_range[0] if opt.x_range else 0.0, opt.x_range[1] if opt.x_range else 1.0, opt.histogram_bins)
                    #bins = opt.histogram_bins
                    #rng = (opt.x_range[0] if opt.x_range else 0.0, opt.x_range[1] if opt.x_range else 1.0)
                    
                    outN, outBins, patches = subplot.hist(dataFunc(series), 
                            bins=bins,
                            color=opt.color_values,
                            label=opt.series_names,
                            normed=opt.normalize_histogram,
                            **allKwargDict['histogramKwargs'])

                else:
                    patches = subplot.bar(dataFunc(series),
                            color=color, label=seriesName, 
                            **allKwargDict['histogramKwargs'])
                            
                    if opt.series_names:
                        subplot.legend(frameon=False, fontsize=24)
                    #this used to be unindented one level, but histogram parts don't seem to have a set_hatch
                    for pat in patches:
                        pat.set_hatch(hatch)
        
        #set the width of the axes box (frame) which I thought I could set on the frame rectangle, but that didn't work
        #kwargs here are <edge><w or s>, i.e. rightw or tops
        dz_axis_kwargs = prepare_plot_kwargs(opt.dz_axis_kwargs, [ ])

        for option, val in dz_axis_kwargs.items():
            which = option[:-1]
            if option[-1] == 'w':
                subplot.spines[which].set_linewidth(int(val))
            elif option[-1] == 's':
                subplot.spines[which].set_linestyle(val)
            else:
                sys.exit('unrecognized dzAxisKwarg: %s=%s' % (option, val))
   
        #this plots x=0 and y=0 axis lines
        if opt.draw_axes_at_zero:
            subplot.axhline(xmin=-sys.maxint - 1, xmax=sys.maxint)
            subplot.axvline(ymin=-sys.maxint - 1, ymax=sys.maxint)

        if opt.draw_vertical_grid:
            subplot.xaxis.grid(color='gray', linestyle='dashed')

        if opt.draw_horizontal_grid:
            subplot.yaxis.grid(color='gray', linestyle='dashed')

    if opt.series_names:
        for ax in axes:
            ax.legend()

    #remove any unused axes
    for ax in axes[numUsedSubplots:]:
        ax.set_axis_off()

    #return the current figure, so more manipulations could potentially be done on it
    return plt.gcf()

    ###################
    if opt.outFile:
        plt.savefig(opt.outFile, transparent=True, bbox_inches='tight')
    else:
        plt.show()

if __name__ == "__main__":
    import doctest
    doctest.testmod()


