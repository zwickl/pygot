#!/usr/bin/env python
import sys
import itertools
import collections
import re
from os import path
from math import ceil

try:
    import matplotlib
except ImportError:
    sys.exit('Sorry, matplotlib is necessary to run this plotting script.\nSee http://matplotlib.org/downloads.html')

#this is important, because the backend can default to an interactive Tk based one even
#on clusters where there is no display
matplotlib.use('pdf')

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D

from pygot.plotutils import PlottingArgumentParser, prepare_plot_kwargs, simpler_path_to_plot_title_nuc_default_cds, filename_to_paren_trees
from pygot.utils import flatten_array, find_shortest_unique_leading_substrings


def filename_to_triple_title(s, oryza=False):
    sp = s.split('/')
    sp = sp[-1]
    sp = sp.split('.')
    if len(sp) == 4:
        nameMap = {'sativaj':'japonica', 'sativai':'indica', 'glaberrimam':'glaberrima', 'glaberrimaf':'glaberrima'}
        first = nameMap.get(sp[0].lower(), sp[0])
        if 'punc' in sp[1].lower() and not 'brach' in s.lower():
            sec = nameMap.get(sp[2].lower(), sp[2])
            third = nameMap.get(sp[3].lower(), sp[3])
        else:
            sec = nameMap.get(sp[1].lower(), sp[1])
            third = nameMap.get(sp[2].lower(), sp[2])
        title = '%s-%s-%s' % (first, sec, third) 
    else:
        #if we can't interpret this string, just return it
        title = s

    return title


def filename_to_full_triple_title(s):
    sp = s.split('/')
    sp = sp[-1]
    sp = sp.split('.')
    if len(sp) == 4:
        nameMap = {'sativaj':'O. sativa subsp. japonica', 'sativai':'O. sativa subsp. indica', 'glaberrimam':'O. glaberrima', 'glaberrimaf':'O. glaberrima', 
                'barthii':'O. barthii', 'nivara':'O. nivara', 'rufipogon':'O. rufipogon', 'brachyantha':'O. brachyantha', 'punctata':'O. punctata', 'officinalis':'O. officinalis', 
                'leersia':'L. perrii', 'perrii':'L. perrii', 'meridionalis':'O. meridionalis', 'glumaepatula':'O. glumaepatula', 'longistaminata':'O. longistaminata'}
        first = nameMap.get(sp[0].lower(), sp[0])
        if 'punc' in sp[1].lower() and not 'brach' in s.lower():
            sec = nameMap.get(sp[2].lower(), sp[2])
            third = nameMap.get(sp[3].lower(), sp[3])
        else:
            sec = nameMap.get(sp[1].lower(), sp[1])
            third = nameMap.get(sp[2].lower(), sp[2])
        title = '%s-%s-%s' % (first, sec, third) 
    else:
        #if we can't interpret this string, just return it
        title = s

    return title


def split_equals_to_tuple(s):
    if len(s) == 2:
        return tuple(s)
    try:
        spl = s.split('=')
    except:
        raise
    if len(spl) != 2:
        exit('problem parsing setting')
    return tuple(spl)


if __name__ == "__main__":
    import doctest
    doctest.testmod()


class OryzaSubplot():
    '''This is a class that holds the quartet data for a single set of trees - i.e. a single pie+cumulative plot
    and plots those data    
    data - a filename with data or a list of lines from such a file
    '''
    def __init__(self, data, columnOrder, reverseSlices=False, title=None):
        
        if isinstance(data, str):
            inp = open(data, 'rbU')
        else:
            inp = data
        self.allData = [ x.split() for x in (line.strip() for line in inp) if len(x) and not 'file' in x]

        if not title:
            #self.plotTitle = path_to_plot_title(filename)
            #self.plotTitle = simpler_path_to_plot_title(filename)
            self.plotTitle = simpler_path_to_plot_title_nuc_default_cds(data)
        else:
            self.plotTitle = title

        #get the support values for each of the triples if they have > 50% support, and sort them
        self.plotColumns = []
        for col in columnOrder:
            self.plotColumns.append(sorted([f for f in (float(d[col]) for d in self.allData) if f > 0.5]))

        self.totCounts = [ len(c) for c in self.plotColumns ]
        self.totLoci = len(self.allData)

        self.fracs = [ c / float(self.totLoci) for c in self.totCounts ]

        #this is a simple hack to reverse the order in which the pie slices appear, but leave all other data alone
        #the colors assigned to the pies need to be flipped in plot() too then
        self.reverseSlices = reverseSlices
        if reverseSlices:
            self.fracs.reverse()
        
        self.fracs.append((float(self.totLoci) - sum(self.totCounts)) / self.totLoci)

        #make a single list with x, y, x, y, x, y coords for 3 lines
        self.xyxyxy_list = flatten_array([ [self.plotColumns[col], xrange(self.totCounts[col])] for col in xrange(len(self.plotColumns)) ])
    
    def plot(self, axis, **kwargs):

        plotKwargs = {}
        dzkwargs = {}
        for k, v in kwargs.items():
            if 'dz' in k:
                dzkwargs[k] = v
                del kwargs[k]
            elif k == 'plotKwargs':
                plotKwargs = v
                del kwargs[k]
       
        cumLines = axis.plot(*self.xyxyxy_list, **plotKwargs)
        
        title = axis.set_title(self.plotTitle, x=0.725, y=0.92, fontsize=dzkwargs['dztitlesize'], horizontalalignment=dzkwargs['dzhorizontalalignment'], verticalalignment=dzkwargs['dzverticalalignment'], weight='bold')
        
        if len(kwargs['labels']) < len(self.fracs):
            kwargs['labels'] = [ '' for _ in self.fracs ]
        
        cumXpos, cumYpos, cumXsize, cumYsize = axis.get_position().bounds
        pieSizeRelativeToCumY = 0.45
        #this is empirically determined, and needed to figure gap to top axis
        proportionOfPieAxisActuallyFilled = 0.9
        
        pieRelativeGapToYaxis = 0.14
        pieRelativeGapToXaxis = 0.10
        
        pieYsize = cumYsize * pieSizeRelativeToCumY
        pieXsize = cumXsize * pieSizeRelativeToCumY
        
        pieXloc = cumXpos + (pieRelativeGapToYaxis - 0.1) / dzkwargs['dzncols']
        pieYloc = cumYpos + cumYsize - cumYsize * pieRelativeGapToXaxis - pieYsize * proportionOfPieAxisActuallyFilled
        
        newAxBounds = [pieXloc, pieYloc, pieXsize, pieYsize]
        pieAx = plt.axes(newAxBounds, aspect='equal')
        
        pieAx.set_xticks([])
        pieAx.set_yticks([])
        
        #this is pretty ridiculous, but the colors need to be reversed for the pie and then restored so that 
        #they don't affect other things
        if 'colors' in kwargs and self.reverseSlices:
            origColors = kwargs['colors'] 
            #flip the first three colors if reversal of pie slices has been requested
            kwargs['colors'] = list(reversed(kwargs['colors'][0:3])) + [kwargs['colors'][3]]

        pats, texts, autotexts = pieAx.pie(self.fracs, **kwargs)
        
        if 'colors' in kwargs and self.reverseSlices:
            #restore the original colors
            kwargs['colors'] = origColors

        for t in autotexts:
            t.set_fontsize(dzkwargs['dzlabelsize'])
            t.set_weight('bold')
        return cumLines, pats


#various defaults
#these are for the pies
defaultHatches = [ ' ', ' ', ' ', ' ' ]
#use my derived Argparse class to parse commandline input
parser = PlottingArgumentParser(description='make a cumulative support plot with embedded pie, possibly as a multiplot', 
        defaultGreys=['0.55', '0.0', '0.25', '1.0'], 
        defaultMarkers="xo*",
        defaultMarkerSizes=[12, 12, 20],
        defaultLineWidth=['4'],
        defaultHatches=defaultHatches, 
        defaultXrange=[0.5, 1.0],
        
        defaultYlabel='Cumulative Count', 
        defaultXlabel='Bootstrap Proportion',
        defaultXlabelLocation='m',
        defaultXlabelKwargs=['x=0.5', 'y=0.01', 'horizontalalignment=center', 'weight=bold'],
        defaultYlabelKwargs=['x=0.35', 'y=0.5', 'verticalalignment=center', 'rotation=90', 'weight=bold'],
        
        defaultSuperTitleKwargs=['x=0.01', 'y=1.15', 'horizontalalignment=left', 'style=italic', 'weight=normal'],
        defaultSuperTitleFontsize=52,

        defaultTickKwargs=['top=False', 'right=False', 'pad=15', 'width=4', 'length=10'],
        defaultXTickKwargs=[],
        defaultYTickKwargs=[],
        defaultAxisTickFontsize=32,
        
        defaultTitleValign='top',
        defaultTitleKwargs=[],
        defaultTitleFontsize=28,

        defaultPlotKwargs=['markerfacecolor=None', 'markeredgewidth=1', 'clip_on=False'],
        defaultSubplotKwargs=[],
        defaultDZaxKwargs=['leftw=4', 'lefts=solid', 'rightw=4', 'rights=dashed', 'topw=0', 'tops=solid', 'bottomw=4', 'bottoms=solid'],
        
        defaultDataColumnFunc='SUPPRESS',
        defaultHistogram='SUPPRESS',
        defaultBarGraph='SUPPRESS',
        defaultNcol=4)

pieArgs = parser.add_argument_group(title='arguments for pie diagram portion of plot')

defaultLabels = ['', '', '', '']
defaultValueString = '%1.0f%%'

#add some more specific arguments to parser
pieArgs.add_argument('-l', '--labels', nargs='*', type=str, default=defaultLabels,
                    help='labels used for pie slices')

pieArgs.add_argument('-ld', '--label-distance', type=float, default=0.6,
                    help='distance from the pie center at which the labels appear')

pieArgs.add_argument('-vd', '--value-distance', type=float, default=1.6,
                    help='distance from the pie center at which the percentage values appear')

pieArgs.add_argument('-vf', '--value-format-string', type=str, default=defaultValueString,
                    help='printf style format string for percentage values')

pieArgs.add_argument('-vs', '--value-size', type=float, default=32,
                    help='font size of percentage values')

pieArgs.add_argument('-plw', '--pie-line-width', type=float, default=3,
                    help='width of lines dividing pie sections')

pieArgs.add_argument('--pie-start-angle', type=float, default=None,
                    help='angle offset of start of first pie slice, from 3 O\'clock, counterclockwise')

pieArgs.add_argument('--reverse-slices', default=False, action='store_true',
                    help='whether pie slices, ordered counterclockwise, should be 21M instead of M12')

'''this is now hard coded below.  It is complicated to keep it consistent with diff numbers of subplots
pieArgs.add_argument('-px', '--pie-x-loc', dest='pieXloc', type=float, default=None,
                    help='relative X location of pie within cumulative plot. Generally automatically determined, but may need to be adjusted when number of subplots change')

pieArgs.add_argument('-py', '--pie-y-loc', dest='pieYloc', type=float, default=None,
                    help='relative Y location of pie within cumulative plot. Generally automatically determined, but may need to be adjusted when number of subplots change')
'''

parser.add_argument('--order', nargs='*', type=int, default=[1, 2, 3], 
                    help='order of quartet columns in input file to plot')

parser.add_argument('--letter-labels', type=str, default='M12U', 
                    help='Letter characters to use for majority, two minority and unresolved resolutions')

parser.add_argument('--paren-labels', nargs='*', type=str, default=[],
                    help='parenthetical notation for three triple resolutions')

parser.add_argument('--oryza-names', action='store_true', default=False,
                    help='Do some name manipulations specific to Oryza data')

parser.add_argument('--no-legend-bars', action='store_true', default=False,
                    help='do not include a legend on the left of the figure indicating the color to triplet correspondence')

#now use the tkinter gui, which may not work, or process the command line
use_tk_gui = False
if use_tk_gui:
    from Tkinter import *
    from tkinterutils import *
    from ttk import *
    root = Tk()
    gui = ArgparseGui(root, parser, height=1000, width=1600)
    
    root.wait_window(gui.frame)
    if gui.cancelled:
        sys.exit('cancelled ...')
    args = gui.make_commandline_list()
    print args

    '''
    args = ['--grey-values', '0.55', '0.0', '0.25', '1.0', '--markers', 'xo*', '--marker-sizes', '6.0', '--line-width', '3.0', '--line-style', '-', '--x-label', 'Bootstrap Proportion', '--x-label-location', 'm', '--y-label', 'Cumulative Count', '--y-label-location', 's', '--axis-label-fontsize', '32.0', '--x-range', '0.5', '1.0', '--axis-tick-fontsize', '16.0', '--draw-axes-at-zero', '--x-tick-label-location', 'm', '--y-tick-label-location', 'm', '--title-horiz-align', 'center', '--title-vert-align', 'top', '--title-fontsize', '18.0', '--super-title-fontsize', '32.0', '--labels', '--label-distance', '0.6', '--value-distance', '0.6', '--value-format-string', '%1.1f%%', '--value-size', '24', '--pie-line-width', '3', '--dz-axis-kwargs', 'leftw=4', 'lefts=solid', 'rightw=4', 'rights=dashed', 'topw=0', 'tops=solid', 'bottomw=4', 'bottoms=solid', '--styles-per-file', '--columns', '2', '--patterns', 'XXX', 'o', '+', '--order', '1', '2', '3']
    args = ['--grey-values', '0.55', '0.0', '0.25', '1.0', '--markers', 'xo*', '--marker-sizes', '6.0', '--line-width', '3.0', '--line-style', '-', '--x-label', 'Bootstrap Proportion', '--x-label-location', 'm', '--y-label', 'Cumulative Count', '--y-label-location', 's', '--axis-label-fontsize', '32.0', '--x-range', '0.5', '1.0', '--axis-tick-fontsize', '16.0', '--draw-axes-at-zero', '--x-tick-label-location', 'm', '--y-tick-label-location', 'm', '--title-horiz-align', 'c', '--title-vert-align', 't', '--title-fontsize', '18.0', '--super-title-fontsize', '32.0', '--labels', '--label-distance', '0.6', '--value-distance', '0.6', '--value-format-string', '%1.1f%%', '--value-size', '24', '--pie-line-width', '3', '--dz-axis-kwargs', 'leftw=4', 'lefts=solid', 'rightw=4', 'rights=dashed', 'topw=0', 'tops=solid', 'bottomw=4', 'bottoms=solid', '--styles-per-file', '--columns', '2', '--patterns', 'XXX', 'o', '+', '--order', '1', '2', '3']
    print args
    '''
    options = parser.parse_args(args)
    print options
else:
    options = parser.parse_args()

outFile = options.outFile or sys.stdout

actualColors = options.color_values if options.color_values else options.grey_values

figureData = []
if not options.inFiles:
    sys.exit('must enter some input files with -i')
for num, f in enumerate(options.inFiles):
    if isinstance(f, str) and not path.exists(f):
        if options.missing_ok:
            sys.stderr.write('WARNING: input file %s missing!\n' % f)
        else:
            sys.exit('ERROR: input file %s missing!\n' % f)
    else:
        sys.stderr.write('Reading file %s ...\n' % f)
        if options.titles:
            figureData.append(OryzaSubplot(f, options.order, reverseSlices=options.reverse_slices, title=options.titles[num]))
        else:
            figureData.append(OryzaSubplot(f, options.order, reverseSlices=options.reverse_slices))

numPlots = len(figureData)

#################
#now the plotting

#don't make obviously empty plots
options.columns = min(options.columns, numPlots)

#set up the individual plots and axes
sys.stderr.write('Plotting ...\n')
nrows = int(ceil(numPlots / float(options.columns)))
ncols = options.columns
sharex = True if not options.x_range or len(options.x_range) <= 2 else False
sharey = True if not options.y_range or len(options.y_range) <= 2 else False
multiplot = True

#DEBUG - trying to change rescaling
individualPlotWidth = 6.35
aspectRatio = 3 / 2.0
individualPlotSize = (individualPlotWidth, individualPlotWidth / aspectRatio)

#figsize here is the axis-to-axis dimensions of the cumulative plots.  Changing it won't have 
#any effect on overall figure size, I think.
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey, figsize=individualPlotSize)

#resize the overall figure so that things don't get squished (otherwise aspect ratio stays same regardless of how many cols and rows)
#default size is apparently 8" x 6"
origInchWidth, origInchHeight = fig.get_size_inches()

if options.no_legend_bars:
    barRegionWidthRelativeToSinglePlotWidth = 0.0
    barRegionInchWidth = 0.0
else:
    barRegionWidthRelativeToSinglePlotWidth = 0.4724409
    barRegionInchWidth = individualPlotWidth * barRegionWidthRelativeToSinglePlotWidth

leftInchSpace = rightInchSpace = bottomInchSpace = topInchSpace = 1.0
interplotInchWidth = interplotInchHeight = origInchWidth + 0.1

#make each actual plot exactly the size of the orig figure, so the padding goes above that
finalInchWidth = origInchWidth * ncols + interplotInchWidth * (ncols - 1) + barRegionInchWidth + leftInchSpace + rightInchSpace
finalInchHeight = origInchHeight * nrows + interplotInchHeight * (nrows - 1) + bottomInchSpace + topInchSpace

#these are proportional to the original dimensions
#now calculate what the borders are proportionally
#these override most of the subplot kwargs that might be specified anywhere else
subplotKwargs = prepare_plot_kwargs(options.subplot_kwargs, [])

allow_plot_kwarg_override = True
if allow_plot_kwarg_override:
    subplotKwargs['bottom'] = subplotKwargs.get('bottom', bottomInchSpace / finalInchHeight)
    subplotKwargs['top'] = subplotKwargs.get('top', (finalInchHeight - topInchSpace) / finalInchHeight)
    subplotKwargs['left'] = subplotKwargs.get('left', (barRegionInchWidth + leftInchSpace) / finalInchWidth)
    subplotKwargs['right'] = subplotKwargs.get('right', (finalInchWidth - rightInchSpace) / finalInchWidth)

    subplotKwargs['wspace'] = subplotKwargs.get('wspace', interplotInchWidth / origInchWidth)
    subplotKwargs['hspace'] = subplotKwargs.get('hspace', interplotInchHeight / origInchHeight)
else:
    subplotKwargs['bottom'] = bottomInchSpace / finalInchHeight
    subplotKwargs['top'] = (finalInchHeight - topInchSpace) / finalInchHeight
    subplotKwargs['left'] = (barRegionInchWidth + leftInchSpace) / finalInchWidth
    subplotKwargs['right'] = (finalInchWidth - rightInchSpace) / finalInchWidth

    subplotKwargs['wspace'] = interplotInchWidth / origInchWidth
    subplotKwargs['hspace'] = interplotInchHeight / origInchHeight

barRegionProportionalWidth = barRegionInchWidth / finalInchWidth
fig.set_size_inches(finalInchWidth, finalInchHeight)
#adjust the location of the plots for the bar region to the left
#this will have PROPORTIONALLY specified the spacing between plots
fig.subplots_adjust(**subplotKwargs)

#if only one subplot ax is a single axes object, otherwise a tuple of them
#rather annoying
#if rows and col > 1, ax will be a 2d array, so flatten
ax = flatten_array(ax) if isinstance(ax, collections.Iterable) else [ax]

#if this isn't explicitly set it is figured out by pyplot and applied
#equally to all subplots
#if a number of ranges are passed in or sharex/sharey are not set then the ranges may be changed below
if options.y_range and len(options.y_range) == 2:
    plt.ylim((options.y_range[0], options.y_range[1]))

if options.x_range and len(options.x_range) == 2:
    plt.xlim((options.x_range[0], options.x_range[1]))
else:
    plt.xlim(0.5, 1.0)

#SUPERTITLE WITH TAXA IN TRIPLE AND NUM LOCI
if options.super_title:
    super_titleKwargs = prepare_plot_kwargs(options.super_title_kwargs, [('fontsize', options.super_title_fontsize)])
    #CAN EITHER DO INITIAL SHORTHAND FOR TAXON NAMES, OR PROPER FULL NAMES
    #supTitle = filename_to_triple_title(options.super_title)
    supTitle = filename_to_full_triple_title(options.super_title)
    includeNumLociInSupertitle = True
    if includeNumLociInSupertitle:
        supTitle += ' (%d loci)' % figureData[0].totLoci
    plt.suptitle(supTitle, **super_titleKwargs)

#just a mapping from various strings acceptable to me to what matplotlib funcs will take
alignments = { loc:loc for loc in [ 'center', 'left', 'right', 'top', 'bottom' ] }
alignments.update({ loc[0]:loc for loc in [ 'center', 'left', 'right', 'top', 'bottom' ] })

xLabelKwargs = prepare_plot_kwargs(options.x_label_kwargs, [('size', options.axis_label_fontsize)])

yLabelKwargs = prepare_plot_kwargs(options.y_label_kwargs, [('size', options.axis_label_fontsize)])

#titleKwargs = prepare_plot_kwargs(options.additionalTitleKwargs, 
#        [('fontsize', options.titleFontsize), ('horizontalalignment', alignments[options.titleHalign]), ('verticalalignment', alignments[options.titleValign])] + [ split_equals_to_tuple(o) for o in defaultTitleKwargs ])

#the tick setup gets very complicated - see below
#first the general tickKwargs are applied to both axes, then any axis-specific ones
tickKwargs = prepare_plot_kwargs(options.tick_kwargs, [('labelsize', options.axis_tick_fontsize)])

yTickKwargs = prepare_plot_kwargs(options.x_tick_kwargs, [])

xTickKwargs = prepare_plot_kwargs(options.y_tick_kwargs, [])

#don't want to get things like options.line_width into plotKwargs here, as can be set for individual lines
plotKwargs = prepare_plot_kwargs(options.plot_kwargs, [])

#make a single centered label
if multiplot:
    if options.x_label and options.x_label_location == 's':
        fig.text(xLabelKwargs.pop('x'), xLabelKwargs.pop('y'), options.x_label, **xLabelKwargs)
    if options.y_label and options.y_label_location == 's':
        fig.text(yLabelKwargs.pop('x') * (origInchWidth / finalInchWidth), yLabelKwargs.pop('y'), options.y_label, **yLabelKwargs)

############
#add bars to identify the three triple resolutions
if not options.no_legend_bars:
    #a dummy invisible axis to allow easy specification of coordinates
    dummyAx = plt.axes([0, 0, 1, 1], axisbg=(1, 1, 1, 0))
    dummyAx.set_axis_off()

    barProportionalHeight = 0.07 / nrows
    barProportionalWidth = 0.3 * barRegionProportionalWidth

    centerLabelsOnBar = True
    labelXoffset = -0.015
    labelYoffset = 0.015

    if options.paren_labels:
        if len(options.paren_labels) not in [3, 4]:
            exit('must pass either zero, three or four paren labels')
        parenLabels = options.paren_labels
    else:
        parenLabels = filename_to_paren_trees(options.inFiles[0], oryza=options.oryza_names)

    parenLabelXoffset = 0.0 / ncols
    parenLabelYoffset = 0.07 / nrows

    if len(options.patterns) == 1:
        newList = []
        for h in options.patterns[0]:
            newList.append(h)
        options.patterns = newList

    barsAtLeft = True
    barForUnres = 1
    if barsAtLeft:
        barStartX = (barRegionProportionalWidth - barProportionalWidth) * 0.25
        #barStartX = 0.02
        if barForUnres:
            barStartY = 0.85
            barVspace = 0.175
            barHspace = 0.0
            barXoffset = 0.0
        else:
            barStartY = 0.75
            barVspace = 0.2
            barHspace = 0.0
            barXoffset = 0.0
        barYoffset = barProportionalHeight + barVspace

    else:
        barStartX = 0.108
        barStartY = 0.98
        barVspace = 0.0
        barHspace = 0.035
        barXoffset = barWidth + barHspace
        barYoffset = 0.0

    for series in range(3 + barForUnres):
        barX = barStartX + barXoffset * series
        barY = barStartY - barYoffset * series
        barXcenter = barX + 0.5 * barProportionalWidth
       
        barLabelFontsize = 28
        
        if centerLabelsOnBar:
            plt.text(barXcenter, barY + barProportionalHeight + labelYoffset, options.letter_labels[series], fontsize=barLabelFontsize, horizontalalignment='center', weight='bold')
        else:
            plt.text(barX - labelXoffset, barY - labelYoffset, options.letter_labels[series], fontsize=barLabelFontsize, horizontalalignment='center', weight='bold')

        plt.text(barXcenter + parenLabelXoffset, barY - parenLabelYoffset, parenLabels[series], fontsize=barLabelFontsize, horizontalalignment='center', weight='bold')
        
        #rec = patches.Rectangle((barX, barY), barProportionalWidth, barProportionalHeight, fc=options.greyValues[series], linewidth=4)
        rec = patches.Rectangle((barX, barY), barProportionalWidth, barProportionalHeight, fc=actualColors[series], linewidth=4)
        rec.set_hatch(options.patterns[series])
        dummyAx.add_patch(rec)
##########

#plot the lines, and set various properties
lines = []
wedges = []
for num, a in enumerate(ax):
    #def plot(self, axis, labels, colors, valueFormat, valueDist, labelDist):
    #if there is an unfilled row, there will be too many plots
    if num < numPlots:
        cumLines, pieWedges = figureData[num].plot(a, 
                labels=options.labels, 
                colors=actualColors,
                startangle=options.pie_start_angle,
                autopct=options.value_format_string, 
                pctdistance=options.value_distance, 
                labeldistance=options.label_distance,
                dzlabelsize=options.value_size,
                dztitlesize=options.title_fontsize,
                dzverticalalignment=alignments[options.title_vert_align],
                dzhorizontalalignment=alignments[options.title_horiz_align],
                dzncols=options.columns,
                plotKwargs=plotKwargs)

        
    #collect the lines and wedges so that we can make changes to them in bulk below.
    lines.append(cumLines)
    wedges.append(pieWedges)

    #this returns a list for some reason
    a.tick_params(**tickKwargs)
    a.tick_params(axis='x', **xTickKwargs)
    a.tick_params(axis='y', **yTickKwargs)
    for label in itertools.chain(a.get_xticklabels(), a.get_yticklabels()):
        label.set_weight('bold')
    #this gets confusing.  If sharex/y, it already deals with making single marginal tick labels (?), so don't touch anything
    #or they will all get cleared.  Only clear them manually if no labels at all were asked for, or marginal labels but not shared axes
    if options.x_tick_label_location == 'n' or (options.x_tick_label_location == 'm' and not sharex and not num >= numPlots - ncols):
        a.get_xaxis().set_ticklabels([])
    
    if options.y_tick_label_location == 'n' or (options.y_tick_label_location == 'm' and not sharey and not num % ncols == 0):
        a.get_yaxis().set_ticklabels([])
    
    #set the width of the axes box (frame) which I thought I could set on the frame rectangle, but that didn't work
    #kwargs here are <edge><w or s>, i.e. rightw or tops
    dzAxisKwargs = prepare_plot_kwargs(options.dz_axis_kwargs, [])

    for opt, val in dzAxisKwargs.items():
        which = opt[:-1]
        if opt[-1] == 'w':
            a.spines[which].set_linewidth(int(val))
        elif opt[-1] == 's':
            a.spines[which].set_linestyle(val)
        else:
            sys.exit('unrecognized dzAxisKwarg: %s=%s' % (opt, val))

    #draw the right side of the axis frame behind everything else
    a.spines['right'].set_zorder(1)
    
    #a list can be passed for ranges, but it is up to the user to make sure that it is the
    #right length, and tick labels will only be shown on the margins
    if not sharey and options.y_range:
        if len(options.y_range) == 2:
            a.set_ylim([options.y_range[0], options.y_range[1]])
        else:
            a.set_ylim([options.y_range[num * 2], options.y_range[num * 2 + 1]])

    if not sharex and options.x_range:
        if len(options.x_range) == 2:
            a.set_xlim([options.x_range[0], options.x_range[1]])
        else:
            a.set_xlim([options.x_range[num * 2], options.x_range[num * 2 + 1]])

    #if location option is m(arginal), add a normal axis label for plots on left or bottom margin
    if options.y_label:
        if options.y_label_location == 'm' or not multiplot:
            if num % ncols == 0:
                a.set_ylabel(options.y_label, **yLabelKwargs)
    if options.x_label:
        if options.x_label_location == 'm' or not multiplot:
            if num >= numPlots - ncols:
                a.set_xlabel(options.x_label, **xLabelKwargs)
   

#now go back and set various properties of the collected plot lines and markers
for lset in lines:
    for l, gv, m, ms, lw in itertools.izip(lset, itertools.cycle(actualColors), itertools.cycle(options.markers), itertools.cycle(options.marker_sizes), itertools.cycle(options.line_width)):
        l.set_marker(m)
        l.set_markersize(ms)
        l.set_color(gv)
        l.set_markeredgecolor(gv)
        l.set_markerfacecolor(gv)
        l.set_linewidth(lw)
        l.set_drawstyle('steps')

#set properties of the wedges
if options.patterns:
    for wedgeSet in wedges:
        for wedge, hatch in itertools.izip(wedgeSet, itertools.cycle(options.patterns)):
            wedge.set_hatch(hatch)

for wedgeSet in wedges:
    for wedge in wedgeSet:
        wedge.set_linewidth(options.pie_line_width)

if options.outFile and options.outFile is not sys.stdout:
    sys.stderr.write('Writing output %s ...\n' % options.outFile)
    plt.savefig(outFile, transparent=True, bbox_inches='tight')
else:
    sys.stderr.write('Displaying to screen ... \n')
    plt.show()

