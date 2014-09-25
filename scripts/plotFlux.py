#!/usr/bin/env python
"""Arrow drawing example for the new fancy_arrow facilities.

Code contributed by: Rob Knight <rob@spot.colorado.edu>

DJZ - this is a fairly confusing sample file that I adapted for my purposes
and stripped of pylab dependency
(see http://matplotlib.org/examples/pylab_examples/arrow_demo.html)

"""

import sys
import itertools
import operator
from math import tan, radians, sin, sqrt, pi 

try:
    import matplotlib
except ImportError:
    sys.exit('Sorry, matplotlib is necessary to run this plotting script.\nSee http://matplotlib.org/downloads.html')
#this is important, because the backend can default to an interactive Tk based one even
#on clusters where there is no display
matplotlib.use('pdf')

import matplotlib.pyplot as plt

from pygot.utils import flatten_array


def add_tuples(t1, t2):
    return tuple(map(operator.add, t1, t2))


def sub_tuples(t1, t2):
    return tuple(map(operator.sub, t1, t2))


def make_arrow_plot(state_and_transition_counts, size=4, display='width', shape='right', 
        max_arrow_width=0.03, arrow_sep=0.02, alpha=0.5, 
        normalize_data=False, ec=None, labelcolor=None, 
        head_starts_at_zero=True, show_arrow_for_scale=False, 
        arrow_function=None, 
        **kwargs):
    """Makes an arrow plot.

    Parameters:

    state_and_transition_counts: dict with probabilities for the bases and pair transitions.
    size: size of the graph in inches.
    display: 'length', 'width', or 'alpha' for arrow property to change.
    shape: 'full', 'left', or 'right' for full or half arrows.
    max_arrow_width: maximum width of an arrow, data coordinates.
    arrow_sep: separation between arrows in a pair, data coordinates.
    alpha: maximum opacity of arrows, default 0.8.

    **kwargs can be anything allowed by a Arrow object, e.g. linewidth and edgecolor.

    DJZ - this is badly hacked up by me based on the example function
    """

    monochrome = True
    
    states = [ 'U', 'S', '1', '2' ]
    #this just allows the internal definition of the states to remain the same, but the output letters to change
    stateOutputMapping = { 'U':'U', 'S':'M', '1':'1', '2':'2' }
    transitions = [ ''.join(t) for t in itertools.product(states, repeat=2) ]

    #Dimensions.  Half the length of one of the triangle edges. This is arbitrary really, 
    #and the locations of things will rescale if it is changed - the size of text will NOT rescale though
    halfEdge = 2.0
    edge = halfEdge * 2.0
    halfVert = halfEdge * tan(radians(30))
   
    #Specify the location of the state characters at the points and center of 
    #the triangle. U is in the middle, S on top, 1 bottom left, 2 bottom right
    state_coords = {}
    coordU = (halfEdge, halfVert)
    state_coords['U'] = coordU

    coordS = (halfEdge, sin(radians(60)) * edge)
    state_coords['S'] = coordS
    
    coord1 = (0, 0)
    state_coords['1'] = coord1

    coord2 = (edge, 0)
    state_coords['2'] = coord2
    
    #spacing around arrows
    arrow_offset_from_point = 0.15 * halfEdge
    diag_arrow_offset_from_point_small = sin(radians(30)) * arrow_offset_from_point
    diag_arrow_offset_from_point_large = sin(radians(60)) * arrow_offset_from_point
    arrow_gap = 0.02 * halfEdge

    #arrow limits
    max_arrow_length = 0.5
    max_head_width = 2.0 * max_arrow_width
    max_head_length = max_head_width / 20.0
    
    #Make the figure and axes
    #The xlim and ylim are the coords of the boundaries of the graph - this doesn't have 
    #anything to do with the actual size of the graph, just defines how the coordinates 
    #used to determine where the text etc. are placed are defined
    #figsize is in inches
    fig = plt.figure(figsize=(size, size))
    #this is the proportion of the figure taken up by the graph, specified
    #as left, bottom, width, height
    axes_rect = (0.125, 0.1, 0.775, 0.8)
    ax = fig.add_axes(axes_rect, frameon=False,
        xlim=(-0.5, edge + 1.0), 
        ylim=(-0.5, (edge + 1.0) * sin(radians(60))), 
        xticks=([]), yticks=([])
        )

    text_params = {'ha':'center', 'va':'center', 'family':'sans-serif', 'fontweight':'bold'}
 
    #this would draw arrows all the way to the states, but other coords would need to be altered first, 
    #since they were specified relative to these origins rather than the states
    #arrow_offset_from_point = 0.0
    
    #tuple of x, y for start position of arrows
    #gotten by taking location of U, 1, 2 or 3 and offsetting by some gap values depending on the direction of the arrow
    arrow_origins = {
        'US': add_tuples(coordU, (-arrow_gap, arrow_offset_from_point)),
        'U1': add_tuples(coordU, (arrow_gap - diag_arrow_offset_from_point_large, -arrow_gap - diag_arrow_offset_from_point_small)),
        'U2': add_tuples(coordU, (arrow_gap + diag_arrow_offset_from_point_large, arrow_gap - diag_arrow_offset_from_point_small)),
        'SU': add_tuples(coordS, (arrow_gap, -arrow_offset_from_point)),
        '1U': add_tuples(coord1, (-arrow_gap + diag_arrow_offset_from_point_large, arrow_gap + diag_arrow_offset_from_point_small)),
        '2U': add_tuples(coord2, (-arrow_gap - diag_arrow_offset_from_point_large, -arrow_gap + diag_arrow_offset_from_point_small)),
        'S2': add_tuples(coordS, (arrow_gap + diag_arrow_offset_from_point_small, arrow_gap - diag_arrow_offset_from_point_large)),
        '2S': add_tuples(coord2, (-arrow_gap - diag_arrow_offset_from_point_small, -arrow_gap + diag_arrow_offset_from_point_large)),
        'S1': add_tuples(coordS, (arrow_gap - diag_arrow_offset_from_point_small, -arrow_gap - diag_arrow_offset_from_point_large)),
        '1S': add_tuples(coord1, (-arrow_gap + diag_arrow_offset_from_point_small, arrow_gap + diag_arrow_offset_from_point_large)),
        '12': add_tuples(coord1, (arrow_offset_from_point, arrow_gap)),
        '21': add_tuples(coord2, (-arrow_offset_from_point, -arrow_gap)),
        }
    
    #sub_tuples just calcs the vector distances between the coords
    #arrow_vectors are the x and y offsets from the start to end of the arrows
    #they aren't data dependent, since lengths are currently fixed and widths vary
    arrow_vectors = { f + t: sub_tuples(state_coords[t], state_coords[f]) for f in states for t in states if f != t }
    
    if monochrome:
        colors = { st:'k' for st in states }
        colors.update({ tran:'k' for tran in transitions })
    else:
        colors = {
            'U':'r',
            'S':'k',
            '1':'g',
            '2':'b',
            'US':'r',
            'SU':'k',
            '1U':'g',
            'U1':'r',
            '2U':'b',
            'U2':'r',
            '1S':'g',
            'S2':'k',
            '2S':'b',
            'S1':'k',
            '12':'g',
            '21':'b',
            }

    #alignment of arrow label text with respect to ???
    label_positions = {
        'US':'center',
        'SU':'center',
        '1U':'center',
        'U1':'center',
        '2U':'left',
        'U2':'left',
        '1S':'left',
        'S2':'left',
        '2S':'center',
        'S1':'center',
        '12':'center',
        '21':'center'
        }

    if show_arrow_for_scale:
        arrow_origins['10%%'] = (0, 1)
        arrow_vectors['10%%'] = (0, halfEdge)
        colors['10%%'] = 'k'
        label_positions['10%%'] = 'center'

    #START OF DATA DEPENDENT STUFF
    #rescale
    #the counts might be used as is, so back them up
    keys = [ x[0] for x in state_and_transition_counts ]
    rawCounts = [ x[1] for x in state_and_transition_counts ]

    #The total number of loci that change resolutions, for text output in the figure (see near the end of this function)
    totChanges = sum([x[1] for x in state_and_transition_counts if len(x[0]) > 1])
    #The number of loci that change resolutions, but ONLY if neither state is ambiguous
    nonAmbigChanges = sum([x[1] for num, x in enumerate(state_and_transition_counts) if len(x[0]) > 1 and num >= 4 and num % 4 != 0])

    #dealing with the fact that the state counts might be tuples indicating the before and after counts, 
    #where in the legacy format they used to only be the after
    if len(rawCounts[0]) > 1:
        #grab the first element of the state tuples into the before states
        beforeStateCounts = [ c[0] for c in rawCounts[::5] ]
        #replace the state tuples with the second element of the tuples
        rawCounts[::5] = [ c[1] for c in rawCounts[::5] ]
    else:
        beforeStateCounts = None
   
    rescaledValues = list(rawCounts)
    upscale = 1.0
    #rescale transition counts with respect to diagonal values (total loci)
    #the '12' in the rescaling below is arbitrary, gotten by experimentation
    totCountOld = sum(rawCounts[x * 4 + x] for x in xrange(4))
    totCount = sum(rawCounts[::5])
    assert totCount == totCountOld
    for i in xrange(4):
        for j in xrange(4):
            if i != j:
                rescaledValues[i * 4 + j] *= upscale / float(totCount)
            else:
                rescaledValues[i * 4 + j] *= 1.0 / (totCount * 12)
    if beforeStateCounts:
        beforeStateCounts = [ bc * 1.0 / (totCount * 12) for bc in beforeStateCounts ]
   
    state_and_transition_counts = dict(zip(keys, rescaledValues))
    countData = dict(zip(keys, rawCounts))
    beforeStateCounts = dict(zip(states, beforeStateCounts))

    #Determine where exactly the resolution characters (e.g., U, S, 1, 2) should be placed so that 
    #they appear to be centered in the circles
    #72 points per inch, need to move coords so that character is centered on point
    #actually, this offset is way too large in that case.  Hmm.
    textSize = 28
    stateOffset = {state:-((textSize / 16.0) / 72.0) for state in ['U', 'S', '2']}
    #the character "1" is odd, and the veritcal bar doesn't fall along the center if there is no bottom bar, 
    #and it looks weird when it is alone, so, different offset for it
    stateOffset['1'] = -((textSize / 9.0) / 72.0)
 
    #zorders determine which elements overlay which others
    zorder_arrow = 1
    zorder_circle = 5
    zorder_text = 10

    for state in states:
        #draw the text indicating the state names, and circles indicating the proportions in each of the states
        #in the left and right treatments
        x, y = state_coords[state]
        lineWeight = 4
        
        ax.text(x + stateOffset[state], y, stateOutputMapping[state], color=colors[state], size=textSize, 
                horizontalalignment='left', zorder=zorder_text, **text_params)
        
        if state_and_transition_counts[state] > beforeStateCounts[state]:
            zorder_circ1 = zorder_circle - 1
        else:
            zorder_circ1 = zorder_circle + 1
        
        #area proportional to count - this seemed to look best
        circ = plt.Circle((x, y), radius=sqrt(state_and_transition_counts[state] / pi) * 5.0, 
                fill=True, ec='k', fc='w', lw=lineWeight, zorder=zorder_circ1)
        circ2 = plt.Circle((x, y), radius=sqrt(beforeStateCounts[state] / pi) * 5.0, 
                fill=True, fc='w', ec='k', lw=lineWeight, ls='dashed', zorder=zorder_circle)

        #radius proportional to count
        #circ = Circle((x, y), radius=state_and_transition_counts[state] * 10, 
        #        fill=False, ec='b', lw=3, zorder=zorder_circ1)
        #circ2 = Circle((x, y), radius=beforeStateCounts[state] * 10, 
        #        fill=False, ec='b', lw=3, ls='dashed', zorder=zorder_circle)
        ax.add_patch(circ)
        ax.add_patch(circ2)

    arrow_params = {'length_includes_head':True, 'shape':shape, 'head_starts_at_zero':head_starts_at_zero}
    
    #This is a closure function defined within the make_arrow_plot function, and is called below
    #what is passed in as the pair here is just a string of length 2, that designates a particular 
    #arrow, i.e. 'U1', '12', etc., which is used to look up details of the arrows as keys in various dicts
    #I would not have implemented it this way, but ...
    def draw_arrow(pair, alpha=alpha, ec=ec, labelcolor=labelcolor):
        #in theory the data can be mapped to the arrow length, transparency or width
        #although I'm not sure that anything besides width is currently working

        #set the length of the arrow
        if display == 'length':
            length = max_head_length + (max_arrow_length - max_head_length) * state_and_transition_counts[pair]
        else:
            length = max_arrow_length
        #set the transparency of the arrow
        if display == 'alph':
            alpha = min(state_and_transition_counts[pair], alpha)
        else:
            alpha = alpha
        #set the width of the arrow
        if display == 'width':
            scale = state_and_transition_counts[pair]
            #this default scaling is totally arbitrary, gotten by experimentation it probably won't work well if there are 
            #proportions > than about 0.20 the addition on the end here just beefs up the smallest arrows a tad, 
            #without affecting the largest much
            #default defined in argparse defaults: lambda scale: (pow(scale, 0.75) / 1.6) * max_arrow_width + 0.04
            widthFunc = eval(arrow_function)
            width = widthFunc(scale)

            #seems like these should scale with the arrow width
            head_width = width * 1.7
            head_length = head_width * 0.5
        else:
            width = max_arrow_width
            head_width = max_head_width
            head_length = max_head_length

        #set the face and (I think) edge colors. If edgecolor isn't defined,
        #default to same as fc
        fc = colors[pair]
        ec = ec or fc

        #the start point of the arrows
        x_pos, y_pos = arrow_origins[pair]
        #the vector indicating the change in x and y coordinate from start to end of arrows
        x_scale, y_scale = arrow_vectors[pair]
        
        #the actual call to matplotlib to make the arrow. Signature is
        #matplotlib.pyplot.arrow(x, y, dx, dy, hold=None, **kwargs)
        #Draws arrow on specified axis from (x, y) to (x + dx, y + dy)
        if countData[pair] > 0:
            plt.arrow(x_pos, y_pos, x_scale * length, y_scale * length, 
                fc=fc, ec=ec, alpha=alpha, width=width, head_width=head_width, 
                head_length=head_length, **arrow_params)

        label_text_size = 30
   
        counts_to_show_text = 10
        #counts_to_show_text = 20

        if countData[pair] >= counts_to_show_text:
            x, y = add_tuples((x_pos, y_pos), (x_scale * length * 0.5, y_scale * length * 0.5))
            ax.text(x, y, '%d' % (countData[pair]), size=label_text_size, ha='center', va='center', weight='bold',
                color=labelcolor or fc)
    #####END OF EMBEDDED FUNCTION####

    for p in arrow_origins.keys():
        draw_arrow(p)

    #originally indicated only number of non-ambiguous changes between loci
    #gca().set_title('%d loci\n%d shifts\n%.1f%%\n%d non-U shifts\n%.1f%%' % 
            #(totCount, totChanges, (100 * totChanges) / float(totCount), nonAmbigChanges, (100 * nonAmbigChanges) / float(totCount)), 
            #x=0.8, y=0.7, fontsize=20)
    
    #indicate the total number of loci that changed states
    #text in the upper right
    #ax.set_title('%d changes\n(%.1f%%)' % (totChanges, (100 * totChanges) / float(totCount)), 
            #x=0.75, y=0.62, fontsize=28, weight='bold')
   
    #total changes centered text below - the triangle isn't actually centered in the figure, so 0.45 is about right
    ax.set_title('%d changes (%.1f%%)' % (totChanges, (100 * totChanges) / float(totCount)), 
            x=0.45, y=-0.025, horizontalalignment='center', fontsize=28, weight='bold')

    return fig

if __name__ == '__main__':
    from argparse import ArgumentParser

    #use argparse module to parse commandline input
    parser = ArgumentParser(description='make gene tree flux figure')

    #add possible arguments
    parser.add_argument('--scale-arrow', action='store_true', default=False,
                        help='show an arrow for scale (default false)')

    parser.add_argument(dest='infile', nargs='?', type=str, default=None, 
                        help='file containing the states and transition counts, or none for stdin')

    parser.add_argument('-o', '--outfile', type=str, default=None,
                        help='File to write figure to. By default determined from infile name. \
                            NOTE: desired file format will be detected from file extension')

    default_arrow_function = "lambda scale: (pow(scale, 0.75) / 1.6) * max_arrow_width + 0.04"
    parser.add_argument('--arrow-function', type=str, default=default_arrow_function,
                        help='lambda function string to convert substitution frequencies to arrow width or length. \
                        (default \'%s\')' % default_arrow_function )
    
    options = parser.parse_args()

    indata = open(options.infile, 'rb') if options.infile else sys.stdin

    vals = flatten_array([ line.split() for line in indata ])
    vals = [float(val) for val in vals]
    if not vals:
        exit('%s is empty?' % sys.argv[1])
    
    #For legacy reasons, passing matrix with total counts of left and right treatments in tuples on 
    #diagonal.  Other matrix elements represent counts of loci changing from resolution i to resolution j 
    #from left to right.  In that case the number of elements passed in will be 20.
    #A more reasonable way of passing the data is a 4x4 matrix in which the diagonal elements
    #are just the counts of loci with the same resolution in both treatments.  In that case the total loci
    #counts are just the sums of the rows and columns.  If 16 values are passed, assume that interpretation
    #but still convert to the 20 element version when passing to make_arrow_plot
    if len(vals) > 16:
        #in this case we're assuming that both the before and after counts were part of the input string
        vals = [ (vals[0], vals[1]), vals[2], vals[3], vals[4], vals[5], (vals[6], vals[7]), 
                vals[8], vals[9], vals[10], vals[11], (vals[12], vals[13]), vals[14], vals[15], 
                vals[16], vals[17], (vals[18], vals[19]) ]

    else:
        leftCounts = [ sum([int(vals[el]) for el in xrange(s * 4, s * 4 + 4)]) for s in xrange(4) ]
        rightCounts = [ sum([int(vals[el]) for el in xrange(s, 16, 4)]) for s in xrange(4) ]

        for state in xrange(4):
            vals[state + state * 4] = ( leftCounts[state], rightCounts[state] )

    order = ['U', 'US', 'U1', 'U2', 'SU', 'S', 'S1', 'S2', '1U', '1S', '1', '12', '2U', '2S', '21', '2']

    #whether to have an arrow indicating the scaling of arrow width to the left of the trianglular figure
    if options.scale_arrow:
        vals.append(0.4)
        order.append('10%%')

    dzData = zip(order, vals)

    display = 'width'
    max_arrow_width = 1.8
   
    fig = make_arrow_plot(dzData, 
            max_arrow_width=max_arrow_width, 
            display=display, 
            linewidth=0.001, 
            edgecolor=None,
            normalize_data=True,
            head_starts_at_zero=True, 
            size=8, 
            show_arrow_for_scale=options.scale_arrow,
            arrow_function=options.arrow_function)
    
    if options.outfile:
        outfile = options.outfile
    elif options.infile:
        outfile = options.infile + '.pdf'
    else:
        #if data came in on stdin, but no outfile specified
        outfile = 'figure.pdf'

    plt.savefig(outfile, transparent=True, bbox_inches='tight')

