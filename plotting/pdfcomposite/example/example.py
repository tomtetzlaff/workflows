'''
Example illustrating how to create a composite pdf figure (vector graphics) from a set of individual pdf figures:

1) a master figure (here: a matplotlib figure) defining the overall layout (size, panel positions, etc), and 

2) two external figures (here: a sketch made using inkscape).

The method is based on LaTeX and TikZ.

Notes:
------

i) On purpose, the method proposed here does not implement any rescaling of the involved figures (though this would be easy). This enforces a clean workflow where all figure components are prepared in exactly the size (figure dimensions, font sizes, etc) needed for the final figure.

ii) The positions of the external figures within the composite figure need to be specified by the user by setting pos_ext_figures in create_composite_figure(). 

iii) For the master figure, a repositionig is not needed (and not implemented here): as the maser figure defines the overall layout of the final composite figure, the composite figure inherits figure dimensions from the master figure. The master figure is therefore positioned at the center of the composite figure.

(Tom Tetzlaff, t.tetzlaff@fz-juelich.de, 2020, 2024)

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

import sys
sys.path.append('..') 
import pdfcomposite

######################################################################
def panel_label(s,pos):
    '''
    Creates a panel label (A,B,C,...) for the current axis object of a matplotlib figure.

    Arguments:
    ----------
    s:   str
         panel label

    pos: tuple
         x-/y- position of panel label (in units relative to the size of the current axis)

    Returns:
    --------
    0

    '''
    
    ax = plt.gca()
    plt.text(pos[0],pos[1],r'\bfseries{}%s' % s,transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',size=12)

    return 0

######################################################################
def some_matplotlib_figure():
    '''
    Creates some matplotlib figure and saves the result as pdf.

    Returns:
    --------
    fname:    str
              Name of matplotlib figure.
    fig_size: tuple
              Figure dimensions (width, height) (in inch).
    '''

    fig_size = (6,4)               ## figure size (width, height) in inches
    font_size = 10                 ## font  size
    dpi = 300                      ## print resolution (need for print versions of this figure, even if the figure is vector graphics)

    panel_label_pos = (-0.18,1.0)  ## x-/y- positions of panel labels (in units of axis size)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size']= font_size
    plt.figure(1,figsize=fig_size,dpi=dpi)
    plt.clf()

    gs = gridspec.GridSpec(2, 2, width_ratios = [1,1], left=0.1, right=0.98, bottom=0.11, top=0.97, wspace=0.3, hspace=0.35)

    ## panel A (placeholder for svg figure to be inserted; see below)
    plt.subplot(gs[0,0])
    plt.axis('off')
    panel_label('A',panel_label_pos)

    ## panel B
    plt.subplot(gs[0,1])
    plt.axis('off')
    panel_label('B',panel_label_pos)
    #plt.plot(np.random.randn(200),lw=1,color='k')
    #plt.xlabel(r'time $t$ (ms)')
    #plt.ylabel(r'activity $x_\mathsf{X}(t)$')
    #panel_label('B',panel_label_pos)

    ## panel C
    plt.subplot(gs[1,0])
    plt.plot(np.random.randn(200),lw=1,color='k')
    plt.xlabel(r'time $t$ (ms)')
    plt.ylabel(r'activity $x_\mathsf{E}(t)$')
    panel_label('C',panel_label_pos)

    ## panel D
    plt.subplot(gs[1,1])
    plt.plot(np.random.randn(200),lw=1,color='k')
    plt.xlabel(r'time $t$ (ms)')
    plt.ylabel(r'activity $x_\mathsf{I}(t)$')
    panel_label('D',panel_label_pos)

    fname = 'matplotlib_figure.pdf'
    plt.savefig("%s" % fname)
    #plt.savefig("%s.eps" % fname) 

    return fname, fig_size

######################################################################

## create a matplotlib master figure and save as pdf (and eps)
# (this figure defines the layout of the final composite figure)
master_file_name, fig_size =  some_matplotlib_figure()

## name (root) of resulting composite figure
composite_figure_name_root = 'composite_example'

## name of external pdf file
ext_file_name = 'inkscape_sketch.pdf'  ## here: created using inkscape
#os.system('inkscape %s.pdf --export-eps=%s.eps' % (ext_file_name,ext_file_name))  ## convert to eps

ext_file_names  = [ ext_file_name, ext_file_name ] ## list of figures to be included
pos_ext_figures = [ (-3.7,2.5)   , (3.7,2.5)     ] ## positions of external figures in composite figure (center = (0,0))

## generate composite figure (with dimensions inherited from the master figure)
pdfcomposite.create(composite_figure_name_root,fig_size,master_file_name,ext_file_names,pos_ext_figures)
