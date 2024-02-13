'''
Example illustrating how to create a composite pdf figure from a maser figure 
generated using matplotlib, and an external pdf figure (created with inkscape).

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
def example():
    
    ##########
    '''
    Step 1: Create a (matplotlib) master figure and save it as pdf.
            This figure defines the layout of the final composite figure,
            including its size.
    '''
    master_figure_name, fig_size =  some_matplotlib_figure()

    ##########
    '''
    Step 2: Specify the names of the external pdf figures to be included, 
            and their positions in the final composite figure.
    '''
    ext_figure_names      = [ 'inkscape_sketch.pdf', 'inkscape_sketch.pdf' ] ## list of figures to be included (here, we use the same figure twice)
    ext_figures_positions = [ (-3.7,2.5)           , (3.7,2.5)     ] ## positions of external figures in composite figure (center = (0,0))

    ##########
    '''
    Step 3: Create the final composite figure using pdfcomposite.
    '''
    pdfcomposite.create('composite_example',
                        fig_size,
                        master_figure_name,
                        ext_figure_names,
                        ext_figures_positions)

######################################################################
if __name__ == '__main__':
    example()

