'''
Example illustrating how to merge a matplotlib figure with an svg figure loaded from file by using the svgutils package.

Workflow:
---------

1) Make sure svgutils is installed (import svgutils). If not, run:

   pip install svgutils

2) Create a matplotlib figure the way you like (see some_matplotlib_figure() in example below). In this step, reserve space (and possibly add a panel label) for the svg figure to be inserted (panel A in the example below). Save this figure as pdf (needed for adjusting the figure style, see step 3).

3) Tune the style of the final figure (size, panel positions, etc) based on the pdf of the matplotlib figure from step 2. To this end, run the figure generation script from the terminal and look at the pdf file. Don't use ipython here!

3) Merge the matplotlib figure with some svg figure loaded from file (e.g. a file generated with inkscape) using the function merge_mpl_with_svg(). The size of the final figure is exactly identical to the size of the matplotlib figure. Nevertheless, this step requires a rescaling of the inserted matplotlib figure (mpl_scale), which needs to be done by hand based on the final pdf. Further, this step involves handtuning of the size (svg_scale) and position (svg_pos) of the inserted svg figure.


(Tom Tetzlaff, t.tetzlaff@fz-juelich.de, 2020)

'''


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import svgutils.transform as sg
import os

#####################################
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

#####################################
def some_matplotlib_figure():
    '''
    Creates some matplotlib figure.

    Returns:
    --------
    fig_mpl:   matplotlib.figure.Figure
               matplotlib figure object
    '''

    fig_size = (6,4)               ## figure size (width, height) in inches
    font_size = 10                 ## font  size
    dpi = 300                      ## print resolution (need for print versions of this figure, even if the figure is vector graphics)

    panel_label_pos = (-0.18,1.0)  ## x-/y- positions of panel labels (in units of axis size)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size']= font_size
    fig_mpl = plt.figure(1,figsize=fig_size,dpi=dpi)
    plt.clf()

    gs = gridspec.GridSpec(2, 2, width_ratios = [1,1], left=0.1, right=0.98, bottom=0.11, top=0.97, wspace=0.3, hspace=0.35)

    ## panel A (placeholder for svg figure to be inserted; see below)
    plt.subplot(gs[0,0])
    plt.axis('off')
    panel_label('A',panel_label_pos)

    ## panel B
    plt.subplot(gs[0,1])
    plt.plot(np.random.randn(200),lw=1,color='k')
    plt.xlabel(r'time $t$ (ms)')
    plt.ylabel(r'activity $x_\mathsf{X}(t)$')
    panel_label('B',panel_label_pos)

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

    plt.savefig('matplotlib_figure.pdf')  ## needed for tuning the core of the figure
    #plt.savefig('matplotlib_figure.svg')

    return fig_mpl

#############

def merge_mpl_with_svg(mpl_figure,svg_file_name,mpl_scale,svg_scale,svg_pos,file_name_root):
    '''
    Merges a matplotlib figure with an svg figure loaded from file, and saves the final file in svg and pdf format.

    Arguments:
    ----------
    mpl_figure:     matplotlib.figure.Figure
                    matplotlib figure object

    svg_file_name:  str
                    file name of inserted svg figure

    mpl_scale:      float
                    factor scaling size of mpl figure

    svg_scale:      float
                    factor scaling size of inserted svg figure

    svg_pos:        tuple
                    x-/y-position of inserted svg plot

    file_name_root: string
                    name root ("xyz") of final svg ("xyz.svg") and pdf file ("xyz.pdf")

    Returns:
    --------
    0


    Comments:
    ---------
    The size of the final (merged) figure is exactly identical to the dimensions of the matplotlib figure. Unfortunately, the size of the matplotlib figure is however not preserved when inserting it into the final (svg) figure. The matplotlib figure therefore needs to be rescaled by hand (mpl_scale). The same holds for the svg figure loaded from file: Its size is not preserved after inclusion, and therefore needs to rescaled (svg_scale). These rescaling problems are nasty, because they prohibit an a priori definition of font sizes (and the styles of other figure elements such as line widths). Not sure what I'm doing wrong here.

    '''
    
    fig_svg_size = mpl_figure.get_size_inches()  ## get size of matplotlib figure

    ## create new svg figure
    fig_svg = sg.SVGFigure("%.4fin" % fig_svg_size[0], "%.4fin" % fig_svg_size[1])

    ## load matplotlib figure
    panel_BCD = sg.from_mpl(mpl_figure)           ## get matplotlib figure object
    #fig1 = sg.fromfile('matplotlib_figure.svg')  ## alternative: load from file

    ## load svg figure (create by inkscape)
    panel_A = sg.fromfile(svg_file_name)

    ## get the figure objects
    panel_BCD_root = panel_BCD.getroot()
    panel_A_root = panel_A.getroot()

    ## reposition and rescale figure objects (this requires hand-tuning)
    panel_BCD_root.moveto(0, 0, scale=mpl_scale)  
    panel_A_root.moveto(svg_pos[0], svg_pos[1], scale=svg_scale)

    ## append objects
    fig_svg.append([panel_BCD_root, panel_A_root])

    ## save final svg figure
    fig_svg.save("%s.svg" % file_name_root)

    ## convert to pdf
    os.system('inkscape --export-pdf=%s.pdf %s.svg' % (file_name_root,file_name_root))

    return 0


#####################################

## create some matplotlib figure
fig_mpl = some_matplotlib_figure()

## merge matplotlib figure with svg loaded from file and save as svg and pdf
merge_mpl_with_svg(fig_mpl,'inkscape_sketch.svg',1.32,2.7,(40,15),'final_figure')


