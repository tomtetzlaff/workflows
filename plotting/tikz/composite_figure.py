'''
Example illustrating how to create a composite figure (vector graphics) from two separate figures (vector graphics):

1) a master figure (here: a matplotlib figure) defining the overall layout (size, panel positions, etc), and 

2) an external figure (here: a sketch made using inkscape).

Each of these two figures is loaded from a pdf file.

The method is based on LaTeX and TikZ.

(Tom Tetzlaff, t.tetzlaff@fz-juelich.de, 2020)

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

    fname = 'matplotlib_figure.pdf'
    plt.savefig(fname) 

    return fname, fig_size

#####################################
def create_composite_figure(composite_figure_name_root,master_file_name,ext_file_name,pos_ext_file = (0,0),draw_grid = False):
    '''
    Creates a composite figure composed of a master figure and an external figure, both availabe as pdf files. The resulting file is saved in pdf format.

    Arguments:
    ----------
    composite_figure_name_root: str
                                File name root of resulting composite figure
    
    master_file_name:              str
                                File name of master figure

    ext_file_name:              str
                                File name of external figure

    pos_ext_file:               tuple
                                Position of external figure within the composite figure (center=(0,0))

    draw_grid:                  bool
                                if True, draws some grid lines to assist positioning of external figure

    Returns:
    --------
    0

    '''
    
    file = open('%s.tex' % composite_figure_name_root , 'w')
    file.write(r"\documentclass{article}")
    file.write("\n")
    file.write(r"\usepackage{geometry}")
    file.write("\n")
    file.write(r"\geometry{paperwidth=%.3fin, paperheight=%.3fin, top=0pt, bottom=0pt, right=0pt, left=0pt}" % (fig_size[0],fig_size[1]))
    file.write("\n")
    file.write(r"\usepackage{tikz}")
    file.write("\n")
    file.write(r"\usepackage{graphicx}")
    file.write("\n")
    file.write(r"\pagestyle{empty}")
    file.write("\n")
    file.write(r"\begin{document}")
    file.write("\n")
    file.write(r"\noindent")
    file.write("\n")
    file.write(r"\resizebox{\paperwidth}{!}{%")
    file.write("\n")
    file.write(r"  \begin{tikzpicture}%")
    file.write("\n")
    file.write(r"    \node[inner sep=-1pt] (matplotlib_figure) at (0,0)")
    file.write("\n")
    file.write(r"    {\includegraphics{%s}};" % (master_file_name))
    file.write("\n")
    file.write(r"    \node[inner sep=-1pt,rectangle] (inkscape_sketch) at (%.4f,%.4f)" % (pos_ext_file[0],pos_ext_file[1]))
    file.write("\n")
    file.write(r"    {\includegraphics{%s}};" % (ext_file_name))
    file.write("\n")
    if draw_grid:
        file.write(r"    \draw[style=help lines] (-6,-4) grid (6,4);")
        file.write("\n")
    file.write(r"  \end{tikzpicture}%")
    file.write("\n")
    file.write(r"}")
    file.write("\n")
    file.write(r"\end{document}")
    file.write("\n")

    file.close()

    ## execute tex script
    os.system('pdflatex %s.tex' % composite_figure_name_root)

    return 0

########################################

## create a matplotlib master figure and save as pdf
# (this figure defines the layout of the final composite figure)
master_file_name, fig_size =  some_matplotlib_figure()

## name (root) of resulting composite figure
composite_figure_name_root = 'composite_example'

## name of external pdf file
ext_file_name = 'inkscape_sketch.pdf'  ## here: created using inkscape

## generate composite figure 
pos_ext_file = (-3.7,2.5)      ## position of external file in composite figure (center = (0,0)
create_composite_figure(composite_figure_name_root,master_file_name,ext_file_name,pos_ext_file)
    
