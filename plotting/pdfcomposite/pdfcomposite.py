'''
Generating a composite figure from

1) a master pdf figure defining the overall layout (size, panel positions, etc), and 

2) a set of external pdf figures.

The method is based on LaTeX and TikZ.

Notes:
------

i)   On purpose, the method proposed here does not implement any rescaling of the 
     involved figures (though this would be easy). This enforces a clean workflow 
     where all figure components are prepared in exactly the size (figure dimensions, 
     font sizes, etc) needed for the final figure.

ii)  The positions of the external figures within the composite figure need to be 
     specified by the user by setting ext_figures_positions in create(). 

iii) For the master figure, a repositionig is not needed (and not implemented here): 
     as the master figure defines the overall layout of the final composite figure, 
     the composite figure inherits figure dimensions from the master figure. The master 
     figure is therefore positioned at the center of the composite figure.

(Tom Tetzlaff, t.tetzlaff@fz-juelich.de, 2020, 2024)

'''

import os

#####################################
def create(composite_figure_name_root, fig_size, master_figure_name, ext_figure_names, ext_figures_positions, draw_grid = False):    
    '''
    Creates a composite figure composed of a master figure and a set of external figures, 
    all of them availabe as pdf files. The resulting file is saved in pdf format.

    Arguments:
    ----------
    composite_figure_name_root: str
                                File name root of resulting composite figure
    
    fig_size:                   tuple
                                Dimensions (width, height) of the composite figure (in inch)

    master_figure_name:         str
                                Name of master figure pdf file (including extension ".pdf")

    ext_figure_names:           list(str)
                                List of file names of external figures (including extensions ".pdf")

    ext_figures_positions:      list(tuple)
                                List of position of external figure within the composite figure 
                                (the master figure is always centered at position (0,0))

    draw_grid:                  bool
                                If True: draws some grid lines to assist positioning of external figures

    Returns:
    --------
    0

    '''

    assert(master_figure_name[-4:]=='.pdf')
    assert(len(ext_figure_names)==len(ext_figures_positions))
    
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
    file.write(r"    {\includegraphics{%s}};" % (master_figure_name))
    file.write("\n")

    ## include external figures
    for i,ext_figure_name in enumerate(ext_figure_names):
        assert(ext_figure_name[-4:]=='.pdf')
        
        file.write(r"    \node[inner sep=-1pt,rectangle] (%s) at (%.4f,%.4f)" % (ext_figure_name, ext_figures_positions[i][0], ext_figures_positions[i][1]))
        file.write("\n")
        file.write(r"    {\includegraphics{%s}};" % (ext_figure_name))
        file.write("\n")

    ## draw grid
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

    ## generate pdf of composite figure using pdflatex
    os.system('pdflatex %s.tex' % composite_figure_name_root)

    ## generate eps of composite figure
    #os.system('latex %s.tex; dvips -o %s.eps %s.dvi' % (composite_figure_name_root,composite_figure_name_root,composite_figure_name_root))

    return 0
