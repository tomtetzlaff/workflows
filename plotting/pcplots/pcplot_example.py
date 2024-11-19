'''
Example of a parallel coordinates (PC) plot showing (some) data read from a csv file.

See also:
https://plotly.com/python/parallel-coordinates-plot/
https://plotly.com/python/reference/parcoords/

(Tom Tetzlaff, t.tetzlaff@fz-juelich.de, 2024)

'''

import pandas as pd
import plotly.graph_objects as go
import os

data_file = 'data.csv'

df = pd.read_csv(data_file)

fig = go.Figure(data=
                go.Parcoords(
                    line = dict(color = df['N'],
                                colorscale = 'Portland',
                                showscale = False,
                                ),
                    dimensions = list([
                        dict(label = '# of neurons N', values = df['N']),
                        dict(label = '# of assemblies A', values = df['A']),
                        dict(label = 'assembly size M', values = df['M']),
                        dict(label = '# of observed neurons N_obs', values = df['Nobs']),
                        dict(label = 'mean pattern size <s>', values = df['s_mean']),
                        dict(label = 'mean pattern count <c>', values = df['c_mean']),
                        dict(label = 'mean multiplicity <m>', values = df['m_mean']),                            
                    ]),
                    unselected=dict(line=dict(color='lightgray',opacity=0.0)),
                )
                )

fig.update_layout(
    font_family="Arial",
    font_color="black",
    font_size = 16,
    title_font_family="Arial",
    title_font_color="black",
    legend_title_font_color="black"
)

os.system('mkdir -p figures')
fig.write_html('figures/pc_plot.html')
