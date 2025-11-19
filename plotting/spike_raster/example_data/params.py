# -*- coding: utf-8 -*-
#
# network.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-2.0-or-later

#####################
'''
Parameters for generation and analysis of reference data.
'''

params = {
    #########################
    # Adapted model and simulation parameters
    #########################
    # scaling factor of the network
    'scaling_factor': 0.2, #1.0,
    # RNG seeds
    'RNG_seeds': ['12345' + str(i) for i in range(0, 10)], # RNG_seeds
    # pre simulation time for network stabilization in ms
    't_presim': 500.0,
    # simulation time for analysis in ms
    't_sim': 1.0e+4,
    #'t_sim': 9.0e+5,
    # local number of threads
    'local_num_threads': 4, #64,
    #########################
    # analysis parameters
    #########################
    # # start of analysis time interval in ms
    # 't_min': 500.0,
    # # seed for neuron subsampling reroducibility (for CC anlysis)
    # 'seed_subsampling': 12345, # seed_subsampling
    # # subsamples for pairwise statistical analysis
    # 'subsample_size': 250, # subsample_size
    # # bin size for generation of spike-count signals (for CC analysis)
    # 'binsize': 2.0, # binsize
    # # limits for rate, CV, and CC histograms
    # 'cc_lim': [-0.015, 0.015],
    # 'rate_lim': [0., 20.],
    # 'cv_lim': [0.5, 1.5],
    # # binsizes for histograms
    # 'cc_binsize': 0.001,
    # 'rate_binsize': 1.0,
    # 'cv_binsize': 0.05,
    # # maximal figure width in inches
    # 'max_fig_width': 7.5, # max_figure_width 
}
