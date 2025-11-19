'''
This example demonstrates how to create high-density a spike raster plots, i.e.,
plots where the density of points is high such that points stronlgy overlap.
The solution proposed here circumvents the use of transparency (alpha blending),
and thereby permits creating high-quality figures suitable for publication in journals
requiring figures in eps format (which does not support transparency). It is based on
coloring spikes according to the instantaneous population rate.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

#####################################################################################
def spike_raster_dense(times, senders, ax, T, markersize = 1, clim = (0.6,0.), dt = 1.):
    """
    Create a scatter plot of spike times vs neuron ids,
    where the color of each spike indicates the instantaneous population rate.
    
    Parameters
    ----------
    times : array-like
        Array of spike times (in ms).
    
    senders : array-like
        Array of neuron ids corresponding to each spike time.

    ax : matplotlib.axes.Axes
        The axes on which to plot the scatter plot.

    T : float
        Total duration of the recording (in ms).
    
    markersize : float, optional
        Size of the markers in the scatter plot. Default is 1.

    clim : tuple of float, optional
        Color limits for the grayscale mapping.
        The 1st value corresponds to the lowest rate, the 2nd to the highest.
        Default is (0.6, 0.).
    
    dt : float, optional
        Bin size (in ms) for rate estimation. Default is 1.0 ms.

    
    Returns:
    --------

    rate : array-like
        The estimated population rate over time.
    
    ts : array-like
        The time bins corresponding to the rate estimation.    
    """

    ## calculate spike densities
    ## (takes very long)
    # from scipy.stats import gaussian_kde
    # xy = np.vstack([times, senders])
    # z = gaussian_kde(xy)(xy)

    #z,_=np.histogram2d(times, senders, bins=10)
    #print(type(z))
    #print(z)
    #stop
    
    ## calculate instantaneous rate
    print("calculating population rate...")
    ts = np.arange(0, T+dt, dt)
    rate = np.histogram(times, bins=ts)[0]/(dt*1e-3) ## population rate (1/s)
    times_idx = np.digitize(times, ts)  ## discretize spike times to the same grid used for the rate estimation
    z = rate[times_idx-1]               ## assign each spike a value correponding to the current rate
    
    z -= z.min()  ## shift to [0, ...]
    z /= z.max()  ## normalize to [0,1]
    clr = (clim[1]-clim[0])*z + clim[0]  ## map z to grayscale in range [clim[0], clim[1]]

    # sort spikes by their rate such that spikes ocurring at high rates are plotted last
    idx = np.argsort(z) #z.argsort()
    times, senders, clr = times[idx], senders[idx], clr[idx]

    plt.scatter(times, senders, marker='o', s=markersize, linewidths = 0, color=clr.astype('str'), rasterized=True, )

    return rate, ts

#####################################################################################
def load_spike_data(path, label, time_interval = None, pop = None, skip_rows = 3):        
    '''
    (adapted from https://github.com/INM-6/microcircuit-PD14-model/blob/main/PyNEST/src/microcircuit/helpers.py)
    
    Load spike data from files.

    Arguments
    ---------
    path:           str
                    Path to folder containing spike files.

    label:          str
                    Spike file label (file name root).

    time_interval:  None (default) or tuple (optional)
                    Start and stop of observation interval (ms). All spikes outside this interval are discarded.
                    If None, all recorded spikes are loaded.

    pop:            None (default), or list (optional)
                    Oberserved neuron population. All spike senders that are not part of this population are discarded.
                    If None, all recorded spikes are loaded.

    skip_rows:      int (optional)
                    Number of rows to be skipped while reading spike files (to remove file headers). The default is 3.

    Returns
    -------
    spikes:   numpy.ndarray
              Lx2 array of spike senders spikes[:,0] and spike times spikes[:,1] (L = number of spikes).

    '''

    if type(time_interval) == tuple:
        print('Loading spike data in interval (%.1f ms, %.1f ms] ...' % (time_interval[0], time_interval[1]) )
    else:
        print('Loading spike data...')        

    files = get_data_file_list(path, label)
    
    ## open spike files and read data
    spikes = []
    for file_name in files:
        try:
            buf = np.loadtxt('%s/%s' % (path,file_name),skiprows=skip_rows) ## load spike file while skipping the header            
            if buf.shape[0]>0:
                if buf.shape == (2,):
                    buf = np.reshape(buf, (1,2))  ## needs to be reshaped to 2-dimensional array in case there is only a single row
                spikes += [buf] 
        except:
            print('Error: %s' % sys.exc_info()[1])
            print('Remove non-numeric entries from file %s (e.g. in file header) by specifying (optional) parameter "skip_rows".\n' % (file_name))

    if len(spikes)>1:
        spikes = np.concatenate(spikes)
    elif len(spikes)==1:
        spikes = np.array(spikes[0])
    elif len(spikes)==0:
        spikes = np.array([])
        
    spike_dict = {}
    
    if spikes.shape == (0,):
        print("WARNING: No spikes contained in %s/%s*." % (path,label))
        spike_dict['senders'] = np.array([])
        spike_dict['times'] = np.array([])        
    else:
        ## extract spikes in specified time interval
        if time_interval != None:
            if type(time_interval) == tuple:
                ind = (spikes[:,1]>=time_interval[0]) * (spikes[:,1]<=time_interval[1]) 
                spikes = spikes[ind,:]
            else:
                print("Warning: time_interval must be a tuple or None. All spikes are loaded.")

        if type(pop) == list:
            spikes_subset = []
            for cn,nid in enumerate(pop):  ## loop over all neurons
                print("Spike extraction from %d/%d (%d%%) neurons completed" % (cn+1, len(pop), 1.*(cn+1)/len(pop)*100), end = '\r')            
                ind = np.where(spikes[:,0] == nid)[0]
                spikes_subset += list(spikes[ind,:])
            spikes = np.array(spikes_subset)
        elif pop == None:
            pass
        else:
            print("Warning: pop must be a list, or None. All spikes are loaded.")
        print()

        spike_dict['senders'] = spikes[:,0]
        spike_dict['times'] = spikes[:,1]

        ind = np.argsort(spike_dict['times'])

        spike_dict['senders'] = spike_dict['senders'][ind]
        spike_dict['times'] = spike_dict['times'][ind]
    
    return spike_dict

#####################################################################################
def get_data_file_list(path, label):
    '''
    (copied from https://github.com/INM-6/microcircuit-PD14-model/blob/main/PyNEST/src/microcircuit/helpers.py

    Searches for files with extension "*.dat" in directory "path" with names starting with "label", 
    and returns list of file names.

    Arguments
    ---------
    path:           str
                    Path to folder containing spike files.

    label:          str
                    Spike file label (file name root).

    Returns
    -------
    files:          list(str)
                    List of file names


    '''
 
    ## get list of files names
    files = []
    for file_name in os.listdir(path):
        if file_name.endswith('.dat') and file_name.startswith(label):
            files += [file_name]
    files.sort()
    
    assert len(files)>0 ,'No files of type "%s*.dat" found in path "%s".' % (label,path)

    return files

#####################################################################################
def spike_raster(rate_coding):
    '''
    Create a figure containing spike raster and instantantous firing rate. 
    
    Parameters
    ----------
    rate_coding : bool
        If True, spikes in raster plot are colored according to the instantaneous population rate.
        If False, all spikes are colored black.
    '''

    ## plotting parameters
    markersize = 1       ## size of the markers in the scatter plot
    clim = (0.6,0.)      ## colors (grayscale) representing minimum and maximum rate
    dt = 1               ## bin size (in ms) for rate estimation

    fig_max_width = 7.5  ## figure width (inch; see PLoS CB; https://journals.plos.org/ploscompbiol/s/figures)
    fig_size = (fig_max_width, 3.0)
    font_size = 8
    font_family = 'sans-serif'
    dpi = 300
    usetex = True

    ######################################
    print('loading spike data...')
    T = 10000. ## (ms)        
    #T = 5000. ## (ms)    
    t_warmup = 500. ## (ms)

    import json
    with open('example_data/data/nodes.json', 'r') as f:
        nodes = json.load(f)

    pops = ['L23E', 'L23I', 'L4E', 'L4I', 'L5E', 'L5I', 'L6E', 'L6I']

    senders = []
    times = []
    N = 0
    for pop in pops:
        N += len(nodes['%s' % pop])
        spikes = load_spike_data('example_data/data', label='spike_recorder-%s' % nodes['spike_recorder_%s' % pop][0])
        senders += spikes['senders'].tolist()
        times += (spikes['times'] - t_warmup).tolist()  ## remove warm-up period

    times = np.array(times)
    senders = np.array(senders)
    
    #print('total number of neurons: N = %d' % N)
    #print(len(np.unique(senders)))

    idx = np.where(times<T)[0]
    times = times[idx]
    senders = senders[idx]
    
    ######################################
    # print('generating mock-up spike data...')
    # N = 100000
    # T = 1000. ## (ms)
    # rate = 1. ## (1/s)
    # n_points =  int(N*rate*T*1e-3)

    # np.random.seed(42)
    # times = np.sort(np.random.uniform(low=0, high=T, size=n_points))
    # senders = np.round(np.random.uniform(low=1, high=N, size=n_points))

    ######################################
    
    if not rate_coding:
        clim = (0.,0.)

    ######################################
    print('creating spike raster plot...')

    from matplotlib import rcParams
    rcParams['figure.figsize']    = fig_size
    rcParams['figure.dpi']        = dpi
    rcParams['font.family']       = font_family
    rcParams['font.size']         = font_size
    rcParams['text.usetex']       = usetex

    plt.figure(1)
    plt.clf()
    gs = gridspec.GridSpec(2, 1,height_ratios = [3,1],left=0.1, right=0.97, bottom=0.17, top=0.95, hspace=0.15)

    ####################
    ax = plt.subplot(gs[0,0])
    rate, ts = spike_raster_dense(times, senders, ax, T, markersize, clim, dt)
    plt.xlim(0,T)
    plt.ylim(1,N)
    #plt.ylim(0,N+1)
    plt.setp(plt.gca(),xticklabels=[])
    plt.ylabel('neuron id')

    ####################
    ax = plt.subplot(gs[1,0])
    dt = ts[1]-ts[0]
    plt.bar(ts[:-1], rate/N, width=dt, align='edge', color='k', edgecolor='k')
    plt.xlim(0,T)
    plt.xlabel('time (ms)')
    plt.ylabel('rate (1/s)')

    ####################
    print('saving figures...')
    os.system('mkdir -p figures')
    if rate_coding:
        fname_root = 'spike_raster_rate_coded'
    else:
        fname_root = 'spike_raster_black'

    plt.savefig('figures/%s.pdf' % fname_root)
    plt.savefig('figures/%s.eps' % fname_root)

######################################
if __name__ == '__main__':

    spike_raster(rate_coding=True)
    
    spike_raster(rate_coding=False)    

    
