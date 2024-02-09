'''
Functions for calculating spike-train power spectra and an illustrating example.

(Tom Tetzlaff; January 2024)

'''

import numpy
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import rcParams

###########################################################################################
def create_poisson_point_proc_realisation(rate,duration):
     '''
     Creates a single realisation of a (temporal) homogeneous Poisson point process.

     Parameters:
     -----------

     rate:         float
                   Rate of the Poisson process (1/s)

     duration:     float
                   Duration of the Poisson process (s)

     Returns:
     --------

     spike_times: ndarray(float)
                  Array of event times (s)

     '''

     n_spikes = numpy.random.poisson(rate*duration)
     spike_times = numpy.sort(numpy.random.uniform(low=0,high=duration,size=n_spikes))
     return spike_times

###########################################################################################
def create_gamma_point_proc_realisation(rate,order,duration):
     '''
     Creates a single realisation of a (temporal) homogeneous Gamma point process with integer order using
     the decimation procedure:

     1) generate a homogeneous Poisson point process of rate k*r
     2) copy each k'th event of 1)

     The intervals of the resulting process are sums of k iid random variables with exponential pdf. 
     The pdf of these sums is a Gamma distribution of integer order k (aka Erlangian distribution; 
     see, e.g., Sec.2.2. eq. (10), in [1]). The rate of the resulting decimated process is r.

     [1] Cox, D. R., & Lewis, P. A. (1966). The statistical analysis of series of events. Springer Dordrecht

     Parameters:
     -----------

     rate:         float
                   Rate of the Gamma process (1/s)

     order:        int
                   Order of the Gamma process

     duration:     float
                   Duration of the Gamma process (s)

     Returns:
     --------

     spike_times:  ndarray(float)
                   Array of event times (s)

     '''

     spike_times = create_poisson_point_proc_realisation(order*rate,duration)
     spike_times = spike_times[::order]  ## decimation method
     return spike_times

###########################################################################################
def create_gamma_point_proc(rate,order,duration,size):
     '''
     Creates an ensemble of independent realisations of a homogeneous Gamma point process with integer order
     (see create_gamma_point_proc_realisation).

     Parameters:
     -----------

     rate:         float
                   Rate of the Gamma process (1/s)

     order:        int
                   Order of the Gamma process

     duration:     float
                   Duration of the Gamma process (s)

     size:         int
                   Ensemble size (number of realisations)

     Returns:
     --------

     spike_times:  list(ndarray(float))
                   List (of size `size`) of arrays of event times (s)

     '''
     
     spike_times = []
     for i in range(size):
          spike_times += [create_gamma_point_proc_realisation(rate,order,duration)]

     return spike_times

###########################################################################################
def spike_times_2_spike_counts(spike_times, t_start, t_stop, binsize):
    '''
    Calculates the spike count signal for a given spike train, within a specified time interval (in disjoint bins).

    Parameters:
    -----------

    spike_times:  ndarray(float)
                  Array of event times (s)

    t_start:      float
                  Start time (lower bound of the time interval; s)

    t_stop:       float
                  Stop time (upper bound of the time interval; s)

    binsize:      float
                  Bin size (s)

    Returns:
    --------

    spike_counts: ndarray(int)
                  Array of spike counts, 
                  i.e., number n_k of events in time interval [t_k,t_k+1) 
                  for k = 0, ... , floor((t_stop-t_start/binsize))-1

    times:        ndarray(float)
                  Array of times t_k (s)

    '''
    
    times = numpy.arange(t_start, t_stop+binsize, binsize)
    spike_counts = numpy.histogram(spike_times,times)[0]

    return spike_counts, times[:-1]

###########################################################################################
def spike_psd(spike_times, t_start, t_stop, binsize, windowlength):
    '''
    Calculates the spike-train power spectral density using the Welch method. 

    Note 1: The length of each data segment is specified by `windowlength`. The overlap between segments
    is given by `windowlength/binsize/2`. As window function the Hann function is used (see scipy.signal.welch).

    Note 2: The spectrum is normalised such that its high-frequency limit is given by the spike rate.

    Parameters:
    -----------

    spike_times:  ndarray(float)
                  Array of event times (s)

    t_start:      float
                  Start time (lower bound of the time interval; s)

    t_stop:       float
                  Stop time (upper bound of the time interval; s)

    binsize:      float
                  Bin size (s)

    windowlength: float
                  Length (duration) of each segment used for PSD estimation (s)

    Returns:
    --------

    P:            ndarray(float)
                  Array of power densities (1/s)

    freqs:        ndarray(float)
                  Array of frequencies (Hz)

    fmin:         float
                  Minimal frequency (frequency resolution; Hz)

    fmax:         float
                  Maximal frequency (Hz)

    '''
    
    spike_times = spike_times[spike_times>=t_start]
    spike_times = spike_times[spike_times<=t_stop]    
     
    spike_counts, times = spike_times_2_spike_counts(spike_times, t_start, t_stop, binsize)    
    samplingfreq = 1./binsize    
    freqs, P = scipy.signal.welch(spike_counts/binsize, fs= samplingfreq, scaling='density', nperseg=windowlength/binsize)
    P *= 0.5  ## to ensure psd converges at `rate` for high frequencies; TODO: Why 1/2?

    fmin = 1/windowlength
    fmax = samplingfreq/2.

    return P, freqs, fmin, fmax

###########################################################################################
def spike_time_list_2_gdf(spike_trains):
     '''
     Transforms spike trains given by list of spike-time lists into gdf format (senders, times).

     Parameters:
     -----------

     spike_trains: list(ndarray(float))
                   List (of size `size`) of arrays of spike times (s)

     Returns:
     --------

     times:   list(float)
              List of spike times

     senders: list(int)
              List of spike-train ids

     '''

     size = len(spike_trains)
     
     senders = []
     times = []
     for i in range(size):
          times += [list(spike_trains[i])]
          senders += [i] * len(spike_trains[i]) 
     times = numpy.concatenate(times)     

     return times, senders

###########################################################################################
def psd_gamma_process_theoretical(freqs,rate,order):
     '''
     Theoretical power spectrum of a homogeneous Gamma process with integer order gamma
     (see, e.g., eqs.(3.11) and (3.12) in [1]).
     
     [1] Tetzlaff et al. (2008). Dependence of neuronal correlations on filter characteristics 
         and marginal spike train statistics. Neural Computation, 20(9), 2133-2184.

     Parameters:
     -----------

     freqs:        ndarray(float)
                   Array of frequencies (Hz)

     rate:         float
                   Rate of the Gamma process (1/s)

     order:        int
                   Order of the Gamma process

     Returns:
     --------
     P:            ndarray(float)
                   Power spectrum (1/s)/

     '''

     ## characteristic function (Fourier transform) of interval density
     P1 = (order * rate / (order*rate + 2.j*numpy.pi*freqs))**order  

     ## power spectrum (of a renewal process)
     P = rate * numpy.real(( (1-P1)**(-1) + (1-numpy.conjugate(P1))**(-1) - 1.0 ))

     return P

###########################################################################################
def example():
     '''
     Example demonstrating how to calculate spike power spectra using the Welch method.

     As test data, we use an ensemble of independent Gamma point-process realisations.

     We calculate and plot the ensemble averaged power spectrum of these processes, 
     as well as the power spectrum of the compound process, i.e., the process resulting 
     from merging all spike trains into a single spike train.

     Parameters:
     -----------
     -

     Returns:
     --------
     -

     '''

     ####################################################
     ## spike-train parameters
     size           = 10000                  ## number of spikes trains
     rate           = 10.0                 ## spike rate (spikes/s)
     order          = 50                   ## order of gamma process (int)
     duration       = 12.                  ## total duration of spike trains (s)
     seed           = 123                  ## rng seed for generation of spike trains

     ## analysis parameters
     t_start        = 1.                   ## observation start time (s)
     t_stop         = 11.                  ## observation stop time (s)
     binsize        = 2e-3                 ## bin size for generating spike counts (s)
     windowlength   = (t_stop-t_start)/5.  ## length of each segment during psd calculation (s)
     windowlength_c = (t_stop-t_start)/5.  ## length of each segment during compound psd calculation (s)

     df_theo = 0.1                         ## frequency resolution for theoretcial spectrum (Hz)
     
     ####################################################
     
     ## creating example spike trains
     numpy.random.seed(seed)
     spike_trains = create_gamma_point_proc(rate,order,duration, size)

     ## population averaged power spectrum
     P = 0
     for i in range(size):
          P_buf,freqs,fmin,fmax = spike_psd(spike_trains[i], t_start, t_stop, binsize, windowlength)
          P += P_buf
     P /= size

     ## theoretical spectrum of the gamma process
     freqs_theo = numpy.arange(fmin,fmax+df_theo,df_theo)     
     P_theo = psd_gamma_process_theoretical(freqs_theo,rate,order)
     
     ## power spectrum of compound process
     compound_spike_train = numpy.sort(numpy.concatenate(spike_trains))     
     P_c,freqs_c,fmin_c,fmax_c = spike_psd(compound_spike_train, t_start, t_stop, binsize, windowlength_c)

     ## theoretical spectrum of superposition of uncorrelated gamma processes
     P_c_theo = size * P_theo
      
     print()
     print("minimal frequency = %.3f Hz" % fmin)
     print("maximal frequency = %.3f Hz" % fmax)
     print()
     print("minimal frequency (compound PSD) = %.3f Hz" % fmin_c)
     print("maximal frequency (compound PSD) = %.3f Hz" % fmax_c)
     print()
     
     ####################################################
     ## plotting

     rcParams['font.family'] = 'serif'
     rcParams['font.size'] = 8
     rcParams['figure.dpi'] = 300
     rcParams['figure.figsize'] = (5,8)
     rcParams['text.usetex'] = True
     plt.rcParams['axes.titley'] = 1.05 

     plt.figure(1)
     plt.clf()

     ### raster plot

     plt.subplot(311)
     times, senders = spike_time_list_2_gdf(spike_trains)
     plt.plot(times, senders, 'k.', ms=1, mfc='k', mew = 0, alpha=0.1, rasterized=True)
     plt.vlines(t_start, 0,size, color = '0.5', ls = '--', lw = 2)
     plt.vlines(t_stop, 0,size, color = '0.5', ls = '--', lw = 2)
     plt.text(t_start,-0.05*size,r'$t_\mathsf{start}$', horizontalalignment='center', verticalalignment='center')
     plt.text(t_stop,-0.05*size,r'$t_\mathsf{stop}$', horizontalalignment='center', verticalalignment='center')
     plt.xlim(0,duration)
     plt.setp(plt.gca(),yticks=[])
     plt.xlim(0,duration)
     plt.ylim(0,size)
     plt.xlabel(r'time (s)')
     plt.ylabel(r'spike train id')
     plt.title(r'\parbox{\linewidth}{\centering Ensemble of Gamma point-process realisations\\[0.5ex] \tiny rate: %.1f\,spikes/s, order: %d, duration: %.1f\,s, size: %d}' % (rate, order, duration, size))

     ### ensemble averaged PSD

     plt.subplot(312)
     plt.plot(freqs,P          ,   'k', lw=2,            label=r'empirical')
     plt.plot(freqs_theo,P_theo,   '0.5', lw=3, alpha=0.6, label=r'theoretical')
     plt.legend(loc=1)
     plt.xlabel(r'frequency $f$ (Hz)')
     plt.ylabel(r'spike-train PSD (1/s)')
     plt.xlim((fmin,10*rate))
     plt.title(r'\parbox{\linewidth}{\centering Ensemble averaged power spectrum\\[0.5ex] \tiny $t_\mathsf{start}= %.1f$\,s, $t_\mathsf{stop}= %.1f$\,s, bin size: %.1f\,ms, window length: %.1f\,s, $f_\mathsf{min}=%.1f$\,Hz, $f_\mathsf{max}=%.1f$\,Hz}' % (t_start, t_stop, binsize*1e3, windowlength,fmin,fmax))

     ### compound PSD

     plt.subplot(313)
     plt.plot(freqs_c,P_c,         'k', lw=2,            label=r'empirical')
     plt.plot(freqs_theo,P_c_theo, '0.5', lw=3, alpha=0.6, label=r'theoretical')
     plt.legend(loc=1)
     plt.xlabel(r'frequency $f$ (Hz)')
     plt.ylabel(r'spike-train PSD (1/s)')
     plt.xlim((fmin,10*rate))
     plt.title(r'\parbox{\linewidth}{\centering Power spectrum of compound process\\[0.5ex] \tiny $t_\mathsf{start}=%.1f$\,s, $t_\mathsf{stop}=%.1f$\,s, bin size: %.1f\,ms, window length: %.1f\,s, $f_\mathsf{min}=%.1f$\,Hz, $f_\mathsf{max}=%.1f$\,Hz}' % \
               (t_start, t_stop, binsize*1e3, windowlength,fmin_c,fmax_c))

     plt.subplots_adjust(left=0.13,right=0.95,bottom=0.05,top=0.93,hspace=0.45)

     plt.savefig("spike_psd.pdf")

###########################################################################################

if __name__ == '__main__':
    example()
