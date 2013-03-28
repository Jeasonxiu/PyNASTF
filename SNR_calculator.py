import numpy as np
import matplotlib.pyplot as plt

def SNR_calculator(tr, t_orig, t_phase, s_tb=-3, s_ta=9, n_tb=-150, n_ta=-30, 
                        method='squared'):
    """
    Signal to Noise ratio calculator
    this function calculates SNR, l1 noise level and l2 noise level
    """
    if method == 'squared':
        phase_signal = tr.slice(t_orig+t_phase+s_tb, t_orig+t_phase+s_ta)
        phase_noise = tr.slice(t_orig+t_phase+n_tb, t_orig+t_phase+n_ta)
        SNR = np.sum(np.square(phase_signal))/np.sum(np.square(phase_noise))
        l1_noise = np.sum(np.abs(phase_noise.data))        
        l2_noise = np.sqrt(np.sum(np.square(phase_noise.data)))
    return SNR, l1_noise, l2_noise, phase_signal

#TRASH
#dt=tr.stats.delta; npts=tr.stats.npts
#t = np.linspace(0., dt*(npts-1), npts)
#plt.plot(t, tr.data, 'b')
#dt=phase_signal.stats.delta; npts=phase_signal.stats.npts
#t_start = phase_signal.stats.starttime - tr.stats.starttime
#t = np.linspace(t_start, t_start+dt*(npts-1), npts)
#plt.plot(t, phase_signal, 'r')
#dt=phase_noise.stats.delta; npts=phase_noise.stats.npts
#t_start = phase_noise.stats.starttime - tr.stats.starttime
#t = np.linspace(t_start, t_start+dt*(npts-1), npts)
#plt.plot(t, phase_noise, 'g')
#plt.show()
