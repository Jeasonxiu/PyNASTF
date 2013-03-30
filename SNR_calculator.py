import matplotlib.pyplot as plt
import math as m
import numpy as np
import os

def SNR_calculator(tr, t_orig, t_phase, s_tb=-3, s_ta=9, n_tb=-150, n_ta=-30, 
                        method='squared', plot_ph_no=True, address=None):
    """
    Signal to Noise ratio calculator
    this function calculates SNR, l1 noise level and l2 noise level
    """
    if method == 'squared':
        try:
            tr_id = '%s.%s.%s.%s' %(tr.stats.network, tr.stats.station, 
                                     tr.stats.location, tr.stats.channel) 
            phase_signal = tr.slice(t_orig+t_phase+s_tb, t_orig+t_phase+s_ta)
            phase_noise = tr.slice(t_orig+t_phase+n_tb, t_orig+t_phase+n_ta)
            if plot_ph_no:
                plt.clf()
                dt=tr.stats.delta; npts=tr.stats.npts
                t = np.linspace(0., dt*(npts-1), npts)
                plt.plot(t, tr.data, 'b')
                dt=phase_signal.stats.delta; npts=phase_signal.stats.npts
                t_start = phase_signal.stats.starttime - tr.stats.starttime
                t = np.linspace(t_start, t_start+dt*(npts-1), npts)
                plt.plot(t, phase_signal, 'r')
                dt=phase_noise.stats.delta; npts=phase_noise.stats.npts
                t_start = phase_noise.stats.starttime - tr.stats.starttime
                t = np.linspace(t_start, t_start+dt*(npts-1), npts)
                plt.plot(t, phase_noise, 'g')
                plt.title(tr_id)
                if not os.path.isdir(os.path.join(address, 'FIGS')): 
                    os.makedirs(os.path.join(address, 'FIGS')) 
                plt.savefig(os.path.join(address, 'FIGS', tr_id + '.png'))
            SNR = np.sum(np.square(phase_signal))/np.sum(np.square(phase_noise))*\
                        (float(phase_noise.stats.npts)/phase_signal.stats.npts)
            l1_noise = np.sum(np.abs(phase_noise.data))/(phase_noise.stats.npts)
            l2_noise = np.sqrt(np.sum(np.square(phase_noise.data)))/(phase_noise.stats.npts)
            if m.isnan(SNR) or m.isnan(l1_noise) or m.isnan(l2_noise):
                raise ValueError('nan')
            return SNR, l1_noise, l2_noise, phase_signal, True
        except Exception, e:
            print 'ERROR %s: %s' %(tr_id, e)
            return False, False, False, False, False
