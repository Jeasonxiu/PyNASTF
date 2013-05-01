import fnmatch
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from obspy.core.util import gps2DistAzimuth
from obspy.signal.rotate import rotate_NE_RT
from obspy.taup.taup import getTravelTimes
import os

#----------------------------- locate -----------------------------
def locate(root = '.', target = 'BH'):
    """
    Locates a subdirectory within a directory.
    """
    matches = []
    for root, dirnames, filenames in os.walk(root):
        for dirnames in fnmatch.filter(dirnames, target):
            matches.append(os.path.join(root, dirnames))
    return matches

#----------------------------- time_window -----------------------------
class time_window:
    """
    time_window class creates an object to calculate:
    - arrival time for the requested phase
    """
    def __init__(self, tr, model='iasp91'):
        self.id='%s.%s.%s.%s' %(tr.stats.network, tr.stats.station,
                                    tr.stats.location, tr.stats.channel)
        self.stats=tr.stats
        self.model=model
    def arr_time(self, epi_dist, req_phase='Pdiff'):
        tt = getTravelTimes(epi_dist, self.stats.sac.evdp, model=self.model)
        t_phase = -12345.0
        for tt_item in tt:
            if tt_item['phase_name'] == req_phase:
                t_phase = tt_item['time']
                break
        return (t_phase)

#----------------------------- preproc -----------------------------
def preproc(tr, filter=True, hfreq=1., lfreq=1./30, 
                  resample=True, sampling_rate=10):
    """
    preprocessing on the traces. Currently:
    - Filter
    - Resample
    """
    if filter:
        # Applying low-pass and high-pass filters
        #tr.filter('highpass', freq=lfreq)
        #tr.filter('lowpass', freq=hfreq)
        
        # XXX maybe better to use SAC instead?!
        tr.filter('bandpass', freqmin=lfreq, freqmax=hfreq, 
                                corners=2, zerophase=True)
    if resample:
        tr.resample(sampling_rate)
    return tr

#----------------------------- azbackaz -----------------------------
def azbackaz(tr):
    """
    calculates azimuth and backazimuth
    """
    return gps2DistAzimuth(tr.stats.sac.evla, tr.stats.sac.evlo,
                               tr.stats.sac.stla, tr.stats.sac.stlo)[1:]

#----------------------------- rotater -----------------------------
def rotater(tr_N, tr_E):
    """
    Rotates horizontal components of a seismogram.
    """
    az, ba = azbackaz(tr_N)
    try:
        return az, rotate_NE_RT(tr_N.data, tr_E.data, ba)[1]
    except Exception, e:
        print 'ERROR %s.%s.%s.%s: %s vs %s %s' %(tr_N.stats.network, 
                tr_N.stats.station, tr_N.stats.location, 
                tr_N.stats.channel, tr_N.stats.npts, tr_E.stats.npts, e)
        return False, tr_N.data

#----------------------------- station_selector -----------------------------
def station_selector(all_ph_data):
    """
    select the stations based on the SNR and azimuth
    """
    # First sort all the stations with their station names
    # check whether there is any stations with the same station ID but
    # with different location ID
    # If there is any, just keep the one with highest SNR
    all_ph_data.sort()
    for i in range(len(all_ph_data)-1, 0, -1):
        if all_ph_data[i][0] == all_ph_data[i-1][0]:
            if all_ph_data[i][2] > all_ph_data[i-1][2]:
                all_ph_data.remove(all_ph_data[i-1])
            else:
                all_ph_data.remove(all_ph_data[i])
    
    # divide the azimuth by 10 degree intervals
    # check each interaval and keep three stations
    azi_ph_all = []
    for iazi in range(0, 360, 10):
        azi_grp = []
        for item in all_ph_data: 
            if iazi <= item[3] < iazi+10:
                azi_grp.append(item)
        if azi_grp:
            if len(azi_grp) > 3:
                azi_grp.sort(key=lambda x: x[2])
                for _i in range(-3, 0):
                    azi_ph_all.append(azi_grp[_i])
            else:
                for _i in range(0, len(azi_grp)):
                    azi_ph_all.append(azi_grp[_i])
    return azi_ph_all

#----------------------------- mapper -----------------------------
def mapper(all_p_data, all_sh_data, address):
    """
    create a map out of all the stations and the event
    it uses the ortho projection
    """
    plt.clf() 
    m = Basemap(projection='ortho', 
                  lon_0=all_p_data[0][-2].stats.sac.evlo, 
                  lat_0=all_p_data[0][-2].stats.sac.evla, 
                  resolution='c')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    m.drawmapboundary()
    x_ev, y_ev = m(all_p_data[0][-2].stats.sac.evlo, 
                    all_p_data[0][-2].stats.sac.evla)
    m.scatter(x_ev, y_ev, color="red", marker="*",
                edgecolor="red", zorder=10)
    
    for _i in xrange(len(all_p_data)):
        x, y = m(all_p_data[_i][-2].stats.sac.stlo, 
                        all_p_data[_i][-2].stats.sac.stla)
        m.scatter(x, y, color="blue", marker="o",
                    edgecolor="black", zorder=10)
 
    for _i in xrange(len(all_sh_data)):
        x, y = m(all_sh_data[_i][-2].stats.sac.stlo, 
                        all_sh_data[_i][-2].stats.sac.stla)
        m.scatter(x, y, color="red", marker="o",
                    edgecolor="black", zorder=10)
            
    plt.savefig(os.path.join(address, 'na.map.png'))

#----------------------------- plot_azi -----------------------------
def plot_azi(all_p_data, all_sh_data, address):
    """
    create two plots (P and SH) for all selected waveforms sorted by azimuth
    """
    plt.clf()
    fig = plt.figure()
    for _i in xrange(len(all_p_data)):
        dt=all_p_data[_i][-2].stats.delta
        npts=all_p_data[_i][-2].stats.npts
        #t = np.linspace(0., dt*(npts-1), npts)
        t = np.linspace(0., 1, npts)
        np_cols = int(np.sqrt(len(all_p_data)))
        np_rows = len(all_p_data)/np_cols+1
        sub_plt = fig.add_subplot(np_rows, np_cols, _i+1)
        sub_plt.spines['right'].set_color('none')
        sub_plt.spines['top'].set_color('none')
        sub_plt.spines['left'].set_color('none')
        sub_plt.spines['bottom'].set_color('none')
        sub_plt.xaxis.set_visible(False)
        sub_plt.yaxis.set_visible(False)
        sub_plt.plot(t, all_p_data[_i][-2].data/abs(all_p_data[_i][-2].max())*10, 'black')
        plt.xticks([])
        plt.yticks([])
        #plt.title('%s\n%s' %(all_p_data[_i][0], all_p_data[_i][3]))
    plt.tight_layout()
    plt.savefig(os.path.join(address, 'na.P.png'))

    plt.clf()
    fig = plt.figure()
    for _i in xrange(len(all_sh_data)):
        dt=all_sh_data[_i][-2].stats.delta
        npts=all_sh_data[_i][-2].stats.npts
        #t = np.linspace(0., dt*(npts-1), npts)
        t = np.linspace(0., 1, npts)
        np_cols = int(np.sqrt(len(all_sh_data)))
        np_rows = len(all_sh_data)/np_cols+1
        sub_plt = fig.add_subplot(np_rows, np_cols, _i+1)
        sub_plt.spines['right'].set_color('none')
        sub_plt.spines['top'].set_color('none')
        sub_plt.spines['left'].set_color('none')
        sub_plt.spines['bottom'].set_color('none')
        sub_plt.xaxis.set_visible(False)
        sub_plt.yaxis.set_visible(False)
        sub_plt.plot(t, all_sh_data[_i][-2].data/abs(all_sh_data[_i][-2].max())*10, 'black')
        plt.xticks([])
        plt.yticks([])
        #plt.title('%s\n%s' %(all_sh_data[_i][0], all_sh_data[_i][3]))
    plt.tight_layout()
    plt.savefig(os.path.join(address, 'na.SH.png'))
