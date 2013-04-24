#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------
#   Filename:  PyNASTF.py
#   Purpose:   PyNASTF (Python Neighbourhood Algorithm STF) main program 
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import ConfigParser
import glob
from ncdump import ncdump
import numpy as np
from obspy import read
from obspy.core.util import locations2degrees
from obspy.core.util import gps2DistAzimuth
from obspy.core.util import kilometer2degrees
import os
import sys
import time

from SNR_calculator import SNR_calculator
from quake_handler import quake_info,pdata_reader
from util_PyNASTF import * 


########################################################################
############################# Main Program #############################
########################################################################

def PyNASTF(**kwargs):
    """
    PyNASTF: Python Neighbourhood Algorithm STF
    """ 
    #----------------------------- input handler -----------------------------
    config = ConfigParser.RawConfigParser()
    inputpath = 'in.na.cfg'

    class input_handler:
        def __init__(self, inputpath):
            self.inpath = inputpath
            self.config = config.read(os.path.join(os.getcwd(), self.inpath))
            self.event_address = config.get('General', 'event_address')
            self.remote_address = config.get('General', 'remote_address')
            self.network = config.get('General', 'network')
            self.station = config.get('General', 'station')
            self.location = config.get('General', 'location')
            self.channel = config.get('General', 'channel')
            self.filter = eval(config.get('General', 'filter'))
            self.lfreq = float(config.get('General', 'lfreq'))
            self.hfreq = float(config.get('General', 'hfreq'))
            self.resample = eval(config.get('General', 'resample'))
            self.sampling_rate = int(config.get('General', 'sampling_rate'))
            self.min_dist = eval(config.get('General', 'min_dist'))
            self.max_dist = eval(config.get('General', 'max_dist'))
            self.bg_model = config.get('General', 'bg_model')
            self.SNR_limit = eval(config.get('General', 'SNR_limit'))
            self.plot_ph_no = eval(config.get('General', 'plot_phase_noise'))
            self.map = eval(config.get('General', 'map'))
            self.plot_azi = eval(config.get('General', 'plot_azi'))

    # create the input class
    inp = input_handler(inputpath)

    # modifying the input objects in a way to be usable in the next steps!
    inp.network = inp.network.split(',')
    for _i in xrange(len(inp.network)):
        inp.network[_i] = inp.network[_i].strip()
    
    inp.channel = inp.channel.split(',')
    for _i in xrange(len(inp.channel)):
        inp.channel[_i] = inp.channel[_i].strip()
    
    # s_tb: Signal time before, s_ta: Signal time after
    # n_tb: Noise Time before, n_ta: Noise time after
    s_tb=-3; s_ta=9
    n_tb=-150; n_ta=-30
    
    ev_name, ev_lat, ev_lon, ev_dp, ev_date = \
            pdata_reader(address = inp.event_address, remote_address=inp.remote_address)
    for _i in range(len(ev_name)):
        ev_name[_i] = os.path.join(ev_name[_i], 'BH')
    for ev_enum in xrange(len(ev_name)):
        e_add = ev_name[ev_enum]
        print '\n==========='
        print 'Event %s/%s: \n%s' %(ev_enum+1, len(ev_name), e_add)
        print '==========='
        if not os.path.isdir(os.path.join(e_add.split('/')[-2], 'infiles')): 
            os.makedirs(os.path.join(e_add.split('/')[-2], 'infiles')) 
        metadata = []
        msg_header = 'Event information; Lat, Lon, Depth\n'
        msg_header += '%.6f %.6f %.6f\n' %(ev_lat[ev_enum],ev_lon[ev_enum], ev_dp[ev_enum])
        msg_p = 'P-wave data ' + 17*'*' + '\n'
        msg_sh = 'SH-wave data ' + 17*'*' + '\n'
        all_p_data = []; all_sh_data = []
        all_sta_add = glob.glob(os.path.join(e_add, '*.*.*.*'))
        all_sta_add.sort()
        print len(all_sta_add)
        for sta_add in all_sta_add:
            print .,
            try:
                tr = read(sta_add)[0]
            except Exception, e:
                print e
                continue
            if not inp.network == ['*']:
                if not tr.stats.network in inp.network: continue
            if not tr.stats.channel in inp.channel: continue
            #epi_dist_prev = locations2degrees(tr.stats.sac.evla, tr.stats.sac.evlo,
            #            tr.stats.sac.stla, tr.stats.sac.stlo)
            epi_km = gps2DistAzimuth(tr.stats.sac.evla, tr.stats.sac.evlo,
                        tr.stats.sac.stla, tr.stats.sac.stlo)[0]
            epi_dist = kilometer2degrees(epi_km/1000.)
            # XXX for testing!
            #epi_dist = tr.stats.sac.gcarc
            if not inp.min_dist<=epi_dist<=inp.max_dist: continue
            if 'Z' in tr.stats.channel:
                tr_tw = time_window(tr, model=inp.bg_model)
                ph_arr = tr_tw.arr_time(epi_dist, req_phase='P')
                # XXX for testing
                #ph_arr = tr.stats.sac.t0
                if ph_arr == -12345.0: continue
                tr = preproc(tr, filter=inp.filter, hfreq=inp.hfreq, lfreq=inp.lfreq, 
                            resample=inp.resample, sampling_rate=inp.sampling_rate)
                SNR, l1_noise, l2_noise, p_data, flag_exist = \
                        SNR_calculator(tr, ev_date[ev_enum], 
                        ph_arr, s_tb=s_tb, s_ta=s_ta, n_tb=n_tb, n_ta=n_ta, method='squared',
                        plot_ph_no=inp.plot_ph_no,
                        address=os.path.join(e_add.split('/')[-2], 'infiles'))
                if not flag_exist: continue
                if SNR < inp.SNR_limit: continue
                innastats_str = '%s\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n' %(tr.stats.station,
                            tr.stats.sac.stla, tr.stats.sac.stlo, ph_arr+s_tb, 
                            SNR, l1_noise, l2_noise)
                az, ba = azbackaz(tr)
                all_p_data.append([tr.stats.station, tr.stats.location, SNR, az, 
                                        p_data, innastats_str]) 
            if 'N' in tr.stats.channel:
                try:
                    tr_E = read(sta_add[:-1] + 'E')[0]
                except Exception, e:
                    print 'Cannot read: \n%s' %(sta_add[:-1] + 'E')
                    continue
                tr = preproc(tr, filter=inp.filter, hfreq=inp.hfreq, lfreq=inp.lfreq, 
                            resample=inp.resample, sampling_rate=inp.sampling_rate)
                tr_E = preproc(tr_E, filter=inp.filter, hfreq=inp.hfreq, lfreq=inp.lfreq, 
                            resample=inp.resample, sampling_rate=inp.sampling_rate)
                tr_sh = tr.copy() 
                az, tr_sh.data = rotater(tr, tr_E)
                if not az: continue 
                tr_tw = time_window(tr_sh, model=inp.bg_model)
                ph_arr = tr_tw.arr_time(epi_dist, req_phase='S')
                if ph_arr == -12345.0: continue
                SNR, l1_noise, l2_noise, sh_data, flag_exist = \
                        SNR_calculator(tr_sh, ev_date[ev_enum], 
                        ph_arr, s_tb=s_tb, s_ta=s_ta, n_tb=n_tb, n_ta=n_ta, method='squared',
                        plot_ph_no=inp.plot_ph_no,
                        address=os.path.join(e_add.split('/')[-2], 'infiles'))
                if not flag_exist: continue
                if SNR < inp.SNR_limit: continue
                innastats_str = '%s\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n' %(tr_sh.stats.station,
                            tr_sh.stats.sac.stla, tr_sh.stats.sac.stlo, ph_arr+s_tb, 
                            SNR, l1_noise, l2_noise)
                all_sh_data.append([tr_sh.stats.station, tr_sh.stats.location, SNR, az, 
                                            sh_data, innastats_str])
        
        
        all_p_data = station_selector(all_p_data)
        all_sh_data = station_selector(all_sh_data)
        
        # sort the stations by azimuth
        all_p_data.sort(key=lambda x: x[3])
        all_sh_data.sort(key=lambda x: x[3])

        msg = msg_header
        print '\n**********************'
        print 'info:'
        print '%s  P-traces\n%s  SH-traces' %(len(all_p_data), len(all_sh_data))
        print '**********************'
        msg += '%s  P-traces\n%s  SH-traces\n' %(len(all_p_data), len(all_sh_data))
        msg += msg_p 
        for _i in xrange(len(all_p_data)): msg += all_p_data[_i][-1]
        msg += msg_sh 
        for _i in xrange(len(all_sh_data)): msg += all_sh_data[_i][-1]
        innastats_open = open(os.path.join(e_add.split('/')[-2], 
                                                'infiles', 'in.na.stats'), 'w')
        innastats_open.write(msg)
        innastats_open.close()
       
        ncdump(os.path.join(e_add.split('/')[-2], 'infiles'), 
                            all_p_data, all_sh_data)
        
        if inp.map: mapper(all_p_data, all_sh_data, 
                            address=os.path.join(e_add.split('/')[-2], 'infiles')) 
        if inp.plot_azi: plot_azi(all_p_data, all_sh_data, 
                            address=os.path.join(e_add.split('/')[-2], 'infiles')) 
########################################################################
########################################################################
########################################################################

def main():
    
    t1_pro = time.time()
    status = PyNASTF()
    print "\n============================="
    print "* Total time %f sec" %(time.time() - t1_pro)
    print "============================="
    # pass the return of main to the command line.
    sys.exit(status)

if __name__ == "__main__":
    main()


#TRASH
#s_arr = tr_tw.arr_time(epi_dist, req_phase='S')
#if s_arr == -12345.0: continue
#SNR_S = SNR_calculator(tr, events[0]['datetime'], s_arr, s_tb=-3, s_ta=9, 
#                    n_tb=-150, n_ta=-30, method='squared')
#if SNR_S < inp.SNR_limit: continue
#m = Basemap(projection='aeqd', 
#              lon_0=all_p_data[0][-1].stats.sac.evlo, 
#              lat_0=all_p_data[0][-1].stats.sac.evla,
#              resolution='c')

#targ_add = locate(root=inp.event_address, target='BH')
 
# XXX the next three lines should be activated for obspyDMT
# for pdata_processed or for psdata this should not work!
#events, e_add_par = \
#    quake_info(address=os.path.join(e_add, os.path.pardir),
#        target = 'info')
