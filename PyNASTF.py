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
import os
import sys
import time

from SNR_calculator import SNR_calculator
from quake_handler import quake_info
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
            self.network = config.get('General', 'network')
            self.station = config.get('General', 'station')
            self.location = config.get('General', 'location')
            self.channel = config.get('General', 'channel')
            self.sampling_rate = int(config.get('General', 'sampling_rate'))
            self.min_dist = eval(config.get('General', 'min_dist'))
            self.max_dist = eval(config.get('General', 'max_dist'))
            self.bg_model = config.get('General', 'bg_model')
            self.SNR_limit = eval(config.get('General', 'SNR_limit'))
            self.map = eval(config.get('General', 'map'))

    # create the input class
    inp = input_handler(inputpath)

    # modifying the input objects in a way to be usable in the next steps!
    inp.network = inp.network.split(',')
    for _i in xrange(len(inp.network)):
        inp.network[_i] = inp.network[_i].strip()

    inp.channel = inp.channel.split(',')
    for _i in xrange(len(inp.channel)):
        inp.channel[_i] = inp.channel[_i].strip()

    targ_add = locate(root=inp.event_address, target='BH')

    for ev_enum in xrange(len(targ_add)):
        e_add = targ_add[ev_enum]
        events, e_add_par = \
            quake_info(address=os.path.join(e_add, os.path.pardir),
                target = 'info')
        print '------------------'
        print '%s/%s' %(ev_enum+1, len(targ_add))
        print 'address: %s' %(e_add)
        metadata = []
        msg_header = 'Event information; Lat, Lon, Depth\n'
        msg_header += '%.6f %.6f %.6f\n' %(events[0]['latitude'],
                                  events[0]['longitude'], events[0]['depth'])
        msg_p = 'P-wave data ' + 17*'*' + '\n'
        msg_sh = 'SH-wave data ' + 17*'*' + '\n'
        p_traces = []; sh_traces = []
        all_p_data = []; all_sh_data = []
        all_sta_add = glob.glob(os.path.join(e_add, '*.*.*.*'))
        for sta_add in all_sta_add:
            tr = read(sta_add)[0]
            if not tr.stats.network in inp.network: continue
            if not tr.stats.channel in inp.channel: continue
            epi_dist = locations2degrees(tr.stats.sac.evla, tr.stats.sac.evlo,
                        tr.stats.sac.stla, tr.stats.sac.stlo)
            if not inp.min_dist<=epi_dist<=inp.max_dist: continue
            if 'Z' in tr.stats.channel:
                tr_tw = time_window(tr, model=inp.bg_model)
                ph_arr = tr_tw.arr_time(epi_dist, req_phase='P')
                if ph_arr == -12345.0: continue
                tr.resample(inp.sampling_rate)
                SNR, l1_noise, l2_noise, p_data = SNR_calculator(tr, events[0]['datetime'], 
                        ph_arr, s_tb=-3, s_ta=9, n_tb=-150, n_ta=-30, method='squared')
                if SNR < inp.SNR_limit: continue
                p_traces.append('%s\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n' %(tr.stats.station,
                            tr.stats.sac.stla, tr.stats.sac.stlo, ph_arr, 
                            SNR, l1_noise, l2_noise))
                az, ba = azbackaz(tr)
                all_p_data.append([tr.stats.station, tr.stats.location, SNR, az, p_data]) 
            if 'N' in tr.stats.channel:
                try:
                    tr_E = read(sta_add[:-1] + 'E')[0]
                except Exception, e:
                    print 'Cannot read: \n%s' %(sta_add[:-1] + 'E')
                    continue
                tr_sh = tr.copy() 
                az, tr_sh.data = rotater(tr, tr_E) 
                tr_tw = time_window(tr_sh, model=inp.bg_model)
                ph_arr = tr_tw.arr_time(epi_dist, req_phase='S')
                if ph_arr == -12345.0: continue
                tr.resample(inp.sampling_rate)
                SNR, l1_noise, l2_noise, sh_data = SNR_calculator(tr_sh, events[0]['datetime'], 
                        ph_arr, s_tb=-3, s_ta=9, n_tb=-150, n_ta=-30, method='squared')
                if SNR < inp.SNR_limit: continue
                sh_traces.append('%s\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n' %(tr_sh.stats.station,
                            tr_sh.stats.sac.stla, tr_sh.stats.sac.stlo, ph_arr, 
                            SNR, l1_noise, l2_noise))
                all_sh_data.append([tr_sh.stats.station, tr_sh.stats.location, SNR, az, sh_data])
        
        
        all_p_data = station_selector(all_p_data)
        all_sh_data = station_selector(all_sh_data)
              
        msg = msg_header
        print '%s  P-traces\n%s  SH-traces\n' %(len(p_traces), len(sh_traces))
        msg += '%s  P-traces\n%s  SH-traces\n' %(len(p_traces), len(sh_traces))
        msg += msg_p 
        for item in p_traces: msg += item
        msg += msg_sh 
        for item in sh_traces: msg += item
        if not os.path.isdir(os.path.join('NASTF-INPUT', e_add.split('/')[-2])): 
            os.makedirs(os.path.join('NASTF-INPUT', e_add.split('/')[-2])) 
        innastats_open = open(os.path.join('NASTF-INPUT', 
                                 e_add.split('/')[-2], 'in.na.stats'), 'w')
        innastats_open.write(msg)
        innastats_open.close()
       
        ncdump(os.path.join('NASTF-INPUT', e_add.split('/')[-2]), 
                            all_p_data, all_sh_data)
        
        if inp.map: mapper(all_p_data, all_sh_data, 
                            address=os.path.join('NASTF-INPUT', e_add.split('/')[-2])) 
        
########################################################################
########################################################################
########################################################################

def main():
    
    t1_pro = time.time()
    status = PyNASTF()
    print "\n=================================================="
    print "* Total time %f sec" %(time.time() - t1_pro)
    print "=================================================="
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


