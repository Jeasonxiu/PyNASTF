from netCDF4 import Dataset
import os

def ncdump(address, all_p_data, all_sh_data, format = 'NETCDF4'):
    """
    Creating and writing a netCDF file out of the data (P and SH)
    the netCDF file has:
    three dimensions (samples, stat_P and stat_SH)
    two variables (data_P and data_SH) in f8 format
    """
    rootgrp = Dataset(os.path.join(address, 'in.na.data'), 'w',
                            format = format)
    rootgrp.createDimension('samples', 512)
    rootgrp.createDimension('stat_P', len(all_p_data))
    rootgrp.createDimension('stat_SH', len(all_sh_data))
    rootgrp.createVariable('data_P', 'f8', ('stat_P','samples'))
    rootgrp.createVariable('data_SH', 'f8', ('stat_SH','samples'))
    for i in xrange(len(all_p_data)):
        rootgrp.variables['data_P'][i, 0:len(all_p_data[i][-1].data)] = \
                                                all_p_data[i][-1].data
    for i in xrange(len(all_sh_data)):
        rootgrp.variables['data_SH'][i, 0:len(all_sh_data[i][-1].data)] = \
                                                all_sh_data[i][-1].data
    rootgrp.close()
