===========================================
PyNASTF: Python Neighbourhood Algorithm STF
===========================================

PyNASTF (Python Neighbourhood Algorithm STF) currently performs the following tasks:

- Filter stations by distance and network
- Load each seismogram, cut out P_noise, P_data, SH_noise, SH_data time windows
- Calculate SNR, L1_noise, L2_noise
- Filter by SNR
- Resample
- Choose three stations with the highes SNR in 10 degree azimuth interval
- Write in.na.stats
- Write data into in.na.data (NetCDF file)
