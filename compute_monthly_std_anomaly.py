# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 21:14:11 2020

@author: guillaume
"""

import xarray as xr 
####https://uoftcoders.github.io/studyGroup/lessons/python/cartography/lesson/

# lecture de la serie complete
file = 'J:/REANALYSES/ERA5/T2m_monthly/Monthly_Mean_T2m_QC_'
for month in ['01','02','03','04']:
    multi_file = [f'{file}{year}{month}.nc' for year in range(1979,2020,1)]
    ds_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
    
    data_all = ds_all.variables['t2m'][:].squeeze()
    
    
    # lecture de la serie pour calculer la clim et la deviation
    multi_file = [f'{file}{year}{month}.nc' for year in range(1990,2020,1)]
    
    ds_clim = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
    
    data_clim = ds_clim.variables['t2m'][:].mean("time")
    data_std = ds_clim.variables['t2m'][:].std("time")
    
    stand_anomalies = xr.apply_ufunc(
        lambda x, m, s: (x - m) / s,
        ds_all.groupby("time"),
        data_clim,
        data_std,
    )
    
    stand_anomalies.to_netcdf('J:/REANALYSES/ERA5/Anomaly_stand/ERA5_T2m_Anomaly_Stand_1979-2019_vs_1990-2019_'+month+'.nc')
    







