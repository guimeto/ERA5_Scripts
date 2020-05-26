# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr
model='era5_t2m_1h'

yi = 2020
yf = 2020
#########################################################
t_in = 'J:/REANALYSES/ERA5/T2m_1h_Outaouais/'
t_out = 'J:/REANALYSES/ERA5/Tmin_daily_Outaouais/' 

for year in range(yi,yf+1):
    for i in range (3,4,1):         
        data = t_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)
        ds = ds - 273.15  # convert from Kelvin to Celcius 
        
        daily_t2m = ds.resample(time = '1D').min()  
        daily_t2m.to_netcdf(t_out + 'ERA5_Outaouais_daily_tmin_'+str(year) +'{:02d}'.format(i)+'.nc')







