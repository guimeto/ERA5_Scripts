# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr

model='ERA5_T2m_1h'

yi = 2020
yf = 2020
#########################################################
t_in = 'J:/REANALYSES/ERA5/T2m_1h/'
max_out = 'J:/REANALYSES/ERA5/Tmax_daily_Outaouais/' 
min_out = 'J:/REANALYSES/ERA5/Tmin_daily_Outaouais/' 
m_out = 'J:/REANALYSES/ERA5/T2m_monthly/' 
for year in range(yi,yf+1):
    for i in range (3,4,1):        
        data = t_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)        
        
        ds = ds - 273.15  # convert from meter to mm 
        ds_min = ds.groupby('time.day').min('time')
        ds_max = ds.groupby('time.day').max('time')
        
        ds_mean = ds.groupby('time.day').mean('time')
        ds_mean = ds_mean.mean('day')
        ds_mean.to_netcdf(m_out + 'Monthly_Mean_T2m_QC_'+str(year) +'{:02d}'.format(i)+'.nc')
        
        ds_min.to_netcdf(min_out + 'ERA5_Outaouais_daily_tmin_'+str(year) +'{:02d}'.format(i)+'.nc')
        ds_max.to_netcdf(max_out + 'ERA5_Outaouais_daily_tmax_'+str(year) +'{:02d}'.format(i)+'.nc')

