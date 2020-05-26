# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:05:43 2020

@author: guillaume
"""
import xarray as xr
model='ERA5_T2m_1h'

yi = 1990
yf = 2020
#########################################################
t_in = 'J:/REANALYSES/ERA5/T2m_1h/'
#daily_out = 'J:/REANALYSES/ERA5/T2m_daily/' 
min_out = 'J:/REANALYSES/ERA5/Month_tasmin/' 
max_out = 'J:/REANALYSES/ERA5/Month_tasmax/' 

for year in range(yi,yf+1):
    print(year)
    for i in range (1,12,1):        
        data = t_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)        
        ds = ds - 273.15  # convert from Kelvin to Celcius 
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')
        
        lat_bnd = [84, 40]
        lon_bnd = [-148, -50]
        ds = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
        
        #daily_mean = ds.groupby('time.day').mean('time')
        #daily_mean.to_netcdf(daily_out + 'Daily_Mean_T2m_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')
        daily_min = ds.groupby('time.day').min('time')
        daily_max = ds.groupby('time.day').max('time')
        
        monthly_min = daily_min.mean('day')
        monthly_max = daily_max.mean('day')
        
        monthly_min.to_netcdf(min_out + 'ERA5_Monthly_Mean_Tasmin_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')
        monthly_max.to_netcdf(max_out + 'ERA5_Monthly_Mean_Tasmax_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')
        
        #monthly_mean = daily_mean.mean('day')
        #monthly_mean.to_netcdf(monthly_out + 'Monthly_Mean_T2m_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')