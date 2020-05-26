# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 14:47:28 2020

@author: guillaume
"""
import xarray as xr
model='era5_pr_1h'

yi = 1979
yf = 2019
#########################################################
t_in = 'J:/REANALYSES/ERA5/PR_1h_Outaouais/'
pr_out = 'J:/REANALYSES/ERA5/PR_daily_Outaouais/' 
mean_out = 'J:/REANALYSES/ERA5/Month_Mean_PR_Outaouais/' 
tot_out = 'J:/REANALYSES/ERA5/Month_PrecTOT_Outaouais/' 
mask = xr.open_mfdataset('Outaouais_ERA5_Grid.nc')

for year in range(yi,yf+1):
    for i in range (12,13,1):         
        data = t_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)
        ds = ds * 1000  # convert from meter to mm  
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')
        daily_mean = ds.groupby('time.day').sum('time')
        daily_mean = daily_mean.where(mask.tp >= 0)
        
        daily_mean.to_netcdf(pr_out + 'ERA5_Outaouais_daily_pr_'+str(year) +'{:02d}'.format(i)+'.nc')
        
        monthly_mean = daily_mean.mean('day')
        monthly_sum = daily_mean.sum('day')
        
        monthly_mean.to_netcdf(mean_out + 'ERA5_Outaouais_Monthly_Mean_PR_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')
        monthly_sum.to_netcdf(tot_out + 'ERA5_Outaouais_Monthly_PrecTOT_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')









