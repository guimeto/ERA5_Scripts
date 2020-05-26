# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr
model='era5_t2m_1h'

yi = 1979
yf = 2019
#########################################################
t_in = 'J:/REANALYSES/ERA5/T2m_1h_Outaouais/'
t_out = 'J:/REANALYSES/ERA5/Tmax_daily_Outaouais/' 
max_out = 'J:/REANALYSES/ERA5/Month_tasmax_Outaouais/' 
mask = xr.open_mfdataset('Outaouais_ERA5_Grid.nc')

for year in range(yi,yf+1):
    for i in range (1,13,1):         
        data = t_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)
        ds = ds - 273.15  # convert from Kelvin to Celcius 
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')
        daily_max = ds.groupby('time.day').max('time')
        daily_max=daily_max.where(mask.tp >= 0)
        
        
        
        if ((year % 4) == 0) & (i == 12): # on supprime le 31 décembre par les années bisextiles
            daily_max.sel(day=slice(0, 30)).to_netcdf(t_out + 'ERA5_Outaouais_daily_tmax_'+str(year) +'{:02d}'.format(i)+'.nc')    
            monthly_max = daily_max.mean('day')
            monthly_max.to_netcdf(max_out + 'ERA5_Outaouais_Monthly_Mean_Tasmax_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')
        else:
            daily_max.to_netcdf(t_out + 'ERA5_Outaouais_daily_tmax_'+str(year) +'{:02d}'.format(i)+'.nc')   
            monthly_max = daily_max.mean('day')
            monthly_max.to_netcdf(max_out + 'ERA5_Outaouais_Monthly_Mean_Tasmax_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')







