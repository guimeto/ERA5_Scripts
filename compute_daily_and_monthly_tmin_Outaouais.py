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
t_out = 'J:/REANALYSES/ERA5/Tmin_daily_Outaouais/' 
min_out = 'J:/REANALYSES/ERA5/Month_tasmin_Outaouais/' 
mask = xr.open_mfdataset('Outaouais_ERA5_Grid.nc')

for year in range(yi,yf+1):
    for i in range (1,13,1):         
        data = t_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)
        ds = ds - 273.15  # convert from Kelvin to Celcius 
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')
        daily_min = ds.groupby('time.day').min('time')
        daily_min=daily_min.where(mask.tp >= 0)
        if ((year % 4) == 0) & (i == 12):
            daily_min.sel(day=slice(0, 30)).to_netcdf(t_out + 'ERA5_Outaouais_daily_tmin_'+str(year) +'{:02d}'.format(i)+'.nc')    
            monthly_min = daily_min.mean('day')
            monthly_min.to_netcdf(min_out + 'ERA5_Outaouais_Monthly_Mean_Tasmin_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')
        else:
            daily_min.to_netcdf(t_out + 'ERA5_Outaouais_daily_tmin_'+str(year) +'{:02d}'.format(i)+'.nc')       
            monthly_min = daily_min.mean('day')
            monthly_min.to_netcdf(min_out + 'ERA5_Outaouais_Monthly_Mean_Tasmin_CAN_'+str(year) +'{:02d}'.format(i)+'.nc')







