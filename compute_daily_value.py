# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr

model='ERA5'
variable = 'PR'

yi = 1979
yf = 1979
#########################################################
rep_data='K:/DATA/REANALYSES/ERA5/'

for year in range(yi,yf+1):
    for i in range (8,9,1):  
        
        data = rep_data + model + '_'+variable+'_'+str(year) +'{:02d}'.format(i)+'.nc'
        ds = xr.open_mfdataset(data)
        
        ds=ds*1000  # convert from meter to mm 
        
        daily_data = ds.resample(time = '1D').mean()       
        daily_data.to_netcdf('Daily_Mean_PR_'+str(year) +'{:02d}'.format(i)+'.nc')
           
       