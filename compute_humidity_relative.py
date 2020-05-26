# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr
import numpy as np
model = 'era5_edna_ea_UU_VV_TD_TT'
modelout = 'ERA5'

yi = 1979
yf = 2019
#########################################################
pr_in = 'J:/REANALYSES/ERA5/UU_VV_TD_TT_17h/'
pr_out = 'J:/REANALYSES/ERA5/HR/' 


for year in range(yi,yf+1):
    for i in range (1,12,1):  
        
        data = pr_in + model + '_'+str(year) +'{:02d}'.format(i)+'_12h.nc'
        ds = xr.open_mfdataset(data)       
        ds2 = 100 - 5 * (ds.t2m - ds.d2m )
        
        ds3 = np.sqrt(np.square(ds.u10) + np.square(ds.v10))
             
        ds['Humidity'] = ds2
        ds['wind'] = 3.6 * ds3
        
        ds.to_netcdf(pr_out + modelout + '_TT_HR_WIND_'+str(year) +'{:02d}'.format(i)+'_12h.nc')



