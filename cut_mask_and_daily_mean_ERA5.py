# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:26:04 2020

@author: guillaume
"""
import gc 
import xarray as xr
import pandas as pd 


m_f=xr.open_dataset('mask_ERA5_Ouganda.nc')
mask = m_f['region'].values

for year in range(1979,2020,1):
    for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:     
        gc.collect()       
        # lecture de la serie complete    
#        file = 'J:/REANALYSES/ERA5/T2m_1h/ERA5_T2m_1h_'+str(year)+month+'_sfc.nc'
#        ds1 = xr.open_mfdataset(file)
#        lat_bnd = [ 10 , -5 ]
#        lon_bnd = [ 25 , 40 ]
#        ds1 = ds1.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
#        daily_mean = ds1.groupby('time.day').mean('time')
#        data1  =  daily_mean['t2m'].where(mask == 1)
#        days = data1.shape[0]
#        basedate = pd.date_range(str(year)+'-'+month+'-01', periods=days)
#        data1['day'] = basedate   
#        
#        data1.to_netcdf('D:/Utilisateurs/guillaume/Desktop/PROJET_AFR44/ERA5/tas/ERA5_tas_ll_'+str(year)+month+'.nc') 
        
        file = 'J:/REANALYSES/ERA5/PR_1h/era5_edna_ea_'+str(year)+month+'_sfc.nc'
        ds2 = xr.open_mfdataset(file)
        lat_bnd = [ 10 , -5 ]
        lon_bnd = [ 25 , 40 ]
        ds2 = ds2.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
        ds2 = ds2 * 1000  # convert from meter to mm
        daily_sum = ds2.groupby('time.day').sum('time')
        data2  =  daily_sum['tp'].where(mask == 1)
        days = data2.shape[0]
        basedate = pd.date_range(str(year)+'-'+month+'-01', periods=days)
        data2['day'] = basedate   
        
        data2.to_netcdf('D:/Utilisateurs/guillaume/Desktop/PROJET_AFR44/ERA5/tp/ERA5_tp_ll_'+str(year)+month+'.nc') 