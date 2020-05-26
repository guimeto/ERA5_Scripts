# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:26:29 2020

@author: guillaume
"""
import gc 
import xarray as xr
#file = ['J:/REANALYSES/ERA5/era5_edna_202001_1.nc']
#file = file + ['J:/REANALYSES/ERA5/era5_edna_202001_2.nc']
#file = file + ['J:/REANALYSES/ERA5/era5_edna_ea_202001_sfc.nc']
#ds_all = xr.concat([xr.open_dataset(f) for f in file], 'time')
#ds_all.to_netcdf('J:/REANALYSES/ERA5/PR_1h/era5_edna_ea_202001.nc')

#month = '01'
#year = 2020

#file = 'J:/REANALYSES/ERA5/PR_1h/era5_edna_ea_'+str(year)+month+'_sfc.nc'
#ds = xr.open_mfdataset(file)
#lat_bnd = [50, 43]
#lon_bnd = [270, 300]
#ds = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
#ds.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/era5_pr_1h_'+str(year)+month+'_sfc.nc') 
 
for year in range(2019,2020,1):
    for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:
      
        gc.collect()
        
        # lecture de la serie complete    
        file = 'J:/REANALYSES/ERA5/SWE_ERA5_Land_daily/ERA5-Land_SDWE_'+str(year)+month+'_daymean.nc4'
        ds = xr.open_mfdataset(file)
        lat_bnd = [50, 43]
        lon_bnd = [270, 300]
        ds = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
        
        ds.to_netcdf('J:/REANALYSES/ERA5/SWE_Outaouais/ERA5-Land_SDWE_'+str(year)+month+'.nc') 


#for year in range(2020,2021,1):
#    for month in ['03']:
#      
#        gc.collect()
#        
#        # lecture de la serie complete    
#        file = 'J:/REANALYSES/ERA5/PR_1h/era5_edna_ea_'+str(year)+month+'_sfc.nc'
#        ds = xr.open_mfdataset(file)
#        lat_bnd = [50, 43]
#        lon_bnd = [270, 300]
#        ds = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
#        
#        ds.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/era5_pr_1h_'+str(year)+month+'_sfc.nc') 
#
#for year in range(2020,2021,1):
#    for month in ['03']:
#        gc.collect()
##        
#        # lecture de la serie complete    
#        file = 'J:/REANALYSES/ERA5/T2m_1h/ERA5_T2m_1h_'+str(year)+month+'_sfc.nc'
#        ds = xr.open_mfdataset(file)
#        lat_bnd = [50, 43]
#        lon_bnd = [270, 300]
#        ds = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
##        
#        ds.to_netcdf('J:/REANALYSES/ERA5/T2m_1h_Outaouais/era5_t2m_1h_'+str(year)+month+'_sfc.nc') 
