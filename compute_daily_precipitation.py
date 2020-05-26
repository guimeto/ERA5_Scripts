# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr
import calendar
import pandas as pd
model='era5_edna_ea'

yi = 1979
yf = 2019
#########################################################
pr_in = 'J:/REANALYSES/ERA5/PR_1h/'
pr_out = 'J:/REANALYSES/ERA5/Prec_daily/' 
range1 = list(range(0,18,1))
range2 = list(range(17,24,1))


for year in range(yi,yf+1):
    for i in range (1,13,1):  
        days = calendar.monthrange(year,i)[1]       
        data = pr_in + model + '_'+str(year) +'{:02d}'.format(i)+'_sfc.nc'
        ds = xr.open_mfdataset(data)        
        ds = ds * 1000  # convert from meter to mm         
        lat_bnd = [62, 43]
        lon_bnd = [276, 306]
        datasets=[]
        for d in range(1,days+1):
            if d == 1 :                
                tmp=ds.sel(time=ds.time.dt.day.isin([d])).sel(time=ds.sel(time=ds.time.dt.day.isin([d])).time.dt.hour.isin(range1)).sum(dim='time')
            elif d>1 and d<=days:
                print('jour' + str(d))
                tmp1=ds.sel(time=ds.time.dt.day.isin([d-1])).sel(time=ds.sel(time=ds.time.dt.day.isin([d-1])).time.dt.hour.isin(range2)).sum(dim='time')
                tmp2=ds.sel(time=ds.time.dt.day.isin([d])).sel(time=ds.sel(time=ds.time.dt.day.isin([d])).time.dt.hour.isin(range1)).sum(dim='time')
                tmp = tmp1+tmp2
            #elif d == days:
            #    tmp=ds.sel(time=ds.time.dt.day.isin([d])).sel(time=ds.sel(time=ds.time.dt.day.isin([1])).time.dt.hour.isin(range2)).sum(dim='time')
            datasets.append(tmp)
        combined = xr.concat(datasets, 'time')
        times =  pd.date_range(str(year)+'-'+str('{:02d}'.format(i))+'-01T17', periods=days, freq='D')     
        combined['time'] = times
        
        daily_pr = combined.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd))   
        daily_pr.to_netcdf(pr_out + 'Daily_Total_PR_17h_QC_'+str(year) +'{:02d}'.format(i)+'.nc')

