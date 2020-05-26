import xarray as xr
from netCDF4 import Dataset
import netCDF4 as nc


for year in range(1981,2020):
    for month in range(1,13):
                
        infile = 'J:/REANALYSES/ERA5/SWE_Outaouais/ERA5-Land_SDWE_'+str(year) +'{:02d}'.format(int(month))+'_BV.nc' 
        outfile = 'J:/REANALYSES/ERA5/SWE_Outaouais/ERA5-Land_SDWE_'+str(year) +'{:02d}'.format(int(month))+'_BV2.nc' 
        nc_Modc=xr.open_dataset(infile)
        if ((year % 4) == 0) & (i == 12): # on supprime le 31 décembre par les années bisextiles
            nc_Modf=Dataset(infile,'r')
            data = nc_Modc.sel(time=~((nc_Modc.time.dt.month == 12) & (nc_Modc.time.dt.day == 31)))
            data.to_netcdf(outfile)
        else:
            nc_Modf=Dataset(infile,'r')
            data = nc_Modc
            data.to_netcdf(outfile)
            
       


  
  
  
  
