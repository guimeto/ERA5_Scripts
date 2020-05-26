import xarray as xr
from netCDF4 import Dataset
import netCDF4 as nc


def create_file_from_source(src_file, trg_file):
    src = nc.Dataset(src_file)
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
   # trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        if name not in tomask:
            trg.variables[name][:] = src.variables[name][:]
            
        else:    
            trg.variables[name][:] = data

    # Save the file
    trg.close()

#create 2d grid mask http://meteo.unican.es/work/xarray_seminar/xArray_seminar.html
tomask = ['sd']

m_f=xr.open_dataset('Outaouais_ERA5_LAND_Grid.nc')
lat2d=m_f.variables['latitude'][:]
lon2d=m_f.variables['longitude'][:]
mask = m_f['sd'].values

for year in range(1981,2020):
    for month in range(1,13):
        
        infile = 'J:/REANALYSES/ERA5/SWE_Outaouais/ERA5-Land_SDWE_'+str(year) +'{:02d}'.format(int(month))+'.nc'           
        outfile = 'J:/REANALYSES/ERA5/SWE_Outaouais/ERA5-Land_SDWE_'+str(year) +'{:02d}'.format(int(month))+'_BV.nc' 
       
        nc_Modc=xr.open_dataset(infile)
        if ((year % 4) == 0) & (i == 12): # on supprime le 31 décembre par les années bisextiles
            nc_Modf=Dataset(infile,'r')
            data = nc_Modc.sel(time=~((nc_Modc.time.dt.month == 12) & (nc_Modc.time.dt.day == 31)))['sd'].where(mask == 1)
            create_file_from_source(infile, outfile)
            
        else:    
            nc_Modf=Dataset(infile,'r')
            data = nc_Modc['sd'].where(mask == 1)
            create_file_from_source(infile, outfile)


  
  
  
  
