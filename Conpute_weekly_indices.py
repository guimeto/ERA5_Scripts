# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr
import gc 
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
import regionmask
import geopandas as gpd
import pandas as pd

model='era5_edna_ea'
#########################################################
t_in = 'J:/REANALYSES/ERA5/T2m_1h/'
t_out = 'J:/REANALYSES/ERA5/T2m_Week/'         
data = t_in + model + '_2018*_sfc.nc'
ds = xr.open_mfdataset(data, chunks = {'time': 10})

# With the function assign_coords the longitude is converted from the 0-360 range to -180,180 
ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')

############################## PART TO MASK THE WHOLE AFRICA CONTINENT 
#http://datapages.com/gis-map-publishing-program/gis-open-files/global-framework/global-heat-flow-database/shapefiles-list

# Load the shapefile
PATH_TO_SHAPEFILE = './continent_shapefile/continent.shp'
nuts = gpd.read_file(PATH_TO_SHAPEFILE)
nuts

# CALCULATE MASK
nuts_mask_poly = regionmask.Regions_cls(name = 'nuts_mask', numbers = list(range(0,8)), names = list(nuts.CONTINENT), abbrevs = list(nuts.CONTINENT), outlines = list(nuts.geometry.values[i] for i in range(0,8)))
print(nuts_mask_poly)

mask = nuts_mask_poly.mask(ds.isel(time = 0).sel(latitude = slice(39, -38), longitude = slice(-25, 55)), lat_name='latitude', lon_name='longitude')
mask
mask.to_netcdf('./mask_all_continent.nc') 

plt.figure(figsize=(12,8))
ax = plt.axes()
mask.plot(ax = ax)
nuts.plot(ax = ax, alpha = 0.8, facecolor = 'none', lw = 1)

ID_COUNTRY = 3
print(nuts.CONTINENT[ID_COUNTRY])
lat = mask.latitude.values
lon = mask.longitude.values
sel_mask = mask.where(mask == ID_COUNTRY).values
sel_mask

id_lon = lon[np.where(~np.all(np.isnan(sel_mask), axis=0))]
id_lat = lat[np.where(~np.all(np.isnan(sel_mask), axis=1))]

out_sel = ds.sel(latitude = slice(id_lat[0], id_lat[-1]), longitude = slice(id_lon[0], id_lon[-1])).compute().where(mask == ID_COUNTRY)
out_sel


plt.figure(figsize=(12,8))
ax = plt.axes()
out_sel.t2m.isel(time = 0).plot(ax = ax)
nuts.plot(ax = ax, alpha = 0.8, facecolor = 'none')



############################## PART TO MASK a SPECIFIC COUNTRY IN AFRICA
# http://www.maplibrary.org/library/stacks/Africa/index.htm

# Load the shapefile
PATH_TO_SHAPEFILE = 'Africa.shp'
nuts = gpd.read_file(PATH_TO_SHAPEFILE)
nuts.head()

# CALCULATE MASK
nuts_mask_poly = regionmask.Regions_cls(name = 'nuts_mask', numbers = list(range(0,762)), names = list(nuts.COUNTRY), abbrevs = list(nuts.COUNTRY), outlines = list(nuts.geometry.values[i] for i in range(0,762)))
print(nuts_mask_poly)

mask = nuts_mask_poly.mask(ds.isel(time = 0).sel(latitude = slice(39, -38), longitude = slice(-25, 55)), lat_name='latitude', lon_name='longitude')
mask
mask.to_netcdf('./mask_africa.nc') 

plt.figure(figsize=(12,8))
ax = plt.axes()
mask.plot(ax = ax)
nuts.plot(ax = ax, alpha = 0.8, facecolor = 'none', lw = 1)



ID_REGION = 1
print(nuts.COUNTRY[ID_REGION])
lat = mask.latitude.values
lon = mask.longitude.values
sel_mask = mask.where(mask == ID_REGION).values
sel_mask

id_lon = lon[np.where(~np.all(np.isnan(sel_mask), axis=0))]
id_lat = lat[np.where(~np.all(np.isnan(sel_mask), axis=1))]

out_sel = ds.sel(latitude = slice(id_lat[0], id_lat[-1]), longitude = slice(id_lon[0], id_lon[-1])).compute().where(mask == ID_REGION)

out_sel


plt.figure(figsize=(12,8))
ax = plt.axes()
out_sel.t2m.isel(time = 0).plot(ax = ax)
nuts.plot(ax = ax, alpha = 0.8, facecolor = 'none')




# Mean weekly nighttime temperature surface data from October 1st 2018 to September 30th 2019 (pixel resolution around 200-300 km^2)
lat_bnd = [39, -38]
lon_bnd = [-20, 55]

DS_date_range = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),time=slice('2018-10-31', '2019-10-31'))
del ds

DS_date_range = DS_date_range - 273.15  # convert from Kelvin to Celcius 
daily_min_t2m = DS_date_range.resample(time = '1D').min()
del DS_date_range
gc.collect()
weekly_t2m = daily_min_t2m.resample(time = '1W').mean()  
del daily_min_t2m

#weekly_t2m.to_netcdf('./weekly_min_t2m.nc') 
# Average weekly temperature data to create one surface of mean annual nighttime temperature at pixel i resolution: i.e. temperature.pixel.i  
mean_weekly_t2m = weekly_t2m.t2m.mean(dim=('time'))

# Standardise and centre the pixel values: s.temperature.pixel.i = (temperature.pixel.i - mean(temperature.all.pixels)) / standard_deviation(temperature.all.pixels)
def standardize(x):
    return x.mean(dim=('time')) - (x.mean(dim=('time')).mean(dim=('latitude','longitude')) / x.mean(dim=('time')).std(dim=('latitude','longitude')))

weekly_t2m_std = weekly_t2m.apply(standardize)

logit_model_pixel_i = -12 + 4.09 * weekly_t2m_std.t2m

probability_pixel_i = np.exp((logit_model_pixel_i))/(1+np.exp((logit_model_pixel_i)))
probability_pixel_i.values
#probability_pixel_i.to_netcdf('./probability_pixel_i.nc') 

def plotMap():
    #Set the projection information
    proj = ccrs.LambertConformal(central_longitude=-97.0,central_latitude=53, standard_parallels=[53])
    #Create a figure with an axes object on which we will plot. Pass the projection to that axes.
    fig, ax = plt.subplots(subplot_kw=dict(projection=proj))
    
    #Zoom in
    #ax.set_extent([-140,-60,10,70])
    
    #Add map features
    ax.add_feature(cfeature.LAND, facecolor='0.9') #Grayscale colors can be set using 0 (black) to 1 (white)
    ax.add_feature(cfeature.LAKES, alpha=0.9)  #Alpha sets transparency (0 is transparent, 1 is solid)
    ax.add_feature(cfeature.BORDERS, zorder=10)
    ax.add_feature(cfeature.COASTLINE, zorder=10)

    #We can use additional features from Natural Earth (http://www.naturalearthdata.com/features/)
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',  name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
    ax.add_feature(states_provinces, edgecolor='gray', zorder=10)
    
    #Add lat/lon gridlines every 20° to the map
    ax.gridlines(xlocs=np.arange(0,361,20), ylocs=np.arange(-80,90,20)) 
    
    return fig, ax


fig, ax = plotMap()
## Choisissons une colormap
cmap0=plt.cm.jet
#cmap0.set_under('w') ## on met en blanc les valeurs inferieures au min de clev
#cmap0.set_over('darkblue')

tt_contour = ax.contourf(probability_pixel_i.longitude.values, 
                         probability_pixel_i.latitude.values, 
                         logit_model_pixel_i.values, zorder=2,  
                          cmap=cmap0, transform = ccrs.PlateCarree())


logit_model_pixel_i.to_netcdf('./logit_model_pixel_i.nc') 

