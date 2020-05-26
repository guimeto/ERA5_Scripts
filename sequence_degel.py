import gc 
import xarray as xr
import pandas as pd
import numpy as np
from netCDF4 import Dataset
# ouverture du masque
mask = xr.open_mfdataset('Outaouais_ERA5_Grid.nc')


# nombre de degré « horaires » de dégel, i.e. cumul de toutes valeurs > 0°C  november to march
year = 2020
file = 'J:/REANALYSES/ERA5/Tmax_daily_Outaouais/ERA5_Outaouais_daily_tmax_'
multi_file = [f'{file}{year-1}{month}.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}.nc' for month in ['01','02','03']])
ds_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
ds_all = ds_all.assign_coords(longitude=(((ds_all.longitude + 180) % 360) - 180)).sortby('longitude') 
count_t2m_0deg = ds_all.where(ds_all.t2m >= 0).count("time") 

ds_all.time[-1].values
x = pd.to_datetime(ds_all.time[-1].values)

last_day = str(x.date().day)

# compute total precipiation from november to march
for yi in range(1981,2020,1):
    gc.collect()
    # lecture de la serie complete    
    file = 'J:/REANALYSES/ERA5/Tmax_daily_Outaouais/ERA5_Outaouais_daily_tmax_'
    multi_file = [f'{file}{yi-1}{month}.nc' for month in ['11','12']]
    multi_file = multi_file + ([f'{file}{yi}{month}.nc' for month in ['01','02','03']])

    ds_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
    ds_all = ds_all.assign_coords(longitude=(((ds_all.longitude + 180) % 360) - 180)).sortby('longitude')
    count_year = ds_all.where(ds_all.t2m >= 0).count("time")  

    count_year.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Count_Tmax_'+str(yi)+'.nc')
    
    
# compute climatologie
# lecture de la serie pour calculer la clim et la deviation
file = 'J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Count_Tmax_'
multi_file = [f'{file}{yi}.nc' for yi in range(1990,2020,1)]

ds_clim = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
data_clim = ds_clim.mean("time")

brute_anomalies = xr.apply_ufunc(
    lambda x, c: x- c,
    count_t2m_0deg,
    data_clim
)



brute_anomalies = brute_anomalies.where(mask.tp >= 0)
data_clim=data_clim.where(mask.tp >= 0)
count_t2m_0deg=count_t2m_0deg.where(mask.tp >= 0)

count_t2m_0deg.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_count_tmax.nc')
brute_anomalies.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_tmax_0deg_brute_Anomaly_'+str(year)+'_vs_1990-2019.nc')
data_clim.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_tmax_0deg_brute_climatology_1990-2019.nc')


# debut des graphiques

import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
## Lecture du fichier 
filename='J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_count_tmax.nc'
nc_fid=Dataset(filename,'r')
data=nc_fid.variables['t2m'][:].squeeze()
lons=nc_fid.variables['longitude'][:].squeeze()
lats=nc_fid.variables['latitude'][:].squeeze()
lon2d, lat2d = np.meshgrid(lons, lats)
# prectot de l annee courant
fig = plt.figure(figsize=(28,16))
ax = plt.subplot(111, projection=ccrs.LambertConformal())
ax.set_extent([-82,-72,45,48])
   # ax.coastlines(resolution='110m');
ax.add_feature(cfeature.OCEAN.with_scale('50m'))      # couche ocean
ax.add_feature(cfeature.LAND.with_scale('50m'))       # couche land
ax.add_feature(cfeature.LAKES.with_scale('50m'))      # couche lac    
ax.add_feature(cfeature.BORDERS.with_scale('50m'))    # couche frontieres
ax.add_feature(cfeature.RIVERS.with_scale('50m'))     # couche rivières 
coast = cfeature.NaturalEarthFeature(category='physical', scale='10m',     # ajout de la couche cotière 
                        facecolor='none', name='coastline')
ax.add_feature(coast, edgecolor='black')

  
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

## Choisissons une colormap
cmap0=plt.cm.jet
cmap0.set_under('darkblue') ## on met en blanc les valeurs inferieures au min de clev
cmap0.set_over('darkred') ## bleu fonce pour les valeurs extremes de pluie

data_levels = np.arange(0, 100, 5.)

mm = ax.contourf(lon2d,\
                   lat2d,\
                   data,\
                   vmin=0,\
                   vmax=100, \
                   transform=ccrs.PlateCarree(),\
                   levels=data_levels,\
                   cmap=cmap0 )

data_contour = ax.contour(lon2d, lat2d, data, 
                          levels = data_levels, 
                          linewidths=2, 
                          colors='k',
                          transform = ccrs.PlateCarree())
#Plot contour labels for the heights, leaving a break in the contours for the text (inline=True)
plt.clabel(data_contour,  data_levels, inline=True, fmt='%1i', fontsize=12)
# Define gridline locations and draw the lines using cartopy's built-in gridliner:
xticks = np.arange(-150.0,-40.0,20)
yticks =np.arange(10,80,10)

fig.canvas.draw()

cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 100.1, 5.), extend='both') # 2018
#cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 400.1, 20.), extend='both') # 2019

cbar.ax.tick_params(labelsize=20) 

string_title=u'Nombre de jours de dégel  Tmax >= 0°C du 1er novembre '+str(year-1)+' au '+ str(last_day) + ' mars '+str(year) 
plt.title(string_title, size='xx-large')
plt.savefig('./figures_update/Degel_BV_'+str(year)+'.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()

import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
## Lecture du fichier 
filename='J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_tmax_0deg_brute_climatology_1990-2019.nc'
nc_fid=Dataset(filename,'r')
data=nc_fid.variables['t2m'][:].squeeze()
lons=nc_fid.variables['longitude'][:].squeeze()
lats=nc_fid.variables['latitude'][:].squeeze()
lon2d, lat2d = np.meshgrid(lons, lats)
# prectot de l annee courant
fig = plt.figure(figsize=(28,16))
ax = plt.subplot(111, projection=ccrs.LambertConformal())
ax.set_extent([-82,-72,45,48])
   # ax.coastlines(resolution='110m');
ax.add_feature(cfeature.OCEAN.with_scale('50m'))      # couche ocean
ax.add_feature(cfeature.LAND.with_scale('50m'))       # couche land
ax.add_feature(cfeature.LAKES.with_scale('50m'))      # couche lac    
ax.add_feature(cfeature.BORDERS.with_scale('50m'))    # couche frontieres
ax.add_feature(cfeature.RIVERS.with_scale('50m'))     # couche rivières 
coast = cfeature.NaturalEarthFeature(category='physical', scale='10m',     # ajout de la couche cotière 
                        facecolor='none', name='coastline')
ax.add_feature(coast, edgecolor='black')

  
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

## Choisissons une colormap
cmap0=plt.cm.jet
cmap0.set_under('darkblue') ## on met en blanc les valeurs inferieures au min de clev
cmap0.set_over('darkred') ## bleu fonce pour les valeurs extremes de pluie
data_levels = np.arange(0, 100, 5.)

mm = ax.contourf(lon2d,\
                   lat2d,\
                   data,\
                   vmin=0,\
                   vmax=100, \
                   transform=ccrs.PlateCarree(),\
                   levels=data_levels,\
                   cmap=cmap0)
data_contour = ax.contour(lon2d, lat2d, data, 
                          levels = data_levels, 
                          linewidths=2, 
                          colors='k',
                          transform = ccrs.PlateCarree())
#Plot contour labels for the heights, leaving a break in the contours for the text (inline=True)
plt.clabel(data_contour,  data_levels, inline=True, fmt='%1i', fontsize=12)

# Define gridline locations and draw the lines using cartopy's built-in gridliner:
xticks = np.arange(-150.0,-40.0,20)
yticks =np.arange(10,80,10)

fig.canvas.draw()

cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 100.1, 5.), extend='both')
cbar.ax.tick_params(labelsize=20) 

string_title=u'Climatologie (1990-2019) du nombre de jours de dégel Tmax > 0°C du 1er novembre au '+ str(last_day) + ' mars' 
plt.title(string_title, size='xx-large')
plt.savefig('./figures_update/Clim_Degel_BV_1990-2019.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()


# Anoamlies brutes
## Lecture du fichier 
filename='J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_tmax_0deg_brute_Anomaly_'+str(year)+'_vs_1990-2019.nc'
nc_fid=Dataset(filename,'r')
data=nc_fid.variables['t2m'][:].squeeze()
lons=nc_fid.variables['longitude'][:].squeeze()
lats=nc_fid.variables['latitude'][:].squeeze()
lon2d, lat2d = np.meshgrid(lons, lats)
# prectot de l annee courant
fig = plt.figure(figsize=(28,16))
ax = plt.subplot(111, projection=ccrs.LambertConformal())
ax.set_extent([-82,-72,45,48])
   # ax.coastlines(resolution='110m');
ax.add_feature(cfeature.OCEAN.with_scale('50m'))      # couche ocean
ax.add_feature(cfeature.LAND.with_scale('50m'))       # couche land
ax.add_feature(cfeature.LAKES.with_scale('50m'))      # couche lac    
ax.add_feature(cfeature.BORDERS.with_scale('50m'))    # couche frontieres
ax.add_feature(cfeature.RIVERS.with_scale('50m'))     # couche rivières 
coast = cfeature.NaturalEarthFeature(category='physical', scale='10m',     # ajout de la couche cotière 
                        facecolor='none', name='coastline')
ax.add_feature(coast, edgecolor='black')

  
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

## Choisissons une colormap
cmap0=plt.cm.jet
cmap0.set_under('darkblue') ## on met en blanc les valeurs inferieures au min de clev
cmap0.set_over('darkred') ## bleu fonce pour les valeurs extremes de pluie
data_levels = np.arange(-20, 20.1, 2.)

mm = ax.contourf(lon2d,\
                   lat2d,\
                   data,\
                   vmin=-20,\
                   vmax=20, \
                   transform=ccrs.PlateCarree(),\
                   levels=data_levels,\
                   cmap=cmap0 )

data_contour = ax.contour(lon2d, lat2d, data, 
                          levels = data_levels, 
                          linewidths=2, 
                          colors='k',
                          transform = ccrs.PlateCarree())

#Plot contour labels for the heights, leaving a break in the contours for the text (inline=True)
plt.clabel(data_contour,  data_levels, inline=True, fmt='%1i', fontsize=12)
# Define gridline locations and draw the lines using cartopy's built-in gridliner:
xticks = np.arange(-150.0,-40.0,20)
yticks =np.arange(10,80,10)

fig.canvas.draw()

cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(-20, 20.1, 2.), extend='both')
cbar.ax.tick_params(labelsize=20) 

string_title=u'Anomalie brute du nombre de jours de dégel avec Tmax > 0°C du 1er novembre au '+ str(last_day) + ' mars '+str(year) 
plt.title(string_title, size='xx-large')
plt.savefig('./figures_update/Anomalie_Brute_Degel_BV_'+str(year)+'_vs_1990-2019.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()

