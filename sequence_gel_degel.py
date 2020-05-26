import gc 
import xarray as xr
import pandas as pd
import numpy as np
from netCDF4 import Dataset
# ouverture du masque
masque = xr.open_mfdataset('Outaouais_ERA5_Grid.nc')


# sequence de gel/degel: (Tmax>0 && Tmin<0)
year = 2020
# ouverture de Tmin
file = 'J:/REANALYSES/ERA5/Tmin_daily_Outaouais/ERA5_Outaouais_daily_tmin_'
multi_file = [f'{file}{year-1}{month}.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}.nc' for month in ['01','02','03']])
ds_min_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
ds_min_all['tmin'] = ds_min_all['t2m']
ds_min_all = ds_min_all.drop(['t2m'])
ds_min_all = ds_min_all.assign_coords(longitude=(((ds_min_all.longitude + 180) % 360) - 180)).sortby('longitude') 

# ouverture de Tmax
file = 'J:/REANALYSES/ERA5/Tmax_daily_Outaouais/ERA5_Outaouais_daily_tmax_'
multi_file = [f'{file}{year-1}{month}.nc' for month in ['11','12']]
multi_file = multi_file + ([f'{file}{year}{month}.nc' for month in ['01','02','03']])
ds_max_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
ds_max_all['tmax'] = ds_max_all['t2m']
ds_max_all = ds_max_all.drop(['t2m'])
ds_max_all = ds_max_all.assign_coords(longitude=(((ds_max_all.longitude + 180) % 360) - 180)).sortby('longitude') 

sequence = np.ones(ds_max_all.tmax.shape)
seqs = xr.Dataset(
        data_vars={'Fr_Th': (('time','latitude','longitude'), sequence)},
        coords={'longitude': ds_max_all.longitude ,
                'latitude': ds_max_all.latitude,
                'time': ds_max_all.time})

mask = (ds_max_all.tmax >= 0) & (ds_min_all.tmin <= 0)
year_seq = seqs.where(mask).sum("time") 

seqs.time[-1].values
x = pd.to_datetime(seqs.time[-1].values)
last_day = str(x.date().day)

# compute total precipiation from november to march
for yi in range(1981,2020,1):
    gc.collect()    
    # ouverture de Tmin
    file = 'J:/REANALYSES/ERA5/Tmin_daily_Outaouais/ERA5_Outaouais_daily_tmin_'
    multi_file = [f'{file}{yi-1}{month}.nc' for month in ['11','12']]
    multi_file = multi_file + ([f'{file}{yi}{month}.nc' for month in ['01','02','03']])
    ds_min_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
    ds_min_all['tmin'] = ds_min_all['t2m']
    ds_min_all = ds_min_all.drop(['t2m'])
    ds_min_all = ds_min_all.assign_coords(longitude=(((ds_min_all.longitude + 180) % 360) - 180)).sortby('longitude') 
    
    # ouverture de Tmax
    file = 'J:/REANALYSES/ERA5/Tmax_daily_Outaouais/ERA5_Outaouais_daily_tmax_'
    multi_file = [f'{file}{yi-1}{month}.nc' for month in ['11','12']]
    multi_file = multi_file + ([f'{file}{yi}{month}.nc' for month in ['01','02','03']])
    ds_max_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
    ds_max_all['tmax'] = ds_max_all['t2m']
    ds_max_all = ds_max_all.drop(['t2m'])
    ds_max_all = ds_max_all.assign_coords(longitude=(((ds_max_all.longitude + 180) % 360) - 180)).sortby('longitude') 
    
    sequence = np.ones(ds_max_all.tmax.shape)
    seqs = xr.Dataset(
            data_vars={'Fr_Th': (('time','latitude','longitude'), sequence)},
            coords={'longitude': ds_max_all.longitude ,
                    'latitude': ds_max_all.latitude,
                    'time': ds_max_all.time})
    
    mask = (ds_max_all.tmax >= 0) & (ds_min_all.tmin <= 0)
    yi_seq = seqs.where(mask).sum("time") 
    yi_seq.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Fr_Th_'+str(yi)+'.nc')



# compute climatologie
# lecture de la serie pour calculer la clim 
file = 'J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Fr_Th_'
multi_file = [f'{file}{yi}.nc' for yi in range(1990,2020,1)]
ds_clim = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')
data_clim = ds_clim.mean("time")

brute_anomalies = xr.apply_ufunc(
    lambda x, c: x- c,
    year_seq,
    data_clim
)


brute_anomalies = brute_anomalies.where(masque.tp >= 0)
data_clim=data_clim.where(masque.tp >= 0)
year_seq=year_seq.where(masque.tp >= 0)

year_seq.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_count_Fr_Th.nc')
brute_anomalies.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Fr_Th_brute_Anomaly_'+str(year)+'_vs_1990-2019.nc')
data_clim.to_netcdf('J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Fr_Th_climatology_1990-2019.nc')


# debut des graphiques

import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
## Lecture du fichier 
filename='J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_count_Fr_Th.nc'
nc_fid=Dataset(filename,'r')
data=nc_fid.variables['Fr_Th'][:].squeeze()
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

data_levels = np.arange(0, 70, 2.)

mm = ax.contourf(lon2d,\
                   lat2d,\
                   data,\
                   vmin=0,\
                   vmax=70, \
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

cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 70.1, 2.), extend='both') # 2018
#cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 400.1, 20.), extend='both') # 2019

cbar.ax.tick_params(labelsize=20) 

string_title=u'Nombre de séquences de gel/dégel Tmax >= 0°C et Tmin <=0°C du 1er novembre '+str(year-1)+' au '+ str(last_day) + ' mars '+str(year) 
plt.title(string_title, size='xx-large')
plt.savefig('./figures_update/Gel_Degel_BV_'+str(year)+'.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()

import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
## Lecture du fichier 
filename='J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Fr_Th_climatology_1990-2019.nc'
nc_fid=Dataset(filename,'r')
data=nc_fid.variables['Fr_Th'][:].squeeze()
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
data_levels = np.arange(0, 70, 2.)

mm = ax.contourf(lon2d,\
                   lat2d,\
                   data,\
                   vmin=0,\
                   vmax=70, \
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

cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 70.1,2.), extend='both')
cbar.ax.tick_params(labelsize=20) 

string_title=u'Climatologie (1990-2019) du nombre de séquences de gel/dégel Tmax >= 0°C et Tmin <=0°C du 1er novembre au '+ str(last_day) + ' mars' 
plt.title(string_title, size='xx-large')
plt.savefig('./figures_update/Clim_gel_Degel_BV_1990-2019.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()


# Anoamlies brutes
## Lecture du fichier 
filename='J:/REANALYSES/ERA5/PR_1h_Outaouais/tmp/ERA5_Fr_Th_brute_Anomaly_'+str(year)+'_vs_1990-2019.nc'
nc_fid=Dataset(filename,'r')
data=nc_fid.variables['Fr_Th'][:].squeeze()
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

string_title=u'Anomalie brute du nombre de séquences de gel/dégel Tmax >= 0°C et Tmin <=0°C du 1er novembre au '+ str(last_day) + ' mars '+str(year) 
plt.title(string_title, size='xx-large')
plt.savefig('./figures_update/Anomalie_Brute_Gel_Degel_BV_'+str(year)+'_vs_1990-2019.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()
    
    
    
    
    