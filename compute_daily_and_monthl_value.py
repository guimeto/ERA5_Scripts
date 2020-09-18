# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr

model='ERA5'

yi = 1985
yf = 2016
#########################################################
pr_in = 'K:/DATA/REANALYSES/ERA5/PR/'
pr_out = 'K:/DATA/REANALYSES/ERA5/PrecTOT/' 

t_in = 'K:/DATA/REANALYSES/ERA5/T2m/'
t_out = 'K:/DATA/REANALYSES/ERA5/Mean_T2m/' 


for year in range(yi,yf+1):
    for i in range (1,13,1):  
        
        data = pr_in + model + '_PR_'+str(year) +'{:02d}'.format(i)+'.nc4'
        ds = xr.open_mfdataset(data)
        
        ds = ds * 1000  # convert from meter to mm 

        monthly_pr = ds.resample(time = '1M').sum()  
        monthly_pr.to_netcdf(pr_out + 'Monthly_Total_PR_'+str(year) +'{:02d}'.format(i)+'.nc')

for year in range(yi,yf+1):
    for i in range (1,13,1):  
        
        data = t_in + model + '_T2m_'+str(year) +'{:02d}'.format(i)+'.nc'
        ds = xr.open_mfdataset(data)
        
        ds = ds - 273.15  # convert from Kelvin to Celcius 
        
        monthly_t2m = ds.resample(time = '1M').mean()  
        monthly_t2m.to_netcdf(t_out + 'Monthly_Mean_t2m_'+str(year) +'{:02d}'.format(i)+'.nc')
           
        

from osgeo import ogr
import numpy as np
import geopandas as gpd
##############################################################
##############################################################
#Fonction créer par Sasha Huziy (Centre ESCER) afin d'aller récupéré sous forme de shapefile le contour d'une région géographique

def get_mask(lons2d, lats2d, shp_path="", polygon_name=None):
    """
    Assumes that the shape file contains polygons in lat lon coordinates
    :param lons2d:
    :param lats2d:
    :param shp_path:
    :rtype : np.ndarray
    The mask is 1 for the points inside of the polygons
    """
    ds = ogr.Open(shp_path)
    """
    :type : ogr.DataSource
    """

    xx = lons2d.copy()
    yy = lats2d

    # set longitudes to be from -180 to 180
    xx[xx > 180] -= 360

    mask = np.zeros(lons2d.shape, dtype=int)
    nx, ny = mask.shape

    pt = ogr.Geometry(ogr.wkbPoint)

    for i in range(ds.GetLayerCount()):
        layer = ds.GetLayer(i)
        """
        :type : ogr.Layer
        """

        for j in range(layer.GetFeatureCount()):
            feat = layer.GetFeature(j)
            """
            :type : ogr.Feature
            """

            # Select polygons by the name property
            if polygon_name is not None:
                if not feat.GetFieldAsString("name") == polygon_name:
                    continue

            g = feat.GetGeometryRef()
            """
            :type : ogr.Geometry
            """

            assert isinstance(g, ogr.Geometry)

            for pi in range(nx):
                for pj in range(ny):
                    pt.SetPoint_2D(0, float(xx[pi, pj]), float(yy[pi, pj]))

                    mask[pi, pj] += int(g.Contains(pt))

    return mask

Imp_Lats =  monthly_t2m['latitude'].values
Imp_Lons =  monthly_t2m['longitude'].values
lon2d, lat2d = np.meshgrid(Imp_Lons, Imp_Lats)

shapes = gpd.read_file("D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp")
list(shapes.columns.values)
for name in shapes['NAME']:
    mask=get_mask(lon2d,lat2d,shp_path="D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp", polygon_name=name)
    np.save('mask_'+str(name.replace(' ','_'))+'.npy',mask)
monthly_t2m_mask = monthly_t2m.where(mask == 1)
monthly_t2m_mask.t2m.mean()
monthly_t2m_mask.t2m.mean(dim=('longitude','latitude')).to_dataframe()


dataDIR = './test.nc'
monthly_t2m_mask.to_netcdf(dataDIR)

# load plotting libraries
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
ax = plt.axes(projection=ccrs.Orthographic(-80, 35))
monthly_t2m_mask.t2m.plot.contourf(ax=ax, transform=ccrs.PlateCarree());
ax.set_global(); ax.coastlines();


# choose a good projection for regional maps
proj=ccrs.LambertConformal(central_longitude=-100)

ax = plt.subplot(111, projection=proj)

monthly_t2m_mask.t2m.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree())

ax.coastlines();


monthly_t2m.plot()
shapes = gpd.read_file("D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp")
list(shapes.columns.values)
for name in shapes['NAME']:
    print(name)
# first features
shapes.head(3)
shapes['NAME'].values[10]




ds = ogr.Open(r"D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp")
lyr = ds.GetLayer()
field_names = [field.name for field in lyr.schema]
print(field_names)
lyr.REGION


ds = ogr.Open(r"D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp")
lyr = ds.GetLayer()

# Apply a SQL WHERE Clause
lyr.SetAttributeFilter("FIELD = 'Value'")

# Loop over features
for feature in lyr:
    # do stuff with features selected by attribute filter
lyr = ds.GetLayer()

# Apply a SQL WHERE Clause
lyr.SetAttributeFilter("FIELD = 'Value'")







