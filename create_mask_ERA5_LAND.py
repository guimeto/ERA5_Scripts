# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:06:25 2019

@author: guillaume
"""
import xarray as xr
import matplotlib.pylab as plt
import warnings; warnings.filterwarnings(action='ignore')
import numpy as np
from osgeo import ogr
import geopandas as gpd

ds = xr.open_mfdataset('J:/REANALYSES/ERA5/SWE_ERA5_Land_daily/ERA5-Land_SDWE_198101_daymean.nc4')
lat_bnd = [50, 43]
lon_bnd = [270, 300]
ds = ds.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd),)
ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')



lat1d =  ds['latitude'].values
lon1d =  ds['longitude'].values
lon2d, lat2d = np.meshgrid(lon1d, lat1d)
shapes = gpd.read_file("./shapefile/BV_Outaouais_masque.shp")

tmpWGS84 = shapes.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})
tmpWGS84.loc[0, 'geometry']
tmpWGS84.to_file('Outaouais_WGS84.shp')

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


mask=get_mask(lon2d,lat2d,shp_path="Outaouais_WGS84.shp")
np.save('Outaouais.npy',mask)
print(u'Terminé')

data = ds['sd'][1].where(mask==1)/ds['sd'][1].where(mask==1)

data.plot()

data.to_netcdf('Outaouais_ERA5_LAND_Grid.nc')


