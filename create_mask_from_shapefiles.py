# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:31:05 2019

@author: guillaume
"""
import xarray as xr
from osgeo import ogr
import numpy as np
import geopandas as gpd


model='ERA5'

#########################################################
file_in = 'K:/DATA/REANALYSES/ERA5/Mean_T2m/'
file_out = 'D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Mask_Old/' 
#file_out = 'D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Mask/' 

data = file_in +'Monthly_Mean_t2m_198501.nc'
ds = xr.open_mfdataset(data)
        
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

Imp_Lats =  ds['latitude'].values
Imp_Lons =  ds['longitude'].values
lon2d, lat2d = np.meshgrid(Imp_Lons, Imp_Lats)

#shapes = gpd.read_file("D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Zones/Masque.shp")
shapes = gpd.read_file("D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp")

list(shapes.columns.values)
for name in shapes['NAME']:
    print(name)
 #   mask=get_mask(lon2d,lat2d,shp_path="D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Zones/Masque.shp", polygon_name=name)
    mask=get_mask(lon2d,lat2d,shp_path="D:/Utilisateurs/guillaume/Documents/GitHub/InSIGHT-PHAC/Countries/Countries_Final-polygon.shp", polygon_name=name)
    
    np.save(file_out + 'mask_'+str(name.replace(' ','_'))+'.npy',mask)
    ds_mask = ds.where(mask == 1)
    ds_mask.to_netcdf(file_out + 'mask_'+str(name.replace(' ','_'))+'.nc')

