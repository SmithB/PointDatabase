# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 22:03:38 2019

@author: ben
"""


import argparse
from osgeo import gdal, gdalconst, osr
import numpy as np

def demBounds(demFile, native=True, EPSG=4326, proj4=None):
    """
        Get the extent of a DEM file
    """   
    ds=gdal.Open(demFile, gdalconst.GA_ReadOnly)
    
    proj=ds.GetProjection()
    band=ds.GetRasterBand(1)
    GT=ds.GetGeoTransform()
    # ii and jj are the pixel center coordinates.  0,0 in GDAL is the upper-left
    # corner of the first pixel.
    ii=np.array([0, band.XSize-1])+0.5
    jj=np.array([0, band.YSize])-0.5
    x=GT[0]+GT[1]*ii
    y=GT[3]+GT[5]*jj    
    if native:
        return np.array([np.min(x), np.max(x)]), np.array([np.min(y), np.max(y)])
    # calculate the projection from latlon to the DEM CS
    outRef = osr.SpatialReference()
    if proj4 is None:
        outRef.ImportFromEPSG(EPSG)
    else:
        try:
            proj4=proj4.decode('utf-8')
        except AttributeError:
            pass
        outRef.ImportFromProj4(proj4)
    demRef = osr.SpatialReference()
    demRef.ImportFromWkt(proj)
    if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
        demRef.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    xform = osr.CoordinateTransformation(demRef, outRef)
    xyOut=np.array(xform.TransformPoints( np.c_[x, y, np.zeros_like(x)]))[:,0:2]
    ds=None
    return np.array([np.min(xyOut[:,0]), np.max(xyOut[:,0])]), np.array([np.min(xyOut[:,1]), np.max(xyOut[:,1])])
