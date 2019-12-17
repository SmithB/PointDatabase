#!/usr/bin/env python3
#
# PURPOSE: Read and reproject an ICESat-2 tiled mask HDF5 file
#    for one northern or southern hemisphere area with all 
#    longitudes (0 - 359.95). Returns logical of land surf_type.
#
# FILES ACCESSED: H5 surf_type file (input)
#
# COMMENTS:
#
#    Usage: surftype surf_type_filename hemisphere_flag 
#
#
# HISTORY:
#
#   YYYY-MM-DD AUID      SCM   Comment
#   ---------- --------- ----- ------------------------------------------------
#   2019-11-19 bjelley   M0265 Initial version for masking ATL11 tiles by surf_type
#  
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import sh
import sys
import os
import time
import datetime
import h5py
import numpy as np
import h5py
from osgeo import osr
import matplotlib
import matplotlib.pyplot as plt
  
PGE_NAME='create_surfmask'
PGE_VERS='Version 1.0'
PGE_INFO='Reads a composite HDF5 file, extracts "land" surf_type.'
#
# Error/Status constants
#
GE_NOERROR=0
GE_NOTICE=1
GE_WARNING=2
GE_FATAL=3
#
# Execution time
#
proc_start=time.time()
proc_end=time.time()
#
#==============================================================================
#
# NAME: msg
#
# PURPOSE: Prints a message in ASAS format
#
# FILES ACCESSED: stdout
#
# COMMENTS:
#
#------------------------------------------------------------------------------
#
def msg(i_res, mod, routine, msg):
#
# Create a timestamp (pulled from asas_common.py)
#
  tstamp=datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
  s_res='{:0>6}'.format(i_res)
  if (i_res==GE_NOERROR):
    mstatus="Status "
  if (i_res==GE_NOTICE):
    mstatus="Notice "
  if (i_res==GE_WARNING):
    mstatus="Warning"
  if (i_res==GE_FATAL):
    mstatus="ERROR  "
  print(tstamp+' | '+mstatus+' | '+s_res+' | '+routine+' | '+msg)

  return
#enddef
#
#==============================================================================
#
# NAME: end_banner
#
# PURPOSE: Prints an ending banner
#
# FILES ACCESSED: stdout
#
# COMMENTS:
#
#------------------------------------------------------------------------------
#
def end_banner(i_res):
#
# Write ending banner
#
  proc_end=time.time()
  msg(GE_NOERROR, PGE_NAME, 'main', '---')
  if (i_res == GE_NOERROR) :
    msg(GE_NOERROR, PGE_NAME, 'main', "Successful execution")
  else:
    msg(GE_NOERROR, PGE_NAME, 'main', "Execution failed")
  msg(GE_NOERROR, PGE_NAME, 'main', "Execution time: "+str(proc_end-proc_start))
  msg(GE_NOERROR, PGE_NAME, 'main', "Result Code : "+str(i_res))
  msg(GE_NOERROR, PGE_NAME, 'main', '---')
  if (i_res == GE_FATAL):
    sys.exit(GE_FATAL)

  return

#enddef
#
#==============================================================================
#
# NAME: ibits
#
# PURPOSE: Return integer value of bits
#
# COMMENTS:
#
#------------------------------------------------------------------------------
def ibits(ival, ipos, ilen):
    """Same usage as Fortran ibits function."""
    ones = ((1 << ilen)-1)
    return (ival & (ones << ipos)) >> ipos
#
#==============================================================================
#
# NAME: read_hdf5
#
# PURPOSE: Read and reproject a HDF5 of surf_type tiles.
#
# FILES ACCESSED: Input- Composite HDF5; Return - Hemisphere specific surf_type
#
# COMMENTS:
#
#------------------------------------------------------------------------------
#
class landmask(object):
  def __init__(self):
    self.x=None
    self.y=None
    self.z=None

  def read_surftype_h5(self,in_file,hemisphere=-1):
#
# Open file and read tile attributes
#
    f_in = h5py.File(in_file, 'r')
    g_in = f_in['TILE_INDEX']
    lon0 = g_in.attrs['LON0']
    lon1 = g_in.attrs['LON1']
    lat0 = g_in.attrs['LAT0']
    lat1 = g_in.attrs['LAT1']
    lon_scale = g_in.attrs['LON_SCALE']
    lat_scale = g_in.attrs['LAT_SCALE']
    tile_name = g_in.attrs['NAME']
    nlon = g_in.attrs['N_LON']
    nlat = g_in.attrs['N_LAT']
#
# Read 1 tile, need dtype for initialization
#
    tile = np.array(f_in[tile_name[0]])
#
# Initialize arrays
#
    xsz = int(np.ceil((np.max(lon1) - np.min(lon0)) / lon_scale[0]))
    ysz = int(np.ceil((np.max(lat1) - np.min(lat0)) / lat_scale[0]))
  
    surf_type = np.empty([ysz,xsz],dtype=tile.dtype)
    lons = np.full([ysz,xsz],np.inf)
    lats = np.full([ysz,xsz],np.inf)
    sum_tile = 0
#
# Loop through all tiles, capturing geolocation and tile data
#
    for lat_tile in range(0,9):
      for lon_tile in range(0,18):
        tile_num = lon_tile * (lat_tile + 1)
        if sum_tile >= np.size(nlon):
          print("breaking loop for beyond array size, size(nlon):",np.size(nlon))
          break
        x = np.linspace(lon0[sum_tile], lon1[sum_tile] - lon_scale[sum_tile], nlon[sum_tile])
        x = np.repeat(x[np.newaxis,:],nlon[sum_tile],0)
        lons[lat_tile*nlat[0]:lat_tile*nlat[0]+nlat[0],lon_tile*nlon[0]:lon_tile*nlon[0]+nlon[0]] = x
        y = np.linspace(lat0[sum_tile], lat1[sum_tile] - lat_scale[sum_tile], nlat[sum_tile])
        y = np.repeat(y[:,np.newaxis],nlat[sum_tile],1)
        lats[lat_tile*nlat[0]:lat_tile*nlat[0]+nlat[0],lon_tile*nlon[0]:lon_tile*nlon[0]+nlon[0]] = y
  
        tile = np.array(f_in[tile_name[sum_tile]])
        surf_type[lat_tile*nlat[0]:lat_tile*nlat[0]+nlat[0],lon_tile*nlon[0]:lon_tile*nlon[0]+nlon[0]] = tile
  
        sum_tile+=1
#
# Resample to every 5th in surf_type (no need for 5m resolution)
#
    subset_size = 3
    lons = lons[0:np.shape(lons)[0]:subset_size,0:np.shape(lons)[1]:subset_size]
    lats = lats[0:np.shape(lats)[0]:subset_size,0:np.shape(lats)[1]:subset_size]
    surf_type = surf_type[0:np.shape(surf_type)[0]:subset_size,0:np.shape(surf_type)[1]:subset_size]
    xsz = int(xsz/subset_size)
    ysz = int(ysz/subset_size)
    surfmask = np.full([ysz,xsz], False)
#
# Set mask True if land bit is set
#
    for j in range(np.shape(surf_type)[1]):
      for i in range(np.shape(surf_type)[0]):
        if int(ibits(surf_type[i,j],0,1)) == 1:
          surfmask[i,j] = True
#
# Subset down to where lats match hemisphere
#  Assumes all longitudes included
#  np.reshape(arr, (-1,int(360/dx))) provides for unknown number of lats
#
    dx=lon_scale[0]
    latlimit=-60.0
    if hemisphere==-1:
      lons = np.reshape(lons[lats <= latlimit], (-1,xsz))
      surfmask = np.reshape(surfmask[lats <= latlimit], (-1,xsz))
      lats = np.reshape(lats[lats <= latlimit], (-1,xsz))
      SRS_proj4='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    else:
      lons = np.reshape(lons[lats >= np.abs(latlimit)], (-1,xsz))
      surfmask = np.reshape(surfmask[lats >= np.abs(latlimit)], (-1,xsz))
      lats = np.reshape(lats[lats >= np.abs(latlimit)], (-1,xsz))
      SRS_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    latsize = np.shape(lats)[0]
    lonsize = np.shape(lons)[1]
  
#
# Reproject to the polar stereographic grid
#
    x, y, z = landmask.reproj(lons, lats, surfmask, SRS_proj4)
    x = np.reshape(x,(-1,lonsize))
    y = np.reshape(y,(-1,lonsize))
    z = np.reshape(z,(-1,lonsize))
    self.x = x
    self.y = y
    self.z = z  
    return self
#
#==============================================================================
#
  def reproj(longitude, latitude, z_in, SRS_proj4):
    out_srs=osr.SpatialReference()
    out_srs.ImportFromProj4(SRS_proj4)
  
    ll_srs=osr.SpatialReference()
    ll_srs.ImportFromEPSG(4326)
    if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
      ll_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct=osr.CoordinateTransformation(ll_srs, out_srs).TransformPoint
  
    x, y, z= list(zip(*[ct(*xy) for xy in zip(np.ravel(longitude), np.ravel(latitude), np.ravel(z_in*1.0))]))
    return np.array(x), np.array(y), np.array(z)
#
#
#-BEGIN_PROLOG-----------------------------------------------------------------
#
# NAME: main
#
# PURPOSE: Creates a composite HDF5 file from NetCDF tiles
#
# FILES ACCESSED: HDF
#
# COMMENTS:
#
#-END_PROLOG-------------------------------------------------------------------
#
if __name__ == '__main__':
  in_file='/Volumes/Data/asas/anc_data/anc12/surfmask_20180608_001_01.h5'
  hemisphere=-1
  surfmask = landmask()
  surfmask.read_surftype_h5(in_file, hemisphere)
  x = surfmask.x
  end_banner(GE_NOERROR)
  sys.exit(GE_NOERROR)

