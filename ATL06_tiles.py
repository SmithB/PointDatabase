# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 09:40:02 2019

@author: ben
"""
import numpy as np
import scipy.interpolate as si
import h5py
from PointDatabase.geo_index import geo_index
from PointDatabase.point_data import point_data
from PointDatabase.point_data import pt_blockmedian
from PointDatabase.xover_search import dilate_bins
from PointDatabase.mapData import mapData
from PointDatabase.ATL06_filters import segDifferenceFilter
from LSsurf.unique_by_rows import unique_by_rows
import sys
from PointDatabase.check_ATL06_blacklist import check_rgt_cycle_blacklist
from glob import glob
from PointDatabase.geo_index import geo_index, index_list_for_files

import os
import re 

def make_tile(args):
    #xy0, SRS_proj4, tile_spacing, pad, bin_W, GI_file, out_dir, field_dict, seg_diff_scale=None, blockmedian_scale=None, blacklist=None):
    xy0=args['xy0']
    SRS_proj4=args['SRS_proj4']
    tile_spacing=args['tile_spacing']
    bin_W=args['bin_W']
    GI_file=args['GI_file']
    out_dir=args['out_dir']
    field_dict=args['field_dict']
    seg_diff_scale=args['seg_diff_scale']
    blockmedian_scale=args['blockmedian_scale']    
    dxb, dyb = np.meshgrid(np.arange(-tile_spacing/2, tile_spacing/2+bin_W, bin_W), 
                           np.arange(-tile_spacing/2, tile_spacing/2+bin_W, bin_W))
    dxb=dxb.ravel()
    dyb=dyb.ravel()

    list_of_fields=[]
    for group in field_dict:
        for ds in field_dict[group]:
            list_of_fields.append(ds)
  
    gI=geo_index().from_file(GI_file, read_file=False)
    out_file=out_dir+('/E%d_N%d.h5' % (xy0[0]/1.e3, xy0[1]/1.e3))
    print("working on "+out_file)
    D=gI.query_xy((xy0[0]+dxb, xy0[1]+dyb), fields=field_dict)
    if D is None:
        return
    file_dict={}
    delete_list=[]
    for file_num, Di in enumerate(D):
        Di.get_xy(SRS_proj4)
        Di.assign({'source_file_num':np.zeros_like(Di.x, dtype=int)+file_num})        
        if seg_diff_scale is not None:  
           Di.h_li[Di.atl06_quality_summary==1]=np.NaN
           segDifferenceFilter(Di, setValid=False, toNaN=True)     
        Di.ravel_fields()
        if blockmedian_scale is not None:
            Di.index(np.isfinite(Di.h_li) & (Di.atl06_quality_summary==0))
            try:
                ind=pt_blockmedian(Di.x, Di.y, Di.h_li, blockmedian_scale, return_index=True)[3]
            except Exception:
                delete_list.append(Di)
                continue
            Di.index(ind[:,0])
        else:
            Di.index(np.isfinite(Di.h_li))
        file_dict[file_num]=Di.filename
    D_all=point_data(list_of_fields=list_of_fields+['x','y','source_file_num']).from_list(D)
    
    y_bin_function=np.round(D_all.y/bin_W)
    x_bin_function=np.round(D_all.x/bin_W)
    x_scale=np.nanmax(x_bin_function)-np.nanmin(x_bin_function)
    t_scale=np.nanmax(D_all.delta_time)-np.nanmin(D_all.delta_time)
    
    xy_bin_function=(y_bin_function-np.nanmin(y_bin_function))*x_scale+(x_bin_function-np.nanmin(x_bin_function))
    xyt_bin_function= xy_bin_function + (D_all.delta_time-np.nanmin(D_all.delta_time))/t_scale
    ind=np.argsort(xyt_bin_function)

    bin_dict={}
    xy_bin_fn_sort=xy_bin_function[ind]
    fn_delta=np.concatenate([[-1], np.where(np.diff(xy_bin_fn_sort))[0], [xy_bin_fn_sort.size]])
    for ii in range(len(fn_delta)-1):
        this_ind=ind[(fn_delta[ii]+1):(fn_delta[ii+1]+1)]
        bin_dict[(x_bin_function[this_ind[0]], y_bin_function[this_ind[0]])]=this_ind
    key_arr=np.array([key for key in bin_dict.keys()])
    key_order=np.argsort(key_arr[:,1]-np.min(key_arr[:,1])*x_scale+(key_arr[:,0]-np.min(key_arr[:,0])))
    key_arr=key_arr[key_order,:]
    for key in key_arr:
        this_group='%dE_%dN' % tuple(key*bin_W)
        D_all.subset(bin_dict[tuple(key)]).to_file(out_file, replace=False, group=this_group)
    
    with h5py.File(out_file,'r+') as h5f:
        grp=h5f.create_group("source_files")
        for key in file_dict:
            grp.attrs['file_%d' % key] = file_dict[key]
        
    
def read_tile(xy0, tile_dir,  W=None):
    tile_file=tile_dir+('/E%d_N%d.h5' % (xy0[0]/1.e3, xy0[1]/1.e3))

    D=dict()
    if W is not None:
        g=re.compile('N(.*)_E(.*).h5').search(os.path.basename(tile_file))
        x0=np.float(g.group(1))*1000
        y0=np.float(g.group(2))*1000
    with h5py.File(tile_file,'r') as h5f:
        # figure out what fields to read:
        # find a group in h5f that contains fields:
        bin_re=re.compile('.*E_.*N')
        for group in h5f:
            if bin_re.match(group) is not None:
                list_of_fields=[field for field in h5f[group]]
                if len(list_of_fields) > 0:
                    break
        # next read all the data
        D={field:[] for field in list_of_fields}
        for group in h5f:
            if bin_re.match(group) is not None:
                for field in list_of_fields:
                    D[field].append(np.array(h5f[group][field]))
        # concatenate the list of arrays in D:
        for field in D:
            D[field]=np.concatenate(D[field])
        # make D into a point_data instance:
        D=point_data(list_of_fields=[key for key in D.keys()]).from_dict(D)
    return reconstruct_tracks(D, x0=x0, y0=y0, W=W)

def reconstruct_tracks(D, x0=None, y0=None, W=None):
    '''
    sort a collection of data (D) into a list of data organized by track
    '''
    if ('cycle' in D.list_of_fields):
        _, bin_dict=unique_by_rows(np.c_[D.cycle, D.rgt, D.BP, D.LR], return_dict=True)
    else:
        _, bin_dict=unique_by_rows(np.c_[D.cycle_number, D.rgt, D.BP, D.LR], return_dict=True)
    D0=[]
    for key in bin_dict:
        if "delta_time" in D.list_of_fields:
            ind=np.argsort(D.delta_time[bin_dict[key]])
        else:
            ind=np.argsort(D.time[bin_dict[key]])
        this_D=D.subset(bin_dict[key][ind])
        if W is not None:
                this_D.index((np.abs(this_D.x-x0)<W/2) &(np.abs(this_D.y-y0)<W/2) )
        D0.append(this_D)   
    return D0

def make_queue(queue_file, hemisphere=-1):

    """
    make a queue of commands to generate the tiles for ATL06
    arguments:
        queue_file: file into which to write the commands
        hemisphere: northern or southern hemisphere, defaults to -1 (southern)
    Edits to adapt to a different machine:
        
        1. Change the locations of the masks
        2. Change the location of the executable (last line in the function)
    """
    
    tile_spacing=1.e5    
    # EDIT HERE TO SET THE MASK LOCATIONS
    if hemisphere==-1:
        mask_G=mapData().from_geotif('/Volumes/ice1/ben/MOA/moa_2009_1km.tif')
        mask_G.z=mask_G.z>100
        XR=[-2800, 2800]
        YR=[-2500, 2500]
    else:
        mask_G=mapData().from_geotif('/Volumes/ice1/ben/GimpMasks_v1.1/GimpIceMask_1km.tif')
        mask_G.z=mask_G.z==1
        XR=[-700, 1000]
        YR=[-3360, -500]
    # tile centers:
    x0, y0=np.meshgrid(np.arange(XR[0], XR[1], tile_spacing/1000)*1000, np.arange(YR[0], YR[1], tile_spacing/1000)*1000)
    mI=si.interp2d(mask_G.x, mask_G.y, mask_G.z.astype(np.float64))
    maski=mI(x0[1,:], y0[:,1]).ravel()>0.5
    x0=x0.ravel()[maski]
    y0=y0.ravel()[maski]
    xyTile=set([xy0 for xy0 in zip(x0, y0)])
    dilate_bins(xyTile, tile_spacing)
    #EDIT BELOW TO SET THE EXECUTABLE LOCATION
    with open(queue_file,'w') as qf:
        for cycle in ['01','02','03']:
            for xy in xyTile:
                qf.write('python3 /home/ben/git_repos/PointDatabase/ATL06_tiles.py %d %d %d %s\n'% (xy[0], xy[1], hemisphere, cycle))
  

def index_tiles(tile_dir_root, hemisphere):
     """
     Generate a geo_index for the tile files
     
     Arguments:
         tile_dir_root: directory under which to search for the tile files
         hemisphere: north(1) or south(-1), used to set the projection
     """
     if hemisphere==1:
         SRS_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
     else:
         SRS_proj4='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '

     files=glob(tile_dir_root+'/E*.h5')
     print("found files:")
     print(files)
     iList = index_list_for_files(files, "indexed_h5", [1.e4, 1.e4], SRS_proj4, dir_root=tile_dir_root)
     geo_index(SRS_proj4=SRS_proj4, delta=[1.e4, 1.e4]).from_list(iList, dir_root=tile_dir_root).to_file(tile_dir_root+'/GeoIndex.h5')

def index_cycle_indices(tile_dir_root, hemisphere):

    if hemisphere==1:
        SRS_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    else:
        SRS_proj4='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

    files=glob(tile_dir_root+'/cycle*/GeoIndex.h5')
    index_list=list()
    for sub_index_file in files:
        temp=geo_index().from_file(sub_index_file)
        xy=temp.bins_as_array()
        index_list.append(geo_index(delta=[1.e4, 1.e4]).from_xy(xy, filename=sub_index_file, file_type='h5_geoindex', fake_offset_val=-1))
        temp=None
    geo_index(delta=[1.e4, 1.e4], SRS_proj4=SRS_proj4).from_list(index_list).to_file(tile_dir_root+'/GeoIndex_v2.h5')
     
    
    #iList=index_list_for_files(files, "h5_geoindex", [1.e4, 1.e4], SRS_proj4)#, dir_root=tile_dir_root)
    #for index in iList:
    #    index.change_root(tile_dir_root)
    #geo_index(SRS_proj4=SRS_proj4, delta=[1.e4, 1.e4]).from_list(iList, dir_root=tile_dir_root).to_file(tile_dir_root+'/GeoIndex.h5')

    
def main():
    #make_queue('/home/ben/temp/GL_tile_queue.txt', hemisphere=1)
    #index_tiles('/Volumes/ice2/ben/scf/GL_06/tiles/001/', 1)
    if len(sys.argv)==3:
        #for dirname in glob(sys.argv[1]+'/cycle*'):
        #    index_tiles(dirname, sys.argv[2])
        index_cycle_indices(sys.argv[1], sys.argv[2])
        sys.exit(0)
    xy=(np.int(sys.argv[1]), np.int(sys.argv[2]))
    hemisphere=np.int(sys.argv[3])
    cycle=sys.argv[4]
    pad=0
    bin_W=1.e4
    tile_spacing=1.e5
    
    blockmedian_scale=None
    #Skip seg_diff_Scale (August 27: changed from 5 to None)
    seg_diff_scale=None
    out_dir='/Volumes/insar10/ben/temp/cycle_%s' % cycle
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
        
    if hemisphere==1:
        SRS_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
        GI_file='/Volumes/ice2/ben/scf/GL_06/001/cycle_%s/index/GeoIndex.h5' % cycle
    else:
        SRS_proj4='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
        #mask_G=mapData().from_geotif('/Volumes/insar5/gmap/OIB_data/AA/masks/MOA_2009_grounded.tif')
        GI_file='/Volumes/ice2/ben/scf/AA_06/001/cycle_%s/index/GeoIndex.h5' % cycle
    field_dict={None:['delta_time','h_li','h_li_sigma','latitude','longitude','atl06_quality_summary','segment_id','sigma_geo_h'], 
            'fit_statistics':['dh_fit_dx'],
            'ground_track':['x_atc', 'sigma_geo_xt','sigma_geo_at'],
            'geophysical' : ['dac','tide_ocean'],
            'orbit_info':['rgt','cycle_number'],
            'derived':['valid','matlab_time', 'n_pixels','LR','BP','spot','rss_along_track_dh']}
    
    arg_dict={'SRS_proj4':SRS_proj4,'tile_spacing':tile_spacing, 'pad':pad, 'bin_W':bin_W, 'GI_file':GI_file,'out_dir':out_dir,'field_dict':field_dict,'cycle':cycle}
    arg_dict.update({'seg_diff_scale':seg_diff_scale, 'blockmedian_scale':blockmedian_scale})
    arg_dict.update({'xy0':xy})
    make_tile(arg_dict)

if __name__=='__main__':
    main()