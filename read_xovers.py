#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 14:05:44 2019

@author: ben
"""
import glob
import re
import h5py
import sys
from PointDatabase.point_data import point_data
from PointDatabase.mapData import mapData
import numpy as np

def read_xovers(xover_base='/Volumes/ice2/ben/scf/AA_06/tiles', release='001', cycle=None, verbose=False, wildcard='*'):
    xover_dir=f'{xover_base}/{release}/cycle_%s/xovers' % cycle
    tiles=glob.glob(xover_dir+'/*.h5')
    with h5py.File(tiles[0],'r') as h5f:
        fields=[key for key in h5f['data_0'].keys()]
    
    D=[]
    meta={'slope_x':[], 'slope_y':[], 'grounded':[]}
    #X=[]    
    for tile in glob.glob(xover_dir+'/'+wildcard+'.h5'):
        try:
            with h5py.File(tile,'r') as h5f:
                for field in ['slope_x', 'slope_y','grounded']:
                    meta[field].append(np.array(h5f['/'+field]))                
        except KeyError:
            if verbose:
                print("failed to read " + tile)
            continue
        D.append([point_data(list_of_fields=fields).from_file(tile, field_dict={gr : fields}) for gr in ['data_0','data_1']])

    for field in meta.keys():
        meta[field]=np.concatenate(meta[field])
    v={}
    for field in fields:
        vi=[]
        for Di in D:
            vi.append(np.r_[[np.sum(getattr(Di[ii], field)*Di[ii].W, axis=0) for ii in [0, 1]]])
        v[field]=np.concatenate(vi, axis=1).T
    delta={field:np.diff(v[field], axis=1) for field in fields}
    bar={field:np.mean(v[field], axis=1) for field in fields}
    return v,  delta,  bar, meta
