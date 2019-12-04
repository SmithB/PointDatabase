# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:40:57 2019

@author: ben
"""

from PointDatabase.ATL06_tiles import read_tile
from PointDatabase.ATL06_filters import segDifferenceFilter
from PointDatabase.xover_search import cross_tracks
from PointDatabase import point_data, mapData

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import glob
import re
import sys

def ATL06_crossovers(xy0, tile_dir, different_cycles=False):
    D=read_tile(xy0, tile_dir)
    for Di in D:
        segDifferenceFilter(Di, tol=2, setValid=True, toNaN=True, subset=False)
        Di.assign({'time':Di.delta_time})
        Di.index(np.isfinite(Di.h_li))
    xover_list=list()
    #plt.clf()
    for ii in np.arange(len(D)):
        for jj in np.arange(len(D)):
            if (D[ii].size <2) or (D[jj].size < 2) or (ii>=jj) or (D[ii].rgt[0]==D[jj].rgt[0]):
                continue
            if different_cycles and D[ii].cycle[0]==D[jj].cycle[0]:
                continue
            xyC, inds, L=cross_tracks([D[ii], D[jj]], delta=20, delta_coarse=1000)
            if xyC is not None:
                try:
                    xover_list.append({'xyC':xyC, 'data_0':D[ii].subset(inds[0]), 'data_1':D[jj].subset(inds[1]), 'L0':L[0], 'L1':L[1]})
                except Exception as e:
                    print("HERE")
    return xover_list



def write_xovers(xover_list, xy0, out_dir):
    tile_name='/E%d_N%d.h5' % (xy0[0]/1.e3, xy0[1]/1.e3)
    out_file=out_dir+'/'+tile_name
    if os.path.isfile(out_file):
        os.remove(out_file)
    with h5py.File(out_file,'w') as h5f:
        for key_D in ['data_0', 'data_1']:
            group='/'+key_D
            h5f.create_group(group)

            key_L=key_D.replace('data_','L')
            L=np.c_[[item[key_L] for item in xover_list]]

            Dtemp=[item[key_D] for item in xover_list]
            Dtemp=point_data(list_of_fields=Dtemp[0].list_of_fields).from_list(Dtemp)
            shape=[np.int(Dtemp.size/2), 2]
            Dtemp.shape=shape
            for key in Dtemp.list_of_fields:
                temp=getattr(Dtemp, key)
                temp.shape=shape
                h5f.create_dataset(group+'/'+key, data=temp)
            h5f.create_dataset(group+'/W', data=np.c_[1-L, L])

        xy=np.c_[[item['xyC'] for item in xover_list]]
        h5f.create_dataset('/x', data=xy[:,0])
        h5f.create_dataset('/y', data=xy[:,1])
        h5f.create_dataset('/slope_x', data=np.array([item['slope_x'] for item in xover_list]))
        h5f.create_dataset('/slope_y', data=np.array([item['slope_y'] for item in xover_list]))
        h5f.create_dataset('/grounded', data=np.array([item['grounded'] for item in xover_list]))
    return #xover_list

def read_xovers(xover_dir):

    tiles=glob.glob(xover_dir+'/*.h5')
    with h5py.File(tiles[0],'r') as h5f:
        fields=[key for key in h5f['D0'].keys()]
    D=[]
    X=[]
    for tile in glob.glob(xover_dir+'/*.h5'):
        D.append([point_data(list_of_fields=fields).from_file(tile, field_dict={gr : fields}) for gr in ['D0','D1']])
        with h5py.open(tile,'r') as h5f:
            X.append(point_data(list_of_fields=['x','y']).from_file(tile,field_dict={None:['x','y']}))
    return D, X

def make_queue(top_dir, cycle, hemisphere, queue_file):
    out_dir=top_dir+'/xovers/cycle_'+cycle
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    with open(queue_file,'w') as qf:
        cycle_dirs=glob.glob(top_dir+'/cycle_'+cycle)
        for cycle_dir in cycle_dirs:
            files=glob.glob(cycle_dir+'/E*.h5')
            for file in files:
                if not os.path.isfile(out_dir+'/'+os.path.basename(file)):
                    qf.write('python3 ~/git_repos/PointDatabase/cross_ATL06_tile.py %s %d %s %s\n' %  ( cycle_dir, hemisphere, out_dir, os.path.basename(file)))

def calc_slope(xovers, hemisphere=-1):

    if hemisphere==-1:
        mask_file='/Volumes/ice1/ben/Quantarctica/Quantarctica2/Scientific/Glaciology/ALBMAP/ALBMAP_Mask.tif'
        dx=1.e4
    else:
        mask_file='/Volumes/ice2/ben/ATL14_test/GimpIceMask_100m_edited.tif'
        dx=200

    xy=np.c_[[item['xyC'] for item in xovers]]

    mask=mapData().from_geotif(mask_file, \
                bounds=[[np.min(xy[:,0].ravel())-dx, np.max(xy[:,0].ravel()+dx)], [np.min(xy[:,1].ravel())-dx, np.max(xy[:,1].ravel()+dx)]])

    grounded=np.abs(mask.interp(xy[:,0], xy[:,1])-1)<.01
    G=np.zeros((4,4))
    G[:,2]=np.array([1, 1, 0, 0])
    G[:,3]=np.array([0, 0, 1, 1])
    for ii, xo in enumerate(xovers):
        x1=np.r_[xo['data_0'].x, xo['data_1'].x]
        y1=np.r_[xo['data_0'].y, xo['data_1'].y]
        G[:,0]=x1-x1.mean()
        G[:,1]=y1-y1.mean()
        m=np.linalg.solve(G, np.r_[xo['data_0'].h_li, xo['data_1'].h_li])
        xo['slope_x'] = m[0]
        xo['slope_y'] = m[1]
        xo['grounded'] = grounded[ii]

def main():
    tile_dir=sys.argv[1]
    hemisphere=int(sys.argv[2])
    out_dir=sys.argv[3]
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if len(sys.argv)<4:
        files=glob.glob(tile_dir+'/*.h5')
    else:
        files=sys.argv[4:]

    for file in files:
        print(file)
        #if os.path.isfile(out_dir+'/'+os.path.basename(file)):
        #    continue
        g=re.compile('E(.*)_N(.*).h5').search(file)
        xy0=(np.int(g.group(1))*1000, np.int(g.group(2))*1000)
        print(xy0)
        xover_list=ATL06_crossovers(xy0, tile_dir)
        if len(xover_list) > 0:
            calc_slope(xover_list, hemisphere=hemisphere)
            write_xovers(xover_list, xy0, out_dir)

if __name__=='__main__':
    main()