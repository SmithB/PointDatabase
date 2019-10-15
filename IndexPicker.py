#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 18:08:00 2019

@author: ben
"""
from geo_index import geo_index
from geo_index import index_list_for_files
from PointDatabase import point_data, mapData
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import scipy.stats as sps

class IndexPicker:
    def __init__(self, fig, index_file, map_file=None, W=4e4):
        self.xy=geo_index().from_file(index_file).bins_as_array()
        self.index_file=index_file
        self.map_file=map_file
        self.fig=fig
        self.fig=plt.figure()
        self.W=4.e4
        if map_file is not None:
            
            backgrnd=mapData().from_geotif(map_file) 
            cr=sps.scoreatpercentile(backgrnd.z[np.isfinite(backgrnd.z)], [16, 84])
            hi=backgrnd.show(cmap='gray', vmin=cr[0]-np.diff(cr), vmax=cr[1]+np.diff(cr))
        
        plt.plot(self.xy[0], self.xy[1],'.')
        self.cid=self.fig.canvas.mpl_connect('button_press_event', self.click)
        plt.show(block=True)
        #self.fig.canvas.draw()
    
    def click(self, event):
        print('click', event)
        if not event.dblclick:
            return
        #xyp=event.xdata+1j*event.ydata
        #best = np.argmin(np.abs(self.xy[0]+1j*self.xy[1] - xyp))
        
        field_dict={None:['delta_time','h_li','h_li_sigma','latitude','longitude','atl06_quality_summary','segment_id','sigma_geo_h'], 
            'fit_statistics':['dh_fit_dx'],
            'ground_track':['x_atc', 'sigma_geo_xt','sigma_geo_at'],
            'geophysical' : ['dac','tide_ocean','r_eff'],
            'orbit_info':['rgt','cycle_number'],
            'derived':['valid','matlab_time','LR','BP','spot','rss_along_track_dh']}
        W=self.W
        self.fig.gca().plot(event.xdata, event.ydata,'m*')
        self.fig.canvas.draw()
        xy0=np.round(np.array([event.xdata, event.ydata])/1.e4)*1.e4
        print(xy0)
        D=geo_index().from_file(self.index_file).query_xy_box(\
                   event.xdata+np.array([-W/2, W/2]), event.ydata+np.array([-W/2, W/2]), fields=field_dict)

        #D=geo_index().from_file(self.index_file).query_xy_box((self.xy[0][best], self.xy[1][best]), fields=field_dict)
        TrackPicker(D, self.map_file, self.index_file)

class TrackPicker:
    def __init__(self, D, map_file, index_file):        
        self.files=[]
        srs_proj4=geo_index().from_file(index_file).attrs['SRS_proj4']
        
        for ii, Di in enumerate(D):
            Di.get_xy(srs_proj4)
            Di.assign({'file_num':np.zeros_like(Di.x)+ii})
            self.files += [Di.filename]
        self.D_all=point_data().from_list(D)
        self.D=D
        XR=[np.nanmin(self.D_all.x), np.nanmax(self.D_all.x)]
        YR=[np.nanmin(self.D_all.y), np.nanmax(self.D_all.y)]
        self.fig=plt.figure()
        if map_file is not None:
            mapData().from_geotif(map_file, bounds=[XR, YR]).show(cmap='gray')
        #for Di in D:           
        #    plt.plot(Di.x, Di.y,'.')
        plt.scatter(self.D_all.x, self.D_all.y, 6, c=self.D_all.r_eff, linewidth=0, vmax=1, vmin=0)
        self.cid=self.fig.canvas.mpl_connect('button_press_event', self)
        plt.show()

    def __call__(self, event):
        xyp=event.xdata+1j*event.ydata
        best = np.argmin(np.abs(self.D_all.x+1j*self.D_all.y - xyp))
        this = int(self.D_all.file_num[best])
        fig=plt.figure()
        #print(self.D[this].list_of_fields)
        plt.subplot(211)
        plt.plot(self.D[this].x_atc, self.D[this].h_li,'.')
        plt.title(self.files[this])
        plt.subplot(212)
        plt.plot(self.D[this].x_atc, self.D[this].r_eff,'.')
        plt.title(self.files[this])
        print(self.files[this])
        fig.canvas.draw()
        plt.show()

def main():
    #fig=plt.figure()
    fig=None
    IndexPicker(fig, sys.argv[1], sys.argv[2])
    #plt.show(block=True)
    
if __name__=='__main__':
    main()
        