#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:32:09 2017

@author: ken.pryor
"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import urllib

def read_BTD_plot(ncf):
    nc_fid = Dataset(ncf, 'r')
    data = nc_fid.variables["data"][:]  # shape lat, lon as shown above
    lats = nc_fid.variables['lat'][:]  # extract/copy the data
    lons = nc_fid.variables['lon'][:]
    np.savez('btd_wvir.npz', btd=data, lats=lats, lons=lons)
    names = nc_fid.variables.keys()
    nc_fid.close()
    return data, lats, lons, names

GOES_E_BTD_file, headers = urllib.request.urlretrieve('ftp://ftp.star.nesdis.noaa.gov/pub/smcd/opdb/mbimg/GOES0100.nc')
#GOES_E_BTD_file = 'GOES0100.nc'
btd_e, lats, lons, names = read_BTD_plot(GOES_E_BTD_file)
btd_e = btd_e[0,:,:]
print(names)
print('btd_e shape', btd_e.shape, btd_e)
btd_e_max = np.amax(btd_e)
btd_e_min = np.amin(btd_e)
print('btd_e min, btd_e max', btd_e_min, btd_e_max)
print('lats shape, lons shape', lats.shape, lons.shape)

# CREATE A MAP

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=25,urcrnrlat=45,llcrnrlon=-90,urcrnrlon=-70) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,btd_e,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
#m.drawparallels(np.arange(41.,46.,1.),labels=[True,False,False,False])
#m.drawmeridians(np.arange(-118,-108.,2.),labels=[False,False,False,True])
m.drawcoastlines(color="white")
m.drawcountries(color="white")
m.drawstates(color="white")
lon = -76.579
lat = 39.267
x, y = m(lon,lat)
print(x,y)
plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color='white')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('GOES-13 WV-IR BTD')
plt.savefig("G13_WVIR_BTD_Map.png",dpi=500,bbox_inches='tight')
plt.show()

GOES_W_BTD_file, headers = urllib.request.urlretrieve('ftp://ftp.star.nesdis.noaa.gov/pub/smcd/opdb/mbimg/GOES0110.nc')
#GOES_W_BTD_file = 'GOES0110.nc'
btd_w, lats, lons, names = read_BTD_plot(GOES_W_BTD_file)
btd_w = btd_w[0,:,:]
print(names)
print('btd_w shape', btd_w.shape, btd_w)
btd_w_max = np.amax(btd_w)
btd_w_min = np.amin(btd_w)
print('btd_w min, btd_w max', btd_w_min, btd_w_max)
print('lats shape, lons shape', lats.shape, lons.shape, lats, lons)

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=25,urcrnrlat=45,llcrnrlon=-125,urcrnrlon=-105) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,btd_w,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
#m.drawparallels(np.arange(41.,46.,1.),labels=[True,False,False,False])
#m.drawmeridians(np.arange(-118,-108.,2.),labels=[False,False,False,True])
m.drawcoastlines(color="white")
m.drawcountries(color="white")
lon1 = -112.65173
lat1 = 43.59413
x1, y1 = m(lon1,lat1)
print(x1,y1)
lon2 = -117.174
lat2 = 32.714
x2, y2 = m(lon2,lat2)
print(x2,y2)
plt.annotate('I', xy=(x1, y1), xycoords='data', xytext=(x1, y1), textcoords='data', color='white')
plt.annotate('SD', xy=(x2, y2), xycoords='data', xytext=(x2, y2), textcoords='data', color='white')
m.drawstates(color="white")
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('GOES-15 WV-IR BTD')
plt.savefig("G15_WVIR_BTD_Map.png",dpi=500,bbox_inches='tight')
plt.show()
